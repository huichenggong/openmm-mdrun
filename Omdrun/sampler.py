from copy import deepcopy

import numpy as np
import math

from openmm import openmm, unit, app
import grand

class CustomGCNCMCSampler(grand.samplers.NonequilibriumGCMCSphereSampler):
    """
    Custom GCNCMC sampler.
    Move can be selected from insertion and deletion
    BoxVectors will also be reset after each move, incase we run NPT or semiisotropic simulations
    """
    def move(self, context, insert):
        """
        Carry out a nonequilibrium GCMC move

        Parameters
        ----------
        context : simtk.openmm.Context
            Current context of the simulation
        insert : bool,
            If True, attempt an insertion move.
            If False, attempt a deletion move.
        """
        # Read in positions
        self.context = context
        state = self.context.getState(getPositions=True, enforcePeriodicBox=True, getVelocities=True)
        self.positions = deepcopy(state.getPositions(asNumpy=True))
        self.velocities = deepcopy(state.getVelocities(asNumpy=True))
        self.box_vectors = state.getPeriodicBoxVectors(asNumpy=True)

        # Update GCMC region based on current state
        self.updateGCMCSphere(state)

        # Set to NCMC integrator
        self.compound_integrator.setCurrentIntegrator(1)

        #  Execute moves
        # Insert or delete a water, based on random choice
        if insert:
            # Attempt to insert a water
            self.logger.info("Insertion")
            success, work = self.insertionMove()
        else:
            # Attempt to delete a water
            self.logger.info("Deletion")
            success, work = self.deletionMove()
        self.n_moves += 1
        self.Ns.append(self.N)

        # Set to MD integrator
        self.compound_integrator.setCurrentIntegrator(0)

        return None

    def insertionMove(self):
        """
        Carry out a nonequilibrium insertion move for a random water molecule
        """
        # Store initial positions
        old_positions = deepcopy(self.positions)

        # Choose a random site in the sphere to insert a water
        new_positions, resid, atom_indices = self.insertRandomWater()

        # Need to update the context positions
        self.context.setPositions(new_positions)

        # Start running perturbation and propagation kernels
        protocol_work = 0.0 * unit.kilocalories_per_mole
        explosion = False
        self.ncmc_integrator.step(self.n_prop_steps_per_pert)
        for i in range(self.n_pert_steps):
            state = self.context.getState(getEnergy=True)
            energy_initial = state.getPotentialEnergy()
            # Adjust interactions of this water
            self.adjustSpecificWater(atom_indices, self.lambdas[i+1])
            state = self.context.getState(getEnergy=True)
            energy_final = state.getPotentialEnergy()
            protocol_work += energy_final - energy_initial
            # Propagate the system
            try:
                self.ncmc_integrator.step(self.n_prop_steps_per_pert)
            except:
                print("Caught explosion!")
                explosion = True
                self.n_explosions += 1
                break

        # Store the protocol work
        self.insert_works.append(protocol_work)

        # Update variables and GCMC sphere
        self.setWaterStatus(resid, 1)
        state = self.context.getState(getPositions=True, enforcePeriodicBox=True)
        self.positions = state.getPositions(asNumpy=True)
        self.updateGCMCSphere(state)

        # Check which waters are in the sphere
        wats_in_sphere = self.getWaterStatusResids(1)

        # Calculate acceptance probability
        if resid not in wats_in_sphere:
            # If the inserted water leaves the sphere, the move cannot be reversed and therefore cannot be accepted
            acc_prob = -1
            self.n_left_sphere += 1
            self.logger.info("Move rejected due to water leaving the GCMC sphere")
        elif explosion:
            acc_prob = -1
            self.logger.info("Move rejected due to an instability during integration")
        else:
            # Calculate acceptance probability based on protocol work
            acc_prob = math.exp(self.B) * math.exp(-protocol_work/self.kT) / self.N  # Here N is the new value

        self.acceptance_probabilities.append(acc_prob)

        # Update or reset the system, depending on whether the move is accepted or rejected
        if acc_prob < np.random.rand() or np.isnan(acc_prob):
            # Need to revert the changes made if the move is to be rejected
            self.adjustSpecificWater(atom_indices, 0.0)

            self.context.setPeriodicBoxVectors(*self.box_vectors)
            self.context.setPositions(old_positions)
            self.context.setVelocities(-self.velocities)  # Reverse velocities on rejection

            self.positions = deepcopy(old_positions)
            self.velocities = -self.velocities
            state = self.context.getState(getPositions=True, enforcePeriodicBox=True)
            self.setWaterStatus(resid, 0)
            self.updateGCMCSphere(state)
            return False, protocol_work
        else:
            # Update some variables if move is accepted
            self.N = len(wats_in_sphere)
            self.n_accepted += 1
            state = self.context.getState(getPositions=True, enforcePeriodicBox=True, getVelocities=True)
            self.positions = deepcopy(state.getPositions(asNumpy=True))
            self.velocities = deepcopy(state.getVelocities(asNumpy=True))
            self.updateGCMCSphere(state)
            return True, protocol_work

    def deletionMove(self):
        """
        Carry out a nonequilibrium deletion move for a random water molecule
        """
        # Store initial positions
        old_positions = deepcopy(self.positions)

        # Choose a random water in the sphere to be deleted
        resid, atom_indices = self.deleteRandomWater()
        # Deletion may not be possible
        if resid is None:
            return None

        # Start running perturbation and propagation kernels
        protocol_work = 0.0 * unit.kilocalories_per_mole
        explosion = False
        self.ncmc_integrator.step(self.n_prop_steps_per_pert)
        for i in range(self.n_pert_steps):
            state = self.context.getState(getEnergy=True)
            energy_initial = state.getPotentialEnergy()
            # Adjust interactions of this water
            self.adjustSpecificWater(atom_indices, self.lambdas[-(2+i)])
            state = self.context.getState(getEnergy=True)
            energy_final = state.getPotentialEnergy()
            protocol_work += energy_final - energy_initial
            # Propagate the system
            try:
                self.ncmc_integrator.step(self.n_prop_steps_per_pert)
            except:
                print("Caught explosion!")
                explosion = True
                self.n_explosions += 1
                break

        # Get the protocol work
        self.delete_works.append(protocol_work)

        # Update variables and GCMC sphere
        # Leaving the water as 'on' here to check that the deleted water doesn't leave
        state = self.context.getState(getPositions=True, enforcePeriodicBox=True)
        self.positions = state.getPositions(asNumpy=True)
        old_N = self.N
        self.updateGCMCSphere(state)

        # Check which waters are in the sphere
        wats_in_sphere = self.getWaterStatusResids(1)

        # Calculate acceptance probability
        if resid not in wats_in_sphere:
            # If the deleted water leaves the sphere, the move cannot be reversed and therefore cannot be accepted
            acc_prob = 0
            self.n_left_sphere += 1
            self.logger.info("Move rejected due to water leaving the GCMC sphere")
        elif explosion:
            acc_prob = 0
            self.logger.info("Move rejected due to an instability during integration")
        else:
            # Calculate acceptance probability based on protocol work
            acc_prob = old_N * math.exp(-self.B) * math.exp(-protocol_work/self.kT)  # N is the old value

        self.acceptance_probabilities.append(acc_prob)

        # Update or reset the system, depending on whether the move is accepted or rejected
        if acc_prob < np.random.rand() or np.isnan(acc_prob):
            # Need to revert the changes made if the move is to be rejected
            self.adjustSpecificWater(atom_indices, 1.0)

            self.context.setPeriodicBoxVectors(*self.box_vectors)
            self.context.setPositions(old_positions)
            self.context.setVelocities(-self.velocities)  # Reverse velocities on rejection

            self.positions = deepcopy(old_positions)
            self.velocities = -self.velocities
            state = self.context.getState(getPositions=True, enforcePeriodicBox=True)
            self.updateGCMCSphere(state)
            return False, protocol_work
        else:
            # Update some variables if move is accepted
            self.setWaterStatus(resid, 0)
            self.N = len(wats_in_sphere) - 1  # Accounting for the deleted water
            self.n_accepted += 1
            state = self.context.getState(getPositions=True, enforcePeriodicBox=True, getVelocities=True)
            self.positions = deepcopy(state.getPositions(asNumpy=True))
            self.velocities = deepcopy(state.getVelocities(asNumpy=True))
            self.updateGCMCSphere(state)
            return True, protocol_work
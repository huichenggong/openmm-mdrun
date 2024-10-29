import argparse
import sys
import os
import time
import shutil
from pathlib import Path
import gzip
import logging

from openmm import openmm, unit
from openmm import app

import Omdrun
from Omdrun.mdp import mdp_parser

def backup_if_exist_gmx(f_name):
    """
    Do gromacs style backup. If md.xtc exist, backup it to #md.xtc.1#
    """
    if os.path.exists(f_name):
        for i in range(1,10000):
            bak = f"#{f_name}.{i}#"
            if not os.path.exists(bak):
                shutil.move(f_name, bak)
                return bak
        raise Exception(f"Cannot backup {f_name}")
    else:
        return None

def load_sys(sys_file):
    """
    Load serialized system file
    :param sys_file: X.xml.gz or X.xml
    :return: openmm.System
    """
    # if sys_file is xml.gz
    if sys_file.endswith(".gz"):
        with gzip.open(sys_file, 'rt') as f:
            system = openmm.XmlSerializer.deserialize(f.read())
    else:
        with open(sys_file, 'r') as f:
            system = openmm.XmlSerializer.deserialize(f.read())
    return system

def set_output_file_name(args):
    """
    Set output file names
    :param args:
    :return:
    """
    # if -deffnm is set, use it as prefix for all non-None output files
    if args.deffnm:
        prefix = args.deffnm
        if args.xtc is None:
            args.xtc = f"{prefix}.xtc"
        if args.cpo is None:
            args.cpo = f"{prefix}.xml"
        if args.e is None:
            args.e = f"{prefix}.csv"
        if args.g is None:
            args.g = f"{prefix}.log"
        if args.pdb is None:
            args.pdb = f"{prefix}.pdb"
    else:
        if args.xtc is None:
            args.xtc = "traj.xtc"
        if args.cpo is None:
            args.cpo = "state.xml"
        if args.e is None:
            args.e = "energy.csv"
        if args.g is None:
            args.g = "md.log"
        if args.pdb is None:
            args.pdb = "final.pdb"

def set_barostat(system, mdp_inputs):
    """
    Add barostat to the system
    :param system:
    :param mdp_inputs:
    :return:
    """
    if mdp_inputs.pcoupltype == "isotropic":
        barostat = openmm.MonteCarloBarostat(mdp_inputs.ref_p[0], mdp_inputs.ref_t, mdp_inputs.nstpcouple)
        system.addForce(barostat)
    elif mdp_inputs.pcoupltype == "membrane":
        barostat = openmm.MonteCarloMembraneBarostat(mdp_inputs.ref_p[0], mdp_inputs.surface_tension,
                                                     mdp_inputs.ref_t,
                                                     openmm.MonteCarloMembraneBarostat.XYIsotropic,
                                                     openmm.MonteCarloMembraneBarostat.ZFree,
                                                     mdp_inputs.nstpcouple)
        system.addForce(barostat)
    elif mdp_inputs.pcoupltype == "anisotropic":
        default_pressure = 1.0 * unit.bar
        px_scale = mdp_inputs.ref_p[0] / default_pressure
        py_scale = mdp_inputs.ref_p[1] / default_pressure
        pz_scale = mdp_inputs.ref_p[2] / default_pressure
        barostat = openmm.MonteCarloAnisotropicBarostat(default_pressure, mdp_inputs.ref_t, px_scale, py_scale, pz_scale, mdp_inputs.nstpcouple)
        system.addForce(barostat)
    else:
        raise ValueError(f"pcoupltype {mdp_inputs.pcoupltype} is not supported")


def set_restraint(system, positions, mdp_inputs, res_index_file):
    """
    Add restraint force to the system
    :param system: openmm.System, will add force to this system
    :param positions: positions of all atoms
    :param mdp_inputs: mdp_parser object, it should have res_fc
    :return: None
    """
    force_con = mdp_inputs.res_fc
    posresPROT = openmm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
    posresPROT.addPerParticleParameter('k')
    posresPROT.addPerParticleParameter('x0')
    posresPROT.addPerParticleParameter('y0')
    posresPROT.addPerParticleParameter('z0')
    with open(res_index_file) as f:
        lines = f.readlines()
    for line in lines:
        segments = line.strip().split()
        atom1 = int(segments[0])
        state = segments[1]
        xpos = positions[atom1].value_in_unit(unit.nanometers)[0]
        ypos = positions[atom1].value_in_unit(unit.nanometers)[1]
        zpos = positions[atom1].value_in_unit(unit.nanometers)[2]
        if state == 'BB' or state == "SC":  # BackBone and SideChain
            posresPROT.addParticle(atom1, [force_con, xpos, ypos, zpos])
    system.addForce(posresPROT)

class StopSimulation(Exception):
    """ Exception raised in order to stop simulation """
    pass

class TimeUpReporter(object):
    """ Exits a Simulation by raising StopSimulation if an exit file exists """
    def __init__(self, t0, maxh, reportInterval):
        """
        Create a new ExitFileReporter.
        :param t0: The time at which the timer starts in seconds
        :param maxh: The maximum number of hours to run
        :param reportInterval: The interval (in steps) at which to report progress
        """
        self.t0 = t0
        self.maxh = maxh
        self._reportInterval = reportInterval


    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval # The number of steps until the next report
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        """
        check time and raise StopSimulation if time is up
        """
        t = time.time()
        time_remaining = self.maxh * 3600 - (t - self.t0)
        step = simulation.currentStep

        if time_remaining > 3600:
            logging.debug(f"Step {step}, Time remaining: {time_remaining // 3600:.0f} h {time_remaining % 3600 // 60} min")
        elif time_remaining > 60:
            logging.debug(f"Step {step}, Time remaining: {time_remaining // 60:.0f} min {time_remaining % 60 :.0f} s")
        else:
            logging.debug(f"Step {step}, Time remaining: {time_remaining:.1f} s")

        if time_remaining < 0:
            logging.info(f"Time is up, stopping simulation")
            raise StopSimulation()




def print_inp(args, mdp_inputs):
    """

    """
    command_line = ""
    for word in sys.argv:
        if " " in word or "%" in word:
            command_line += f'"{word}" '
        else:
            command_line += word + " "
    logging.info(f"Command line: {command_line}")
    logging.info(f"openmm_mdrun version {Omdrun.__version__}")
    for key, value in args.__dict__.items():
        logging.info(f"-{key:7}: {value}")

    for key, value in mdp_inputs.__dict__.items():
        logging.info(f"{key:18} = {value}")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sys", metavar="   sys.xml.gz  ", default="sys.xml.gz",
                        help="Serialized system file. It will be loaded as a `openmm.System` object. This system should "
                             "include all bonded/non-bonded, constraints, but not pressure coupling. It can be xml or xml.gz.")
    parser.add_argument("-p", metavar="     top.psf/top.parm7", default="top.psf",
                        help="Input topology file")
    parser.add_argument("-t", metavar="     restart.xml ",
                        help="Restart coordinates/velocities/box from this file and reset the simulation time")
    parser.add_argument("-cpt", metavar="   state.cpt   ",
                        help="Restart coordinates/velocities/box from this file and continue the simulation")
    parser.add_argument("-mdp", metavar="   md.mdp      ", default="md.mdp",
                        help="Input file with MD parameters. Only limitied gmx mdp keywords are supported.")
    parser.add_argument("-r", metavar="     restrain.xml",
                        help="Restrain target")
    parser.add_argument("-rind", metavar="  index.dat   ",
                        help="Restrain index")
    parser.add_argument("-maxh", metavar="   23.8 hour  ", default=23.8, type=float,
                        help="Maximal number of hours to run")
    parser.add_argument("-deffnm", metavar="md          ",
                        help="The default filename for all output")
    parser.add_argument("-xtc", metavar="     traj.xtc    ",
                        help="xtc trajectory")
    parser.add_argument("-cpo", metavar="   state.xml   ",
                        help="Checkpoint file. Serialized `openmm.State` object")
    parser.add_argument("-e", metavar="     energy.csv  ",
                        help="csv log file")
    parser.add_argument("-g", metavar="     md.log      ",
                        help="Log file")
    parser.add_argument("-pdb", metavar="   final.pdb   ",
                        help="pdb output for the final frame. After the simulation properly finished, the final frame will be saved to this file")
    parser.add_argument("--debug", action="store_true",
                        help="Print debug information")

    time_start = time.time()

    args = parser.parse_args()
    set_output_file_name(args)

    if args.cpt and (args.t is None):    restart_xml, continuation = args.cpt, True
    elif args.t and (args.cpt is None):  restart_xml, continuation = args.t,   False
    else:
        raise ValueError("Only one of -t or -cpt should be set")
    info_list = []
    if not continuation:
        for f_name in [args.xtc, args.e, args.g, args.cpo]:
            f_bak = backup_if_exist_gmx(f_name)
            if f_bak:
                info_list.append(f"Backup {f_name} to {f_bak}")

    # set INFO level logging to args.g
    if args.debug:
        logging.basicConfig(filename=args.g, level=logging.DEBUG, filemode = 'a',
                            format='%(asctime)s - %(levelname)s - %(message)s'
                            )
    else:
        logging.basicConfig(filename=args.g, level=logging.INFO, filemode = 'a',
                            format='%(message)s'
                            )
    for info in info_list:
        logging.info(info)


    # Load mdp file
    mdp_inputs = mdp_parser().read(args.mdp)
    print_inp(args, mdp_inputs)


    # Load system
    system = load_sys(args.sys)
    logging.info(f"System loaded from {args.sys}")
    logging.info(f"Number of particles: {system.getNumParticles()}")

    # Load state
    with open(restart_xml, 'r') as f:
        s_restart = openmm.XmlSerializer.deserialize(f.read())


    # Load topology
    if args.p.endswith(".psf"):
        psf = app.CharmmPsfFile(args.p, periodicBoxVectors=s_restart.getPeriodicBoxVectors())
        topology = psf.topology
    elif args.p.endswith(".parm7"):
        prmtop = app.AmberPrmtopFile(args.p)
        topology = prmtop.topology
    else:
        raise ValueError(f"Topology file {args.p} is not supported")

    # Set restraints
    if mdp_inputs.restraint:
        with open(args.r, 'r') as f:
            s_restraints = openmm.XmlSerializer.deserialize(f.read())
        set_restraint(system, s_restraints.getPositions(), mdp_inputs, args.rind)

    # Set barostat
    if mdp_inputs.pcoupltype:
        set_barostat(system, mdp_inputs)

    # Set context
    integrator = mdp_inputs.integrator(mdp_inputs.ref_t, mdp_inputs.tau_t, mdp_inputs.dt)
    sim = app.Simulation(topology, system, integrator)
    sim.context.setState(s_restart)
    if not continuation:
        sim.context.setStepCount(0)
        sim.context.setTime(0)
        # Set gen_vel
        if mdp_inputs.gen_vel:
            logging.debug(f"Set random velocities to temperature {mdp_inputs.gen_temp}")
            sim.context.setVelocitiesToTemperature(mdp_inputs.gen_temp)

    # Set reporters, if continuation, xtc and csv will be appended
    if mdp_inputs.nstxout_compressed > 0:
        logging.debug(f"Set xtc reporter to {args.xtc}")
        sim.reporters.append(app.XTCReporter(args.xtc, mdp_inputs.nstxout_compressed, continuation))
    logging.debug(f"Set csv reporter to {args.e}")
    sim.reporters.append(
        app.StateDataReporter(
            args.e, mdp_inputs.nstlog, step=True, time=True, potentialEnergy=True, temperature=True,
            volume=True, density=True,
            progress=True,
            remainingTime=True, speed=True, totalSteps=mdp_inputs.nsteps, separator='\t,', append=continuation))

    logging.debug(f"Set log reporter to {args.g}")
    sim.reporters.append(
        app.CheckpointReporter(args.cpo, mdp_inputs.nstlog, writeState=True))
    logging.debug(f"Set time up reporter")
    sim.reporters.append(
        TimeUpReporter(time_start, args.maxh, mdp_inputs.nstmaxh))

    # Run simulation
    if continuation:
        nsteps = mdp_inputs.nsteps - sim.context.getState().getStepCount()
    else:
        nsteps = mdp_inputs.nsteps

    # print DEBUG info

    for i, f in enumerate(system.getForces()):
        logging.debug(f"Force {i}: {f}")
    logging.debug(f"Time {sim.context.getTime()}")
    logging.debug(f"Step {sim.context.getStepCount()}")
    logging.debug(f"Run further simulation by {nsteps} steps")

    if nsteps > 0:
        try:
            sim.step(nsteps)
            sim.saveState(args.cpo)
            state = sim.context.getState(getPositions=True)
            with open(args.pdb, 'w') as f:
                app.PDBFile.writeFile(sim.topology, state.getPositions(), f)
            logging.info("Simulation finished")
        except StopSimulation:
            logging.info(f"Running time exceeded maxh {args.maxh}. The last step is {sim.currentStep}")
            sim.saveState(args.cpo)
            logging.info(f"State saved to {args.cpo}")
    else:
        logging.info(f"No more steps to run")


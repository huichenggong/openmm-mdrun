import argparse
import time
import sys
from pathlib import Path
import logging
import json

from openmm import openmm, unit, app
from openmmtools.integrators import BAOABIntegrator, NonequilibriumLangevinIntegrator

try:
    import grand
except ImportError:
    raise ImportError("This script requires grand. Please install it from : https://github.com/essex-lab/grand")

import Omdrun
from Omdrun.mdp import mdp_parser
from Omdrun.utils import backup_if_exist_gmx, load_sys, load_top, print_inp

def read_moves(filename):
    n_completed = 0
    n_accepted = 0
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines[-1::-1]:
        if " move(s) completed (" in line:
            n_completed = int(lines[-1].split()[4])
            n_accepted = int(lines[-1].split()[7].strip('('))
            break
    return n_completed, n_accepted

def check_continuation(args):
    """
    Check the input arguments and determine if the simulation is a continuation or a fresh start.
    if rst is provided,
        and it exists, then it is a continuation.
        and it does not exist, then it is a fresh start. state file should exist.
    if rst is not provided, then it is a fresh start. state file should exist.
    """
    if args.rst is not None:
        if Path(args.rst).is_file():
            logging.info(f"Load/continue simulation from {args.rst}")
            return True

    if args.t is not None and Path(args.t).is_file():
        logging.info(f"Arg -rst is not provided. Start a new simulation from {args.t}. Random velocities will be generated.")
        return False
    else:
        raise FileNotFoundError(f"Arg -rst is not provided and -t is not provided or not exist.")


def main():
    time_start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-sys", metavar="     sys.xml.gz  ", default="sys.xml.gz",
                        help="Serialized system file. It will be loaded as a `openmm.System` object. This system should "
                             "include all bonded/non-bonded, constraints, but not pressure coupling. It can be xml or xml.gz.")
    parser.add_argument("-p", metavar="       top.psf/top.parm7", default="top.psf",
                        help="Input topology file")
    parser.add_argument("-mdp", metavar="     md.mdp      ", default="md.mdp",
                        help="Input file with MD parameters. Only limitied gmx mdp keywords are supported.")
    parser.add_argument("-t", metavar="       start.rst7  ",
                        help="Initial coordinates/velocities/box file")
    parser.add_argument("-rst", metavar="     md_rst.rst7 ",
                        help="Restart input file, amber rst7 format. If rst is provided, the simulation will continue from the rst file.")
    parser.add_argument("-ighosts", metavar=" ghosts-0.txt", default="ghosts-0.txt",
                        help="Input ghost water file for starting or continuing")
    parser.add_argument("-ilog", metavar="    md_gcmc.log     ",
                        help="gcmc log file from the previous run")
    parser.add_argument("-atom", metavar="    atom.json  ", default="atom.json",
                        help="Input atom json file for atom selection")
    parser.add_argument("-maxh", metavar="    23.8 hour   ", default=23.8, type=float,
                        help="Maximal number of hours to run. Time will only be checked after each cycle.")
    parser.add_argument("-deffnm", metavar="  md          ", default="md",
                        help="The default filename for all output.")
    parser.add_argument("-zfree", action="store_true",
                        help="If provided, the z dimension is free in the membrane barostat")
    parser.add_argument("--debug", action="store_true",
                        help="Print debug information.")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.info(f"openmm_mdrun, version {Omdrun.__version__}")
    logging.info(f"Command line: {' '.join(sys.argv)}")

    # determine continue or fresh start
    continuation = check_continuation(args)
    if continuation:
        rst_in = app.AmberInpcrdFile(args.rst)
        n_moves, n_accepted = read_moves(args.ilog)
    else:
        rst_in = app.AmberInpcrdFile(args.rst)
        n_moves, n_accepted = 0, 0

    mdp_inputs = mdp_parser().read(args.mdp)
    print_inp(args, mdp_inputs)

    output_name_dict = {"gcmc_log": f"{args.deffnm}.log",
                        "dcd":      f"{args.deffnm}.dcd",
                        "rst":      f"{args.deffnm}.rst7",
                        "ghosts":   f"{args.deffnm}.dat",
                        "energy":   f"{args.deffnm}.csv"
                        }

    # use deffnm for all output files
    for f_key in ["dcd", "gcmc_log", "ghosts", "energy", "rst"]:
        fname = backup_if_exist_gmx(output_name_dict[f_key])
        if fname is not None:
            logging.info(f"Backup {output_name_dict[f_key]} to {fname}")
    # print output filenames
    for f_key, f_name in output_name_dict.items():
        logging.info(f"Output {f_key:8} : {f_name}")

    # load atom selection
    with open(args.atom, 'r') as f:
        atoms_region = json.load(f)
    atoms_region['radius'] *= unit.nanometer

    # load topology
    topology = load_top(args.p)
    topology.setPeriodicBoxVectors(rst_in.getBoxVectors())

    # load system
    system = load_sys(args.sys)

    # read ghost water
    ghosts = grand.utils.read_ghosts_from_file(args.ighosts)[-1]
    logging.info(f"Read {len(ghosts)} ghost waters from {args.ighosts}")
    logging.info(f"Ghost water res index: {ghosts}")
    for res in topology.residues():
        if res.index in ghosts:
            logging.info(f"Residue {res.index}:{res.name} has {len(list(res.atoms()))} atoms")

    integrator = BAOABIntegrator(mdp_inputs.ref_t, 1 / mdp_inputs.tau_t, mdp_inputs.dt)
    # Define the NCMC Sampler
    gcncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(
        system=system,
        topology=topology,
        temperature=mdp_inputs.ref_t,
        timeStep=mdp_inputs.dt,
        integrator=integrator,
        nPertSteps=399,  # number of perturbation steps (Hamiltonian switching)
        nPropStepsPerPert=50,  # number of propagation steps per perturbation step (constant Hamiltonian, relaxation)
        excessChemicalPotential=mdp_inputs.ex_potential,
        standardVolume=mdp_inputs.standard_volume,
        referenceAtoms=atoms_region["ref_atoms"],
        sphereRadius=atoms_region["radius"],
        log=output_name_dict["gcmc_log"],
        dcd=output_name_dict["dcd"],
        rst=output_name_dict["rst"],
        ghostFile=output_name_dict["ghosts"],
        overwrite=False
    )
    if mdp_inputs.pcoupltype == "membrane":
        if args.zfree:
            zmode = openmm.MonteCarloMembraneBarostat.ZFree
        else:
            zmode = openmm.MonteCarloMembraneBarostat.ZFixed
        barostat = openmm.MonteCarloMembraneBarostat(mdp_inputs.ref_p[0], mdp_inputs.surface_tension,
                                                     mdp_inputs.ref_t,
                                                     openmm.MonteCarloMembraneBarostat.XYIsotropic,
                                                     zmode,
                                                     mdp_inputs.nstpcouple)
        system.addForce(barostat)

    # set up the simulation(context), all information will be sent to GPU from here
    sim = app.Simulation(topology, system, gcncmc_mover.compound_integrator)

    sim.reporters.append(
        app.StateDataReporter(
            output_name_dict["energy"],
            1000,
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            density=True,  # this density is not so useful because it cannot ignore ghost water
            volume=True,
            separator='\t,', append=False))

    sim.context.setPositions(rst_in.getPositions())
    sim.context.setVelocities(rst_in.getVelocities())
    sim.context.setPeriodicBoxVectors(*rst_in.getBoxVectors())

    gcncmc_mover.initialise(sim.context, ghosts)
    gcncmc_mover.n_moves = n_moves
    gcncmc_mover.n_accepted = n_accepted

    logging.info(f"Number of HOH residues   : {len(gcncmc_mover.getWaterResids("HOH"))}")
    logging.info(f"Number of GCMC  waters   : {len([i for i, s in gcncmc_mover.water_status.items() if s == 1])}")
    logging.info(f"Number of outside  waters: {len([i for i, s in gcncmc_mover.water_status.items() if s == 2])}")
    logging.info(f"Number of ghost waters   : {len([i for i, s in gcncmc_mover.water_status.items() if s == 0])}")
    logging.info(f"{[i for i, s in gcncmc_mover.water_status.items() if s == 0]}")

    while gcncmc_mover.n_moves < mdp_inputs.ncycle:
        time_left = args.maxh*3600 - (time.time() - time_start)
        logging.info(f"### {time_left // 3600 :.0f} h {time_left/60 % 60 :.1f} min left ################################")
        logging.info(f"Cycle {gcncmc_mover.n_moves}, MD")
        sim.step(mdp_inputs.nsteps) # 500=1ps, 1000=2ps, 5000=10ps, 25000=50ps
        logging.info(f"Cycle {gcncmc_mover.n_moves}, MC")
        gcncmc_mover.move(sim.context, 1) # insertion and deletion will be randomlly chosen
        gcncmc_mover.report(sim)

        # check time/maxh
        time_passed = (time.time() - time_start)
        if time_passed > args.maxh*3600:
            logging.info(f"Time limit reached: {time_passed / 3600 :.2f} hours")
            break
    logging.info(f"Insertion work:")
    for w in gcncmc_mover.insert_works:
        logging.info(f"{w}")
    logging.info(f"Deletion  work:")
    for w in gcncmc_mover.delete_works:
        logging.info(f"{w}")
    time_passed = (time.time() - time_start)
    logging.info(f"GCNCMC script finished in {time_passed / 3600 :.0f} h {time_passed/60 % 60 :.0f} min {time_passed% 60 :.1f} s")

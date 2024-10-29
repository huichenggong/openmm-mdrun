import argparse
import sys
import os
import shutil
from pathlib import Path
import gzip
import logging

from openmm import openmm, unit
from openmm import app
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
        with gzip.open(sys_file, 'rb') as f:
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
        if args.o is None:
            args.o = f"{prefix}.dcd"
        if args.cpo is None:
            args.cpo = f"{prefix}.xml"
        if args.e is None:
            args.e = f"{prefix}.csv"
        if args.g is None:
            args.g = f"{prefix}.log"
    else:
        if args.o is None:
            args.o = "traj.dcd"
        if args.cpo is None:
            args.cpo = "state.xml"
        if args.e is None:
            args.e = "energy.csv"
        if args.g is None:
            args.g = "md.log"

def set_barostat(system, mdp_inputs):
    """
    Add barostat to the system
    :param system:
    :param mdp_inputs:
    :return:
    """
    if mdp_inputs.pcoup_type == "isotropic":
        barostat = openmm.MonteCarloBarostat(mdp_inputs.ref_p[0], mdp_inputs.ref_t, mdp_inputs.tau_p)
        system.addForce(barostat)
    elif mdp_inputs.pcoup_type == "membrane":
        barostat = openmm.MonteCarloMembraneBarostat(mdp_inputs.ref_p[0], mdp_inputs.surface_tension,
                                                     mdp_inputs.ref_t,
                                                     openmm.MonteCarloMembraneBarostat.XYIsotropic,
                                                     openmm.MonteCarloMembraneBarostat.ZFree,
                                                     mdp_inputs.tau_p)
    elif mdp_inputs.pcoup_type == "anisotropic":
        raise NotImplementedError("Anisotropic barostat is not implemented yet")
        system.addForce(barostat)

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



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sys", metavar="   sys.xml.gz  ", default="sys.xml.gz",
                        help="Serialized system file. It will be loaded as a `openmm.System` object. This system should "
                             "include all bonded/non-bonded, constraints, but not pressure coupling. It can be xml or xml.gz.")
    parser.add_argument("-p", metavar="     top.psf     ", default="top.psf",
                        help="Input topology file")
    parser.add_argument("-t", metavar="     restart.xml ", default="restart.xml",
                        help="Restart coordinates/velocities/box from this file")
    parser.add_argument("-mdp", metavar="   md.mdp      ", default="md.mdp",
                        help="Input file with MD parameters. Only limitied gmx mdp keywords are supported.")
    parser.add_argument("-r", metavar="     restrain.xml",
                        help="Restrain target")
    parser.add_argument("-rind", metavar="  index.dat   ",
                        help="Restrain index")
    parser.add_argument("-deffnm", metavar="md          ",
                        help="The default filename for all output")
    parser.add_argument("-o", metavar="     traj.dcd    ",
                        help="DCD trajectory")
    parser.add_argument("-cpo", metavar="   state.xml   ",
                        help="Checkpoint file. Serialized `openmm.State` object")
    parser.add_argument("-e", metavar="     energy.csv  ",
                        help="csv log file")
    parser.add_argument("-g", metavar="     md.log      ",
                        help="Log file")

    args = parser.parse_args()
    set_output_file_name(args)
    for f_name in [args.o, args.cpo, args.e, args.g]:
        backup_if_exist_gmx(f_name)

    # set up logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)  # Set the logger to the lowest level (DEBUG) to capture all levels

    file_handler = logging.FileHandler(args.g)
    file_handler.setLevel(logging.INFO)  # Log INFO and above to the file
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.ERROR)  # Log ERROR and above to the console

    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # Load system
    system = load_sys(args.sys)

    # Load topology and state
    with open(args.t, 'r') as f:
        s_restart = openmm.XmlSerializer.deserialize(f.read())
    psf = app.CharmmPsfFile(args.p, periodicBoxVectors=s_restart.getPeriodicBoxVectors())

    # Load mdp file
    mdp_inputs = mdp_parser().read(args.mdp)

    # Set restraints
    if args.r:
        with open(args.r, 'r') as f:
            s_restraints = openmm.XmlSerializer.deserialize(f.read())
        set_restraint(system, s_restraints.getPositions(), mdp_inputs, args.rind)

    # Set barostat
    if mdp_inputs.pcoupltype:
        set_barostat(system, mdp_inputs)

    # Set context
    integrator = mdp_inputs.integrator(mdp_inputs.ref_t, mdp_inputs.tau_t, mdp_inputs.dt)
    sim = app.Simulation(psf.topology, system, integrator)
    sim.context.setState(s_restart)
    sim.context.setStepCount(0)
    sim.context.setTime(0)

    # Set gen_vel
    if mdp_inputs.gen_vel:
        sim.context.setVelocitiesToTemperature(mdp_inputs.gen_temp)

    # Set reporters
    sim.reporters.append(app.DCDReporter(args.o, mdp_inputs.nstdcd))
    sim.reporters.append(
        app.StateDataReporter(
            args.e, mdp_inputs.nstlog, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=mdp_inputs.nsteps, separator='\t,'))

    # Run simulation
    sim.step(mdp_inputs.nsteps)

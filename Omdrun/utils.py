import os
import sys
import logging
import gzip
import shutil

from openmm import openmm, unit, app
import Omdrun


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

def load_top(top_file):
    """
    Load topology file. It can be psf(Charmm) or parm7/prmtop(Amber).
    return: openmm.Topology
    remember to run topology.setPeriodicBoxVectors if you provide a charmm psf file. If you don't set
    the box for topology, later trajectory will not have box information.
    """
    if top_file.endswith(".psf"):
        psf = app.CharmmPsfFile(top_file)
        topology = psf.topology
    elif top_file.endswith(".parm7") or top_file.endswith(".prmtop"):
        prmtop = app.AmberPrmtopFile(top_file)
        topology = prmtop.topology
    else:
        raise ValueError(f"Topology file {top_file} is not supported")
    return topology

def print_inp(args, mdp_inputs):
    """
    :param args      : argparse.ArgumentParser().parse_args()
    :param mdp_inputs: mdp_parser().read("md.mdp")
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


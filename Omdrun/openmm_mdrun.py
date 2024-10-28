import argparse
import sys
from pathlib import Path
import gzip

from openmm import openmm, unit
from openmm import app

# TODO: mdp reader/parser


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sys", metavar="   sys.xml.gz  ", default="sys.xml.gz",
                        help="Serialized system file. It will be loaded as a `openmm.System` object. This system should "
                             "include all bonded/non-bonded, but not pressure coupling. It can be xml or xml.gz.")
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
    parser.add_argument("-o", metavar="     traj.dcd    ", default="traj.dcd",
                        help="DCD trajectory")
    parser.add_argument("-cpo", metavar="   state.xml   ", default="state.xml",
                        help="Checkpoint file. Serialized `openmm.State` object")
    parser.add_argument("-e", metavar="     energy.csv  ", default="energy.csv",
                        help="csv log file")
    parser.add_argument("-g", metavar="     md.log      ", default="md.log",
                        help="Log file")

    args = parser.parse_args()
    print(args.p)
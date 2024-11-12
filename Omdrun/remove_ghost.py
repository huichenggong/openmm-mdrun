import sys
import argparse
from copy import deepcopy
import json

import numpy as np
from tqdm import tqdm

import grand
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import Omdrun
from Omdrun.utils import load_top, backup_if_exist_gmx

def main():
    parser = argparse.ArgumentParser(description=
                 f"""Omdrun {Omdrun.__version__}. This script can shift the ghost water out side the simulation box, and
                     save the processed trajectory to a xtc. Automatic centering is still under development.""", )
    parser.add_argument("-p", metavar="       top.psf/top.parm7", default="top.psf",
                        help="Input topology file")
    parser.add_argument("-idcd", metavar="    mdX.dcd", nargs='+',
                        help="Input dcd files")
    parser.add_argument("-idat", metavar="    mdX.dat", nargs='+',
                        help="Input ghosts.dat files")
    parser.add_argument("-atom", metavar="    atom.json  ",
                        help="Input atom json file for atom selection. If provided, we will try to center the channel.")
    parser.add_argument("-dframe", metavar="  1", type=int, default=1,
                        help="Down sample the trajectory by dframe")
    parser.add_argument("-oxtc", metavar="    md.xtc", default="md.xtc",
                        help="Output xtc files, the ghost waters will be shifted outside the simulation cell.")
    parser.add_argument("-pbc", metavar="     atom", default="atom", choices=["atom", "mol"],
                        help="wrap with atom or mol. mol is very slow")
    args = parser.parse_args()


    print(f"Omdrun Version: {Omdrun.__version__}")
    print(f"Command: {" ".join(sys.argv)}")

    top = load_top(args.p)

    ghost_list = []
    for ghost_f, dcd in zip(args.idat, args.idcd):
        print(f"Load ghost {ghost_f}, and trajectory {dcd}")
        ghosts = grand.utils.read_ghosts_from_file(ghost_f)
        ghost_list.extend(ghosts)


    print(f"Length of the ghosts list :{len(ghost_list)}")

    f_name = backup_if_exist_gmx(args.oxtc)
    if f_name:
        print(f"Backup {args.oxtc} to {f_name}")

    u = mda.Universe(top, args.idcd)
    if args.atom:
        with open(args.atom, 'r') as f:
            atom_dict = json.load(f)["ref_atoms"][0]
        selection_str = ""
        if "chain" in atom_dict:
            selection_str += f"segid {atom_dict['chain']} and "
        if "resname" in atom_dict:
            selection_str += f"resname {atom_dict['resname']} and "
        if "resid" in atom_dict:
            res = [r for r in top.residues() if r.id == atom_dict['resid'] ] [0]
            selection_str += f"resnum {res.index+1} and "
        if "name" in atom_dict:
            selection_str += f"name {atom_dict['name']} and "
        selection_str = selection_str[:-4]
        center = u.select_atoms(selection_str)
        print(f"Centering the channel with {center}")
    not_protein = u.select_atoms('not protein')
    protein = u.select_atoms('protein')
    all_atoms = u.select_atoms("name *")

    res_list = list(top.residues())
    with mda.Writer(args.oxtc, u.atoms.n_atoms) as xtc_writer:
        for ts, ghosts in tqdm(zip(u.trajectory[::args.dframe], ghost_list[::args.dframe])):
            if args.atom:
                dim = ts.triclinic_dimensions
                box_center = np.sum(dim, axis=0) / 2
                u.atoms.translate(box_center - center.center_of_mass())
                if args.pbc == "mol":
                    protein.wrap(compound='segments')
                    not_protein.wrap(compound='residues')
                elif args.pbc == "atom":
                    all_atoms.wrap(compound="atoms")
            z = -20
            for resid in ghosts:
                res = res_list[resid]
                at_list = list(res.atoms())
                shift_array = np.array([-5, -5, z*2]) - ts.positions[at_list[0].index]
                for atom in res.atoms():
                    ts.positions[atom.index] += shift_array
                z += 1
            xtc_writer.write(u.atoms)

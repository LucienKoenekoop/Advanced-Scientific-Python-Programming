import residue
import argparse
import numpy as np
import pandas as pd

"""
Module that makes hybrid .pdb, .lib and .prm files of two given residues with a maximum common substructure
"""

parser = argparse.ArgumentParser(description="""
Command-line tool for creating a maximum common substructure hybrid residue 
""")
reqarg = parser.add_argument_group('Required')
reqarg.add_argument('res1', help='residue 1')
reqarg.add_argument('res2', help='residue 2')

args = parser.parse_args()


class Hybrid(object):
    def __init__(self, res1, res2):
        self.res1 = residue.Residue(res1)
        self.res2 = residue.Residue(res2)
        self.res1.get_atoms(self.res1.lib)
        self.res2.get_atoms(self.res2.lib)
        self.res1.get_pdb(self.res1.pdb_file)
        self.res2.get_pdb(self.res2.pdb_file)

    def get_hybrid_atoms(self):
        check = self.res1.atoms.isin({'atoms': self.res2.atoms['atoms']})['atoms']
        mcs = self.res1.atoms[check == True]
        res1_rest = self.res1.atoms[check == False]
        res2_rest = self.res2.atoms[check == False]
        res2_rest['atoms'] = res2_rest['atoms'].str.lower()
        self.hybrid_atoms = mcs.append([res1_rest, res2_rest], ignore_index=True)
        # print(self.hybrid_atoms)
        return self.hybrid_atoms

    def make_hybrid_pdb(self, pdb1, pdb2, hybrid_atoms):
        res1_check = pdb1.isin({'atom_name': hybrid_atoms['atoms']})['atom_name']
        swap_atoms = hybrid_atoms['atoms'].str.swapcase()
        compare_atoms = pdb2['atom_name']
        res2_check = compare_atoms.isin(swap_atoms)
        atoms_res1 = pdb1[res1_check == True]
        atoms_res2 = pdb2[res2_check == True]
        atoms_res2['atom_name'] = atoms_res2['atom_name'].str.lower()
        hybrid_pdb = atoms_res1.append(atoms_res2, ignore_index=True)
        hybrid_pdb = hybrid_pdb.assign(atom_num=[x for x in range(1, len(hybrid_pdb)+1)], res_name='HYB')

        # Loop over rows of hybrid_pdb and parse out each line to pdb format
        # save result as hybrid.pdb

test = Hybrid(res1=args.res1, res2=args.res2)
test.get_hybrid_atoms()
test.make_hybrid_pdb(test.res1.pdb, test.res2.pdb, test.hybrid_atoms)
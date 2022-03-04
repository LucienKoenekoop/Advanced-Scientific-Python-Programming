import residue
import argparse
from biopandas.pdb import PandasPdb

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
        self.res1.get_atom_types(self.res1.prm)
        self.res2.get_atoms(self.res2.lib)
        self.res2.get_atom_types(self.res2.prm)
        self.res1.get_pdb(self.res1.pdb_file)
        self.res2.get_pdb(self.res2.pdb_file)

    def get_hybrid_atoms(self):
        """
        Function to find the atoms of the hybrid residue, returns these as a DataFrame
        """
        check = self.res1.atoms.isin({'atoms': self.res2.atoms['atoms']})['atoms']
        mcs = self.res1.atoms[check == True]
        res1_rest = self.res1.atoms[check == False]
        check, other = check.align(self.res2.atoms, axis=0)
        res2_rest = self.res2.atoms[check == False]
        res2_rest['atoms'] = res2_rest['atoms'].str.lower()
        self.hybrid_atoms = mcs.append([res1_rest, res2_rest], ignore_index=True)
        return self.hybrid_atoms

    def make_hybrid_pdb(self, pdb1, pdb2, hybrid_atoms):
        """
        Function to make a .pdb file of the hybrid residue
        """
        res1_check = pdb1.df['ATOM'].isin({'atom_name': hybrid_atoms['atoms']})['atom_name']
        swap_atoms = hybrid_atoms['atoms'].str.swapcase()
        compare_atoms = pdb2.df['ATOM'].atom_name
        res2_check = compare_atoms.isin(swap_atoms)
        atoms_res1 = pdb1.df['ATOM'][res1_check == True]
        atoms_res2 = pdb2.df['ATOM'][res2_check == True]
        atoms_res2['atom_name'] = atoms_res2['atom_name'].str.lower()
        hybrid_pdb = atoms_res1.append(atoms_res2, ignore_index=True)
        hybrid_pdb = hybrid_pdb.assign(atom_number=[x for x in range(1, len(hybrid_pdb)+1)], residue_name='HYB')
        hybrid = PandasPdb()
        hybrid.df['ATOM'] = hybrid_pdb
        hybrid.to_pdb('tmp/hybrid.pdb')

    def make_hybrid_lib(self, lib1, lib2):
        """
        Function to make a .lib file of the hybrid residue
        """
        # {HYB}
        # print('{HYB}')

        # [atoms]
        # for index, row in self.hybrid_atoms.iterrows():
        #     print(row.atoms, row.type, row.charge)

        # [bonds]
        # self.res1.get_bonds(self.res1.lib)
        # self.res2.get_bonds(self.res2.lib)
        # check = ~self.res1.bonds.isin(self.res2.bonds)
        # bonds = self.res1.bonds.append(self.res2.bonds[check.atom2 == True])
        # a = 0
        # for i, row in bonds.iterrows():
        #     if int(i) == a:
        #         print(row.atom1, row.atom2)
        #     else:
        #         if int(i) == (a-1):
        #             print(row.atom1.lower(), row.atom2.lower())
        #         else:
        #             print(row.atom1, row.atom2.lower())
        #             a = int(i)+1
        #     a += 1

        # [impropers]
        # if self.res1.get_impropers(self.res1.prm).empty and self.res2.get_impropers(self.res2.prm).empty:
        #     pass
        # else:
        #     self.res1.get_impropers(self.res1.prm)
        #     if not self.res1.impropers.empty:
        #         for i, row in self.res1.impropers.iterrows():
        #             print(row.type1, row.type2, row.typw3, row.type4, row.kx, row.x0)
        #     self.res2.get_impropers(self.res2.prm)
        #     if not self.res2.impropers.empty:
        #         for i, row in self.res2.impropers.iterrows():
        #             print(row.type1, row.type2, row.type3, row.type4, row.kx, row.x0)

        # [charge_groups]



# test = Hybrid(res1=args.res1, res2=args.res2)
# test.get_hybrid_atoms()
# test.make_hybrid_pdb(test.res1.pdb, test.res2.pdb, test.hybrid_atoms)
# test.make_hybrid_lib(test.res1.atoms, test.res2.atoms)
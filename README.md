# Advanced-Scientific-Python-Programming

The goal of the project is to create a package for the setup of dual topology molecular dynamics (MD) - free energy perturbation (FEP) of amino acid side chains.
The package should contain several modules that will generate the input files necessary for the requested single point mutation of any given amino acid into another.

The first stage of the project is to aim at the mutation of the so-called side chain mimic, i.e. only the molecule that makes up the side chain of the amino acid excluding the backbone atoms (e.g. methane for alanine, methanol for serine, etc).
The module should be able to read in the atoms of the WT side chain / mimic and the atoms of the mutant atoms, assess approriate library and parameter files readable for the force field used in the MD simulations (MD software used will be Q, force field for initial design will be OPLS-AA), create the corresponding files for the dual topology hybrid side chain / molecule, and generate the input FEP and MD files. 
Herein, it is important to store the molecular and atomic properties describing the bonding and non-bonding parameters such as atomic charges, bonding, angle bending, torsional and dihedral constants, and van der Waals and electrostatic interaction constants.

Ultimately, in a later (challenging) stage the modules should also work for the full amino acid mutations, for both in aqueous environment (estimate of the unfolded protein state) as well as a solvated protein environment (the folded state). This will allow for the calculation of the free energy differnce of the single point mutation for either protein stability or substrate binding purposes. 

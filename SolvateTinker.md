# Solvate a protein from tinker xyz format

This assumes that you have a TINKER XYZ file of your complete PROTONATED protein chain (and ligand if applicable). If you are missing residues or protons go back and follow the instructions in [ProteinPrep](./ProteinPrep.md) or go to [PDBtoTinkerXYZ](./PDBtoTinkerXYZ.md) to convert your file type.

This is the alternative method to using the gromacs solvate to create your waterbox and ions around your protein. This method is best for small molecules, perhaps paramterized from Poltype2, or proteins already in tinker xyz format that you do not wish to convert back to pdb

If you have a pdb file of your protein handy, we will still use pacmol to get some general information about our system that will be useful for the solvation protocol. Otherwise you will need to figure out the size of the box needed, and charge of your protein manually. 

### Pacmol - system charge and box size
The charge of the protein is computed by adding the charge of all residues. Only HIS, ARG, LYS, GLU and ASP are charged so the problem is equivalent to counting how many of these residues you have in your protein. 

There is a very simple tool to help you do that, which you can download with the [Packmol software](http://m3g.iqm.unicamp.br/packmol/download.shtml). Once you have downloaded the binaries, look for the executable `solvate.tcl`

Then navigate to the directory where you have your protein chain PDB and execute the script:

```sh
cd $PATH_to_directory
~/Downloads/packmol/solvate.tcl 1kjl_complete.pdb

```
Note that it may be a good idea to move the binaries from the `Downloads` folder. Don't forget to update the path to `solvate.tcl` if that is the case. 

This should return the total structure charge and other information:

```sh
  ###########################################################################
 solvate.tcl script: Solvates a (protein) structure with water and ions,
                     using Packmol. 
 ###########################################################################
 Input pdb file to be solvated: 1kjl_complete.pdb 
 User set structure charge: none  - will compute from structure. 
 User set shell size: none  - using default: 15.0 Angs. 
 User set density: none  - using default: 1.0 g/ml 
 -------------------------------------------------------
  Minimum and maximum coordinates of structure atoms 
  X_min = -12.695     X_max = 26.753 
  Y_min = -22.049     Y_max = 23.335 
  Z_min = -24.602     Z_max = 13.634 
 -------------------------------------------------------
  Box side length in each direction: 
  x: 69.44800000000001 
  y: 75.384 
  z: 68.236 
 -------------------------------------------------------
 -------------------------------------------------------
  HIS = 0 (associated charge = 0) 
  ARG = 9 (associated charge = +9) 
  LYS = 8 (associated charge = +8) 
  GLU = 6 (associated charge = -6) 
  ASP = 7 (associated charge = -7) 
 -------------------------------------------------------
  Total structure charge = 4
 -------------------------------------------------------
 Unsure about mass of element HG3: 1.00800. Is this correct? (y/n)
```


Type 'y' and hit enter for any hydrogens it cannot identify


```sh
 -------------------------------------------------------
  Molar mass of structure: 16772.934399999976
 -------------------------------------------------------
  Number of water molecules to be put:   11016 
  Total volume: 357233.75 A^3
  Volume occupied by water: 329372.54 A^3 
  Number of Sodium ions to be put: 28
  Number of Chloride ions to be put: 32
  Wrote packmol input. 
 -------------------------------------------------------
  Use these lengths for periodic boundary conditions: 
  x: 70.44800000000001
  y: 76.384
  z: 69.236
 -------------------------------------------------------
 -------------------------------------------------------
  Now, run packmol with: packmol < packmol_input.inp       
 -------------------------------------------------------
  The solvated file will be: solvated.pdb 
 -------------------------------------------------------
```

We use packmol for the charge information, number of ions to add, and to get a box size estimate. 

Our system has a positive charge of 4 meaning to neutralize it we will need to add 4 negative ions, probably chloride ions. 

We also get a nice measurement of the protein dimensions.You can manually go into VMD to measure your system size as well. We will need to know how big the protein is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. I will usually round up the largest side that pacmol recommends and use that for all sides of my box. 


For the galectin-3 system I use an 80x80x80 box, which is slightly overestimating the box size but you want to make sure that the protein is never interacting with its mirror images across periodic boundaries. ***This is especially tricky for flexible molecules!*** The downside to a larger volume box be the longer simulation time that is required. 

### 







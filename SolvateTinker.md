# Solvate a protein from tinker xyz format

This assumes that you have a PDB file of your complete PROTONATED protein chain (and ligand if applicable). If you have multiple chains, select one now. If you are missing residues or protons go back and follow the instructions in [ProteinPrep](./ProteinPrep.md)

This is the alternative method to using the gromacs solvate to create your waterbox and ions around your protein. This method is best for small molecules, perhaps paramterized from Poltype2, or proteins already in tinker xyz format that you do not wish to convert back to pdb

If you have a pdb file of your protein handy, we will still use pacmol to get some general information about our system that will be useful for the solvation protocol. Otherwise you will need to figure out the size of the box needed, and charge of your protein manually. 

### pacmol
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
  HIS = 4 (associated charge = 0) 
  ARG = 9 (associated charge = +9) 
  LYS = 8 (associated charge = +8) 
  GLU = 6 (associated charge = -6) 
  ASP = 7 (associated charge = -7) 
 -------------------------------------------------------
  Total structure charge = 4
 -------------------------------------------------------
 Unsure about mass of element HG3: 1.00800. Is this correct? (y/n)

```

Note that you obtain additional information from the routine that you may find useful, such as

```sh
 -------------------------------------------------------
  Molar mass of structure: 17224.708600000013
 -------------------------------------------------------
  Number of water molecules to be put:   10991 
  Total volume: 357233.75 A^3
  Volume occupied by water: 328622.34 A^3 
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

You can manually go into VMD to measure your system size as well. We will need to know how big the protein is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. I will usually round up the largest side that pacmol recommends and use that for all sides of my box. For the galectin-3 system I use an 80x80x80 box, which is slightly overestimating the box size but you want to make sure that the protein is never interacting with its mirror images across periodic boundaries. The downside is that a larger volume is going to require longer simulation times so it is a balance. 

### PDB to Tinker XYZ

If coming from the protein prep steps, the first step in solvation will be converting your pdb to a tinker xyz using the Tinker8 executable ```xyzedit```

First download a set of tinker [linux executables](https://dasher.wustl.edu/tinker/) if you are carrying out the solvation on ARC, or a set of mac executables if carrying these steps out on your local macbook.

For a quick overview of how Tinker works and what Tinker XYZ files are, you can look at these external tutorials: [Tinker Explained](http://chembytes.wikidot.com/tinker-s-wiki) and [Tinker Preparation](https://tinker-hp.org/wp-content/uploads/2022/10/Tinker_preparation_tutorial.pdf)



Adding water and ions to the system was probably your last step so your starting PDB here will be the output of Gromacs solvate (see `SolvateProtein.md`).


###Edit the PDB file
1) Change the name of histidine residues from HIS to HIE so that the protonation state remains when converting to Tinker XYZ. In vim you can do it easily by typing:
`:%s/HIS/HIE/g`

2) Rename the water "residues" as `HOH` rather than the default Gromacs `SOL`. You can also do this in vim by typing `:%s/SOL/HOH/g`

3) If a ligand (or any complexed small molecule) is present, make sure it is at the end of the file. In general, it will be included after the protein coordinates and will need to be cut and paste after the water+ions. Open the PDB in VMD and save it again to update the atom numbering (check the result by making sure the ligand atoms appear at the end with its atom numbers in sequential order). 

###Convert to Tinker XYZ
1) Use the pdbxyz tool by typing `~/Tinker/pdbxyz filename.pdb `. This will create the file `filename.xyz` where all atoms from the pdb were converted to Tinker XYZ and assigned `0` as atom type. CHECK (please) that you indeed have the same number of atoms as in the PDB. 

2) Make a copy of the pdb file `filename.pdb` and change `ATOM` for `HETATM` for the water, ions and ligand atoms. Do this by typing:

```sh
cp filename.pdb filename_hetatm.pdb
vi hetatm.pdb 

```
and look for the line numbers (type `:set number` in vim) corresponding to the water, ion and ligand coordinates. For example if the protein coordinates go from line 2 to line 4500 and the rest of the coordinates go from line 4501 to 65000, type `:4501,65000 s/ATOM  /HETATM/g` (make sure there are two spaces after `ATOM` to conserve the spacing between columns). This will change `ATOM` for `HETATM` between lines 4501 and 65000, therefore leaving the protein coordinates as it were. 

3) Convert `filename_hetatm.pdb` to Tinker XYZ with: `~/Tinker/pdbxyz filename_hetatm.pdb `. This creates a second XYZ file (`filename_hetatm.xyz`) with actual atom types assigned to the protein, water and ions (check). Note however that the ligand coordinates are missing here (check that the number of atoms in this new XYZ file = total number of atoms - number of ligand atoms). 

4) Cut the ligand coordinates from the first XYZ file you converted (`filename.xyz`). Since you put the ligand coordinates at the end, this is trivially done by typing `tail -n 81 filename.xyz >> ligand.xyz` for example. Replace `81` by the actual number of atoms in your ligand. To assign atom types to the ligand, you will need a script that pairs atom names to atom types (make sure all atoms have a unique name for this, for example C1, C2, C3, ... N1, N2, ... H1, H2, ...H30). Such a script (`atomtype.sh` looks like:

```sh
#!/bin/bash

awk '{
split($0, a, FS, seps);
if (a[2]=="N6B") a[6]="253";
if (a[2]=="H63") a[6]="260";
if (a[2]=="H64") a[6]="260";
if (a[2]=="C6B") a[6]="249";
if (a[2]=="N1B") a[6]="257";
if (a[2]=="C2B") a[6]="254";
if (a[2]=="H2B") a[6]="259";
if (a[2]=="N3B") a[6]="256";
...
if (a[2]=="H62") a[6]="260";
for (i=1;i<=NF;i++) printf("%s%s", a[i], seps[i]); print ""}' $1 >> $2

```

and can be used on `ligand.xyz` as follows `./atomtype.sh ligand.xyz ligand_atomtype.xyz`. This will create the file `ligand_atomtype.xyz` that contains the Tinker XYZ formatted ligand coordinates and corresponding atom types. 

5) Add the ligand coordinates at the end of `filename_hetatm.xyz` for a Tinker XYZ file that contains ALL atoms. This can be done by typing `cat filename_hetatm.xyz ligand_atomtype.xyz >> Final_input.xyz` making sure to delete any blank likes (if there are any) between the two files. 

6) Update the total number of atoms on the first line (to include the ligand atoms) and CHECK that your `Final_input.xyz` file has (i) an atom type assigned to all atoms (ii) the correct number of atoms. 







Using the path to your tinker executables, call xyzedit program with the PDB file as the target. 

```

```














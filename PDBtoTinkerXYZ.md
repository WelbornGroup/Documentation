

# Convert PDB to Tinker XYZ

The [Tinker software](https://dasher.wustl.edu/tinker/) reads input files with extension `.xyz` that are formatted differently from conventional XYZ files. 

For a quick overview of how Tinker works and what are Tinker XYZ files, please refer to [Tinker Explained](http://chembytes.wikidot.com/tinker-s-wiki) and [Workshop Exercise](https://sites.google.com/site/amoebaworkshop/exercise-1#TOC-Restarting-a-simulation).

The below assumes that you have Tinker downloaded and installed on Cascades (binaries available) as well as a complete PDB file of your single protein chain in water (no residue is missing, there are protons and water + ions surround the protein and complexed ligand, filling a box of appropriate dimension). 

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
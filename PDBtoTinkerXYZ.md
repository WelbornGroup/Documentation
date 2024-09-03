

### PDB to Tinker XYZ

If solvating in Tinker, the first step will be converting your PDB to a tinker xyz file

If you have already solvated your system using gromacs (protein + water + ions), or have a ligand present, there will be a few extra steps to consider 

First download a set of tinker [executables](https://dasher.wustl.edu/tinker/) - linux version if you are carrying out the solvation on ARC, or a set of mac executables if carrying these steps out on your local macbook. Also download the amoebabio18.prm file, which contains parameters to simulate canonical proteins in AMOEBA. 

For a quick overview of how Tinker works and what Tinker XYZ files are, you can look at these external tutorials:[Tinker Preparation](https://tinker-hp.org/wp-content/uploads/2022/10/Tinker_preparation_tutorial.pdf) and [Tinker Explained](http://chembytes.wikidot.com/tinker-s-wiki)


### Prepare the PDB file
1) Change the name of histidine residues from HIS to HIE *or* HID if necessary, so that the protonation state remains when converting to Tinker XYZ. These labels state whether the histidine is protonated in the epsilon or delta position, and the HIS stands for the doubly protonated state where the histidine is positively charged.

From our 1kjl galectin-3 example we can look at the pacmol output we see that we have 4 histidines that are contributing no charge, therefore they are all in HIE or HID states. You can look back at the output PDB from the reduce program to see which of the hydrogens was omitted to determine which state that HIS should be changed to. It may be easier to open the 1kjl_complete.pdb in VMD and visually assess the histidines. The selection `all resname HIS` should select the residues for you. You will see that the four histidines in our system, from lower to higher atom index, are in states HID HIE HIE HID.

Alternatively, you can make it simple and make all of the histidines a single protonation state by changing all HIS to HIE.

You can accomplish this with find and replace in a text editor, or an example with terminal is:
```sh
sed 's/HIS/HIE/g' 1kjl_complete.pdb > 1kjl_prep1.pdb
```

If there is a histidine in an active site or one that is positively charged you may not want to take the all HIE approach.


2) If coming from gromacs solvate, rename the water "residues" as `HOH` rather than the default Gromacs `SOL`. You can do this in terminal with `sed 's/SOL/HOH/g' 1kjl_prep1.pdb > 1kjl_prep2.pdb` 

3) If a ligand (or any complexed small molecule) is present, make sure it is at the end of the file. In general, it will be included after the protein coordinates and will need to be cut and paste after any water+ions that are present. Open the PDB in VMD and save it again to update the atom numbering (check the result by making sure the ligand atoms appear at the end with its atom numbers in sequential order). 



### Convert to Tinker XYZ
1) Use the pdbxyz tool by typing the following command in terminal: `~/Tinker/pdbxyz filename.pdb`

  It will prompt you for a parameter file, type `amoebabio18.prm` (you should have it in your working directory alongside the PDB) and hit enter

This will create the file `filename.xyz` where all atoms from the pdb were converted to Tinker XYZ. CHECK (please) that you indeed have the same number of atoms as in the PDB. For our case we can see that we actually gained 3 atoms in this process. This is not a problem and it usually happens at the terminal amino acids when converting to a tinker xyz. For this specific example our N-terminus proline gained 2 hydrogens giving it a positive charge, and our C-terminus gained an oxygen to complete a carboxylic acid group on the isoleucine. Both of these changes are fine since we normally want charged terminal ends for our simulations, so our new atom total is 2230. (If we had instead left our histidine groups named HIS we would see our atom number at 2234.)

If you had only a protein present you are done with the PDB to Tinker xyz conversion! You can move on to the minimization and dynamics simulation in [RunningTinkerBasics](./RunningTinkerBasics.md)

2) If you have water, ions, or ligands present we need to continue:

Make a copy of the pdb file `filename_hetatm.pdb` and change `ATOM` for `HETATM` for the water, ions and ligand atoms using a command like:

```sh
sed '2232,51076 s/ATOM  /HETATM/g' filename.pdb > filename_hetatm.pdb
```

This will change `ATOM` for `HETATM` for lines 2232 to 51076, range inclusive (make sure there are two spaces after `ATOM` to conserve the spacing between columns). These lines should be corresponding to the water, ion and ligand coordinates in our file and not the protein lines.

3) Convert `filename_hetatm.pdb` to Tinker XYZ with: `~/Tinker/pdbxyz filename_hetatm.pdb `. This creates a second XYZ file (`filename_hetatm.xyz`) with actual atom types assigned to the protein, water and ions (check). Note however that the ligand coordinates are missing here since they were not present in the parameter file. Check that the number of atoms in this new XYZ file = total number of atoms - number of ligand atoms + any terminal residue additions. 

4) Cut the ligand coordinates from the first XYZ file you converted (`filename.xyz`). Since you put the ligand coordinates at the end, this can be done in terminal with `tail -n 81 filename.xyz >> ligand.xyz` for example. Replace `81` by the actual number of atoms in your ligand. To assign atom types to the ligand, you will need a script that pairs atom names to atom types (make sure all atoms have a unique name for this, for example C1, C2, C3, ... N1, N2, ... H1, H2, ...H30). A bash script `atomtype.sh` that accomplishes this looks like:

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

Alternatively a python script can be used to similar effect. ChatGPT is a useful tool for small scripts like this but keep backup files while testing scripts so you do not lose your progress!

5) Add the ligand coordinates at the end of `filename_hetatm.xyz` for a Tinker XYZ file that contains ALL atoms. This can be done by typing `cat filename_hetatm.xyz ligand_atomtype.xyz >> final_input.xyz` making sure to delete any blank likes (if there are any) between the two files. 

6) Update the total number of atoms on the first line (to include the ligand atoms) and CHECK that your `final_input.xyz` file has (i) an atom type assigned to all atoms (ii) the correct number of atoms. 








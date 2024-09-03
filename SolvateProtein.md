# Solvate a protein from PDB

This assumes that you have a PDB file of your complete PROTONATED protein chain (and ligand if applicable). If you have multiple chains, select one now. If you are missing residues or protons go back and follow the instructions in [ProteinPrep](./ProteinPrep.md)


### What is the charge of the system ?
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


Type 'y' and hit enter for any hydrogens it cannot identify


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

You can manually go into VMD to measure your system size as well. We will need to know how big the protein is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. I will usually round up the largest side that pacmol recommends and use that for all sides of my box. For the galectin-3 system I use an 80x80x80 box, which is slightly overestimating the box size but you want to make sure that the protein is never interacting with its mirror images across periodic boundaries. ***This is especially tricky for flexible molecules!*** The downside to a larger volume box be the longer simulation time that is required. 


### Solvate the protein on ARC
There are a number of routines available to solvate a protein. Instad of using pacmol, we will use the Gromacs `solvate` routine as it allows water to be placed within as well as outside the protein. Gromacs is already installed on the tinkercliffs CPU on ARC so we continue there. 

Copy the pdb file over to ARC and log in:

```sh
scp Name_of_file.pdb username@tinkercliffs.arc.vt.edu:/path_to_working_directory/
ssh username@tinkercliffs.arc.vt.edu
cd path_to_working_directory

```

***OR*** use OpenOnDemand and open a tinkercliffs cpu terminal. I prefer this and would suggest it to beginners over using purely a terminal to access ARC.

Write the script `solvate.sh` that contains:

```sh
#!/bin/bash -l
#SBATCH -p dev_q
#SBATCH -J Solvate
#SBATCH -N 1
#SBATCH --ntasks-per-node 32
#SBATCH -t 00:30:00 
#SBATCH -A welbornlab
 
module load GROMACS

# Center the protein in a box. Change the name of the input file (protein_input.pdb and the size of the box. Watch out the units are in nm here!
gmx editconf -f $1 -o box.pdb -box 8.000 8.000 8.000 -c 

# Fill the box with water. You can change the name of the output file (Box_Water.pdb)
gmx solvate -cp box.pdb -o box_water.pdb 

# Insert ions, making sure the system is overall neutral, here we will add both ions to create 0.1M conditions. You will need to change how many ions you want and the name of the output file. chloride.pdb and sodium.pdb are input files you will need to copy from below.

gmx insert-molecules -f box_water.pdb -nmol 29 -ci sodium.pdb -o box_water_Na.pdb 

gmx insert-molecules -f box_water_Na.pdb -nmol 33 -ci chloride.pdb -o box_water_NaCl.pdb

```

`chloride.pdb`

```sh
COMPND    Chloride
AUTHOR    GENERATED BY OPEN BABEL 2.3.90
HETATM    1 CL    CL     1      -1.694  -0.065   0.000  1.00  0.00          Cl  
MASTER        0    0    0    0    0    0    0    0    1    0    1    0
END
```

`sodium.pdb`

```sh
COMPND    Sodium
AUTHOR    GENERATED BY OPEN BABEL 2.3.90
HETATM    1 NA    NA     1      -1.694  -0.065   0.000  1.00  0.00          Na  
MASTER        0    0    0    0    0    0    0    0    1    0    1    0
END
```

Note that since the solvation process is very short, you have two options to execute the script.

Option 1 (preferred): run it interactively on the login node by typing `./solvate.sh`.

Option 2: use the `dev_q` to submit the job on the queue by typing `sbatch solvate.sh`

The development queues on ARC have a smaller limit on the amount of jobs you can submit and also a max time limit of 4 hours per job, they can be useful for short and intensive jobs that may be too much to run on the login node.



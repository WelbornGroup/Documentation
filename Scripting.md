Start with the tutorials described in `READMEFIRST_workflow.md `. You can only script efficiently once you know well the sequence of tasks to perform. 

###Solvation and PDB to XYZ with a script

The solvation script described in `SolvateProtein.md` already does all the steps automatically, provided you edit the name of the file prior to running it. To make it more efficient to run multiple times, it would be better if we didn't have to update the name each time. To do that, we write the script such that it takes the name of the input file from an argument in the command line as opposed to hardcopy it within. The new script, `solvate_script.sh` is then

```sh
#!/bin/bash -l
##SBATCH -p normal_q
#SBATCH -p dev_q
#SBATCH -J Solvate
#SBATCH -N 1
#SBATCH --ntasks-per-node 32
#SBATCH -t 00:01:00 
#SBATCH -A welbornlab
 
module load gcc/5.2.0  mvapich2/2.2
module load gcc/5.2.0  openmpi/3.0.0
module load gcc/5.2.0  openmpi/3.1.2
module load gcc/6.1.0  openmpi/3.0.0
module load gcc/6.1.0  openmpi/3.1.2
module load intel/15.3  mvapich2/2.2
module load intel/15.3  openmpi/3.0.0
module load intel/15.3  openmpi/3.1.2
module load intel/18.2  openmpi/3.0.0
module load fftw/3.3.8
module load gromacs/5.1.2 

gmx editconf -f $1 -o Box.pdb -box 8.600 9.700 7.900 -c 
gmx solvate -cp Box.pdb -o Box_Water.pdb 
gmx insert-molecules -f Box_Water.pdb -nmol 27 -ci Sodium.pdb -o Box_Water_Na.pdb 
#gmx insert-molecules -f Box_Water_Na.pdb -nmol 49 -ci Chloride.pdb -o Box_Water_Na_Cl.pdb
```

that is now executed as `./solvate_script.sh inputfile.pdb`. Note that here, you still need to edit the box dimensions and the number of sodium and chloride ions you need and it will write the output as `Box_Water_Na_Cl.pdb` or `Box_Water_Na.pdb`. The whole script could be changed such that all of these options are read from the command line but here, given how we are planning to use it, it is enough to only edit the input file command.   


Now, assuming the output of `solvate_script.sh` was `Box_Water_Na_Cl.pdb` you can script all the steps to convert the PDB into a Tinker XYZ file, including centering in the box. The script, `PDBtoXYZ.sh` is given below and ASSUMES that you have a `tinker.key` file (do not change the name) in your directory. Without this `tinker.key` file, a few steps in this script will not work. Note also that this is an example script for a protein that has a ligand named `X77` containg 67 atoms - you will need to edit the corresponding lines for your case as well as the paths to the executables prior to execution.

```sh
sed -i 's/SOL/HOH/g' Box_Water_Na_Cl.pdb 
sed -i 's/HIS/HIE/g' 
grep "X77" Box_Water_Na_Cl.pdb >> ligand.pdb
grep -v "X77" Box_Water_Na_Cl.pdb >> noligand.pdb
cat noligand.pdb ligand.pdb >> Reorder.pdb
egrep -v "TER|ENDMDL" noligand.pdb
~/Tinker/pdbxyz Reorder.pdb
cp Reorder.pdb Reorder_copy.pdb
sed -i '/X77/ s/ATOM  /HETATM/' Reorder_copy.pdb
sed -i '/CL/ s/ATOM  /HETATM/' Reorder_copy.pdb
sed -i '/NA/ s/ATOM  /HETATM/' Reorder_copy.pdb
sed -i '/HOH/ s/ATOM  /HETATM/' Reorder_copy.pdb
~/Tinker/pdbxyz Reorder_copy.pdb 
tail -n 67 Reorder.xyz >> ligand.xyz
~/X77_atomtypes.sh ligand.xyz ligand_at.xyz
cat Reorder_copy.xyz ligand_at.xyz >> Input.xyz
```

where `X77_atomtypes.sh` is the script to assign ligand atom types:

```sh
#!/bin/bash

awk '{
split($0, a, FS, seps);
if (a[2]=="C1") a[6]="901";
if (a[2]=="C2") a[6]="902"; 
if (a[2]=="C3") a[6]="903";
if (a[2]=="C4") a[6]="904";
if (a[2]=="C5") a[6]="905";
if (a[2]=="C6") a[6]="906";
if (a[2]=="C7") a[6]="907";
if (a[2]=="C8") a[6]="908";
if (a[2]=="C9") a[6]="909";
if (a[2]=="C10") a[6]="910";
if (a[2]=="C11") a[6]="911";
if (a[2]=="C12") a[6]="912";
if (a[2]=="C13") a[6]="913";
if (a[2]=="C14") a[6]="914";
if (a[2]=="C15") a[6]="915";
if (a[2]=="C16") a[6]="916";
if (a[2]=="C17") a[6]="917";
if (a[2]=="C18") a[6]="918";
if (a[2]=="C19") a[6]="919";
if (a[2]=="C20") a[6]="920";
if (a[2]=="C21") a[6]="920";
if (a[2]=="C22") a[6]="920";
if (a[2]=="C23") a[6]="921";
if (a[2]=="C24") a[6]="922";
if (a[2]=="C25") a[6]="923";
if (a[2]=="C26") a[6]="924";
if (a[2]=="C27") a[6]="925";
if (a[2]=="N1") a[6]="926";
if (a[2]=="N2") a[6]="927";
if (a[2]=="N3") a[6]="928";
if (a[2]=="N4") a[6]="929";
if (a[2]=="N5") a[6]="930";
if (a[2]=="O1") a[6]="931";
if (a[2]=="O2") a[6]="932";
if (a[2]=="H1") a[6]="933";
if (a[2]=="H2") a[6]="934";
if (a[2]=="H3") a[6]="935";
if (a[2]=="H4") a[6]="935";
if (a[2]=="H5") a[6]="936";
if (a[2]=="H6") a[6]="936";
if (a[2]=="H7") a[6]="937";
if (a[2]=="H8") a[6]="937";
if (a[2]=="H9") a[6]="938";
if (a[2]=="H10") a[6]="938";
if (a[2]=="H11") a[6]="939";
if (a[2]=="H12") a[6]="939";
if (a[2]=="H13") a[6]="940";
if (a[2]=="H14") a[6]="941";
if (a[2]=="H15") a[6]="942";
if (a[2]=="H16") a[6]="943";
if (a[2]=="H17") a[6]="944";
if (a[2]=="H18") a[6]="945";
if (a[2]=="H19") a[6]="946";
if (a[2]=="H20") a[6]="946";
if (a[2]=="H21") a[6]="946";
if (a[2]=="H22") a[6]="946";
if (a[2]=="H23") a[6]="946";
if (a[2]=="H24") a[6]="946";
if (a[2]=="H25") a[6]="946";
if (a[2]=="H26") a[6]="946";
if (a[2]=="H27") a[6]="946";
if (a[2]=="H28") a[6]="947";
if (a[2]=="H29") a[6]="948";
if (a[2]=="H30") a[6]="949";
if (a[2]=="H31") a[6]="950";
if (a[2]=="H32") a[6]="951";
if (a[2]=="H33") a[6]="952";
for (i=1;i<=NF;i++) printf("%s%s", a[i], seps[i]); print ""}' $1 >> $2
```
This script will also be specific to your ligand and will need to be activated the first time with `chmod +x X77_atomtypes.sh`. 


###Looping over a number of structures
In the group, we often have to run the same simulations for a number of starting protein conformations, which improves our sampling. This doesn't represent much more work if it is set up properly. In the section above, we saw how to write a script to do all the file preparation work, from adding solvent to convert the PDB into an XYZ file. Now, we can run these as many times as we want, in a single command line. 

Go in the directory where you have, say 25 protein conformation PDB files, named `myprotein_1.pdb` through `myprotein_25.pdb`. We will create one directory for each simulation, numbering them 1 through 25. You will need a number of files that are common to all 25 cases, such as tinker.key, solvate.sh ... Start by creating a directory where you will keep all these common files and edit them for your specific project before copying them over to where you will run the simulations. 

```sh
cd MyDirectory		#where you have your 25 starting PDB files
mkdir CommonFiles
cd CommonFiles
cp or scp (depending on where you have them) all necessary files.
cd ../
```


The necessary files are:

- `solvate_script.sh` where you will need to edit the box dimension, number of Na+ or Cl- ions needed,
- `sodium.pdb` and `chloride.pdb` to be able to add the ions in the box, 
- `tinker.key` where you will need to edit the box dimension and the path to the parameter file
- `PDBtoXYZ.sh` to convert your PDB into a Tinker XYZ file
- `Name_atomtypes.sh` script to assign atom type to ligand or small molecules. This script will need to be fully updated to your particular ligand. Remember to `chmod +x Name_atomtypes.sh` the first time you create it so that it is recognized as an executable.


It is easy to make a mistake while editing these files (you can forget a file, forget to change a path, ...) so the first thing is to test it on one structure (assuming you are now in `MyDirectory`):

```sh
mkdir Test
cp myprotein_1.pdb Test/
cp CommonFiles/* Test/
cd Test
./solvate_script.sh myprotein_1.pdb
./PDBtoXYZ.sh
```

This should run without error and create, among other files, the Tinker input file `Input_final.xyz_2`. Fix any errors that arise and once it works go back to `MyDirectory`. IF you had to edit the scripts to correct a mistake, DO NOT FORGET to copy them back to the `CommonFiles` directory. Say you had a typo in `Name_atomtypes.sh` do `cp Test/Name_atomtypes.sh CommonFiles`. 

Now, you are ready to loop over all the conformations (assuming you are now in `MyDirectory`). You can type:



```sh
for i in `seq 1 25`; do mkdir Simulation_${i}; done
for i in `seq 1 25`; do cp CommonFiles/* Simulation_${i}; done
for i in `seq 1 25`; do cd Simulation_${i}; ./solvate_script.sh myprotein_${i}.pdb; cd ../; done
for i in `seq 1 25`; do cd Simulation_${i};./PDBtoXYZ.sh; cd ../ ; done
```
or, combining them into a single line (note the `;` between each steps):

```sh
for i in `seq 1 25`; do mkdir Simulation_${i}; cp CommonFiles/* Simulation_${i}; cd Simulation_${i}; ./solvate_script.sh myprotein_${i}.pdb; ./PDBtoXYZ.sh; cd ../ ; done
```

This should get you ready to start the energy minimization with a `Input_final.xyz_2` in each of the 25 directory. 

Do check on a couple of directories (take 2 at random from the 25), that the input files are there (not a empty file) and different from one another. 

Say now that you have the script to launch the energy minimization (see `RunningTinkerBasics.sh`) and that you edited the name of the input file in that script for `Input_final.xyz_2`, you can start all energy minimization at once (assuming you are now in `MyDirectory`):

```sh
for i in `seq 1 25`; do cp launch_minimization.sh Simulation_${i}; done
for i in `seq 1 25`; do cd Simulation_${i}; sbatch launch_minimization.sh; cd ../; done
```

And similarly for the dynamics once the energy minimization is done. 

This file is detailed







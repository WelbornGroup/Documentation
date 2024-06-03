# Workflow for REU students - Summer 2024

## Create Carbohydrate Structures

Divide up the monosaccharides between each other, just make sure you both are not doing the same structures twice!

Use the glycam carbohydrate builder to generate the initial carbohydrate structure. Make sure that the stereochemistry and substituent groups are correct for each structure. (double check please) Then save the minimized PDB file. 

Open the `.pdb` file in Avogadro, and save the molecule again as a `.mol` file

## Poltype

Poltype2 is the program we will use to generate the parameters for these carbohydrate molecules

This is the first major step in your project this summer, we will work on getting most of the carbohydrates parameterized first before working on the following steps.

First you will need to install Poltype2 on ARC, see [Poltype_Install.md](https://github.com/WelbornGroup/Documentation/blob/REU_update/Poltype_Install.md) for directions.

Once installed, see the [Poltype_Usage](https://github.com/WelbornGroup/Documentation/blob/REU_update/Poltype_Usage.md) for instructions on running Poltype2 for a molecule on ARC 

## Solvation 

Now we will use Tinker to solvate the carbohydrate molecules 

Copy the `final.xyz` from the poltype result into a new directory and rename the file (i.e. galactose.xyz)

Combine the resulting parameters in `final.key` with `amoebabio18.prm` naming it the same as the structure (i.e. galactose.prm) and place this new paramter file in the same directory

Copy the pre-made water box located in `projects/welborn/Newman/Path` to the directory 

This example will use the Tinker8 executables located `/projects/welbornlab/Poltype2/TinkerEx`
(You can copy the TinkerEx directory to your own space and use them from there if you would like)

Now in your directory, open a terminal and run `/projects/welbornlab/Poltype2/TinkerEx/xyzedit galactose.xyz`

You should see the pop-up menu for xyzedit, read through the various ways this code can be used to manipulate your structure file.

Select `(13) Translate Center of Mass to the Origin` by typing 13 and hitting enter
This translates our carbohydrate to the origin, (centering it on x=0,y=0,z=0) and it is good practice to center your structures before carrying out other manipulations

Next select `(24) Soak Current Molecule in Box of Solvent`

This will prompt you for your solvent box name `waterbox.xyz`

There will be a new structure file generated (`galactose.xyz_2`) which will be your centered carbohydrate in the minimized waterbox

You should download and open this file in VMD to make sure it looks okay with no obvious errors

## Minimize
First create a key file named the same as your carbohydrate and parameter files (i.e. `galactose.key`)

```
integrator nose-hoover

a-axis 31.00 
b-axis 31.00
c-axis 31.00

neighbor-list
polar-eps 0.00001
vdw-cutoff 12.0
vdw-correction

ewald
ewald-cutoff 7.0

polar-predict
polarization mutual

verbose
```
This key file houses the settings for our dynamics simulation such as the box dimensions, integrator, barostat, thermostat, etc.
This keyfile will be universal for all the carbohydrates we will do, but you can look up some of the settings and their meanings in the tinker handbook


Next we will create a script to minimize our solvated carbohydrates

In the same directory before, create a file `minimization.sh` which will have the following contents:
```
#!/bin/bash
#SBATCH --account=welbornlab
#SBATCH --partition=v100_normal_q
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=0-00:05:00

# Each node type has different modules avilable. Resetting makes the appropriate stack available
module reset
module load infer-skylake_v100/tinker9/1.4.0-nvhpc-21.11

# Run the example
echo "-------- Starting tinker9 minimization: `date` -------"

tinker9 minimize galactose.xyz_2 0.1 > min.log

echo "------- tinker9 minimization has exited: `date` --------"
```

Make sure to change the structure file name between molecules, and that you are referencing the solvated molecule (*.xyz_2*) and not the original (*.xyz*)

`sbatch minimization.sh` to submit the job to the queue 

Once finished check the min.log file to see if the minimization was completed successfully, and look for the resulting minimized structure file (`galactose.xyz_3`) 

## Molecular Dynamics
Finally we will run a molecular dynamics simulation with Tinker9 on our system
First we will center our molecule again: `/projects/welbornlab/Poltype2/TinkerEx/xyzedit galactose.xyz_3` 
Select `(13) Translate Center of Mass to the Origin`
The result will be `galactose.xyz_4`


Create the dynamics submission script `dynamic.sh`: 
```
#!/bin/bash
#SBATCH --account=welbornlab
#SBATCH --partition=v100_normal_q
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00

# Each node type has different modules avilable. Resetting makes the appropriate stack available
module reset
module load infer-skylake_v100/tinker9/1.4.0-nvhpc-21.11

# Run the example
echo "-------- Starting tinker9: `date` -------"

tinker9 dynamic galactose.xyz_4 30000000 1 10 4 300.00 1.0 > dynamics.log

echo "------- tinker9 has exited: `date` --------"

```

This is running a dynamics simulation using Tinker 9 the settings are the following: 
- 30000000 is the number of steps in femtoseconds which converts to 30 nanoseconds of simulation time
- 1 is the time step in femtoseconds
- 10 is the number of picoseconds between frames that will get printed out to our trajectory
- 4 determines that this is running in the NPT ensemble keeping pressure and temperature to target values
- 300.0 is the degrees in Kelvin
- 1.0 is the pressure in atm

We will keep all of these setting consistent for the carbohydrates, however make sure you are referencing your carbohydrate structure file correctly (*i.e. change galactose.xyz_4 >> example-carbohydrate.xyz_4*)

Make sure your structure file, key file, and parameter files all have the same names (galactose.xyz_4, galactose.key, galactose.prm)

`sbatch dynamic.sh` to submit the dynamics job to the queue 

Check the dynamics.log file to see the progress of your simulation and when it is completed the `.arc` file will contain the trajectory of the system

Next we will load this trajectory file into VMD and analyze our simulation!

## Analysis with VMD


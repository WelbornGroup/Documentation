# Tinker simulations

Log onto ARC and create/go to the directory where you have your input Tinker XYZ file (see `PDBtoTinkerXYZ.md` if you don't have such input). 


To run a tinker molecular dynamics simulation you will always need the following file types: 

.xyz - this is your structure file that we have prepared and solvated

.key - this is the tinker key file that will contain settings for the simulation

.prm - this is the parameter file, we have been using amoebabio18.prm for general proteins

.sh - these are bash scripts and we will use these files to submit minimization and dynamic jobs to the supercomputer queue



First lets create a Tinker key file, you can name this `tinker.key` or my preference is to name it `filename.key` where filename matches my .xyz and .prm file names 

The following is a basic key file with our normal settings for a protein system:

```sh
parameters amoebabio18.prm 

integrator nose-hoover

a-axis 80.00 
b-axis 80.00
c-axis 80.00

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

You will need to edit (i) the path the the parameter file if it is not in the directory with the key file, you can get rid of this parameter line in your keyfile if you name your .prm file the same as your .xyz and .key (ii) the size of the box since it will change per system and (iii) the integrator line if you are not performing a simulation in the NPT ensemble. In this case the integrator is velocity verlet, the thermostat and barostat are coupled nose-hoover. The rest of the keywords are transferable to most simulations we will run. Look into them to discover what they mean.

We should have all of our files in a directory now, except for the bash submission scripts. 

#### Center the system

We want to make sure our system is centered

`/Tinker/xyzedit Final_input.xyz`


Choose option #12 `(12) Translate Center of Mass to the Origin` and type enter to exit. You will find a new file `filename.xyz_2` that contains the centered system. 

If you saved the `.key` file under a different name, you will first be asked to provide the path to the parameter file. 


#### Perform an energy minimization

To launch an energy minimization, write a submission script `launch_minimize.sh` that contains:

```sh
#!/bin/bash
#SBATCH --account=welbornlab
#SBATCH --partition=v100_dev_q
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=0-01:00:00

# Each node type has different modules avilable. Resetting makes the appropriate stack available
module reset
# This loads tinker9 on the gpu
module load infer-skylake_v100/tinker9/1.4.0-nvhpc-21.11

# This lets us know when the job starts and stops 
echo "-------- Starting tinker9 minimization: `date` -------"

# This line runs the minimization using tinker9
tinker9 minimize filename.xyz_2 0.1 > min.log

echo "------- tinker9 minimization has exited: `date` --------"
```

Run this script by opening a terminal on Infer, the GPU, and typing `sbatch launch_minimize.sh`

In addition to the output file `min.log`, output system coordinates will be saved as `Final_input.xyz_3`. It is good practice to center the system before any simulation, so do it again before running the molecular dynamics:
`~/Tinker/xyzedit filename.xyz_3`, which will create a new file `filename.xyz_4`.

####  Start the molecular dynamics simulation
Using the latest system coordinate file, write the submission script `launch_dynamic.sh` containing:

```sh
#!/bin/bash
#SBATCH --account=welbornlab
#SBATCH --partition=v100_normal_q
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=3-12:00:00

# Each node type has different modules avilable. Resetting makes the appropriate stack available
module reset
module load infer-skylake_v100/tinker9/1.4.0-nvhpc-21.11

# Run the example
echo "-------- Starting tinker9: `date` -------"

tinker9 dynamic filename.xyz_4 30000000 1 10 4 300.00 1.0 > dynamics.log

echo "------- tinker9 has exited: `date` --------"
```

This will launch a 30,000,000 step simulation of 1 fs timestep, saving frames every 10 ps in the ensemble 4 (NPT) with temperature 300 K and pressure 1 atm. 

Run this script by opening a terminal on Infer, the GPU, and typing `sbatch launch_dynamic.sh`

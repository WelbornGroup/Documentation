# Tinker MD simulations

Log onto ARC and create/go to the directory where you have your input TinkerXYZ file (the one we made after parameterization). 

To run a Tinker MD simulation you will always need the following file types: 

`polymer_filename.xyz` - this is your structure file that we have prepared and solvated

`tinker.key` - this is the tinker key file that will contain settings for the simulation

`beta_glucose_params.prm` - this is the parameter file, we have been using beta_glucose_params.prm for parameterizing the beta-glucose

`minimize.sh` and `dynamics.sh` - these are bash scripts and we will use these files to submit minimization and dynamic jobs respectively to ARC queue


First lets create the `tinker.key` file - you can name this `tinker.key` (my preference) or `filename.key` where filename matches the `.xyz` and `.prm` file names.

The following is a basic key file with our normal settings for a polymer system (example parameters from the `beta_glucose_params.prm`, change filename as required):

````sh
parameters beta_glucose_params.prm

a-axis 50.00 
b-axis 50.00
c-axis 50.00

integrator verlet
thermostat bussi
barostat montecarlo
polar-eps 0.000010
vdw-cutoff 10.0
ewald
neighbor-list
polar-predict
polarization mutual
````

You will need to edit (i) the path the the parameter file if it is not in the directory with the key file, you can get rid of this parameter line in your keyfile if you name your .prm file the same as your .xyz and .key (ii) the size of the box since it will change per system and (iii) the integrator line if you are not performing a simulation in the NPT ensemble. In this case the integrator is velocity verlet, the thermostat and barostat are bussi and montecarlo respectively (you can also use coupled nose-hoover, in which case the three lines: `integrator verlet thermostat bussi barostat montecarlo` will be replaced by a single line `integrator nose-hoover`). The rest of the keywords are transferable to most simulations we will run. 

We should have all of our files in a directory now, except for the bash submission scripts. 

#### Center the system

We want to make sure our system is centered

`/Tinker/xyzedit polymer_filename.xyz`

Choose option #12 `(12) Translate Center of Mass to the Origin` and type enter to exit. You will find a new file `polymer_filename.xyz_2` that contains the centered system. 

If you saved the `tinker.key` file under a different name, you will first be asked to provide the path to the parameter file. 


#### Perform an energy minimization

To launch an energy minimization, write a submission script `minimize.sh` that contains:

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
tinker9 minimize polymer_filename.xyz_2 0.1 > min.log

echo "------- tinker9 minimization has exited: `date` --------"
```
NOTE - If you're using `#SBATCH --partition=v100_dev_q`, then the module to be loaded will be `module load infer-skylake/tinker9/1.4.0-nvhpc-21.11`.

Run this script by opening a terminal on Infer, the GPU, and typing `sbatch minimize.sh`.

In addition to the output file `min.log`, output system coordinates will be saved as `polymer_filename.xyz_3`. 

NOTE - you should look for `Normal Termination` at the end of the `min.log` file. If that is not present, you need to repeat some or all of the previous steps till you get it.

It is good practice to center the system before any simulation, so do it again before running the molecular dynamics:
`~/Tinker/xyzedit polymer_filename.xyz_3`, which will create a new file `polymer_filename.xyz_4`.

####  Start the molecular dynamics simulation
Using the latest system coordinate file, write the submission script `dynamics.sh` containing:

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

tinker9 dynamic polymer_filename.xyz_4 10000000 1 10 4 300.00 1.0 > dyn.log

echo "------- tinker9 has exited: `date` --------"
```

This will launch a 10,000,000 step simulation of 1 fs timestep (total 10ns) , saving frames every 10 ps in the ensemble 4 (NPT) with temperature 300 K and pressure 1 atm. 

Run this script by opening a terminal on Infer, the GPU, and typing `sbatch dynamics.sh`.

When the dynamics has finished (or the time is up), you should see a `polymer_filename.arc` (the trajectory file which you will load onto VMD), `polymer_filename.dyn` (the velocities/forces update at the last frame), `slurm-xxx.out` where xxx is the JOBID on ARC, and a `dyn.log` file (which will have some basic system information over the trajectory, and errors, if any).

If the MD didn't start or stopped with a error, or any abnormal ending - look for the error in the `dyn.log` file (if it is there) and the `slurm-xxx.out` file.


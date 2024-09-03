# Tinker simulations

Log onto ARC and create/go to the directory where you have your input Tinker XYZ file (see `PDBtoTinkerXYZ.md` if you don't have such input). 

Create a Tinker key file `tinker.key` that contains the following:


```sh
parameters amoebabio18.prm 

integrator nose-hoover

a-axis 80.00 
b-axis 80.00
c-axis 80.00

polar-eps 0.000010
cutoff 10.0
ewald

neighbor-list
polar-predict
polarization mutual

```

Edit (i) the path the the parameter file, (ii) the size of the box and (iii) remove the barostat line if NOT performing an NPT simulation. The rest of the keywords are transferable to most simulations we will run. Make sure you know what they mean. 

###Center the system
Each software has a different convention to center the system, it is important to be in line with Tinker's convention. Use the Tinker tools by typing:

`~/Tinker/xyzedit Final_input.xyz`

Each Tinker routine reads files named `tinker.key` by default. If you have saved the file above under that name then the command will prompt you with:

```sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                          Version 8.7  June 2019                          ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2019               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 The Tinker XYZ File Editing Utility Can :

    (1) Offset the Numbers of the Current Atoms
    (2) Deletion of Individual Specified Atoms
    (3) Deletion of Specified Types of Atoms
    (4) Deletion of Atoms Outside Cutoff Range
    (5) Insertion of Individual Specified Atoms
    (6) Replace Old Atom Type with a New Type
    (7) Assign Connectivities for Linear Chain
    (8) Assign Connectivities Based on Distance
    (9) Convert Units from Bohrs to Angstroms
   (10) Invert thru Origin to Give Mirror Image
   (11) Translate All Atoms by an X,Y,Z-Vector
   (12) Translate Center of Mass to the Origin
   (13) Translate a Specified Atom to the Origin
   (14) Translate and Rotate to Inertial Frame
   (15) Move to Specified Rigid Body Coordinates
   (16) Move Stray Molecules into Periodic Box
   (17) Delete Molecules Outside of Periodic Box
   (18) Append a Second XYZ File to Current One
   (19) Create and Fill a Periodic Boundary Box
   (20) Soak Current Molecule in Box of Solvent
   (21) Place Monoatomic Ions around a Solute

 Number of the Desired Choice [<Enter>=Exit] :  

```
Choose option #12 and type enter to exit. You will find a new file `Final_input.xyz_2` that contains the centered system. 
If you saved the `.key` file under a different name, you will first be asked to provide the path to the parameter file. 


###Perform an energy minimization

To launch an energy minimization, write a submission script `launch_minimization.sh` that contains:

```sh
#!/bin/bash -l
#SBATCH -p normal_q
#SBATCH -J EnergyMinimization
#SBATCH -N 1
#SBATCH --ntasks-per-node 32 
#SBATCH -t 04:59:00 
#SBATCH -A welbornlab
​
~/Tinker/minimize Final_input.xyz_2 -k tinker.key 0.1 > Minimization.log

```
In addition to the output file `Minimization.log`, output system coordinates will be saved as `Final_input.xyz_3`. It is good practice to center the system before any simulation, so do it again before running the molecular dynamics:
`~/Tinker/xyzedit Final_input.xyz_3`, which will create a new file `Final_input.xyz_4`.

###Start the molecular dynamics simulation
Using the latest system coordinate file, write the submission script `launch_dynamics.sh` containing:

```sh
#!/bin/bash -l
#SBATCH -p normal_q
#SBATCH -J MolecularDynamics
#SBATCH -N 1
#SBATCH --ntasks-per-node 32 
#SBATCH -t 04:59:00 
#SBATCH -A welbornlab
​
~/Tinker/dynamic Final_input.xyz_4 -k tinker.key 1000000 1.0 1.0 4 300.00 1.0 > Dynamics.log
```

This will launch a 1,000,000 step simulation of 1fs timestep, saving frames every 1ps in the ensemble 4 (NPT) with temperature 300K and pressure 1atm. Note that NVT is ensemble 2, in which case you will have `1000000 1.0 1.0 2 300.00 ` at the end of the command line. 

#Quantum Espresso simulations

###Set up


1) Log onto TinkerCliffs by typing in Terminal

```sh
ssh username@tinkercliffs1.arc.vt.edu
```

 and type your password when prompted. Go to your WORK directory: `cd $WORK`. 
 
2) Create a directory where you will want to run your QuantumEspresso simulations. For example, create a directory named "ExampleBenzene" by typing `mkdir ExampleBenzene`. Go to the directory by typing `cd ExampleBenzene`. This is where you will put the Quantum Espresso input file that you will need to run your calculation. 

3) Below I provide a working example that you can copy and paste (we will see the meaning of the different section of the input file in a later document). To do that, open a new file (named input.in here) by typing `vi input.in`. Then press the key `i` to enter the "insert mode" of vi (it should say so at the bottom left of the terminal window). Once you are in the insert mode, you can copy the text below and paste with Cmd+c and Cmd+v. Press `esc` to exit the insert mode if needed and `:w` (write) to save your changes. 
 
 ```sh
  &CONTROL
    prefix='benzene',
    pseudo_dir = '/home/vwelborn/Espresso/pseudo', 
 /

 &SYSTEM
    assume_isolated      =  'martyna-tuckerman',
    
    ibrav = 6
    A = 11.0
    C = 7.0

    ecutwfc =  20.0,
    ecutrho =  200.0,

    nat =  12,
    ntyp =  2,

    nbnd = 16
 /

 &ELECTRONS
 /

ATOMIC_SPECIES 
   C   1.0   C.pbe-rrkjus.UPF 
   H   1.0   H.pbe-rrkjus.UPF 

ATOMIC_POSITIONS angstrom
   H   5.5000000   7.98563953   3.5
   C   5.5000000   6.89520922   3.5
   C   6.7089386   6.19812524   3.5
   H   7.6529918   6.74309454   3.5
   C   6.7089386   4.80187470   3.5
   H   7.6529918   4.25690561   3.5
   C   5.5000000   4.10479062   3.5
   H   5.5000000   3.01436043   3.5
   C   4.2910612   4.80187468   3.5
   H   3.3470081   4.25690556   3.5
   C   4.2910613   6.19812528   3.5
   H   3.3470082   6.74309458   3.5

K_POINTS gamma
 
 ```

Two things on the input file: 1) at the very top, there is a line that specify the path to the directory containing pseudopotentials. In this example, it refers to a directory named `pseudo` that is in another directory named `Espresso` in my home directory. You will need to change that for the path to the directory where you keep all the `.UPF` files. 2) You see that you need `.UPF` files in the `ATOMIC_SPECIES` section where `C.pbe-rrkjus.UPF` and `H.pbe-rrkjus.UPF` are listed. See below on how to do that. For now you can quit the input file by typing `:wq` in vi (write and quit). 

###Secure copy files or directory from your laptop to ARC
1) To get the `pseudo` directory and its content, the easiest is to copy over the directory I emailed you. Say you downloaded the `pseudo` directory from your emails. You should be able to access it via the Terminal by opening a new Terminal window and typing :
`cd Downloads/pseudo`. Type `ls` to list the content of the directory. You should see a list of `.UPF` files, including the two whose name are in the input file above. 

2) You want to copy this directory over to ARC. Type `cd ../` to be in the directory that contains `pseudo` (i.e. `Downloads` here) and type 

```sh
scp -r pseudo username@tinkercliffs1.arc.vt.edu:~/
```

3) After this is complete, you can log onto TInkerCliffs (`ssh username@tinkercliffs1.arc.vt.edu`) and type `ls`: you should see the directory `pseudo` listed. To know the path to write in the input file, go in the directory (`cd pseudo`) and type `pwd`. Copy the path printed on the Terminal window and paste it in place of the example on the line `pseudo_dir` in the input file. 


###Run a simulation

You will need a couple more files to run a simulation.

1) The first one, `environment_variables`, will be the same for all simulations: you can secure copy from your emails by going back to your `Downloads` folder and typing:

```sh
scp environment_variables username@tinkercliffs1.arc.vt.edu:/work/username/ExampleBenzene
```

2) Let's go back to `ExampleBenzene` (go back to your first Terminal window or open a new one, log onto Tinkercliffs and navigate to the directory). If you type `ls`, you should see the input file you created earlier and the newly copied `environment_variables` file.

3) The last thing you need is a script (with the extension `.sh`) to launch your simulation. That tells ARC how many cores or nodes to use for your simulation and how you run your simulation. You can create this script using `vi` just like we did for the input file. For example, if we call the script `launch_QE.sh`, we will type `vi launch_QE.sh`, go in insert mode and copy/paste the following:

```sh
#! /bin/bash
#
#SBATCH -t 00:10:00
#SBATCH -N1 --ntasks-per-node=4
#SBATCH -p dev_q
#SBATCH -J Espresso
#SBATCH -A personal 


#
module reset
module load QuantumESPRESSO

source ./environment_variables


# how to run executables
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"
$ECHO
$ECHO "  running pw.x as: $PW_COMMAND"
$ECHO


$PW_COMMAND < input.in > output.out
```

4) Note that the last line call the input file `input.in`: if you called it differently above, you will need to edit the name here. 

5) You will need to change the name of the allocation (after the `#SBATCH -A` on the 7th line of the file). The name should be available on ColdFront. If you are not sure of what it is, write a ticket on the ARC website and ask.

6) Once all is done you can launch the simulation by running
`sbatch launch_QE.sh`
This is an example and it is set up to run on the `debug` queue so it should be very fast (less than one minute). You will see that on time you will simulations that last hours and we will need to update the fields in the header of the script to account for that. 


###Running geometry optimization of periodic materials

1) Creat a directory for the calculation and copy the `enviroment_variables` and `launch_QE.sh` into the directory.

2) Copy and paste the below input file as discribed above.

```sh
&CONTROL
  prefix='Ni_bulk',
  pseudo_dir = '/home/vwelborn/Espresso/pseudo', 
/

&SYSTEM
   ibrav = 2
   A = 3.52

   ecutwfc =  50.0,
   ecutrho =  400.0,

   nat =  1,
   ntyp =  1,

   occupations = 'smearing',
   smearing = 'mv',
   degauss = 0.01,
   
   nspin = 2
   starting_magnetization(1)= 0.8
/

&ELECTRONS
/

ATOMIC_SPECIES 
  Ni  58.693  Ni.pbe-spn-kjpaw_psl.1.0.0.UPF  

ATOMIC_POSITIONS crystal
   Ni 0.0  0.0  0.0

K_POINTS automatic
   15 15 15   1 1 1

```

Again, change the path to the pseudopotentials directory according to your setting. 

3) Lanuch the simulation by running `sbatch launch_QE.sh`.

4) This calculation will compute a total energy for the unit cell specified, but will not change the unit cell parameters. This is what we call a single point calculation. To optimize the lattice parameters under the chosen theory, you can run several single point calculations with varying unit cell parameters (in our case, changing the value of 'A'), and fit the total energy to the lattice parameter to find the optimal geometry.


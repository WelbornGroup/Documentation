### Energy barrier calculation with Quantum Espresso

The nudged elastic band (NEB) is the method for finding saddle points and minimum energy paths between two conformations. 

In the test calculation, we will run NEB calculation on the inversion of NH3 molecule. The geometry of NH3 is from a geometry optimization calculation. In real case scenario, you should complete the geometry optimization calculations of your structures as instructed in previous documents. 

Create a directory `NH3_NEB` for this example calculation. In this directory, create the following files:

`input.in`:
```sh
BEGIN
BEGIN_PATH_INPUT
&PATH
  restart_mode = 'from_scratch'
  string_method = 'neb'
  nstep_path = 50
  ds = 1.D0
  opt_scheme = 'broyden'
  num_of_images = 25
  CI_scheme = 'no-CI'
  path_thr = 0.5D0
/
END_PATH_INPUT

BEGIN_ENGINE_INPUT
&CONTROL
  prefix = 'NH3'
  pseudo_dir = '/home/tychen/bin/EspressoPseudo'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 10
  
  nat = 4
  ntyp = 2
  
  ecutwfc = 50
  ecutrho = 400
  
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.01
/

&ELECTRONS
/

ATOMIC_SPECIES 
  N   14.007  N.pbe-rrkjus.UPF  
  H    1.008  H.pbe-rrkjus.UPF

BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS angstrom
   N    0.000000000    0.000000000    0.393329392     0   0   1
   H    0.440363454    0.828540737    0.000000000     1   1   0
   H    0.497678529   -0.795539773    0.000000000     1   1   0
   H   -0.937240770   -0.033064356    0.000000000     1   1   0
LAST_IMAGE
ATOMIC_POSITIONS angstrom
   N    0.000000000    0.000000000   -0.393329392     0   0   1
   H    0.440363454    0.828540737    0.000000000     1   1   0
   H    0.497678529   -0.795539773    0.000000000     1   1   0
   H   -0.937240770   -0.033064356    0.000000000     1   1   0
END_POSITIONS

K_POINTS automatic
   4 4 4   0 0 0   
END_ENGINE_INPUT
END   
```

`lanch_NEB.sh`:
```sh
#! /bin/bash
#
#SBATCH -t 4:00:00
#SBATCH --output="%j.%N.out"
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH --ntasks-per-node 24 
#SBATCH -J Espresso
#SBATCH -A welborn

module reset
module load QuantumESPRESSO

source ./environment_variables

NEB_COMMAND="$PARA_PREFIX $BIN_DIR/neb.x $PARA_POSTFIX"
$ECHO
$ECHO "  running neb.x as: $NEB_COMMAND"
$ECHO


    $NEB_COMMAND -inp input.in > output.out
```

Remember to change the directory of potential files and allocation name to your setting.

Launch the calculation by running `sbatch launch_NEB.sh`.

If the calculation finished succefully, the program should generate `output.out`, `NH3.axsf` and some other files. 

In `output.out`, similar to previous calculations, it contains the results for several iterations. In the section of the very last iteration, look for the line that contains `activation energy`, the energy barrier for the conformation change between the two images in `input.in` is listed. 

In `NH3.axsf`, the geometries of the structure's images along the path are listed. 
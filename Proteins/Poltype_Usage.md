# This is a Quick Guide to using Poltype2
(Please see the Poltype2 [install guide](./Poltype_Install.md) first to get it set up)

To utilize this software you will require 4 input files:

1. A structure file for your molecule (`.mol`)
2. An initiation file which contains the settings or commands for Poltype (`poltype.ini`)
3. A file that contains the paths to certain software that poltype needs (`paths.sh`)
4. The bash script to submit the poltype job to the queue on ARC (`run-poltype.sh`)


### Structure File

Generate a structure file for your molecule using Avogadro or glycam's carbohydrate builder

Poltype requires a `.mol` or `.sdf` structure file. 

### Ini File

Here is a basic example for the `poltype.ini` file:

```
structure=galactose.mol
atmidx=400
dontfrag=True
```

The documentation for settings in this file can be found on the [Poltype2 Github](https://github.com/TinkerTools/poltype2/blob/master/README/README_HELP.MD)

In this example the structure file is named, the starting atom index for the resulting parameters is set, and the fragmenter is turned off which leads to more memory being used in the QM calculations but avoiding an error that was causing poltype to crash during the fragmentation process.

### Path File

The `paths.sh` should contain the following: 

```
export PATH=/projects/welbornlab/Poltype2/TinkerEx/:$PATH
export GDMADIR=/projects/welbornlab/Poltype2/bin/
export PATH=/projects/welbornlab/Poltype2/bin/:$PATH
export PATH=/home/lnewman2/anaconda3/envs/xtbenv/bin/:$PATH
export PSI_SCRATCH=/localscratch/
```

The only change will be the fourth line, where the path should reflect the location within your own personal directory on ARC

### Poltype Job Submission Script

The `run-poltype.sh` is a bash script to submit the poltype job to the queue on ARC. The only changes should be made to the header for example if you want to submit to another CPU queue or request more time/resources

```
#!/bin/bash
#SBATCH -J Poltype2
#SBATCH -A welbornlab
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --time=3-00:00:00

# Run the example
echo "-------- Starting Poltype2: `date` -------"

source activate poltype_psi417

source paths.sh

python /projects/welbornlab/Poltype2/master/PoltypeModules/poltype.py

echo "------- Poltype2 has exited: `date` --------"
```

### Running Poltype

Make a new directory on ARC (`mkdir GAL_Poltype`) 

Copy the four input files to this directory (`cp poltype.ini GAL_Poltype`)

Move into the directory (`cd GAL_Poltype`) and submit the job to the queue (`sbatch run-poltype.sh`)



A successful poltype run will result in an OPEN_ME directory, the minimized structure (`final.xyz`), and the parameters (`final.key`)



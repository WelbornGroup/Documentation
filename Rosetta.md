#Generating a conformational ensemble with the Rosetta Software Package

### Download and install the software

Follow the instructions [here](https://www.rosettacommons.org/software/license-and-download) to register for a license (free of charge) and download/copy the corresponding binaries on Cascades.

### Generate a parameter file for the ligand (if applicable)

If a ligand (or other small molecule) is complexed with the protein, you will need a specific parameter file for Rosetta. 

1) Start by renaming the ligand atoms in your PDB file (necessary when converting to .mol). In order of appearance, name carbons C1 through CX, hydrogens H1 though HY, ... where X and Y is the total number of carbons and hydrogens, respectively. This also ensures that all ligand atoms are uniquely named, which you will likely need when derving force field parameters for Tinker. It is also good practice to draw the ligand structure with corresponding atom names to know what they are. 

2) Convert the renamed PDB file into a `.mol` file. This can be done with [OpenBabel](http://openbabel.org/wiki/Main_Page) or opening the PDB in [Avogadro](https://avogadro.cc/) and saving it as `.mol` file. Note that VMD for example does not support this format (do not use `.mol2`). 

3) Use the Python script from the Rosetta package to create the parameter file from the `.mol` file by typing:

`~/rosetta/main/source/scripts/python/public/molfile_to_params.py -n X77 -p X77 ligand.mol`

Here, `X77` is the 3 letter name of the ligand that is used as a residue name in the PDB file. Replace X77 by the appropriate 3 letter code for your molecule. The output looks something like:

```sh
Centering ligands at ( -20.096,   18.844,  -27.394)
Atom names contain duplications -- renaming all atoms.
WARNING:  structure contains double bonds but no aromatic bonds
  Aromatic bonds must be identified explicitly --
  alternating single/double bonds (Kekule structure) won't cut it.
  This warning does not apply to you if your molecule really isn't aromatic.
Total naive charge -2.255, desired charge 0.000, offsetting all atoms by 0.034
WARNING: fragment 1 has 67 total atoms including H; protein residues have 7 - 24 (DNA: 33)
WARNING: fragment 1 has 34 non-H atoms; protein residues have 4 - 14 (DNA: 22)
WARNING: fragment 1 has 9 rotatable bonds; protein residues have 0 - 4
Average 67.0 atoms (34.0 non-H atoms) per fragment
(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)
WARNING:  no root atom specified, using NBR atom instead.
Wrote PDB file X77_0001.pdb
Wrote params file X77.params
```

and creates the `X77.params` file that you will need later. 

### Backbone conformational ensemble

We use the [Backrub application](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/backrub) of Rosetta to sample backbone conformations. We will usually create between 25 and 50 different backbone conformations for each starting PDB. Parameters for the simulation can be written to a file with the `write_flags_backrub.py` Python script below.

```sh
#!/usr/bin/python
from os import listdir
import sys
import operator
import math
import os
if __name__=='__main__':
	a = os.getcwd()
	pdb=str(sys.argv[1])
	lig=str(sys.argv[2])
	b=listdir(a)
	ntrial=10000
	PivotResidues=range(1,215)
	nstruct=25
	print len(PivotResidues)	
#	Write out output file
	fin=open('flags_backrub','w')
	fin.write(str('-s ')+str("%s" %pdb)+'\n')
	fin.write(str('-extra_res_fa ')+str("%s" %lig)+'\n')
	fin.write(str('-ignore_unrecognized_res')+'\n')
	fin.write('\n')
	fin.write(str('-backrub:ntrials ')+str("%d" %ntrial)+'\n')
	fin.write('\n')
	fin.write(str('-pivot_residues '))
	for i in range(0,len(PivotResidues)):
		fin.write(str("%d" %PivotResidues[i])+' ')
	fin.write('\n'+'-pivot_atoms CA'+'\n')
	fin.write('-out:nooutput'+'\n')
	fin.write('-nstruct '+str("%d" %nstruct)+'\n')

```
Edit (i) the number of residues in the line setting the range for `PivotResidues` (going from 1 to #residue+1). If you just created this file, you probably need to make it an executable by typing `chmod +x write_flags_backrub.py`.

Execute the script by typing `./write_flags_backrub.py name_of_file.pdb X77.params`. This will create the file `flags_backrub` that contains the parameter for the simulation (CHECK IT).

Finally, launch the simulation on Cascades by typing `sbatch launch_backrub.sh` where `launch_backrub.sh` is:



```sh
#!/bin/bash -l
#SBATCH -p normal_q
#SBATCH -J Backrub
#SBATCH -N 1
#SBATCH --ntasks-per-node 32 
#SBATCH -t 06:59:00 
#SBATCH -A welbornlab

~/rosetta/main/source/bin/backrub.static.linuxgccrelease @flags_backrub
```


### Side chain conformational ensemble
We use the [Fixbb application](https://www.rosettacommons.org/docs/latest/application_documentation/design/fixbb) to repack and minimize the side chains around the new backbone conformations. This simulation will have to be run on each of the 25 (or however many) backrub output structures. Similarly to the backrub application, we write the parameters into a file, `flags_fixbb`, writen via the script `write_flags_fixbb.py`:

```sh
#!/usr/bin/python
from os import listdir
import sys
import operator
import math
import os
if __name__=='__main__':
	a = os.getcwd()
	pdb=str(sys.argv[1])
	lig=str(sys.argv[2])
	b=listdir(a)
	nstruct=1
#	Write out output file
	fin=open('flags_fixbb','w')
	fin.write(str('-s ')+str("%s" %pdb)+'\n')
	fin.write(str('-extra_res_fa ')+str("%s" %lig)+'\n')
	fin.write(str('-ignore_unrecognized_res')+'\n')
	fin.write(str('-in:file:fullatom')+'\n')
	fin.write('\n')
	fin.write(str('-resfile resfile.txt')+'\n')
	fin.write(str('-min_pack')+'\n')
	fin.write(str('-ex1')+'\n')
	fin.write(str('-ex2')+'\n')
	fin.write(str('-overwrite')+'\n')
	fin.write('-nstruct '+str("%d" %nstruct)+'\n')
```

Notice that you will need the file `resfile.txt` that contains two lines:

```sh
NATAA 
start

```
Here, `NATAA` tells `fixbb` to leave every amino acid at their position and only allows the rotamers to change (Dunbrack library). 

`NATRO` would instead leave the natural rotamer as well as the amino acid. 
Other options are available but they involve design (i.e. change) of amino acids: `ALLAA` will allow the full design of any amino acid while `PIKAA` followed by a list of single letter code restricts the design to just these amino acids (for example `1 A PIKAA NT`).

Launch the simulation by typing `sbatch launch_fixbb.sh` where `launch_fixbb.sh` is:


```sh
#!/bin/bash -l
#SBATCH -p normal_q
#SBATCH -J Fixbb
#SBATCH -N 1
#SBATCH --ntasks-per-node 32 
#SBATCH -t 4:59:00 
#SBATCH -A welbornlab

~/rosetta/main/source/bin/fixbb.static.linuxgccrelease @flags_fixbb > output.txt

```






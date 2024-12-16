# Steps to follow to parameterize a polymer.

#### Poltype2
For the parameterization, we will be using the Poltype2.

To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

## Tutorial
Here we will parameterize a beta-blucose dimer, to be able to parameterize a polymer with beta-glucose subunits. It looks like this:
![bgdimer](https://github.com/user-attachments/assets/e2cfb451-6813-42a1-8ef6-5f6e8d051ecf)

### Making the molecule

We will first make the molecule using PyMOL builder. You can download PyMOL [here](https://pymol.org/), and get an academic license if you need it.

The first step is to make a beta-blucose monomer and paramterizing it with Poltype2. Open PyMOL, click `Builder` on the top-right side, and build the glucose ring.

** HINT - You can use the 6-membered ring to start with, that way, you can get the hexagonal ring in one go. **

NOTE - If you think the bonda/angles does not look the way they is supposed to look, click the `Sculpt` button on the builder window to get an approximation of the lowest energy structure.
The beta-glucose monomer looks like this:
![monbuilder](https://github.com/user-attachments/assets/adc46ff0-94bc-482e-aef0-9cbd6ca0741a)
Once you're satisfied with the structure, Go to `File` -> `Export` and export it as a PDB.

NOTE - Since we will use `.sdf` files as inputs to Poltype2, export another copy as a SDF file for later.

However, you see that the `OH` group on the `C1` is not how it occurs in the dimer. This oxygen atom in the dimer (and similarly a polymer, if need be) is connected to a carbon atom, and not a hydrogen.
Since we would need to model the behaviour of that oxygen as it bonds to the next ring. Ideally, we would run poltype for a beta-glucose dimer, and get the parameters for the linkage from there. But, running poltype for the dimer requires too much time and memory to run correctly on ARC, so we improvise. We add a `CH3` group to the oxygen, and since that carbon is supposed to be tertiary, we add one `CH3` group on either side, to account for that environment. The molecule we end up parameterizing is this (check the functional group attached to the OH:
![bgm_oc](https://github.com/user-attachments/assets/8d04e930-14f5-40e4-9d70-686d3e4cca0c)


Save it as a PDB and an SDF file. (Here I've named them `BGM_OC`, and will be referring to this molecule as the "functionalized monomer".)

NOTE - Also make the beta-glucose dimer, and save it as a PDB file. We will use it later, when we check out parameters.

### Structure input files
So far, we have the following SDF files we will use as inputs:
The monomer file `BGM.sdf` looks like:
````sh
BGM
  PyMOL3.0          3D                             0

 24 24  0  0  0  0  0  0  0  0999 V2000
   -0.5700    1.1470   -2.6280 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9670    0.1070   -0.9140 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3920    0.2200   -1.4120 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6560   -0.6770    0.6530 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2830    0.6960   -0.2530 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5580   -1.0810    2.5410 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0860   -0.1860    0.9900 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8900    0.3020    2.0670 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3920   -0.2200    1.4120 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6560    0.6770   -0.6530 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2830   -0.6960    0.2530 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6910    1.7230   -2.9790 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0100   -1.5690    0.6260 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8800   -1.9360    2.2480 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9960    1.2530    1.9840 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090   -0.2090   -0.5440 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9550    0.5720   -3.4710 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9440    2.3680   -2.3150 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2750    1.9410   -2.3790 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7110   -0.8010   -1.7310 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0320    1.7550   -0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4340   -1.2240    0.7690 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7110    0.8010    1.7310 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0320   -1.7550    0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0
  1 12  1  0  0  0  0
  1 17  1  0  0  0  0
  1 19  1  0  0  0  0
  1  3  1  0  0  0  0
  2  3  1  0  0  0  0
  3  5  1  0  0  0  0
  3 20  1  0  0  0  0
  4 13  1  0  0  0  0
  5  7  1  0  0  0  0
  5 10  1  0  0  0  0
  5 21  1  0  0  0  0
  6 14  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 22  1  0  0  0  0
  8 15  1  0  0  0  0
  6  9  1  0  0  0  0
  9 11  1  0  0  0  0
  9 23  1  0  0  0  0
 10 16  1  0  0  0  0
  2 11  1  0  0  0  0
  4 11  1  0  0  0  0
 11 24  1  0  0  0  0
 12 18  1  0  0  0  0
M  END
>  <pdb_header>


$$$$

````

The new file for the functionalized monomer `BGM_OC.sdf` looks like:
````sh
bgm_oc
  PyMOL3.0          3D                             0

 33 33  0  0  0  0  0  0  0  0999 V2000
   -0.1991   -0.5323   -0.8443 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0904   -0.0563    0.3141 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8937   -0.9381    1.5577 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5841   -0.9724    1.9791 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4755   -1.4484    0.8208 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3761   -2.7577    1.1795 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3778    0.3947   -2.0610 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7375   -2.7909    1.8987 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4014   -3.7204    1.8829 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7622   -0.1805   -2.9033 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1366    1.6161   -1.7479 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0823    1.1883   -1.8120 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0726   -2.6885    2.8155 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8031    0.5004    2.5510 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8167   -0.9611    0.0237 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5063   -3.0654    0.1420 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2165   -3.3728    2.8993 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1166   -3.8127    1.9133 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6178   -2.4325    2.9211 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4428   -2.1498    1.3698 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8355   -4.7198    1.9134 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4608   -3.7501    1.3329 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1600   -0.6452   -0.3466 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8832    0.9706   -2.4118 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7501   -1.8332    3.1089 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6977   -0.4501    2.6349 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4631   -0.0754   -0.0861 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8483   -1.4293    1.2209 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5184   -1.5537   -1.1634 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8392    1.0023    0.5655 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2412   -1.9759    1.3361 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9034    0.0489    2.2982 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2243   -2.5070    0.5693 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 29  1  0  0  0  0
  2  3  1  0  0  0  0
  2 27  1  0  0  0  0
  2 30  1  0  0  0  0
  3  4  1  0  0  0  0
  3 26  1  0  0  0  0
  3 31  1  0  0  0  0
  4  5  1  0  0  0  0
  4 25  1  0  0  0  0
  4 32  1  0  0  0  0
  5 28  1  0  0  0  0
  5 33  1  0  0  0  0
  6  8  1  0  0  0  0
  6  9  1  0  0  0  0
  6 16  1  0  0  0  0
  7 10  1  0  0  0  0
  7 12  1  0  0  0  0
  7 24  1  0  0  0  0
  8 18  1  0  0  0  0
  8 19  1  0  0  0  0
  8 20  1  0  0  0  0
  9 17  1  0  0  0  0
  9 21  1  0  0  0  0
  9 22  1  0  0  0  0
  1 23  1  0  0  0  0
  5 23  1  0  0  0  0
 11 24  1  0  0  0  0
 13 25  1  0  0  0  0
 14 26  1  0  0  0  0
 15 27  1  0  0  0  0
  6 28  1  0  0  0  0
M  END
$$$$

````

### Poltype input files
Now we move on to using Poltype.

REMINDER - To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

From the Usage guidelines, we know that we need three more files other than the input sdf file. To parameterize the `BGM_OC.sdf` structure, I have the following three input files:

#### paths.sh
The first is the `paths.sh` file, which looks like this for me:
````sh
export PATH=/projects/welbornlab/Poltype2/TinkerEx/:$PATH
export GDMADIR=/projects/welbornlab/Poltype2/bin/
export PATH=/projects/welbornlab/Poltype2/bin/:$PATH
export PATH=/home/mondal/miniconda3/envs/xtbenv/bin/:$PATH
export PSI_SCRATCH=/localscratch/   

````
The only line you need to change is the `export PATH=/home/mondal/miniconda3/envs/xtbenv/bin/:$PATH` line, which will be the path to your `xtbenv/bin`.

#### poltype.ini
The second is the `poltype.ini` file, which looks like this for me:
````sh
structure=BGM_OC.sdf
atmidx=400
new_gdma=True
gdmacommand_Radius_S=0.80
prmmodfile=dma4_hfe2023
````
Check the first line - it is to have the name of the input `sdf` file. The second line is the first value of the atom types. Poltype will assign atom types here starting from 400. I have used 400 for functionalized beta-glucose, and 500 for the beta-glucose monomer but you are free to choose differently, if you have an existing parameter file that overlaps with this numbering.

#### run-poltype.sh
Finally, we have a submit script for Tinkercliffs on ARC, named `run-poltype.sh` which looks like this:
````sh
#!/bin/bash
#SBATCH -J bgm_oc
#SBATCH -A welbornlab
#SBATCH -p welborn_q
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --time=5-23:00:00
#SBATCH --mem=128G
# Run the example
echo "-------- Starting Poltype2: `date` -------"
source activate poltype_psi417
source paths.sh
python /projects/welbornlab/Poltype2/master/PoltypeModules/poltype.py
echo "------- Poltype2 has exited: `date` --------"
````
The details you should know from the Poltype Usage file.

#### Submitting the Job
Upload the four files on ARC (to a folder) and submit the job using the command `sbatch run-poltype.sh`. This will run Poltype for the dimer.
Assuming everything goes as planned, you should see a `final.xyz` and a `final.key` file in that folder. Download the contents of the folder.

TASK - Now repeat this step in a different folder for the monomer, and run Poltype on that one as well.

### Identifying the atom types
Now that we have the parameters for both the monomer and the functionalized monomer. Lets talk about the monomer first. Open the `final.xyz` file of the monomer, and look at the contents.
My `final.xyz` file for the monomer looks like this:
````sh
    24
     1  C      2.591438   -0.511873   -0.335508      500    3    12    17    19
     2  O      0.376997   -1.307845   -0.080451      501    3    11
     3  C      1.131526   -0.162896   -0.476801      502    1     2     5    20
     4  O     -1.657686   -2.326217    0.028653      503   11    13
     5  C      0.726455    1.026047    0.383543      504    3     7    10    21
     6  O     -2.915971    0.253632    0.176162      505    9    14
     7  C     -0.769743    1.256992    0.293519      506    5     8     9    22
     8  O     -1.106868    2.240439    1.273274      507    7    15
     9  C     -1.556303   -0.016723    0.515279      508    6     7    11    23
    10  O      1.421236    2.183361   -0.077005     509    5    16
    11  C     -1.010998   -1.142317   -0.343434     510    2     4     9    24
    12  O      2.900673   -1.467836   -1.344478     511    1    18
    13  H     -1.848650   -2.836579   -0.774144     512    4
    14  H     -3.468190   -0.282456    0.768040     513    6
    15  H     -1.843396    2.760216    0.912316     514    8
    16  H      1.473920    2.791107    0.679336     515   10
    17  H      3.180285    0.406509   -0.447465     516    1
    18  H      3.812661   -1.766067   -1.194504     517   12
    19  H      2.750134   -0.919583    0.672821     516    1
    20  H      0.935553    0.080759   -1.534662     518    3
    21  H      0.953252    0.820114    1.440731     519    5
    22  H     -1.013476    1.615998   -0.717499     520    7
    23  H     -1.454259   -0.354853    1.555096     521    9
    24  H     -1.143683   -0.927498   -1.417403     522   11
````
This is the format of a Tinker xyz file. Column 1 is the atom index, Column 2 is the atom name, Columns 3-5 are the x, y and z coordinates. Column 6 is the atom type, which we need, and the rest of the columns are connectivities. These are the atom indices of the atoms this one is connected to. Which means, atom 1 which is carbon of type 500, is connected to atoms 3, 12, 17 and 19.
Draw the structure of the monomer on paper, that will help you assign the atom types.

Check the bonding of the different atoms in the structure, and the xyz file, and find out which atom is assigned which number as the atom type.

HINT - There is only one carbon that is bonded to two hydrogens (which will be equivalent, i.e. having the same atom type), and there is only one oxygen bonded to two carbons.
This is what my monomer looks like:

![monomer](https://github.com/user-attachments/assets/4c0c7188-8b56-494c-843b-d8edff136147)


TASK - Do the same for the functionalized monomer. For that molecule, my `final.xyz` file looks like this:
````sh
    33 
     1  C     -0.960980    1.206418   -0.263471      400    2     7    23    29
     2  C     -1.970961    0.073744   -0.382255      401    1     3    27    30
     3  C     -1.306140   -1.231243   -0.823281      402    2     4    26    31
     4  C     -0.070436   -1.546024   -0.008426      403    3     5    25    32
     5  C      0.822548   -0.319565    0.062707      404    4    23    28    33
     6  C      2.948836    0.330620    0.907843      405    8     9    16    28
     7  C     -1.483568    2.481025    0.373418      406    1    10    12    24
     8  C      3.917535   -0.017958   -0.208521      407    6    18    19    20
     9  C      3.584998    0.267916    2.279976      407    6    17    21    22
    10  H     -2.047367    2.224056    1.282091     414    7
    11  H     -2.657055    3.948035   -0.126940     415   24
    12  H     -0.608261    3.071664    0.677573     414    7
    13  H      1.141022   -3.057654    0.082234     416   25
    14  H     -2.238731   -2.692437    0.073054     417   26
    15  H     -3.757624   -0.080162   -1.115076     418   27
    16  H      2.523042    1.327392    0.744315     419    6
    17  H      2.834970    0.461742    3.049275     420    9
    18  H      3.430215    0.012603   -1.186066     420    8
    19  H      4.752357    0.689294   -0.225868     420    8
    20  H      4.313389   -1.024579   -0.050850     420    8
    21  H      4.005001   -0.727406    2.449673     420    9
    22  H      4.388166    1.004609    2.367676     420    9
    23  O      0.097664    0.779383    0.604809     408    1     5
    24  O     -2.282116    3.172382   -0.575499     409    7    11
    25  O      0.623721   -2.625116   -0.617532     410    4    13
    26  O     -2.272409   -2.282879   -0.808904     411    3    14
    27  O     -2.958624    0.424489   -1.339580     412    2    15
    28  O      1.863732   -0.625926    0.930784     413    5     6
    29  H     -0.572032    1.448889   -1.265810     421    1
    30  H     -2.417946   -0.117889    0.607715     422    2
    31  H     -1.011139   -1.117505   -1.873023     423    3
    32  H     -0.325822   -1.763749    1.042349     424    4
    33  H      1.188047   -0.057090   -0.946331     425    5

````
and the atom types look like:
![bgm_oc](https://github.com/user-attachments/assets/a46f4343-f058-459e-bacc-8606057dae88)


### Building a parameter file
The file `final.key` contains the parameters for the regular and functionalized monomer. We need to copy them to a `.prm` file, which can be read by Tinker when we run and MD using the AMOEBA Force Field. But the `final.key` file also has a lot of comments, which start with a `#` symbol, which I prefer removing. To do so, I have a python script `copyfile.py`.
````sh
def copy_file_excluding_comments(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
input_file = 'final.key'
output_file = 'beta_glucose_test.prm' #change this line if you need a different filename (like for the BGM_OC)
copy_file_excluding_comments(input_file, output_file)
print(f"Contents of {input_file} have been copied to {output_file} excluding lines that start with #.")
````
You can create this file in the folder that you downloaded for the `BGM` and `BGM_OC`. You can do so by using `cd` to go to that folder, then running the command `vi copyfile.py`, then pressing `I` to insert, and then pasting these lines. Hit `Esc`, then `:wq` to save this file. You can then run it in the same folder using the command `python copyfile.py`.

Now I have a `beta_glucose_test.prm` file in the folder for the monomer that I will be using and adding to, to parameterize the dimer.
This file should have the atom types described in the beginning which look something like this (only showing the first few lines):
````sh
atom          500    500    C     "BGM                 "         6    12.011    4
atom          501    501    O     "BGM                 "         8    15.999    2
atom          502    502    C     "BGM                 "         6    12.011    4
atom          503    503    O     "BGM                 "         8    15.999    2
atom          504    504    C     "BGM                 "         6    12.011    4
````
Copy all the atom definitions, and add them to the list of atom definitions in the `amoebabio18.prm` file (after the atom types for amino acids, NAs, Ions etc, and before the `SOLUTE` parameters). 
For the other parameters in the `beta_glucose_test.prm` file (torsion, SOLUTE, polarize, mutipole, bond, angle, strbend, vdw - the order that I have in my test file), simply add all those lines at the bottom of the `amoebabio18.prm` file. Save this file, and rename it to `beta_glucose_params.prm`.
We will use this file to add new parameters, and run MD in Tinker.

### Adding Parameters for the Linkage
Save the pdb file for the dimer in a new folder, let's call it `BGD` and the file name for the dimer is `BGD.pdb`. Copy the `beta_glucose_params.prm` file to that folder as well. You'll need another file to run Tinker commands, called `tinker.key` which will let Tinker know about the parameter file and general system information. You can learn more about running Tinker commands [here](https://github.com/WelbornGroup/Documentation/blob/master/RunningTinkerBasics.md). For now, my `tinker.key` file looks like this:
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

It shows a cubic box of side length 50 A, and that the parameters can be found in `beta_glucose_params.prm`. We will keep the discussion on the other terms for some other time.

Now convert the pdb structure of the beta-glucose dimer to Tinker XYZ format. You can do that using Tinker's `pdbxyz` executable, with the input as `BGD.pdb` (the pdb file for the dimer), and the parameters as `beta_glucose_params.prm`. To do so, run the command `~/path-to-tinker/analyze BGD.pdb -k tinker.key`. Change the path to Tinker as needed. This will give you a `BGD.xyz` file, which will be the Tinker XYZ file for the dimer. Open the file in a text editor, and change the atom types of the atoms as per the atom types you have in the monomer.

Note that the dimer will have `2*n-3` atoms, where `n` is the number of atoms in the monomer. In our case, we have assumed that atoms with atom types 509 (O), 515 (H) and 512 (H) leave as a water molecule, so the bridging oxygen will have the atom type of 503.

(How do you know that is the assumption? - Because we wanted to parameterize what parameters the 503 oxygen will have, which is why we functionalized the glucose at that point.)

My 'BGD.xyz` file looks like this:
````sh
    45
    50.000000   50.000000   50.000000   90.000000   90.000000   90.000000
     1  C      6.152000   -3.680000   -2.161000   502     2     6    15    19
     2  C      5.175000   -3.272000   -1.031000   504     1     3    11    35
     3  C      5.551000   -3.985000    0.293000   506     2     4    17    20
     4  C      6.963000   -3.513000    0.716000   508     3     5    18    21
     5  C      7.978000   -3.809000   -0.430000   510     4    14    15    22
     6  C      5.741000   -3.227000   -3.571000   500     1     7     9    16
     7  H      4.760000   -3.638000   -3.812000   516     6
     8  H      6.466000   -3.346000   -5.394000   517    16
     9  H      5.684000   -2.134000   -3.593000   516     6
    10  H      9.582000   -3.900000    0.729000   512    14
    11  H      4.243000   -3.771000   -1.113000   519     2
    12  H      4.514000   -2.840000    1.535000   514    17
    13  H      7.404000   -5.128000    1.804000   513    18
    14  O      9.278000   -3.406000   -0.038000   503     5    10
    15  O      7.467000   -3.211000   -1.684000   501     1     5
    16  O      6.711000   -3.663000   -4.521000   511     6     8
    17  O      4.596000   -3.777000    1.326000   507     3    12
    18  O      7.367000   -4.173000    1.915000   505     4    13
    19  H      6.185000   -4.797000   -2.222000   518     1
    20  H      5.548000   -5.089000    0.130000   520     3
    21  H      6.954000   -2.425000    0.929000   521     4
    22  H      7.994000   -4.904000   -0.616000   522     5
    23  C      1.926000   -0.410000   -3.048000   502    24    28    36    41
    24  C      1.009000   -0.072000   -1.850000   504    23    25    38    42
    25  C      1.064000   -1.226000   -0.827000   506    24    26    39    43
    26  C      2.502000   -1.540000   -0.373000   508    25    27    40    44
    27  C      3.436000   -1.767000   -1.616000   510    26    35    36    45
    28  C      1.793000    0.609000   -4.184000   500    23    30    31    37
    29  H      2.537000    0.820000   -5.998000   517    37
    30  H      2.069000    1.597000   -3.824000   516    28
    31  H      0.767000    0.639000   -4.556000   516    28
    32  H     -0.713000   -0.667000   -2.672000   515    38
    33  H      0.588000   -0.123000    0.769000   514    39
    34  H      2.114000   -3.457000    0.030000   513    40
    35  O      4.699000   -2.226000   -1.146000   503     2    27
    36  O      3.284000   -0.606000   -2.496000   501    23    27
    37  O      2.643000    0.213000   -5.263000   511    28    29
    38  O     -0.355000    0.115000   -2.242000   509    24    32
    39  O      0.265000   -0.908000    0.315000   507    25    33
    40  O      2.413000   -2.667000    0.487000   505    26    34
    41  H      1.612000   -1.386000   -3.486000   518    23
    42  H      1.348000    0.888000   -1.404000   519    24
    43  H      0.636000   -2.139000   -1.312000   520    25
    44  H      2.950000   -0.702000    0.199000   521    26
    45  H      2.989000   -2.481000   -2.347000   522    27
````

Run Tinker `analyze`. We do that by running the command `~/path-to-tinker/analyze BGD.xyz -k tinker.key M > analyze.log`. Change the path to Tinker as needed. Option M generates the principal moments, which we will also give us an idea of the net charge on the molecule (which will help us later - since here it has to be zero).

Running this gives us an `analyze.log` file, which will probably read something like this:
````sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                       Version 8.10.1  October 2021                       ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2021               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 Undefined Bond Stretching Parameters :

 Type             Atom Names           Atom Classes

 Bond            2-C      35-O           504  503

 Undefined Angle Bending Parameters :

 Type                  Atom Names                   Atom Classes

 Angle           1-C       2-C      35-O           502  504  503
 Angle           3-C       2-C      35-O           506  504  503
 Angle          11-H       2-C      35-O           519  504  503
 Angle           2-C      35-O      27-C           504  503  510

 Undefined Torsional Parameters :

 Type                        Atom Names                        Atom Classes

 Torsion         6-C       1-C       2-C      35-O         500  502  504  503
 Torsion        15-O       1-C       2-C      35-O         501  502  504  503
 Torsion        19-H       1-C       2-C      35-O         518  502  504  503
 Torsion        35-O       2-C       3-C       4-C         503  504  506  508
 Torsion        35-O       2-C       3-C      17-O         503  504  506  507
 Torsion        35-O       2-C       3-C      20-H         503  504  506  520
 Torsion         1-C       2-C      35-O      27-C         502  504  503  510
 Torsion         3-C       2-C      35-O      27-C         506  504  503  510
 Torsion        11-H       2-C      35-O      27-C         519  504  503  510
 Torsion        26-C      27-C      35-O       2-C         508  510  503  504
 Torsion        36-O      27-C      35-O       2-C         501  510  503  504
 Torsion        45-H      27-C      35-O       2-C         522  510  503  504

 MECHANIC  --  Some Required Potential Energy Parameters are Undefined

 Tinker is Unable to Continue; Terminating the Current Calculation

````

Notice how all the errors we see here, are because of undefined parameters around the linkage atom. Check to see if you're missing any other parameters - if you do, you might have made a mistake while assigning the atom types.

The goal here is to look at each of the missing parameters, and try to find it soemwhere in either the parameters of the monomer, or the functionalized monomer. We will not go through all of them here, but I will show you one example of each.

Starting with the ` Bond            2-C      35-O           504  503` parameter. This line means that the parameter for the bond between atom of type 504 and atom of type 503 is missing. Looking at the structure of the monomer, this bond should be the same as the bond between atom of type 504 and atom of type 509. So we locate the line `bond   504   509   204.1189   1.4263` from `beta_glucose_test.prm` and copy-paste this line at the bottom of the file `beta_glucose_params.prm`. But we need to edit this line, so that it applies to the bond in question. So the new line becomes `bond   504   503   204.1189   1.4263`. Save this new version of the `beta_glucose_params.prm` file, and run `~/path-to-tinker/analyze BGD.xyz -k tinker.key M > analyze.log` again. Now you will see the `analyze.log` file no longer has the missing parameter, and looks like this:
````sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                       Version 8.10.1  October 2021                       ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2021               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 Undefined Angle Bending Parameters :

 Type                  Atom Names                   Atom Classes

 Angle           1-C       2-C      35-O           502  504  503
 Angle           3-C       2-C      35-O           506  504  503
 Angle          11-H       2-C      35-O           519  504  503
 Angle           2-C      35-O      27-C           504  503  510

 Undefined Torsional Parameters :

 Type                        Atom Names                        Atom Classes

 Torsion         6-C       1-C       2-C      35-O         500  502  504  503
 Torsion        15-O       1-C       2-C      35-O         501  502  504  503
 Torsion        19-H       1-C       2-C      35-O         518  502  504  503
 Torsion        35-O       2-C       3-C       4-C         503  504  506  508
 Torsion        35-O       2-C       3-C      17-O         503  504  506  507
 Torsion        35-O       2-C       3-C      20-H         503  504  506  520
 Torsion         1-C       2-C      35-O      27-C         502  504  503  510
 Torsion         3-C       2-C      35-O      27-C         506  504  503  510
 Torsion        11-H       2-C      35-O      27-C         519  504  503  510
 Torsion        26-C      27-C      35-O       2-C         508  510  503  504
 Torsion        36-O      27-C      35-O       2-C         501  510  503  504
 Torsion        45-H      27-C      35-O       2-C         522  510  503  504

 MECHANIC  --  Some Required Potential Energy Parameters are Undefined

 Tinker is Unable to Continue; Terminating the Current Calculation
````

Now you get the picture. We need to do the same for each line that shows up as an error. For angle where the missing parameter is for `A  B  C` we either try to find the relevant `A  B  C` line, or `C  B  A` line (same with torsions, for `ABCD`, we need `ABCD` or `DCBA`) any other option will not be for the relevant angle/torsion.

For `Angle           1-C       2-C      35-O           502  504  503`, I use the line `angle   502   504   509   91.7411   108.7595` and change it to `angle   502   504   503   91.7411   108.7595`.

For `Torsion         6-C       1-C       2-C      35-O         500  502  504  503`, I use the line `torsion 500 502 504 509 0.854 0.0 1 -0.374 180.0 2 0.108 0.0 3` and change it to `torsion 500 502 504 503 0.854 0.0 1 -0.374 180.0 2 0.108 0.0 3`.

NOTE: Make sure to COPY the old line and PASTE it at the BOTTOM of the `beta_glucose_params.prm` file. DO NOT EDIT THE PARAMETERS IN THEIR ORIGINAL LOCATION. Pasting it at the bottom allows you to quickly find the new/added parameters, should a different error arise.

Now, after I run the `analyze` command again, my `analyze.log` file looks like this:
````sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                       Version 8.10.1  October 2021                       ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2021               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 Undefined Angle Bending Parameters :

 Type                  Atom Names                   Atom Classes

 Angle           3-C       2-C      35-O           506  504  503
 Angle          11-H       2-C      35-O           519  504  503
 Angle           2-C      35-O      27-C           504  503  510

 Undefined Torsional Parameters :

 Type                        Atom Names                        Atom Classes

 Torsion         6-C       1-C       2-C      35-O         500  502  504  503
 Torsion        15-O       1-C       2-C      35-O         501  502  504  503
 Torsion        19-H       1-C       2-C      35-O         518  502  504  503
 Torsion        35-O       2-C       3-C       4-C         503  504  506  508
 Torsion        35-O       2-C       3-C      17-O         503  504  506  507
 Torsion        35-O       2-C       3-C      20-H         503  504  506  520
 Torsion         1-C       2-C      35-O      27-C         502  504  503  510
 Torsion         3-C       2-C      35-O      27-C         506  504  503  510
 Torsion        11-H       2-C      35-O      27-C         519  504  503  510
 Torsion        26-C      27-C      35-O       2-C         508  510  503  504
 Torsion        36-O      27-C      35-O       2-C         501  510  503  504
 Torsion        45-H      27-C      35-O       2-C         522  510  503  504

 MECHANIC  --  Some Required Potential Energy Parameters are Undefined

 Tinker is Unable to Continue; Terminating the Current Calculation
````

Now perform this for every one of the missing parameters. 

NOTE: You can also look for the parameter line in the file we have for the functionalized monomer. Check to see the relevant parameter in the functionalized monomer has similar values for the one you have in the monomer. 

ADDITIONAL NOTE: When you copy a line from the parameters for the functionalized monomer, make sure to change ALL atom types in that line to the one that is in your error.

Once you're done with all the lines, your `analyze.log` file should look like this:
````sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                       Version 8.10.1  October 2021                       ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2021               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 Total Electric Charge :                  0.11950 Electrons

 Dipole Moment Magnitude :                  4.897 Debye

 Dipole X,Y,Z-Components :                 -3.717        0.179       -3.183

 Quadrupole Moment Tensor :               -27.520       -9.243       26.610
      (Buckinghams)                        -9.243       25.035       -6.328
                                           26.610       -6.328        2.485

 Principal Axes Quadrupole :              -43.409       10.959       32.450

 Radius of Gyration :                       3.847 Angstroms

 Center of Mass Coordinates :            4.244217    -2.188986    -1.365000
 Euler Angles (Phi/Theta/Psi) :           -29.673      -16.755      -73.823

 Moments of Inertia and Principal Axes :

             Moments (amu Ang^2)            X-, Y- and Z-Components of Axes

                   1319.359              0.831981     0.378482     0.405659
                   3458.028             -0.474028     0.105013     0.874225
                   4487.027              0.288279    -0.919632     0.266780
````

Now we find out why we used the `analyze` function of Tinker. Finding out the missing parameters could also be achieved with the `minimize` function, but that wouldn't give us the total charge of the molecule. Here we see that the total charge on the dimer is `0.11950 Electrons`, whereas it should actually be zero. This is because of the multipole parameters, which contain the partial charge for the atom. These partial charges are added up to get the total charge of the system.

Take for example the line:
````sh
multipole   500  511  502               0.03477
                                        0.22616    0.00000    0.15847
                                        0.00263
                                        0.00000   -0.12909
                                       -0.16273    0.00000    0.12646
````
Look at the first line. It tells us that these are the multipole parameters for atom of type 500, and the partial charge on this atom is `0.03477` (electrons). Right now, we will not discuss the rest of the numbers and how the dipole/multipole axes are chosen, but the only thing we will talk about is that the point for atom 500 is taken into account using atoms of types 511 and 502. Most commonly, these are atoms which are bonded to the atom 500.

That tells us that there is one type of connection that hasn't been taken into account when the multipoles based on the electronic environment of that atom were calculated. We know that we created the bond (by adding the parameters) between atoms of type 504 and 503. Looking at the multipoles for those two atoms we find (and the other atom connected to that oxygen):
````sh
multipole   503  510  512              -0.35743
                                        0.11505    0.00000    0.25337
                                        0.13440
                                        0.00000   -0.52984
                                       -0.30603    0.00000    0.39544
multipole   504  509  519               0.03073
                                       -0.11875    0.00000    0.00335
                                       -0.32144
                                        0.00000    0.06846
                                        0.19475    0.00000    0.25298
multipole   510  501  503               0.15806
                                        0.12817    0.00000    0.30494
                                        0.01311
                                        0.00000   -0.16794
                                        0.03297    0.00000    0.15483
````

We see that the multipoles for atom 503 depends on 512, and that of 504 depends on 509. Looking at the equivalent parameters in the functionalized monomer, we see:
````sh
multipole   413  404  405              -0.43986
                                        0.16851    0.00000    0.19264
                                        0.02689
                                        0.00000   -0.66297
                                       -0.49595    0.00000    0.63608
multipole   405  413  419               0.00122
                                       -0.02026    0.00000    0.36625
                                       -0.19893
                                        0.00000   -0.02011
                                        0.07995    0.00000    0.21904
multipole   404  408  413               0.16018
                                        0.26937    0.00000    0.18525
                                        0.01466
                                        0.00000   -0.30277
                                       -0.02986    0.00000    0.28811
````
So we see that we do not have any parameters for the bridging oxygen when it is connected to two carbons (like the `multipole   413  404  405` line). To remedy that, we copy that line, and append it to the bottom of the `beta_glucose_params.prm` file, and we will edit the partial charge to be:
````sh
multipole   503  510  504              -0.11950
                                        0.16851    0.00000    0.19264
                                        0.02689
                                        0.00000   -0.66297
                                       -0.49595    0.00000    0.63608
````

Note that the partial charge here balances out the extra charge that we saw in the analyze file. Run `Tinker/analyze` again, and now we can see that the charge has been balanced and is zero. The new `analyze.log` file looks like this:
````sh

     ######################################################################
   ##########################################################################
  ###                                                                      ###
 ###            Tinker  ---  Software Tools for Molecular Design            ###
 ##                                                                          ##
 ##                       Version 8.10.1  October 2021                       ##
 ##                                                                          ##
 ##               Copyright (c)  Jay William Ponder  1990-2021               ##
 ###                           All Rights Reserved                          ###
  ###                                                                      ###
   ##########################################################################
     ######################################################################


 Total Electric Charge :                 -0.00000 Electrons

 Dipole Moment Magnitude :                  6.057 Debye

 Dipole X,Y,Z-Components :                 -4.931        0.687       -3.449

 Quadrupole Moment Tensor :               -28.384       -7.732       26.213
      (Buckinghams)                        -7.732       25.927       -6.097
                                           26.213       -6.097        2.457

 Principal Axes Quadrupole :              -43.575       11.691       31.884

 Radius of Gyration :                       3.847 Angstroms

 Center of Mass Coordinates :            4.244217    -2.188986    -1.365000
 Euler Angles (Phi/Theta/Psi) :           -29.673      -16.755      -73.823

 Moments of Inertia and Principal Axes :

             Moments (amu Ang^2)            X-, Y- and Z-Components of Axes

                   1319.359              0.831981     0.378482     0.405659
                   3458.028             -0.474028     0.105013     0.874225
                   4487.027              0.288279    -0.919632     0.266780
````

TASK - Make a trimer, convert to Tinker XYZ, and edit the atom types. Run `Tinker/analyze` for the trimer without the added multipole line. You should see a charge of `0.23900 Electrons`, which is double the charge we saw for the dimer. Adding the new multipole parameter for atom type 503 (the same thing we just added), will bring the total charge to zero. You can confirm this by using this latest version of the `beta_glucose_params.prm` file for a polymer of any `n` number of beta-glucose monomers, and the charge should still come out to be zero.

### Running Minimization and Dynamics

Now that you have a parameter file, you can run `Tinker/minimize` and `Tinker/dynamic` on the molecule (or you can solvate it first) by following the steps [here](https://github.com/WelbornGroup/Documentation/blob/master/RunningTinkerBasics.md).

## Parameterizing Polymers
You should be able to parameterize any polymer you want. The steps to follow are:

1. Run poltype on the monomers, and cap the ends like we did with the functionalization. The caps should reflect the type of bonding in the polymer.

2. Make the polymer step by step, and parameterize ONE linkage at a time. So if your monomers are X,Y and Z and the sequence for the polymer is XYZYZX then make the molecules in this order: X, then XY, then XYZ, then XYZY then XYZYZ and finally XYZYZX.

3. At each step, make sure the missing parameters are added, and that the charge adds up to the total charge of the polymer at that step.

4. Run a 1ns MD for each step of the polymer, so that you know the MD runs without any issues.

Happy Parameterizing!

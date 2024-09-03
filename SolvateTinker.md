# Solvate a protein from tinker xyz format

This assumes that you have a TINKER XYZ file of your complete PROTONATED protein chain (and ligand if applicable). If you are missing residues or protons go back and follow the instructions in [ProteinPrep](./ProteinPrep.md) or go to [PDBtoTinkerXYZ](./PDBtoTinkerXYZ.md) to convert your file type.

This is the alternative method to using the gromacs solvate to create your waterbox and ions around your protein. This method is best for small molecules, perhaps paramterized from Poltype2, or proteins already in tinker xyz format that you do not wish to convert back to pdb

If you have a pdb file of your protein handy, we will still use pacmol to get some general information about our system that will be useful for the solvation protocol. Otherwise you will need to figure out the size of the box needed, and charge of your protein manually. 

### Pacmol - system charge and box size
The charge of the protein is computed by adding the charge of all residues. Only HIS, ARG, LYS, GLU and ASP are charged so the problem is equivalent to counting how many of these residues you have in your protein. 

There is a very simple tool to help you do that, which you can download with the [Packmol software](http://m3g.iqm.unicamp.br/packmol/download.shtml). Once you have downloaded the binaries, look for the executable `solvate.tcl`

Then navigate to the directory where you have your protein chain PDB and execute the script:

```sh
cd $PATH_to_directory
~/Downloads/packmol/solvate.tcl 1kjl_complete.pdb

```
Note that it may be a good idea to move the binaries from the `Downloads` folder. Don't forget to update the path to `solvate.tcl` if that is the case. 

This should return the total structure charge and other information:

```
  ###########################################################################
 solvate.tcl script: Solvates a (protein) structure with water and ions,
                     using Packmol. 
 ###########################################################################
 Input pdb file to be solvated: 1kjl_complete.pdb 
 User set structure charge: none  - will compute from structure. 
 User set shell size: none  - using default: 15.0 Angs. 
 User set density: none  - using default: 1.0 g/ml 
 -------------------------------------------------------
  Minimum and maximum coordinates of structure atoms 
  X_min = -12.695     X_max = 26.753 
  Y_min = -22.049     Y_max = 23.335 
  Z_min = -24.602     Z_max = 13.634 
 -------------------------------------------------------
  Box side length in each direction: 
  x: 69.44800000000001 
  y: 75.384 
  z: 68.236 
 -------------------------------------------------------
 -------------------------------------------------------
  HIS = 0 (associated charge = 0) 
  ARG = 9 (associated charge = +9) 
  LYS = 8 (associated charge = +8) 
  GLU = 6 (associated charge = -6) 
  ASP = 7 (associated charge = -7) 
 -------------------------------------------------------
  Total structure charge = 4
 -------------------------------------------------------
 Unsure about mass of element HG3: 1.00800. Is this correct? (y/n)
```


Type 'y' and hit enter for any hydrogens it cannot identify


```
 -------------------------------------------------------
  Molar mass of structure: 16772.934399999976
 -------------------------------------------------------
  Number of water molecules to be put:   11016 
  Total volume: 357233.75 A^3
  Volume occupied by water: 329372.54 A^3 
  Number of Sodium ions to be put: 28
  Number of Chloride ions to be put: 32
  Wrote packmol input. 
 -------------------------------------------------------
  Use these lengths for periodic boundary conditions: 
  x: 70.44800000000001
  y: 76.384
  z: 69.236
 -------------------------------------------------------
 -------------------------------------------------------
  Now, run packmol with: packmol < packmol_input.inp       
 -------------------------------------------------------
  The solvated file will be: solvated.pdb 
 -------------------------------------------------------
```

We use packmol for the charge information, number of ions to add, and to get a box size estimate. 

Our system has a positive charge of 4 meaning to neutralize it we will need to add 4 negative ions, probably chloride ions. 

We also get a nice measurement of the protein dimensions.You can manually go into VMD to measure your system size as well. We will need to know how big the protein is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. I will usually round up the largest side that pacmol recommends and use that for all sides of my box. 


For the galectin-3 system I use an 80x80x80 box, which is slightly overestimating the box size but you want to make sure that the protein is never interacting with its mirror images across periodic boundaries. ***This is especially tricky for flexible molecules!*** The downside to a larger volume box be the longer simulation time that is required. 

### Make a waterbox

Make a file named `water.xyz` in your directory with your protein file and amoeba parameters
```
     3
     1  O     -0.054822    0.026553   -0.013027   349     2     3
     2  H      0.199460   -0.864828   -0.246637   350     1
     3  H      0.670683    0.443372    0.453400   350     1
```

We will use this single water to populate our waterbox with many waters. First we need to calculate the number of waters to add for a specific volume based on the density of water 1000 kg/m^3. An example of this calculation is done on page 21 of this [tinker tutorial](https://tinker-hp.org/wp-content/uploads/2022/10/Tinker_preparation_tutorial.pdf).

For our 80x80x80 box we will need 17,129 water molecules. 


Now in your directory with your tinker xyz protein file, open a terminal and run `/Tinker/xyzedit water.xyz` and name the parameter file when prompted `amoebabio18.prm`

You should see the pop-up menu for xyzedit, read through the various ways this tool can be used to manipulate your structure file. (Note: different versions of tinker may have different numbering for the commands in xyzedit)

Select `(24) Create and Fill a Periodic Boundary Box`

`17129` when prompted for the number of copies

`80,80,80` to indicate the box size

`Y` to refine the box and minimize the waters

`amoebabio18.prm` provide the parameters once more for the minimization

This will result in a water.xyz_2 file which can be renamed to waterbox.xyz or waterbox_80.xyz for future use. This box can also be trimmed down using xyzedit if a smaller box is needed later. 

To check that the calculations were correct we could use `/Tinker/analyze waterbox.xyz`  and select `General System and Force Field Information [G]`

*The system density should be close to 1* 


### Insert Protein

Now in your directory with your tinker xyz protein file, open a terminal and run `/Tinker/xyzedit protein_filename.xyz`

Select `(13) Translate Center of Mass to the Origin`

This translates our carbohydrate to the origin, (centering it on x=0,y=0,z=0) and it is good practice to center your structures before carrying out other manipulations

Next select `(24) Soak Current Molecule in Box of Solvent`

This will prompt you for your solvent box name: waterbox.xyz

There will be a new structure file generated (protein_filename.xyz_2) which will be your centered carbohydrate in the minimized waterbox

You should download and open this file in VMD to make sure it looks okay with no obvious errors


### Adding Ions

To add ions `/Tinker/xyzedit protein_filename.xyz_2` then `(26) Place Monoatomic Ions around a Solute`

This will prompt you with  *Enter Atom Numbers in Solute Molecules :*

This is looking for any atom lines that are NOT water. This includes your protein, ligand, and other ions!

To answer this prompt you can use tinker's range syntax, if your protein has 2230 atoms then you can input `-1 2230`

Next prompt will be *Enter Ion Atom Type Number & Copies to Add :*

Input should look like `352 29` to add 29 Na+ ions to the system, resulting in the protein_filename.xyz_3

Repeat this process now to add chloride ions:

`/Tinker/xyzedit protein_filename.xyz_3` then `(26) Place Monoatomic Ions around a Solute`

Solute range now needs to include the new ions! `-1 2230 -51005 51033`

And the atom type will be different for Chloride `363 33`

With that you can pull up the final protein_filename.xyz_4 in VMD to visually check the box! Perhaps duplicate and rename this file to `protein_solv.xyz`









# Solvate a polymer from TinkerXYZ

This assumes that you have a TinkerXYZ file for your polymer. If you don't, export the structure from PyMOL as a PDB, and then use Tinker's `pdbxyz` executable to convert it to TinkerXYZ file format. (You need to have a `.key` file for that, which can be the same as that mentioned in [parameterization-for-amoeba.md](./parameterization-for-amoeba.md) - you'll have to change the box size (mentioned below).

## System charge and box size
The charge of the polymer is computed by adding the charge of all monomers. Since you've already parameterized the monomers (or have the parameters from somewhere), you know the charge of each monomer (***do not forget to add the charges at the ends of the polymer, if they exist!!***). Add them up to get a total charge of the system.

For the box size, open the structure of the polymer you made in PyMOL, and get an approximate idea of the size of the polymer by calculating the distance between the farthest atoms across each of the X,Y,Z axes. For that, there is a nice PyMOL tutorial [here](https://www.compchems.com/how-to-measure-distances-and-display-interactions-in-pymol/#distances).

We will need to know how big the polymer is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. It is generally better to slightly overestimate the box size but you want to make sure that the polymer is never interacting with its mirror images across periodic boundaries (unless that is specifically something you're trying to simulate). ***This is especially tricky for flexible molecules!*** The downside to a larger volume box be the longer simulation time that is required. 

## Solvate the polymer on ARC/Local
There are a number of routines available to solvate a protein (Pacmol, GROMACS etc). Here we will use Tinker's `xyzedit` executable.

### Make a waterbox

Make a file named `water.xyz` in your directory with your polymer file and amoeba parameters.
```
     3
     1  O     -0.054822    0.026553   -0.013027   349     2     3
     2  H      0.199460   -0.864828   -0.246637   350     1
     3  H      0.670683    0.443372    0.453400   350     1
```

We will use this single water to populate our waterbox with many waters. First we need to calculate the number of waters to add for a specific volume based on the density of water `1000 kg/m^3`. An example of this calculation is done on page 21 of this [Tinker tutorial](https://tinker-hp.org/wp-content/uploads/2022/10/Tinker_preparation_tutorial.pdf).

Let's assume for our `50x50x50` box we will need `11,000` water molecules ***you need to change this number to what you calculate***. 

I prefer to do this in Tinker8 on my local machine, so the commands here are in Tinker8. For information on these commands in Tinker9, look into [running-tinker9-basics.md](./running-tinker9-basics.md)

Now in your directory with your TinkerXYZ polymer file, open a terminal and run `/Tinker/xyzedit water.xyz` and name the parameter file when prompted `beta_glucose_params.prm`.

You should see the pop-up menu for `xyzedit`, read through the various ways this tool can be used to manipulate your structure file. (Note: different versions of tinker may have different numbering for the commands in `xyzedit`).

Select `(24) Create and Fill a Periodic Boundary Box`.

`11000` when prompted for the number of copies.

`50,50,50` to indicate the box size.

`Y` to refine the box and minimize the waters.

`beta_glucose_params.prm` provide the parameters once more for the minimization.

This will result in a `water.xyz_2` file which can be renamed to `waterbox.xyz` or `waterbox_50.xyz` for future use. This box can also be trimmed down using `xyzedit` if a smaller box is needed later. 

To check that the calculations were correct we could use `/Tinker/analyze waterbox.xyz` and select `G`, you can also choose option `M` like we did for parameterization, and check to see that the charge on the waters should be EXACTLY 0.

*The system density should be close to 1*.

### Insert Polymer

Now in your directory with your TinkerXYZ polymer file, open a terminal and run `/Tinker/xyzedit polymer_filename.xyz`.

Select `(13) Translate Center of Mass to the Origin`.

This translates our polymer to the origin, (centering it on x=0,y=0,z=0) and it is good practice to center your structures before carrying out other manipulations.

Next select `(24) Soak Current Molecule in Box of Solvent`.

This will prompt you for your solvent box name: `waterbox_50.xyz`

There will be a new structure file generated (`polymer_filename.xyz_2`) which will be your centered polymer in the minimized waterbox.

You should download and open this file in VMD to make sure it looks okay with no obvious errors.


### Adding Ions

Use Tinker `analyze` on the new file (`polymer_filename.xyz_2`) to confirm the charge on the system.

Based on the charge of the system, you'll have to add positive or negative counterions in the system. You might also need some concentration of salts on a case-to-case basis, each ion species will have to be added based on the following steps. ***Make sure to do that calculation for number of ions of each species BEFORE you do any of the steps mentioned below.***

To add ions `/Tinker/xyzedit polymer_filename.xyz_2` then `(26) Place Monoatomic Ions around a Solute`.

This will prompt you with  *Enter Atom Numbers in Solute Molecules :*

This is looking for any atom lines that are NOT water. This includes your polymer, protein, small molecules, ligand, and other ions!

To answer this prompt you can use Tinker's range syntax, if your SOLUTE (polymer + anything else not water) has 2230 atoms then you can input `-1 2230`.

Next prompt will be *Enter Ion Atom Type Number & Copies to Add :*.

Input should look like `352 29` to add 29 Sodium (Na+) ions to the system, resulting in the `polymer_filename.xyz_3`.

Repeat this process now to add chloride ions:

`/Tinker/xyzedit protein_filename.xyz_3` then `(26) Place Monoatomic Ions around a Solute`.

Solute range now needs to include the new ions! `-1 2230 -51005 51033`.

And the atom type will be different for Chloride (Cl-) `363 33` (assuming I'm adding `33` Chloride ions).

With that you can pull up the final (`polymer_filename.xyz_4`) in VMD to visually check the box! Duplicate and rename this file to `polymer_solv.xyz`. This will be your input file for whatever MD you might run. We usually minimize this system before we run the MD, but it depends on the project you're working on.

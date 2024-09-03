 ## Protein File
To begin, head over to the [RSCB Protein Data Bank](https://www.rcsb.org) to find your protein of interest. This is where many protein and DNA structures have been deposited by experimentalists and where many of our simulation projects begin! 

Search for the protein by name (**galectin-3**) or better by its unique PDB ID (**1KJL**). 

It is usually helpful to read the information about the protein entry on the RSCB website page and use the built in visualization tools to view the system. Look out for the number of modeled residues, you may need to add/remove amino acids, and any complexed molecules. For our example we see that there are 138 residues and a LacNAc molecule is present in the binding site of the protein. 

Download the file in a `.pdb` format. 

### Clean the Protein
The protein file will usually include crytallography artifacts such as water, ions, solvent molecules, etc. There may also be multiple chains of the protein that are unnecessary. We will want to remove these extra molecules to obtain a 'clean' protein. The bound ligand (LacNAc) can also be removed at this step to make a protein-only system. 

***Ligands make things more complicated:*** If you wish to include the ligand in the bound position from the crystal structure you can keep it in the pdb file however there will be some complications when converting to a tinker xyz (.txyz) file. Most solutions involve manually labeling the individual atoms and transferring the xyz coordinates to the correct atoms/atom types. If you have recently parameterized the molecule in Poltype2, you can take the poltype output txyz file, transfer the coordinates from the pdb to the txyz, using atom labels to match up the correct atom lines. You can then delete the ligand for the rest of the protein preparation and only add it after you have your protein in a txyz file. Paste the ligand at the end of the file and just renumber the index. Python scripts can help automate this. Note this wont work for every situation, for example if you are solvating with gromacs you may want the ligand in the pbd file when the water is added. 

You can **clean the protein** file in multiple ways:
1. Open the `.pdb` file in a software such as ChimeraX. Here you can use the selection tools to specify the unwanted molecules and remove them. This is a beginner friendly way to clean a file and is the most intuitive since you can visualize what is being removed.

2. Alternatively, open the `.pdb` file in a text editor and manually delete the unwanted atom lines from the file. This is a more hands-on approach but will allow you to familiarize yourself with the internal structure of these file types. The top of the PDB contain a header with information about the protein. You can remove this portion for the working file after you have read through the important sections. The `ATOM` lines signify the start of the actual protein structure. Normally the protein atoms will be listed first, followed by the water and other extraneous molecules. The file will sometimes end with a connectivity section, especially for the non-protein molecules. You will want to keep the proper 'ATOM' lines and the 'END' line for your cleaned protein, while removing the rest. You can load the pdb file into VMD to check and save it again to renumber it if needed. It is good practice to save new files when changes are made, this will allow you to compare the files to see the specific chagnes and have a backup to revert to if needed.


```
ATOM      1  N   PRO A 113      25.131  -4.504 -10.006  1.00 31.14           N  
ATOM      2  CA  PRO A 113      23.749  -4.861  -9.555  1.00 30.32           C  
ATOM      3  C   PRO A 113      23.547  -4.400  -8.112  1.00 29.57           C  
ATOM      4  O   PRO A 113      24.520  -4.262  -7.348  1.00 30.41           O  
ATOM      5  CB  PRO A 113      23.669  -6.314  -9.512  1.00 31.16           C  
...
ATOM   1130  C   ILE A 250      10.599  18.442 -13.748  1.00 26.64           C  
ATOM   1131  O   ILE A 250      10.593  19.149 -12.722  1.00 26.31           O  
ATOM   1132  CB  ILE A 250      10.536  16.065 -14.493  1.00 25.19           C  
ATOM   1133  CG1 ILE A 250       9.183  15.838 -13.813  1.00 25.15           C  
ATOM   1134  CG2 ILE A 250      11.338  14.758 -14.571  1.00 25.43           C  
ATOM   1135  CD1 ILE A 250       8.149  15.133 -14.690  1.00 26.48           C  
END                                                                             
```

Note the protein should have 1135 atoms once cleaned, for our galectin-3 example.

Save the new pdb separate under a name like ```ikjl_clean.pdb```



### Add Hydrogens

Most pdb files will need polar hydrogens added, if not all hydrogens added. This can be done in software like ChimeraX, however our group uses the '[reduce](https://github.com/rlabduke/reduce/blob/master/README.md)' program from the Richardson lab for predicting protonation states. There is also a [reduce2](https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce) out now which will continue to have support, however it needs to be installed as part of the [cctbx](https://github.com/cctbx/cctbx_project/tree/master) package. The commands may be slightly different using the reduce2 version. 

Download a version of reduce, open a terminal, and change to the directory where your pdb file is located.

```
reduce -build 1kjl_clean.pdb > 1kjl_prot.pdb
```

The build argument is to make sure the histidines are correctly protonated. The histidines are the most important to pay attention to during the protonation step and will decide the overall charge of your protein, along with your charged residues (Arg, Lys, Asp, Glu). 

Open the output.pdb and you should see new lines corresponding to the hydrogens atoms that were added. Note however that these atoms are not numbered sequentially. It is good practice to open the file in VMD (or similar software) to CHECK the structure and save it again as ```1kjl_complete.pdb```, this will also fix the numbering issue. ***Always check your pdb files regularly in both VMD and text editors to make sure problems are found early in the process.***

Note that our galectin-3 protein now has 2289 atoms.

You will want to keep an eye on the total number of atoms in your protein system at this point, record it, and keep track of it in the following steps (solvation & tinker xyz conversion). Sometimes the number of atoms can change during these steps due to protonation changes and you will want to know when this happens.


### Structure Modification
If the protein is missing residues you can add them using software like ([Modeller](AddMissingResidues.md)). Tinker also has an amino acid builder under the 'protein' command. If a structured domain needs to be built then homology modelling like Swiss-Model or AlphaFold can be useful.





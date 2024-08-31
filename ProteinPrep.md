 ## Protein File
To begin, head over to the [RSCB Protein Data Bank](https://www.rcsb.org) to find your protein of interest. This is where many protein and DNA structures have been deposited by experimentalists and where many of our simulation projects begin! 

Search for the protein by name (**galectin-3**) or better by its unique PDB ID (**1KJL**). 

It is usually helpful to read the information about the protein entry on the RSCB website page and use the built in visualization tools to view the system. Look out for the number of modeled residues, you may need to add/remove amino acids, and any complexed molecules. For our example we see that there are 138 residues and a LacNAc molecule is present in the binding site of the protein. 

Download the file in a `.pdb` format. 

### Clean the Protein
The protein file will usually include crytallography artifacts such as water, ions, solvent molecules, etc. There may also be multiple chains of the protein that are unnecessary. We will want to remove these extra molecules to obtain a 'clean' protein. The bound ligand (LacNAc) can also be removed at this step to make a protein-only system. 

You can **clean the protein** file in multiple ways:
1. Open the `.pdb` file in a software such as ChimeraX. Here you can use the selection tools to specify the unwanted molecules and remove them. This is a beginner friendly way to clean a file and is the most intuitive since you can visualize what is being removed.

2. Alternatively, open the `.pdb` file in a text editor and manually delete the unwanted atom lines from the file. This is a more hands-on approach but will allow you to familiarize yourself with the internal structure of these file types. The top of the PDB contain a header with information about the protein. You can remove this portion for the working file after you have read through the important sections. The `ATOM` lines signify the start of the actual protein structure. Normally the protein atoms will be listed first, followed by the water and other extraneous molecules. The file will sometimes end with a connectivity section, especially for the non-protein molecules. You will want to keep the proper 'ATOM' lines and the 'END' line for your cleaned protein. You can load the pdb file into VMD and save it again to reformat it nicely when you are done.


***Ligands make things more complicated:*** If you wish to include the ligand in the bound position from the crystal structure you can keep it in the pdb file however there will be some complications when converting to a tinker xyz (.txyz) file. Most solutions involve manually labeling the individual atoms and transferring the xyz coordinates to the correct atoms/atom types. If you have recently parameterized the molecule in Poltype2, you can take the poltype output txyz file, transfer the coordinates from the pdb to the txyz, using atom labels to match up the correct atom lines. You can then delete the ligand for the rest of the protein preparation and only add it after you have your protein in a txyz file. Paste the ligand at the end of the file and just renumber the index. Python scripts can help automate this. Note this wont work for every situation, for example if you are solvating with gromacs you may want the ligand in the pbd file when the water is added. 


### Add Hydrogens

Most pdb files will need polar hydrogens added. This can be done in software like ChimeraX, however our group uses the 'Reduce' program to try and predict which residues should be protonated. Histidine 






### Structure Modification
If any residues are missing you will need to add them. If the missing residues will not be involved in a structured domain, the 'protein' Tinker executable can be useful to add them later on when the protein is already in a tinker xyz file format. Structure editors such as Pymol and Avagodro are also options. Otherwise homology modelling software such as Swiss-Model or Alpha Fold can be used to predict the new protein structure.


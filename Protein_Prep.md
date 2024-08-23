 ## Protein File
Head over the [RSCB Protein Data Bank](https://www.rcsb.org) to find the protein of interest. This is where many protein and DNA structures have been deposited by experimentalists and where many of our simulation projects begin! 

Search for the protein by name (**galectin-3**) or better by its unique PDB ID (**1KJL**). 

It is usually helpful to read the information about the protein file on the RSCB website page and use the built in visualization tools to view the system. For our example we see that there are 138 residues and a LacNAc molecule is present in the binding site of the protein. 

Download the file in a `.pdb` format. 


### Clean the Protein
The protein file will usually include crytallography artifacts such as water, ions, solvent molecules, etc. We will want to remove these extra molecules to obtain a 'clean' file.
The bound ligand (LacNAc) can also be removed at this step if a protein-only system is desired. 

You can clean the protein file in multiple ways:
1. Open the `.pdb` file in a software such as ChimeraX. Here you can use the selection tools to specify the unwanted molecules and remove them. This is a beginner friendly way to clean a file and is the most intuitive since you can visualize what is being removed.
2. Alternatively, open the `.pdb` file in a text editor and manually delete the unwanted atom lines from the file. This is a more hands-on approach but will allow you to familiarize yourself with the internal structure of these file types. The top of the PDB will be a header that contains information about the protein, you can remove this portion for the working file after you have read through the important sections. The `ATOM` lines signify the start of the actual protein structure. Normally the protein atoms will be listed first, followed by the water and other extraneous molecules. The file will sometimes end with a connectivity section, especially for the non-protein molecules.

### Add Hydrogens

### Structure Modification
If any residues are missing you will need to add them. If the missing residues will not be involved in a structured domain, the 'protein' Tinker executable can be useful to add them later on when the protein is already in a tinker xyz file format. Structure editors such as Pymol and Avagodro are also options. Otherwise homology modelling software such as Swiss-Model or Alpha Fold can be used to predict the new protein structure.


#Workflow for protein simulations with Tinker

1) Get a PDB file of your protein complexed with all necessary small molecules (if any). This is usually downloaded from the PDB databank. 

2) The structure in PDB files often contain multiple chain with ligand geometries of varying quality and incomplete protein sequence. Clean up the file and add all missing residues by following the instructions in `AddMissingResidues.md`

3) Hydrogen atoms are also often missing, add them following the instructions in `AddHydrogensPDB.md`

4) [Optional depending on the project] Create a conformational ensemble by parsing through the backbone and side chain degrees of freedom with Rosetta (see instructions in `Rosetta.md`)

5) Add explicit water following the instructions in `SolvateProtein.md`

6) Convert the PDB(s) into Tinker readable coordinate files following the instructions in `PDBtoTinkerXYZ.md`

7) Run Tinker simulations! (`RunningTinkerBasics.md`)


NB: the documentation is written for someone using Cascades. This is not a requirement as other machines are available on ARC. However, if you use a different machine, make sure you make all necessary changes to the submission scripts. 

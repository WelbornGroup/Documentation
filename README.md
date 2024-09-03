# Welborn Lab Documentation

This is where we keep documentation for the tasks and software used in the group!


Below you will find a general tutorial and walkthrough on preparing a protein for Tinker MD simulations.

An older workflow for simulating proteins is [here](previous_workflow.md) and for simple carbohydrates [here](REU_Workflow.md).

Note, not all of the documentation files are linked through the workflows.


## Introduction
This workflow below will go through the process of preparing a protein system for a molecular dynamics (MD) simulation using the Tinker MD engine. This is a general overview and can be adapted to fit your needs and specific systems. The galectin-3 protein will be used as an example throughout this tutorial. There will always be multiple ways to accomplish some steps and these will be mentioned when possible! 


Workflow:

**1.** Obtain and prepare your protein structure file. [Protein Prep](./ProteinPrep.md)

**2.** [Optional] Parameterize small molecules with Poltype2. [Poltype Usage](./Poltype_Usage.md)

**3.** [Optional] Use Rosetta to create multiple starting conformations of your protein. [Rosetta](./Rosetta.md)

**4a.** [Convert PDB to Tinker XYZ](PDBtoTinkerXYZ.md) ***then*** [Solvate with Tinker](./SolvateTinker.md) 
(*Best for small molecules and files already in .txyz format*)



***OR***

**4b.** [Solvate with Gromacs](./SolvateProtein.md) ***then*** [Convert PDB to Tinker XYZ](PDBtoTinkerXYZ.md) 
(*Best for proteins with internal cavities that need hydration*)

 
     
**5.** Run the dynamics simulation with Tinker9. 

**6.** [Optional] Tinker analysis and helpful tools.





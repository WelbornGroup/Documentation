# Welborn Lab Documentation

This is where we keep documentation for the tasks and software used in the group!


Below you will find a general tutorial and walkthrough on preparing a protein for Tinker MD simulations.

An older workflow for simulating proteins is [here](previous_workflow.md) and for simple carbohydrates [here](REU_Workflow.md)

Note that not all of the documentation is linked through the workflows.


## Introduction
The workflow below will go through the process of preparing a protein system for a molecular dynamics (MD) simulation using the Tinker MD engine. This is a general overview and can be adapted to fit your needs and specific systems. The galectin-3 protein will be used as an example throughout this tutorial. There will always be multiple ways to accomplish some steps and these will be mentioned when possible! 

Workflow:

Obtain and prepare your protein structure file. [Protein_Prep](./Protein_Prep.md)

[Optional] Use Rosetta to create multiple starting conformations of your protein. [Rosetta](./Rosetta.md)

Solvate your protein system with Gromacs. (Best for proteins with internal cavities that need hydration) [SolvateProtein](./SolvateProtein.md)
OR
Solvate your system with Tinker. (Best for small molecules and files already in .txyz format) [SolvateTinker](./SolvateTinker.md)

Run the dynamics simulation with Tinker9. 

[Optional] Tinker analysis and helpful tools.





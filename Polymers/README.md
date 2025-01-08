# INTRODUCTION

This workflow below will go through the process of preparing a Polymer system for a molecular dynamics (MD) simulation using the Tinker MD engine. This is a general overview and can be adapted to fit your needs and specific systems. A beta-glucose dimer will be used as an example throughout this tutorial. There will always be multiple ways to accomplish some steps and these will be mentioned when possible!

## Workflow:

All steps are optional, and will depend of the system you're working on and the files already available to you.

1. Create the polymer input structure. [create-polymer-with-pymol.md](./create-polymer-with-pymol.md)
2. Parameterize the monomers with Poltype2 to be used for Tinker MD simulation with AMOEBA Force Field. [parameterization-for-amoeba.md](./parameterization-for-amoeba.md)
3. Solvate the polymer and add ions. [solvate-polymer.md](./solvate-polymer.md)
4. Basics of running Tinker MD. [running-tinker-basics.md](./running-tinker-basics.md)
5. Trajectory analysis and helpful tools. (Under Development)

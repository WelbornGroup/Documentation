#Add hydrogens to PDB files

Download and install the [Reduce program](http://kinemage.biochem.duke.edu/software/reduce.php).

In the directory where you have your pdb, type the command

```sh
reduce -build input.pdb > output.pdb
```

The `build` argument is to make sure the Histdines are correctly protonated. 

Open `output.pdb`. You should see new lines corresponding to the hydrogens atoms that were added. 
Note however that these atoms are not numbered sequentially. It is good practice to open the file in VMD (or similar software) to CHECK the structure and save it again as a pdb. This will fix the numbering issue. 

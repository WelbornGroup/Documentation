# This is a Quick Guide to using Poltype2
### (Please see the Poltype2 install guide first to get it set up on ARC)


Poltype2 is an automated software that generates AMOEBA MD parameters for small molecules

To utilize this software you will require 4 input files:

1. A structure file for your molecule (`.mol`)
2. An initiation file which contains the settings or commands for Poltype (`poltype.ini`)
3. A file that contains the paths to certain software that poltype needs (`paths.sh`)
4. The bash script to submit the poltype job to the queue on ARC (`run-poltype.sh`)


**Structure File** 
Generate a structure file for your molecule using Avogadro or glycam's carbohydrate builder

Poltype requires a `.mol` or `.sdf` structure file. If using glycam download the minimized carbohydrate `.pdb` file, open it in Avogadro, and save it as a new `.mol` file

**Ini File**
Here is a basic example for the ini file:
```
structure=galactose.mol
dontfrag=True
atmidx=400
```

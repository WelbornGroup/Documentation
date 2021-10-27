### VASP simulations

Vienna Ab initio Simulation Package (VASP) is a software for ab initiao quantum mechanical calculations using pseudopotentials and plane wave basis set.

For a single point calculation, VASP requires four input files to run, `INCAR`, `KPOINTS`, `POSCAR`, and `POTCAR`.


The `INCAR` specify the calculation parameters for the calculation.

```sh
SYSTEM = density generation

#### SCF Routine
ENCUT = 750.00 eV # basis set energy cutoff
PREC = Accurate
LREAL = .FALSE.
NELMIN = 6; # minimum number of SCF cycles
NELM = 120
NELMDL = -4
EDIFF = 0.00001

ALGO = Fast; # uses Davidson initially then switches to DIIS
TIME = 0.4 # 0.4 is the default, controls the size of the time step for these algorithms

NSW = 0
ISPIN = 2 # spin polarized system
LCHARG = .TRUE.
LAECHG = .TRUE.
LWAVE=.FALSE.

# U parameters setting
LDAU = .TRUE.
LDAUL = -1 2 2 -1 -1
LDAUU = 0.00 4.00 3.50 0.00 0.00
LDAUJ = 0.00 0.00 0.00 0.00 0.00
LMAXMIX = 4

### Symmetry determination
ISYM = 0;

### Parameters related to fractional occupations of the orbitals
ISMEAR = -1; # Fermi semearing
SIGMA = 0.05; # in eV for smearing width

###Parallelization options
LPLANE = .TRUE.
NCORE = 16
LSCALU = .FALSE.
NSIM = 1
```

For charged system, you need to set the number of electrons in the unit cell with `NELECT` tag.

For geometry optimization, you need to set the optimization parameter with `ISIF` and `IBRION` tages. See VASPwiki for different options.

`POSCAR` specifies the unit cell information, including the unit cell vectors, elements in the unit cell and atom positions. 

```sh
1.0
6.970873 0.000000 4.024635
2.323624 6.572202 4.024635
0.000000 0.000000 12.073905
Li Mn Cr O F
8 3 1 8 4
direct
0.000000 0.000000 0.000000 
0.000000 0.000000 0.333333 
0.000000 0.000000 0.666667 
1.000000 0.500000 0.000000 
0.500000 0.000000 0.333333 
0.500000 0.500000 1.000000 
0.500000 0.500000 0.333333 
0.500000 0.500000 0.666667 
1.000000 0.500000 0.333333 
0.500000 0.000000 1.000000 
0.500000 0.000000 0.666667 
1.000000 0.500000 0.666667 
0.250000 0.250000 0.166667 
0.250000 0.250000 0.500000 
0.250000 0.250000 0.833333 
0.250000 0.750000 0.166667 
0.750000 0.250000 0.500000 
0.750000 0.250000 0.833333 
0.750000 0.750000 0.500000 
0.750000 0.750000 0.833333 
0.250000 0.750000 0.500000 
0.250000 0.750000 0.833333 
0.750000 0.250000 0.166667 
0.750000 0.750000 0.166667 
```

The first line is the scaling factor of the unit cell, followed by the three unit cell vectors. Then the fifth line specify the elements in the system, followed by the number of atoms for each element in the unit cell. The order of this two line must match, e.g., if you switch position of Li and Mn in this example, you need to change the sixth line to `3 8 1 8 4` as well.

The seventh line specifies how the atom positions are defined. In this example, `direct` means using fractional coordinates. You can also use `cartesian`, which defines the atom position in cartesian coordinates in the unit of angstrom.

The rest of the file are the coordinates of each atom. The order of the atoms must match the order in line five, meaning in this particular example, you must specify the coordinates for all the Li atoms before Mn atoms, and so on and so forth. 

`POTCAR` contains the pseudopotentials for each elements for the calculations. The order of elements in this file must match the order of elements in line 5 of the `POSCAR` as well. 

You need to copy the content in `VASP-redommended-pseudopotentials` on GitHub to your work directory. In the directory for the calculation, type: `cat path/to/pseudopotential/element1/POTCAR path/to/pseudopotential/element2/POTCAR ... >POTCAR`. This will create a `POTCAR` file in the current directory with contents from each element's POTCARs.

`KPOINTS` contains the specification of the kpoints settings of the calculation. 

```sh
Automatic mesh
0
Monkhorst-Pack
3 3 2
0 0 0
```

Line 4 is the kpoints on each vector and line 5 is the offset. Note that if you have `1 1 1` for kpoints, you need to change line 3 to `Gamma`.
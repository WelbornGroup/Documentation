#Types of calculations

Two major types of calculations can be carried out: electronic structure calculation and geometry optimization. 

Electronic structure calculation yields the ground state energy of the system under the level of theory you used. It also yields the electron and spin density of the system.

Geometry optimization calculation starts from a specified inital system geometry and finds the system geometry that has the lowest energy under the level of theory you used. Note that you should start geometry optimization from a geometry that closes to the optimized geometry. Usually experimental geometry data reported in the literature is a good starting point (experimental geometry is not necessarily the optimized geometry, depending on the level of theory used, there may be slight differences).


#Quantum Espresso input file

First We will use the previous benzene calculation as an example to explain the meaning of the different section of the input file.

The input file consists of several namelists, followed by other sections introduced by keywords. The required namelists are:

    &CONTROL: general variables controlling the calculation
    &SYSTEM: strictiral information on the system
    &ELECTRONS: electronic variables

In our example calculations, namelist &CONTROL defines the calculation name and the directory containing pseudopotentials. 

In namelist &SYSTEM:
    `assume_isolated` instructs the caluclation being performed assuming the system to be isolated, `martyna-tuckerman` (or `m-t`, `mt) means to calculate both total energy and scf potential of the system. Details of this method can be found in http://doi.org/10.1063/1.477923. 
    The default value of this tag is `None`, which treats the system as periodic. This should be used (or just remove this tag) when running calculation on solids. 

    `ibrav` specifies the the Bravais-lattice index of the unit cell. You should also specify the lattice parameters for the lattice you chose. In our example, `ibrav = 6` defines a tetragonal unit cell. For this particular choice, we need to specify the length of vector a and c in angstrom.

    `ecutwfc` specifies the kinetic energy cutoff in Ry for wavefunctions.
    `ntyp` specifies the kinetic energy cutoff in Ry for charge density and potential.
    In general, higher the cutoff, higher the accuracy of the calculation, the computational cost will increase as well.

    `nat` specifies the number of atoms in the unit cell.
    `ntyp` specifies the number of types of atoms in the unit cell. 

    `nbnd` specifies the number of electronic states (bands) to be calculated.

ATOMIC_SPECIES section instructs the program to read corresponding pseudopotential files for each element in the system. 
The second field (i.e., 1.0) is the atomic mass. This field is not used except for MD calculations.

ATOMIC_POSITIONS section specifies the coordinates of each atom in angstrom.

`KPOINTS` are sampling points in reciprocal space for plane-wave calculation. This should be specified for calculation of periodic materials. In some cases, the optimization of Kpoints should be performed. 

Let`s now take a look at the input file for magnatic solid system using our Ni calculation. 

In namelist &SYSTEM, we have two new tags: `nspin` and `starting_magnetization(i)`. The former specifies the magnetization of the system, `1` for non-polarized (default), `2` for collinear system, and `4`for noncollinear system. The latter specifies the initial guess of the magnetization magnitude for each atomic species (each one in ATOMIC_SPECIES section). Also, we removed the tag `assume_isolated` because we want to run the calculation as a periodic system. Another difference is, instead of `nbnd`, we used `occupations` instead. By setting the value to `smearing`, it allows fractional occupatition of the eigenstates, avoiding the electron occupation changes coused by every small pertubations when there are several degenerate eigenstates.

In K_POINTS section, we use `automatic` instead of `gamma`. This means generate (nk1, nk2, nk3) grid (i.e., the first group of three numbers) with (sk1, sk2, sk3) offset (i.e., the second group of three numbers). For periodic system, you may need to optimize the grid, i.e., changing the first group of numbers to find the group that gives the lowest energy. A good starting point of Kpoints is nk1 * x_axis >= 32.

Note that you need to optimize the Kpoints (preferably using the experimental geometry) first before optimizing the lattice parameters.


#Quantum Espresso output file

The submission script should generate an output file of the program `output.out`. This file contains the stepwise details of the calculation.

1) Go to your calculation directory: `cd path/to/calculation/directory`.

2) Open the file by typing `cat output.out`. This command will display the contant of the file in your terminal.

3) In the beginning of the file there will be a section specifying the settings of the calculation. Double check to make sure these are what you specified.

```sh
     bravais-lattice index     =            6
     lattice parameter (alat)  =      20.7870  a.u.
     unit-cell volume          =    5715.8394 (a.u.)^3
     number of atoms/cell      =           12
     number of atomic types    =            2
     number of electrons       =        30.00
     number of Kohn-Sham states=           16
     kinetic-energy cutoff     =      20.0000  Ry
     charge density cutoff     =     200.0000  Ry
     convergence threshold     =      1.0E-06
     ⋮
     ⋮
```

4) Look for the line that says `Self-consistent Calculation`. Below that line, are the self-consistent field cycles. Each cycle will have a total energy, and a total magnetization for polarized system. These energies (and magnetizations) are only intermediate results, do not mistake them as the final energy/magnetization of the system.

5) Look for the line that says `End of self-consistent calculation`. The total energy and total magnetization below that line are the final results of the system.
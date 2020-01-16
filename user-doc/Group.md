\page Group Groups and Virtual Atoms 

\section atomSpecs Specifying Atoms

The vast majority of the CVs implemented in PLUMED are calculated from a list of atom positions.  Within PLUMED
atoms are specified using their numerical indices in the molecular dynamics input file. 

In PLUMED lists of atoms can be either provided directly inside the definition of each collective variable, or
predefined as a \subpage GROUP that can be reused multiple times. Lists of atoms can be written as:

- comma separated lists of numbers (`GROUP ATOMS=10,11,15,20 LABEL=g1`)
- numerical ranges.  So `GROUP ATOMS=10-20 LABEL=g2` is equivalent to `GROUP ATOMS=10,11,12,13,14,15,16,17,18,19,20 LABEL=g2`
- numerical ranges with a stride. So `GROUP ATOMS=10-100:10 LABEL=g3 is equivalent to `GROUP ATOMS=10,20,30,40,50,60,70,80,90,100 LABEL=g3`
- atoms ranges with a negative stride. So `GROUP ATOMS=100-10:-10 LABEL=g4 is equivalent to `GROUP ATOMS=100,90,80,70,60,50,40,30,20,10 LABEL=g4`

In addition, there are a few shortcuts that can be used:

- `@mdatoms` indicate all the physical atoms present in the MD engine (e.g. `DUMPATOMS ATOMS=@mdatoms`).
- `@allatoms` indicates all atoms, including \ref vatoms "those defined only in PLUMED" (e.g. `DUMPATOMS ATOMS=@allatoms`).

The list of the virtual atoms defined in PLUMED can be obtained by using the command `GROUP ATOMS=@allatoms REMOVE=@mdatoms`.

Other shortcuts are available if you loaded the structure of the molecule using the \ref MOLINFO command.

All the above methods can be combined just putting one name after the other separated by a comma:
\plumedfile
DUMPATOMS ATOMS=1,2,10-20,40-60:5,100-70:-2 LABEL=g5 FILE=test.xyz
\endplumedfile

Some collective variable must accept a fixed number of atoms, for example a \ref DISTANCE is calculated
using two atoms only, an \ref ANGLE is calculated using either 3 or 4 atoms and \ref TORSION is calculated using 4 atoms.

Additional material and examples can be also found in the tutorial \ref belfast-1. 

\subsection mols Molecules

In addition, for certain colvars, pdb files can be read in using the following keywords and used to select ATOMS:

@TOPOLOGY@

\subsection pbc Broken Molecules and PBC 

PLUMED is designed so that for the majority of the CVs implemented the periodic boundary conditions are treated 
in the same manner as they would be treated in the host code.  In some codes this can be problematic when the colvars
you are using involve some property of a molecule.  These codes allow the atoms in the molecules to become separated by 
periodic boundaries, a fact which PLUMED could only deal with were the topology passed from the MD code to PLUMED.  Making this
work would involve a lot laborious programming and goes against our original aim of having a general patch that can be implemented 
in a wide variety of MD codes.  Consequentially, we have implemented a more pragmatic solution to this problem - the user specifies
in input any molecules (or parts of molecules) that must be kept in tact throughout the simulation run.  In PLUMED 1 this was done
using the ALIGN_ATOMS keyword.  In PLUMED 2 the same effect can be achieved using the \subpage WHOLEMOLECULES command.

The following input computes the end-to-end distance for a polymer of 100 atoms and keeps it at a value around 5.

\plumedfile
WHOLEMOLECULES ENTITY0=1-100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
\endplumedfile

Notice that NOPBC is used to be sure in \ref DISTANCE that if the end-to-end distance is larger than half the simulation box the distance 
is compute properly. Also notice that, since many MD codes break molecules across cell boundary, it might be necessary to use the 
\ref WHOLEMOLECULES keyword (also notice that it should be before distance).

Notice that most expressions are invariant with respect to a change in the order of the atoms,
but some of them depend on that order. E.g., with \ref WHOLEMOLECULES it could be useful to
specify atom lists in a reversed order.

\plumedfile
# to see the effect, one could dump the atoms as they were before molecule reconstruction:
# DUMPATOMS FILE=dump-broken.xyz ATOMS=1-20
WHOLEMOLECULES STRIDE=1 ENTITY0=1-20
DUMPATOMS FILE=dump.xyz ATOMS=1-20
\endplumedfile

Notice that there are other ways to manipulate the coordinates stored within PLUMED:
- Using the \subpage FIT_TO_TEMPLATE they can be aligned to a template structure.
- Using \subpage WRAPAROUND you can bring a set of atom as close as possible to another set of
  atoms.
- Using \subpage RESET_CELL you can rotate the periodic cell.

\section vatoms Virtual Atoms

Sometimes, when calculating a colvar, you may not want to use the positions of a number of atoms directly.  Instead
 you may wish to use the position of a virtual atom whose position is generated based on the positions of a collection 
of other atoms.  For example you might want to use the center of mass of a group of atoms.  Plumed has a number of routines
for calculating the positions of these virtual atoms from lists of atoms:

@VATOM@

To specify to a colvar that you want to use the position of a virtual atom to calculate a colvar rather than one of the atoms
in your system you simply use the label for your virtual atom in place of the usual numerical index.  Virtual
atoms and normal atoms can be mixed together in the input to colvars as shown below:

\plumedfile
COM ATOMS=1,10 LABEL=com1
DISTANCE ATOMS=11,com1
\endplumedfile

If you don't want to calculate CVs from the virtual atom.  That is to say you just want to monitor the position of a virtual atom 
(or any set of atoms) over the course of your trajectory you can do this using \ref DUMPATOMS.


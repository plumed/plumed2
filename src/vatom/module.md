Sometimes, when calculating a colvar, you may not want to use the positions of a number of atoms directly.  Instead
 you may wish to use the position of a virtual atom whose position is generated based on the positions of a collection
of other atoms.  For example you might want to use the center of mass of a group of atoms.  This module contains 
actions that calculatie the positions of these virtual atoms from lists of input atoms.

To specify to a colvar that you want to use the position of a virtual atom to calculate a colvar rather than one of the atoms
in your system you simply use the label for your virtual atom in place of the usual numerical index.  Furthermore, virtual
atoms and normal atoms can be mixed together in the input to colvars as shown below:

```plumed
com1: COM ATOMS=1,10
d: DISTANCE ATOMS=11,com1
```

If you don't want to calculate CVs from the virtual atom.  That is to say you just want to monitor the position of a virtual atom
(or any set of atoms) over the course of your trajectory you can do this using [DUMPATOMS](DUMPATOMS.md).

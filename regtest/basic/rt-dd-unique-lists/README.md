# Sharing of positions and forces through DOMAIN_DECOMPOSITION in serial

`main.cpp` runs the same calculation with three layouts of the arrays passed by the
MD code: all the atoms in order, all the atoms in reversed order (`setAtomsGatindex`),
and only the atoms needed by PLUMED, as NAMD does (`createFullList`/`getFullList`).
The input mixes actions that are active only at some steps, so that the atoms needed
change from step to step and the union of the atoms needed is sometimes larger and
sometimes smaller than half of the atoms, a neighbor list, so that an action requests
different atoms during a step, and actions with and without forces.
The forces on the atoms and the bias of the first layout are printed, and the other
layouts must agree with it.

The same `main.cpp` is used by `rt-dd-unique-lists-yes` and `rt-dd-unique-lists-no`
with `PLUMED_FORCE_UNIQUE=yes` and `no`; the three tests share the same reference,
since the choice of the algorithm must not change the results.
The reference is identical to the output of the code before the lists of atoms were
cached (with `./exe skip`, which skips the error checks below).

Finally, the test checks that PLUMED reports an error, instead of crashing or reading
out of bounds, when forces have to be applied but the MD code did not pass the array
for the forces (on atoms and on a value created with PUT), and when the positions of
shuffled atoms are passed in an array that is too short.

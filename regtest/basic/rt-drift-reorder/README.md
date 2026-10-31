# EFFECTIVE_ENERGY_DRIFT with atoms that change order

`main.cpp` passes the same trajectory twice, in serial: once with the atoms in order,
once with a different order at every step (`setAtomsGatindex`), as LAMMPS does when it sorts
its atoms. The effective energy drift must be the same. Before the fix, the stored positions
and forces of the previous step were only remapped with more than one process, so the second
run gave a wrong drift.

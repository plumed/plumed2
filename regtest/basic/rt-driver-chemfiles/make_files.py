# Writes the trajectories of this test with the chemfiles python package (pip install chemfiles==0.10.4):
# the same 5 frames of 12 atoms in a triclinic cell, in the formats that the native readers also read
# (xtc, trr, dcd, pdb) and in formats read only with chemfiles (Amber NetCDF, LAMMPS dump, compressed xyz and pdb).
# The coordinates are multiples of 0.25 A, exact in every format (the xtc precision is 0.001 nm).
# Every frame stores its step (100*i, in the formats that can) and its time (0.2*i ps).
import chemfiles

natoms = 12
nframes = 5


def frame(i):
    f = chemfiles.Frame()
    f.cell = chemfiles.UnitCell([32.0, 30.0, 28.0], [80.0, 95.0, 105.0])
    for j in range(natoms):
        name = ["O", "H", "H"][j % 3]
        x = 0.25 * ((7 * j + 3 * i) % 100)
        y = 0.25 * ((11 * j + 5 * i) % 90)
        z = 0.25 * ((13 * j + 2 * i + 4) % 80)
        f.add_atom(chemfiles.Atom(name), [x, y, z])
    f.step = 100 * i
    f["time"] = 0.2 * i
    return f


for name, fmt in [
    ("traj.xtc", ""),
    ("traj.trr", ""),
    ("traj.dcd", ""),
    ("traj.pdb", ""),
    ("traj.nc", ""),
    ("traj.lammpstrj", ""),
    ("traj.xyz.gz", ""),
    ("traj.pdb.gz", ""),
    ("traj.dump", "LAMMPS"),
]:
    with chemfiles.Trajectory(name, "w", fmt) as t:
        for i in range(nframes):
            t.write(frame(i))

# two frames with a different number of atoms
with chemfiles.Trajectory("natoms.xyz", "w") as t:
    t.write(frame(0))
    f = frame(1)
    f.remove(0)
    t.write(f)

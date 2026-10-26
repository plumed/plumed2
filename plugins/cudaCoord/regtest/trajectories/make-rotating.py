#! /usr/bin/env python3
# Generates traj-rotating.xyz and template-rotating.pdb (deterministic).
#
# A rigid 20-atoms "host" and 200 "water" atoms in a cubic box. In every frame
# the whole system is rigidly rotated by a different, large, rotation about the
# centre of the box and wrapped back into the box; a small noise is added to
# the host so that the fit is not trivially exact.
#
# With FIT_TO_TEMPLATE TYPE=OPTIMAL on the host, plumed brings each frame back
# to the frame of the template AND rotates the box accordingly: the actions
# downstream receive a box that is orthorhombic but NOT diagonal. This is what
# a funnel-like setup (FIT_TO_TEMPLATE + coordination with the solvent) does.
# The host sits off-centre so that part of its solvation shell is across the
# periodic boundary.
import numpy as np

rng = np.random.default_rng(20261026)
L = 6.0
nwat = 200
nframes = 6

# host: a helix, radius 0.5, 1.5 long, centred in (1,1,1)
nhost = 20
t = np.linspace(0.0, 4.0 * np.pi, nhost)
host0 = np.stack(
    [0.5 * np.cos(t), 0.5 * np.sin(t), np.linspace(-0.75, 0.75, nhost)], axis=1
)
host0 += np.array([1.0, 1.0, 1.0])
water0 = rng.uniform(0.0, L, size=(nwat, 3))
centre = np.array([L / 2, L / 2, L / 2])


def random_rotation():
    # uniform random rotation (Shoemake)
    u1, u2, u3 = rng.uniform(size=3)
    q = np.array(
        [
            np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
            np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
            np.sqrt(u1) * np.sin(2 * np.pi * u3),
            np.sqrt(u1) * np.cos(2 * np.pi * u3),
        ]
    )
    x, y, z, w = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


with open("traj-rotating.xyz", "w") as f:
    for frame in range(nframes):
        R = np.eye(3) if frame == 0 else random_rotation()
        host = host0 + rng.normal(scale=0.02, size=host0.shape)
        pos = np.concatenate([host, water0 + rng.normal(scale=0.05, size=water0.shape)])
        pos = (pos - centre) @ R.T + centre
        pos = np.mod(pos, L)
        f.write(f"{nhost + nwat}\n")
        f.write(f"{L} 0 0 0 {L} 0 0 0 {L}\n")
        for i, p in enumerate(pos):
            name = "C" if i < nhost else "O"
            f.write(f"{name} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n")

with open("template-rotating.pdb", "w") as f:
    for i, p in enumerate(host0):
        # occupancy and beta are the weights for the alignment and for the displacement
        f.write(
            "ATOM  %5d  C   HST     1    %8.3f%8.3f%8.3f  1.00  1.00\n"
            % (i + 1, 10 * p[0], 10 * p[1], 10 * p[2])
        )
    f.write("END\n")

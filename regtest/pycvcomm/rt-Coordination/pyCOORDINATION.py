# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
from sys import stderr as log

# import plumedUtilities
print("Imported pyCoord.", file=log)

N: int = 6
M: int = 12
D0: float = 0.0
R0: float = 2.0
INVR0: float = 1.0 / R0
DMAX = D0 + R0 * (0.00001 ** (1.0 / (N - M)))
STRETCH = 1.0
SHIFT = 0.0


def switch(d: np.ndarray) -> np.ndarray:
    ret = np.zeros_like(d)
    dfunc = np.zeros_like(d)
    WhereToCalc =  d < DMAX# & d > D0
    #print(f"{dfunc=}", file=log)
    #not doinf d-D0 for now, so no need for rdist<=0
    rdist = d[WhereToCalc] * INVR0
    rNdist = (rdist) ** (N - 1)
    ret[WhereToCalc] = (1.0 / (1 + rdist * rNdist)) * STRETCH
    ret[WhereToCalc] += SHIFT
    dfunc[WhereToCalc] = -N * rNdist * ret[WhereToCalc] * ret[WhereToCalc] / d[WhereToCalc]
    dfunc *= STRETCH * INVR0
    return ret, dfunc


def stretchSwitch(mydmax):
    s, _ = switch(np.array([0.0, mydmax]))
    stretch = 1 / (s[0] - s[1])
    return stretch, -s[1] * stretch


# some juggling for calculationg the stretch
mydmax = DMAX
DMAX += 0.1
STRETCH, SHIFT = stretchSwitch(mydmax)
DMAX = mydmax
print(f"{N=} {M=} {D0=} {R0=} {INVR0=} {DMAX=} {STRETCH=} {SHIFT=}", file=log)
      
def pyCoord(action: plumedCommunications.PythonCVInterface):
    
    atoms = action.getPositions()
    nat = atoms.shape[0]
    nl = action.getNeighbourList()

    assert nl.size() == ((nat - 1) * nat) // 2
    pbc = action.getPbc()
    couples = nl.getClosePairs()
    #from here we are in "pure python"
    d = atoms[couples[:, 0]] - atoms[couples[:, 1]]
    #print(f"before {d}",file=log)
    d = pbc.apply(d)
    #print(f"after {d}",file=log)
    dist = np.linalg.norm(d, axis=1)
    sw, dfunc = switch(dist)
    dev = np.zeros_like(atoms)

    predev = d * dfunc.reshape((-1,1))
    for atomID in range(nat):
        # wherePlus =couples[:, 0]==atomID
        # whereMinus=couples[:, 1]==atomID
        dev[atomID] = np.sum(predev[couples[:, 0] == atomID], axis=0) - np.sum(
            predev[couples[:, 1] == atomID], axis=0
        )
    virial=np.zeros((3,3))
    for i in range(predev.shape[0]):
        virial-=np.outer(predev[i],d[i])

    return np.sum(sw), dev, virial

# This is a not optimzed implementation of the COORDINATION

import numpy as np
import plumedCommunications
from sys import stderr as log
from jax import jit

# import plumedUtilities
print("Imported pyCoord.", file=log)

N: int = 6
M: int = 12
D0: float = 0.0
R0: float = 0.4
INVR0: float = 1.0 / R0
DMAX = D0 + R0 * (0.00001 ** (1.0 / (N - M)))
STRETCH = 1.0
SHIFT = 0.0


@jit
def jaxSwitch(d):
    rdist = d * INVR0
    rNdist = (rdist) ** (N - 1)
    ret = (1.0 / (1 + rdist * rNdist)) * STRETCH
    ret += SHIFT
    dfunc = -N * rNdist * ret * ret / d
    dfunc *= STRETCH * INVR0
    return ret, dfunc


def switch(d: np.ndarray) -> np.ndarray:
    ret = np.zeros_like(d)
    dfunc = np.zeros_like(d)
    WhereToCalc = d < DMAX  # & d > D0
    # print(f"{dfunc=}", file=log)
    # not doinf d-D0 for now, so no need for rdist<=0
    ret[WhereToCalc], dfunc[WhereToCalc] = jaxSwitch(d[WhereToCalc])
    return ret, dfunc


def stretchSwitch():
    s, _ = jaxSwitch(np.array([0.0, DMAX]))
    stretch = 1 / (s[0] - s[1])
    return stretch, -s[1] * stretch


# some juggling for calculationg the stretch
STRETCH, SHIFT = stretchSwitch()
print(f"{N=} {M=} {D0=} {R0=} {INVR0=} {DMAX=} {STRETCH=} {SHIFT=}", file=log)


def pyCoord(action: plumedCommunications.PythonCVInterface):
    atoms = action.getPositions()
    nat = atoms.shape[0]
    nl = action.getNeighbourList()

    assert nl.size() == ((nat - 1) * nat) // 2
    pbc = action.getPbc()
    couples = nl.getClosePairs()
    absoluteIndexes=[]
    #not so fast, but speed here is not important
    for i in action.absoluteIndexes:
        absoluteIndexes.append(i.index)
    absoluteIndexes=np.array(absoluteIndexes)
    #sameIndex = np.where(absoluteIndexes[couples[:, 0]]==absoluteIndexes[couples[:, 1]])
    sameIndex = absoluteIndexes[couples[:, 0]]==absoluteIndexes[couples[:, 1]]
    d = atoms[couples[:, 0]] - atoms[couples[:, 1]]
    d = pbc.apply(d)
    # from here we are in "pure python"
    dist = np.linalg.norm(d, axis=1)
    
    dist[sameIndex] += DMAX*2.0
    
    sw, dfunc = switch(dist)
    dev = np.zeros_like(atoms)

    predev = d * dfunc.reshape((-1, 1))
    for atomID in range(nat):
        # wherePlus =couples[:, 0]==atomID
        # whereMinus=couples[:, 1]==atomID
        dev[atomID] = np.sum(predev[couples[:, 0] == atomID], axis=0) - np.sum(
            predev[couples[:, 1] == atomID], axis=0
        )
    virial = np.zeros((3, 3))
    for i in range(predev.shape[0]):
        virial -= np.outer(predev[i], d[i])

    return np.sum(sw), dev, virial

#this tests also the concatenating special keywords
plumedInit = {
    "Value": plumedCommunications.defaults.COMPONENT,
    "GROUPA": "@mdatoms,@mdatoms",
}

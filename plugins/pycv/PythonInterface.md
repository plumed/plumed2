# The Python interface

## Getting the manual
To obtain a <i>**VERY BASIC**</i> reference/manual with all the function and method definitions, run the test rt-doc or call `PYFUNCTION` in with: 

**plumed.dat**
```
LOAD GLOBAL FILE=path/to/PythonCVInterface.so
PYFUNCTION IMPORT=pycv
```

**pycv.py**
```python
import plumedCommunications
import pydoc

def plumedInit(_):
    pydoc.writedoc(plumedCommunications)
    pydoc.writedoc(plumedCommunications.defaults)
    return {"Value":plumedCommunications.defaults.COMPONENT_NODEV}

def plumedCalculate(_):
    return 0.0
```

## Specifics
In the following section
### Common interface
Both `PYFUNCTION` and `PYCVINTERFACE` have the following attribute:
 - `label` (readonly) returns the label

And the following functions:
 - `log(s:object)` puts a string in the PLUMED output
 - `lognl(s:object)` puts a string in the PLUMED output (and appends a newline)
 - `getStep()` Returns the current step
 - `getTime()` Return the present time
 - `getTimeStep()` Return the timestep
 - `isExchangeStep()` Check if we are on an exchange step
 - `isRestart()` Return true if we are doing a restart

### PYFUNCTION: arguments

`PYFUNCTION` accepts arguments in the plumed file with the `ARG` keyword, the arguments are then accessible in python with the functions:
 - `PythonFunction.argument(argID:int)` Get value of the of argID-th argument (as `float`)
 - `PythonFunction.arguments()` Retuns a ndarray with the values of the arguments (as `numpy.ndarray[numpy.float64]`)

and the argument:
- `PythonFunction.nargs` (readonly) Get the number of arguments

`PYFUNCTION` also has a few functions that can be used to interact with the periodicity of the arguments:
 - `PythonFunction.bringBackInPbc(argID:int,x:float)` Takes one value and brings it back into the pbc of argument argID
 - `PythonFunction.difference(argID:int,x:float,y:float)` Takes the difference taking into account pbc for argument argID

### PYCVINTERFACE: atoms, pbc, neigbourlist, makeWhole
`PYCVINTERFACE` works with atoms:

Has the following attributes:
 - `PythonCVInterface.nat` (readonly) Return the number of atoms

For getting the atomic data you can call the following functions:
 - `PythonCVInterface.getPosition(atomID:int)` Returns an ndarray with the position of the atomID-th atom 
 - `PythonCVInterface.getPositions()` Returns a numpy.array that contains the atomic positions of the atoms 
 - `PythonCVInterface.mass(atomID:int)` Get mass of atomID-th atom
 - `PythonCVInterface.masses()` Returns and ndarray with the masses 
 - `PythonCVInterface.charge(atomID:int)` Get charge of atomID-th atom
 - `PythonCVInterface.charges()` Returns and ndarray with the charges 
 - `PythonCVInterface.absoluteIndexes()` Get the npArray of the absolute indexes (like in AtomNumber::index()).

 And in general you can use some support-function from plumed:
 - `PythonCVInterface.getNeighbourList()` returns an interface to the current Neighborlist
 - `PythonCVInterface.getPbc()` returns an interface to the current pbcs
 - `PythonCVInterface.makeWhole()` Make atoms whole, assuming they are in the proper order

`plumedCommunications.Pbc` and `plumedCommunications.NeighborList` are simple methods container to help with the calculations:

`plumedCommunications.Pbc` has the following functions:
 - `Pbc.apply(d:ndarray)` pply PBC to a set of positions or distance vectors (`d.shape==(n,3)`)
 - `Pbc.getBox()` Get a numpy array of shape (3,3) with the box vectors
 - `Pbc.getInvBox()` Get a numpy array of shape (3,3) with the inverted box vectors

`plumedCommunications.NeighborList` has the following functions and attributes:
- `NeigbourList.getClosePairs()` get the (nl.size,2) nd array with the list of couple indexes
- `NeigbourList.size` and `len(nl:NeigbourList.size)` give the number of couples
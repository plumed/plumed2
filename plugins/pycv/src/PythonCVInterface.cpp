/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019-2023 of Toni Giorgino, Daniele Rapetti

The pycv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The pycv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "PythonCVInterface.h"

#include "plumed/core/ActionRegister.h"
#include "plumed/core/PlumedMain.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

//+PLUMEDOC COLVAR PYCVINTERFACE
/*
Define collective variables in the Python language.

Plumed will import the module chosen wiht the `IMPORT` keyword.
The module can be a simple py file or a directory with the standard module
configuration (a directory that contains at least an __init__.py, see the
rt-persistentData test for an example).

In the module there must be a function that will be used in caclulation and the
definition of the component(s) (meaning period and if derivatives will be
returned) of your cv.

PYCVINTERFACE will call two or more elements from the module.
Each element or function can be selected with its dedicated keyword:
 - `CALCULATE` will select the function to be called at each step (defaults to
    `"plumedCalculate"`)
 - `INIT` will select the function (or the dict, see down) to be called during
    the initilization (defaults to `"plumedInit"`)
 - `UPDATE` will select the function to be called at each step, during the
    update() invocation
 - `PREPARE` will select the function to be called at each step, during the
    prepare() invocation
All the function called will need a single argument of type
  `plumedCommunications.PythonCVInterface`


\par Getting started

The minimal plumed input is:
\plumedfile
cv1: PYTHONCV IMPORT=hello ATOMS=1
PRINT FILE=colvar.out ARG=*
\endplumedfile

this should be paired with the `hello.py` file:
@code{.py}
import plumedCommunications as PLMD

plumedInit={"Value":PLMD.defaults.COMPONENT_NODEV}

def plumedCompute(action: PLMD.PythonCVInterface):
    action.log("Hello, world, from calculate!")
    return 0.0
@endcode

If `INIT` is not specified, plumed will search for an object "plumedInit",
that can be either a function that returns a dict or a dict
This dict MUST contain at least the informations about the presence of the
derivatives and on the periodicity of the variable.
We will refer to this dict as "the init dict" from now.

If `CALCULATE` is not specified, plumed will search for a function named
"plumedCalculate" plumed will read the variable returned accordingly to what it
was specified in the initialization dict.

The init dict will tell plumed how many components the calculate function will
return and how they shall behave.
Along this the dict can contain all the keyword that are compatible with
PYCVINTERFACE.
Mind that if the same keyword is specified both in the init dict and in the
plumed file the calculation will be aborted to avoid unwated settings confict.
In case of flags the dict entry must be a bool, differently from the standard
plumed input.

The only keyword that can only be specified in python is `COMPONENTS`.
The `COMPONENTS` key must point to a dict that has as keys the names of the
componensts.
Each component dictionary must have two keys:
 - `"period"`: `None` of a list of two values, min and max (like `[0,1]` or also
 strings like `["0.5*pi","2*pi"]`)
 - `"derivative"`: `True` or `False`
If you want to use a single component you can create the `"COMPONENTS"` dict
with as single key, the name will be ignored.
In the previous example the key `"Value"` is used instead of `"COMPONENTS"`:
it is a shorter form for `"COMPONENTS":{"any":{...}}`.
To avoid confusion you cannot specify both `"COMPONENTS"` and `"Value"` in the
 same dict.

To speed up the declarations of the components the `plumedCommunications` module
contains a submodule `defaults` with the default dictionaries already set up:
 - plumedCommunications.defaults.COMPONENT={"period":None, "derivative":True}
 - plumedCommunications.defaults.COMPONENT_NODEV={"period":None, "derivative":False}

\par The calculate function

The calculate funtion must, as all the other functions accept a
PLMD.PythonCVInterface object as only input.

The calculate function must either return a float or a tuple or, in the case of
multiple components, a dict whose keys are the name of the components, whose
elements are either float or tuple.

Plumed will assign automatically assign the result to the CV (to the key named
element), if the name of the component is missing the calculation will be
interrupted with an error message.
If derivatives are disabled it will expect a float(or a double).
In case of activated derivatives it will interrupt the calculation if the
return value would not be a tuple.
The tuple should be (float, ndArray(nat,3),ndArray(3,3)) with the first
elements the value, the second the atomic derivatives and the third the box
derivative (that can also have shape(9), with format (x_x, x_y, x_z, y_x, y_y,
y_z, z_x, z_y, z_z)), if the box derivative are not present a WARNING will be
raised, but the calculation won't be interrupted.

\par The prepare and update functions and the "data" attribute

If the `PREPARE` keyword is used, the defined function will be called at
prepare time, before calculate.
The prepare dictionary can contain a `"setAtomRequest"` key with a parseable
ATOM string, like in the input (or a list of indexes, 0 based).
@code{.py}
#this , with "PREPARE=changeAtom" in the plumed file will select a new atom at each new step
def changeAtom(plmdAction: plumedCommunications.PythonCVInterface):
    toret = {"setAtomRequest": f"1, {int(plmdAction.getStep()) + 2}"}
    if plmdAction.getStep() == 3:
        toret["setAtomRequest"] = "1,2"
    return toret
@endcode

If the `UPDATE` keyword is used, the defined function will be called at update
time, after calculate. As now plumed will ignore the return of this function
(but it stills need to return a dict) and it is intended to accumulate things
or post process data afer calculate

In the example `plmdAction.data["pycv"]=0` is intialized in `pyinit` and its
value is updated in calculate.

\plumedfile
cv1: PYCVINTERFACE  ...
  ATOMS=@mdatoms
  IMPORT=pycvPersistentData
  CALCULATE=pydist
  INIT=pyinit
...

PRINT FILE=colvar.out ARG=*
\endplumedfile

@code{.py}
import plumedCommunications as PLMD
from plumedCommunications.defaults import COMPONENT_NODEV

def pyinit(plmdAction: PLMD.PythonCVInterface):
    plmdAction.data["pycv"]=0
    print(f"{plmdAction.data=}", file=log)
    return {"Value":COMPONENT_NODEV}

def pydist(plmdAction: PLMD.PythonCVInterface):
    plmdAction.data["pycv"]+=plmdAction.getStep()
    d=plmdAction.data["pycv"]
    return d
@endcode

The plumedCommunications.PythonCVInterface has a `data` attribute that is a
dictionary and can be used to store data during the calculations


\par Getting the manual
You can obtain the manual of the module (and of its components)
by running plumed driver with this plumed file
\plumedfile
LOAD GLOBAL FILE=PythonCVInterface.so
cvdist: PYCVINTERFACE IMPORT=pyhelp
PRINT FILE=colvar.out ARG=*
\endplumedfile
and this py module
@code{.py}
import plumedCommunications
import pydoc

def plumedInit(_):
    with open('PythonCVInterface.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications.PythonCVInterface)
    with open('plumedCommunications.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications)
    with open('plumedCommunications.defaults.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications.defaults)
    return {"Value":plumedCommunications.defaults.COMPONENT_NODEV, "ATOMS":"1"}

def plumedCalculate(_):
    return 0
@endcode


\par Tips and tricks and examples

Automatic differentiation and transparent compilation (including to
GPU) can be performed via Google's [JAX
library](https://github.com/google/jax): see the example below.


The following input tells PLUMED to print the distance between atoms 1
and 4.

\plumedfile
cv1: PYTHONCV ATOMS=1,4 IMPORT=distcv CALCULATE=cv
PRINT FILE=colvar.out ARG=*
\endplumedfile

The file `distcv.py` should contain something as follows.

@code{.py}
import numpy as np
import plumedCommunications as PLMD

plumedInit={"Value":PLMD.defaults.COMPONENT}
# Define the distance function
def dist_f(x):
    r = x[0,:]-x[1,:]
    d2 = np.dot(r,r)
    return np.sqrt(d2)

def grad_dist(x):
    d = dist_f(x)
    r = x[0,:]-x[1,:]
    g = r/d
    return np.array([g,-g])

# The CV function actually called
def cv(action:PLMD.PythonCVInterface):
    return dist_f(action.getPositions()), grad_dist(action.getPositions())

@endcode


\par JAX for automatic differentiation and compilation

Automatic differentiation and transparent compilation (including to
GPU) can be performed via Google's [JAX
library](https://github.com/google/jax). In a nutshell, it's sufficient
to replace `numpy` with `jax.numpy`. See the following example.


\plumedfile
cv1: PYTHONCV ATOMS=1,2,3 IMPORT=jaxcv CALCULATE=angle
PRINT FILE=colvar.out ARG=*

\endplumedfile


And, in `jaxcv.py`...

@code{.py}
# Import the JAX library
import jax.numpy as np
from jax import grad, jit, vmap
import plumedCommunications as PLMD

plumedInit={"Value":PLMD.defaults.COMPONENT}

# Implementation of the angle function
def angle_f(x):
    r1 = x[0,:]-x[1,:]
    r2 = x[2,:]-x[1,:]

    costheta = np.dot(r1,r2) / np.linalg.norm(r1) / np.linalg.norm(r2)
    theta = np.arccos(costheta)
    return theta

# Use JAX to auto-gradient it
angle_grad = grad(angle_f)

def cv(action:PLMD.PythonCVInterface):
    return angle_f(action.getPositions()), angle_grad(action.getPositions())
@endcode


There are however
[limitations](https://github.com/google/jax#current-gotchas), the most
notable of which is that indexed assignments such as `x[i]=y` are not
allowed, and should instead be replaced by functional equivalents such
as `x=jax.ops.index_update(x, jax.ops.index[i], y)`.


\par Multiple components

It is possible to return multiple components at a time. This may be
useful e.g. if they reuse part of the computation. To do so, pass the
declare the `"COMPONENTS"` key in the init dict, and assign to each key a dict
with the same rules of the key `"Value"`. In this case, the function must
return a dict with the component names as keys, see the "calculate" paragraph.
Inside PLUMED, component names
will be prefixed by `py-`.

Note that you can use JAX's Jacobian function `jax.jacrev()` to
conveniently compute the gradients all at once (see regtests). For
example:

\plumedfile
cv1:  PYTHONCV ATOMS=1,3,4 IMPORT=distcv FUNCTION=cv COMPONENTS=d12,d13
\endplumedfile

@code{.py}
import jax.numpy as np
from jax import jacrev, jit
import plumedCommunications as PLMD

plumedInit = dict(
    COMPONENTS=dict(d12=PLMD.defaults.COMPONENT, d13=PLMD.defaults.COMPONENT)
)

# Define the distance function
@jit
def dist_f(X):
    return {
     'd12': np.linalg.norm( X[0,:]-X[1,:] ),
     'd13': np.linalg.norm( X[0,:]-X[2,:] )
     }

dist_grad=jacrev(dist_f)

def cv(action: PLMD.PythonCVInterface):
    toret = {}
    d=dist_f(action.getPositions())
    g=dist_grad(action.getPositions())
    for key in ["d12","d13"]:
        toret[key] = (d[key],g[key])
    return  toret
@endcode

\par Installation

A use of an virtual environment or equivalent is recomended.
To compile pycv you just need numpy and pybind11, jax is not necessary for
compilation and installation.

To compile the shared object library you need also plumed in your path.
You need to export the following  environmental variables:
\verbatim
export PLUMED_MKLIB_CFLAGS="$(python3-config --cflags --embed) $(python -m pybind11 --includes)"
export PLUMED_MKLIB_LDFLAGS="$(python3-config --ldflags --embed)"
\endverbatim
and then compile the shared object:
\verbatim
plumed mklib PythonCVInterface.cpp ActionWithPython.cpp PlumedPythonEmbeddedModule.cpp
\endverbatim

If you are on linux you can use pycv only with a plumed version that is
compatible with the `GLOBAL` keyword for the action `LOAD`

*/
//+ENDPLUMEDOC


#define vdbg(...)                                                              \
  std::cerr << std::boolalpha<<std::setw(4) << __LINE__ << ":" << std::setw(20)                \
            << #__VA_ARGS__ << " " << (__VA_ARGS__) << '\n'

namespace py = pybind11;

using std::string;
using std::vector;

namespace PLMD {
namespace pycv {

PLUMED_REGISTER_ACTION(PythonCVInterface, "PYCVINTERFACE")

void PythonCVInterface::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of atoms to be passed to the function");
  //NL
  keys.add("atoms","GROUPA","First list of atoms for the neighbourlist");
  keys.add("atoms","GROUPB","Second list of atoms for the neighbourlist (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  //python components
  keys.add("hidden","COMPONENTS","if provided, the function will return multiple components, with the names given");
  keys.addOutputComponent(PYCV_COMPONENTPREFIX.data(),"COMPONENTS","Each of the components output py the Python code, prefixed by py-");
  //python calling
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","CALCULATE",PYCV_DEFAULTCALCULATE,"the function to call as calculate method of a CV");
  keys.add("compulsory","INIT",PYCV_DEFAULTINIT,"the function to call during the construction method of the CV");
  // python: add other callable methods
  keys.add("compulsory","PREPARE",PYCV_NOTIMPLEMENTED,"the function to call as prepare method of the CV");
  keys.add("compulsory","UPDATE", PYCV_NOTIMPLEMENTED,"the function to call as update() method of the CV");

  // NOPBC is in Colvar!
}

//TODO: add callable checks!!!

PythonCVInterface::PythonCVInterface(const ActionOptions&ao) ://the catch only applies to pybind11 things
  PLUMED_COLVAR_INIT(ao),
  ActionWithPython(ao) {
  try {
    py::gil_scoped_acquire gil;
    //Loading the python module
    std::string import;
    parse("IMPORT",import);
    //setting up the calculate function
    std::string calculateFunName;
    parse("CALCULATE",calculateFunName);
    log.printf("  will import %s and call function %s\n", import.c_str(),
               calculateFunName.c_str());
    // Initialize the module and function pointers
    pyModule = py::module::import(import.c_str());
    if (!py::hasattr(pyModule,calculateFunName.c_str())) {
      error("the function " + calculateFunName + " is not present in "+ import);
    }

    pyCalculate = pyModule.attr(calculateFunName.c_str());
    std::string initFunName;
    parse("INIT",initFunName);
    py::dict initDict;
    if(py::hasattr(pyModule,initFunName.c_str())) {
      log.printf("  will use %s during the initialization\n", initFunName.c_str());
      auto initFcn = pyModule.attr(initFunName.c_str());
      if (py::isinstance<py::dict>(initFcn)) {
        initDict = initFcn;
      } else {
        initDict = initFcn(this);
      }
    } else if(initFunName!=PYCV_DEFAULTINIT) {
      //If the default INIT is not preset, is not a problem
      error("the function "+ initFunName + " is not present in "+ import);
    }

    std::string prepareFunName;
    parse("PREPARE",prepareFunName);
    if (prepareFunName!=PYCV_NOTIMPLEMENTED) {
      if (!py::hasattr(pyModule,prepareFunName.c_str())) {
        error("the function " + prepareFunName + " is not present in "+ import);
      }
      hasPrepare=true;
      pyPrepare=pyModule.attr(prepareFunName.c_str());
      log.printf("  will use %s while calling prepare() before calculate()\n", prepareFunName.c_str());
    }

    std::string updateFunName;
    parse("UPDATE",updateFunName);
    if (updateFunName!=PYCV_NOTIMPLEMENTED) {
      if (!py::hasattr(pyModule,updateFunName.c_str())) {
        error("the function " + updateFunName + " is not present in " + import);
      }
      pyUpdate=pyModule.attr(updateFunName.c_str());
      hasUpdate=true;
      log.printf("  will use %s while calling update() after calculate()\n", updateFunName.c_str());
    }

    {
      std::vector<std::string> components;
      parseVector("COMPONENTS", components);
      if (components.size()>1) {
        error("Please define multiple COMPONENTS from INIT in python.");
      }
    }

    if(initDict.contains("COMPONENTS")) {
      if(initDict.contains("Value")) {
        error("The initialize dict cannot contain both \"Value\" and \"COMPONENTS\"");
      }
      if(!py::isinstance<py::dict>(initDict["COMPONENTS"])) {
        error("COMPONENTS must be a dictionary using with the name of the components as keys");
      }
      py::dict components=initDict["COMPONENTS"];
      for(auto comp: components) {
        auto settings = py::cast<py::dict>(comp.second);
        if(components.size()==1) { //a single component
          initializeValue(dynamic_cast<::PLMD::ActionWithValue&>(*this), settings);
          valueSettings(settings,getPntrToValue());
        } else {
          auto name=std::string(PYCV_COMPONENTPREFIX)
                    +"-"+py::cast<std::string>(comp.first);
          initializeComponent(dynamic_cast<::PLMD::ActionWithValue&>(*this),
                              name,
                              settings);
          valueSettings(settings,getPntrToComponent(name));
        }
      }

    } else if(initDict.contains("Value")) {
      py::dict settingsDict=initDict["Value"];
      initializeValue(dynamic_cast<::PLMD::ActionWithValue&>(*this),settingsDict);
      valueSettings(settingsDict,getPntrToValue());
    } else {
      warning("  WARNING: by defaults components periodicity is not set and component is added without derivatives - see manual\n");
      //this will crash with an error, beacuse periodicity is not explicitly set
      addValue();
    }

    std::vector<AtomNumber> atoms;
    pyParseAtomList("ATOMS",initDict,atoms);
    std::vector<AtomNumber> groupA;
    pyParseAtomList("GROUPA",initDict,groupA);
    std::vector<AtomNumber> groupB;
    pyParseAtomList("GROUPB",initDict,groupB);

    if(atoms.size() !=0 && groupA.size()!=0) {
      error("you can choose only between using the neigbourlist OR the atoms");
    }

    if(atoms.size()==0&& groupA.size()==0 && groupB.size()==0) {
      error("At least one atom is required");
    }

    if (atoms.size() != 0 && groupA.size() != 0) {
      error("you can choose only between using the neigbourlist OR the atoms");
    }

    if (atoms.size() == 0 && groupA.size() == 0 && groupB.size() == 0) {
      error("At least one atom is required");
    }

    bool nopbc;
    pyParseFlag("NOPBC",initDict, nopbc);
    pbc = !nopbc;

    if (groupA.size() > 0) {
      // parse the NL things only in the NL case
      bool dopair;
      pyParseFlag("PAIR",initDict, dopair);
      // this is a WIP
      bool serial = false;
      bool doneigh;
      pyParseFlag("NLIST",initDict,doneigh);
      double nl_cut = 0.0;
      int nl_st = 0;
      if (doneigh) {
        pyParse("NL_CUTOFF", initDict, nl_cut);
        if (nl_cut <= 0.0) {
          error("NL_CUTOFF should be explicitly specified and positive");
        }
        pyParse("NL_STRIDE",initDict, nl_st);
        if (nl_st <= 0) {
          error("NL_STRIDE should be explicitly specified and positive");
        }
      }
      // endof WIP
      if (groupB.size() > 0) {
        if (doneigh)
          nl = Tools::make_unique<NeighborList>(
                 groupA, groupB, serial, dopair, pbc, getPbc(), comm, nl_cut, nl_st);
        else
          nl = Tools::make_unique<NeighborList>(groupA, groupB, serial, dopair,
                                                pbc, getPbc(), comm);
      } else {
        if (doneigh)
          nl = Tools::make_unique<NeighborList>(groupA, serial, pbc, getPbc(),
                                                comm, nl_cut, nl_st);
        else
          nl = Tools::make_unique<NeighborList>(groupA, serial, pbc, getPbc(),
                                                comm);
      }
      requestAtoms(nl->getFullAtomList());
    } else {
      requestAtoms(atoms);
    }

    if (getNumberOfComponents()>1) {
      log.printf("  it is expected to return dictionaries with %d components\n",
                 getNumberOfComponents());
    }

    log << "  Bibliography " << plumed.cite(PYTHONCV_CITATION) << "\n";
    // NB: the NL kewywords will be counted as error when using ATOMS
    checkRead();
  } catch (const py::error_already_set &e) {
    error(e.what());
    //vdbg(e.what());
  }
}

void PythonCVInterface::prepare() {
  try {
    if (nl) {
      if (nl->getStride() > 0) {
        if (firsttime || (getStep() % nl->getStride() == 0)) {
          requestAtoms(nl->getFullAtomList());
          invalidateList = true;
          firsttime = false;
        } else {
          requestAtoms(nl->getReducedAtomList());
          invalidateList = false;
          if (getExchangeStep())
            error("Neighbor lists should be updated on exchange steps - choose a "
                  "NL_STRIDE which divides the exchange stride!");
        }
        if (getExchangeStep()) {
          firsttime = true;
        }
      }
    }
    if (hasPrepare) {
      py::gil_scoped_acquire gil;
      py::dict prepareDict = pyPrepare(this);
      if (prepareDict.contains("setAtomRequest")) {
        //should I use "interpretAtomList"?
        std::vector<PLMD::AtomNumber> myatoms;
        if(py::isinstance<py::tuple>(prepareDict["setAtomRequest"])||
            py::isinstance<py::list>(prepareDict["setAtomRequest"])) {
          py::tuple t = prepareDict["setAtomRequest"];

          for (const auto &i : t) {
            auto at = PLMD::AtomNumber::index(i.cast<unsigned>());
            myatoms.push_back(at);
          }
        } else {
          auto atomlist=PLMD::Tools::getWords(
                          py::str(prepareDict["setAtomRequest"]).cast<std::string>(),
                          "\t\n ,");
          interpretAtomList( atomlist, myatoms );
        }
        requestAtoms(myatoms);
      }
    }
  } catch (const py::error_already_set &e) {
    plumed_merror(e.what());
  }
}

void PythonCVInterface::update() {
  try {
    if(hasUpdate) {
      py::gil_scoped_acquire gil;
      py::dict updateDict=pyUpdate(this);
      //See what to do here
    }
  } catch (const py::error_already_set &e) {
    plumed_merror(e.what());
  }
}

// calculator
void PythonCVInterface::calculate() {
  try {
    if (nl) {
      if (nl->getStride() > 0 && invalidateList) {
        nl->update(getPositions());
      }
    }
    py::gil_scoped_acquire gil;
    // Call the function
    py::object r = pyCalculate(this);
    if(getNumberOfComponents()>1) {		// MULTIPLE NAMED COMPONENTS
      calculateMultiComponent(r);
    } else { // SINGLE COMPONENT
      readReturn(r, getPntrToValue());
    }

  } catch (const py::error_already_set &e) {
    plumed_error_nested()<<"("<<getLabel()<<") caught a python exception:\n"<<e.what();
  }
}

void PythonCVInterface::readReturn(const py::object &r, Value* valPtr) {
  // Is there more than 1 return value?
  if (py::isinstance<py::tuple>(r)||py::isinstance<py::list>(r)) {
    // 1st return value: CV
    py::list rl=r.cast<py::list>();
    pycvComm_t value = rl[0].cast<pycvComm_t>();
    valPtr->set(value);
    //shape returns long int
    auto natoms = static_cast<long int > (getPositions().size());
    if (rl.size() > 1) {
      if(!valPtr->hasDerivatives()) {
        error(valPtr->getName()+" was declared without derivatives, but python returned with derivatives");
      }
      // 2nd return value: gradient: numpy array of (natoms, 3)
      py::array_t<pycvComm_t> grad(rl[1]);
      // Assert correct gradient shape
      if (grad.ndim() != 2 || grad.shape(0) != natoms || grad.shape(1) != 3) {
        log.printf("Error: wrong shape for the gradient return argument: should be "
                   "(natoms=%d,3), received %ld x %ld\n",
                   natoms, grad.shape(0), grad.shape(1));
        error("Python CV returned wrong gradient shape error");
      }
      // To optimize, see "direct access"
      // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
      for (unsigned i = 0; i < natoms; i++) {
        Vector3d gi(grad.at(i, 0), grad.at(i, 1), grad.at(i, 2));
        setAtomsDerivatives(valPtr, i, gi);
      }
    } else if (valPtr->hasDerivatives()) {
      error(valPtr->getName()+" was declared with derivatives, but python returned none");
    }

    if (rl.size() > 2) {
      if(!valPtr->hasDerivatives()) {
        plumed_merror(valPtr->getName()+" was declared without derivatives, but python returned with box derivatives");
      }
      py::array_t<pycvComm_t> pyBoxDev(rl[2]);
      // expecting the box derivatives
      Tensor boxDev;
      if (pyBoxDev.ndim() == 2 &&
          (pyBoxDev.shape(0) == 3 && pyBoxDev.shape(1) == 3)) { // boxDev is 3x3
        boxDev =
          Tensor({pyBoxDev.at(0, 0), pyBoxDev.at(0, 1), pyBoxDev.at(0, 2),
                  pyBoxDev.at(1, 0), pyBoxDev.at(1, 1), pyBoxDev.at(1, 2),
                  pyBoxDev.at(2, 0), pyBoxDev.at(2, 1), pyBoxDev.at(2, 2)});
      } else if (pyBoxDev.ndim() == 1 && pyBoxDev.shape(0) == 9) {
        boxDev = Tensor({pyBoxDev.at(0), pyBoxDev.at(1), pyBoxDev.at(2),
                         pyBoxDev.at(3), pyBoxDev.at(4), pyBoxDev.at(5),
                         pyBoxDev.at(6), pyBoxDev.at(7), pyBoxDev.at(8)});
      } else {
        log.printf(
          "Error: wrong shape for the box derivatives return argument: "
          "should be (size 3,3 or 9), received %ld x %ld\n",
          natoms, pyBoxDev.shape(0), pyBoxDev.shape(1));
        error("Python CV returned wrong box derivatives shape error");
      }
      setBoxDerivatives(valPtr, boxDev);
    } else if (valPtr->hasDerivatives()) {
      warning(valPtr->getName()+" was declared with derivatives, but python returned no box derivatives");
    }
  } else {
    // Only value returned. Might be an error as well.
    if (valPtr->hasDerivatives()) {
      warning(BIASING_DISABLED);
    }
    pycvComm_t value = r.cast<pycvComm_t>();
    valPtr->set(value);
  }
  //TODO: is this ok?
  if (!pbc) {
    setBoxDerivativesNoPbc(valPtr);
  }
}


void PythonCVInterface::calculateMultiComponent(py::object &r) {

  const auto nc = getNumberOfComponents();
  if (py::isinstance<py::dict>(r)) {
    py::dict dataDict = r.cast<py::dict>(); // values
    for(int i=0; i < nc; ++i) {
      auto component=getPntrToComponent(i);
      //get the without "label.prefix-"
      std::string key=component->getName().substr(
                        2 + getLabel().size()
                        +PYCV_COMPONENTPREFIX.size());
      if (dataDict.contains(key.c_str())) {
        readReturn(dataDict[key.c_str()], component);
      } else {
        error( "python did not returned " + key );
      }
    }
  } else {
    // In principle one could handle a "list" return case.
    error("Multi-components pyCVs need to return dictionaries");
  }
}

void PythonCVInterface::pyParseAtomList(const char* key, const ::pybind11::dict &initDict, std::vector<AtomNumber> &myatoms) {
  parseAtomList(key,myatoms);

  if(initDict.contains(key)) {
    if (myatoms.size()>0) {
      error(std::string("you specified the same keyword ").append(key)+ " both in python and in the settings file");
    }
    auto atomlist=PLMD::Tools::getWords(
                    py::str(initDict[key]).cast<std::string>(),
                    "\t\n ,");
    interpretAtomList( atomlist, myatoms );
  }
}

NeighborList &PythonCVInterface::getNL() {
  return *nl;
}
} // namespace pycvs
} // namespace PLMD

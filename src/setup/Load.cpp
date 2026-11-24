/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionAnyorder.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"

namespace PLMD {
namespace setup {

//+PLUMEDOC GENERIC LOAD
/*
Loads a library, possibly defining new actions.

The LOAD action is only available on systems that allow dynamic loading.
This action allows you load new funcionality into PLUMED at runtime.
This new functionality can be in a .so file or a .cpp file. If the
new functionality is in a cpp file then the code is compiled and the
the resulting object is loaded.

Extra header search directories can be supplied when compiling a `.cpp` file:

````
LOAD FILE=x.cpp INCLUDE_PATH={headers:other headers}
````

`INCLUDE_PATH` is a colon-separated, ordered list. Enclose the whole value in
braces when a directory contains spaces. Relative directories are resolved from
the working directory, not from the source or input file directory. Empty values
and empty entries are errors; directory names containing colons are not supported.
The option is rejected for precompiled libraries. Existing compiler flags and
include environment variables remain in effect; these additional directories
apply only to this LOAD invocation.

With `INCLUDE_PATH`, PLUMED compiles on every LOAD into a separate temporary
location next to the source and removes the temporary library after loading.
This avoids reusing a library compiled with different include paths, including
when multiple PLUMED objects load the same source concurrently, and reflects
header changes on the next LOAD. The source directory must be writable.
Without `INCLUDE_PATH`, the existing source-timestamp cache is unchanged (it does
not track header dependencies).

Using the LOAD action is useful for making quick tests while developing your own
actions. Of course, once your implementation is ready you should
add it to the PLUMED source tree and recompile PLUMED.

One way to use the LOAD action is to directly load a cpp file as illustrated in the input below.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt15/Distance2.cpp
# load the new definition
# this is a cpp file so it will be compiled
LOAD FILE=regtest/basic/rt15/Distance2.cpp
# compute standard distance
d: DISTANCE ATOMS=1,10
# compute modified distance
d2: DISTANCE2 ATOMS=1,10
# print them on a file
PRINT ARG=d,d2 FILE=compare-them
```

When PLUMED reads the input above it first compiles the code in `Distance2.cpp`. The resulting object file
is then loaded into PLUMED. If you look at the cpp file that is input in the command above you can see that
a new action called `DISTANCE2` is defined within it.  We can thus use this new action in the input above.

Instead of compiling the code at runtime you can construct a shared object from the Distance2.cpp file and LOAD the shared library instead of
the cpp file.  To construct the shared object you would use the following command:

````
> plumed mklib Distance2.cpp
````

When you run this command a new file called `Distance2.so` (or `Distance2.dylib` on a mac) is created.
You can then load this shared object in PLUMED by using the following LOAD command

````
LOAD FILE=Distance2.so
````

The new functionality within the Distance2.cpp file can then be used in the remainder of the input in the same way it was
used in the previous example.
(Notice that the only reason for not using the usual PLUMED syntax highlightling in the example above is that we have no
file called Distance2.so.)

From PLUMED 2.10 onwards the LOAD action can be placed at any point of the input
file. The loaded functionality will then only affect commands that are placed after the LOAD action.
These features are illustrated in the input below:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt15/Distance3.cpp
# compute standard distance
d: DISTANCE ATOMS=1,10
# load the new definition
LOAD FILE=regtest/basic/rt15/Distance3.cpp
# compute modified distance
d2: DISTANCE ATOMS=1,10
# print them on a file
PRINT ARG=d,d2 FILE=compare-them
```

Notice that `Distance3.cpp` replaces DISTANCE2 in the command that was used in `Distance2.cpp`; namely:

```c++
PLUMED_REGISTER_ACTION(Distance,"DISTANCE2")
```

with DISTANCE. Consequently, when we load the `Distance3.cpp` file here we redefine the DISTANCE command.
The functions that compute `d` and `d2` in the above input are thus different.

A final point to note is that, starting with PLUMED 2.10, the LOAD action can be used in contexts where
multiple Plumed objects exist. A possible example is multithreading: loading an action
from a Plumed object used in one thread will not affect other threads.
Another example is if multiple Plumed objects are created in the C/C++ or Python interface.
If a LOAD command is used in one of these objects, the loaded action will not affect
the other objects.


## Code written for older versions of PLUMED

The C++ interface that is used by the files you LOAD changes between versions of PLUMED.
If a file that you wrote for an older version does not compile, the following changes are usually enough:

- `#include "core/Atoms.h"` and `plumed.getAtoms()`: the Atoms class was removed in PLUMED 2.10.
  The include is removed automatically by `plumed mklib` (with a warning). Use the methods of the action instead:
  `getKBoltzmann()` instead of `plumed.getAtoms().getKBoltzmann()`, `getkBT()` instead of `plumed.getAtoms().getKbT()`
  (`getkBT()` also reads the TEMP keyword), `getUnits()`, `usingNaturalUnits()` and `getTimeStep()`.
- `error: invalid use of incomplete type 'class PLMD::Communicator'`: since PLUMED 2.10 the headers of the actions
  do not include the definition of the communicators any more. If your code uses `comm` or `multi_sim_comm`, add

  ```c++
  #include "tools/Communicator.h"
  ```

- `componentsAreNotOptional(keys)` was removed. It only changed the text of the manual, so you can delete the line.
- `keys.addOutputComponent(name,flag,description)` is deprecated: give the type of the component as third argument,
  e.g. `keys.addOutputComponent("x","COMPONENTS","scalar","the x component")`.
- `keys.use("ARG")` after the `registerKeywords` of the base class is not needed any more (ARG is added by the base class),
  but it is still accepted.
- A flag cannot be on by default: `keys.addFlag("COMPONENTS",true,...)` stops with "the second argument to addFlag must be false".
  Use `false` and, if needed, a flag with the opposite meaning.
- Include paths such as `#include "colvar/ActionRegister.h"` are replaced by `#include "core/ActionRegister.h"`
  automatically by `plumed mklib` (with a warning).
- The classes of the old multicolvar module (`multicolvar/MultiColvarBase.h`) were removed in PLUMED 2.10.
  A collective variable that is computed for many groups of atoms is now written for one group of atoms and then
  used with the `MultiColvarTemplate` class: see how `src/colvar/Distance.cpp` defines DISTANCE_SCALAR and DISTANCE_VECTOR.
- Some helper functions of the tools that PLUMED itself did not use were removed in PLUMED 2.11:
  in `PDB` the setters (`setAtomNumbers`, `setAtomPositions`, `setArgumentNames`, `setArgumentValue`, `addBlockEnd`), `print`, `hasFlag`, `getMtype`,
  `getAtomsInChain`, `checkForResidue`, `getChainID` of a residue number (use the one of an atom), and the versions without a chain
  of `getResidueName` and `getNamedAtomFromResidue` (use the ones that take a chain, with `"*"` for any chain);
  in `Grid` `addKernel`, `writeCubeFile`, `integrate`, `findSetOfPointsOnContour`, `logAllValuesAndDerivatives` and `mpiSumValuesAndDerivatives`;
  in `RMSD` `calc_SOMA`, `getAlign` and `getDisplace`; in `HistogramBead` `generateBins` and `setKernelType`
  (give the kernel type to the constructor); in `KernelFunctions` the constructors from a string and from a pointer, and `normalize`.

!!! note ""

    The example inputs on this page appear as not working because you cannot use the dynamic loading on GitHub Actions. On a machine where this functionality is available these inputs should work.

*/
//+ENDPLUMEDOC

class Load :
  public virtual ActionAnyorder {
public:
  static void registerKeywords( Keywords& keys );
  explicit Load(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Load,"LOAD")

void Load::registerKeywords( Keywords& keys ) {
  ActionAnyorder::registerKeywords(keys);
  keys.add("optional","FILE","file to be loaded");
  keys.add("optional","INCLUDE_PATH","colon-separated header search directories for a .cpp file; relative to the working directory");
  keys.add("optional","MODULE","installed module to be loaded, together with its hard dependencies");
}

Load::Load(const ActionOptions&ao):
  Action(ao),
  ActionAnyorder(ao) {
  std::string f, module;
  parse("FILE",f);
  parse("MODULE",module);
  const bool hasIncludePath=!getKeyword("INCLUDE_PATH").empty();
  std::string includePath;
  parse("INCLUDE_PATH",includePath);
  // An explicitly empty keyword must not silently select the ordinary cache.
  if(hasIncludePath && includePath.empty()) {
    error("INCLUDE_PATH must contain at least one directory");
  }
  std::vector<std::string> includeDirs;
  if(hasIncludePath) {
    size_t start=0;
    while(true) {
      const auto end=includePath.find(':',start);
      const auto dir=includePath.substr(start,end==std::string::npos ? end : end-start);
      if(dir.empty()) {
        error("INCLUDE_PATH cannot contain empty directories");
      }
      includeDirs.push_back(dir);
      if(end==std::string::npos) {
        break;
      }
      start=end+1;
    }
    if(f.size()<4 || f.substr(f.size()-4)!=".cpp") {
      error("INCLUDE_PATH is only supported for FILE ending in .cpp");
    }
  }
  checkRead();

  if(f.empty()==module.empty()) {
    error("Specify exactly one of FILE and MODULE");
  }
  if(!module.empty()) {
    plumed.loadModule(module);
  } else {
    plumed.load(f,includeDirs);
  }
}

}
}


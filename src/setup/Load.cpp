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

It is available only
on systems allowing for dynamic loading. It can also be fed with a cpp file,
in which case the file is compiled first.

\par Examples

If you have a shared object named extensions.so and want to
use the functions implemented within it within PLUMED you can
load it with the following syntax

\plumedfile
LOAD FILE=extensions.so
\endplumedfile

As a more practical example, imagine that you want to make a
small change to one collective variable that is already implemented
in PLUMED, say \ref DISTANCE . Copy the file `src/colvar/Distance.cpp`
into your work directory, rename it as `Distance2.cpp`
and  edit it as you wish. It might be better
to also replace any occurrence of the string DISTANCE within the file
with DISTANCE2, so that both old and new implementation will be available
with different names. Then you can compile it into a shared object using
\verbatim
> plumed mklib Distance2.cpp
\endverbatim
This will generate a file `Distance2.so` (or `Distance2.dylib` on a mac)
that can be loaded.
Now you can use your new implementation with the following input
\plumedfile
# load the new library
LOAD FILE=Distance2.so
# compute standard distance
d: DISTANCE ATOMS=1,10
# compute modified distance
d2: DISTANCE2 ATOMS=1,10
# print them on a file
PRINT ARG=d,d2 FILE=compare-them
\endplumedfile

You can even skip the initial step and directly feed PLUMED
with the `Distance2.cpp` file: it will be compiled on the fly.
\plumedfile
# load the new definition
# this is a cpp file so it will be compiled
LOAD FILE=Distance2.cpp
# compute standard distance
d: DISTANCE ATOMS=1,10
# compute modified distance
d2: DISTANCE2 ATOMS=1,10
# print them on a file
PRINT ARG=d,d2 FILE=compare-them
\endplumedfile

This will allow to make quick tests while developing your own
variables. Of course, after your implementation is ready you might
want to add it to the PLUMED source tree and recompile
the whole PLUMED.

Starting with PLUMED 2.10, the LOAD action can be placed in any point of the input
file, and will only affect commands that are placed after the LOAD action.
In other words, you can create a file named `Distance.cpp` and that reimplement
the \ref DISTANCE action and use it like this:

\plumedfile
# compute standard distance
d: DISTANCE ATOMS=1,10
# load the new definition
LOAD FILE=Distance.so
# compute modified distance
d2: DISTANCE ATOMS=1,10
# print them on a file
PRINT ARG=d,d2 FILE=compare-them
\endplumedfile

In addition, starting with PLUMED 2.10 the LOAD action can be used in contexts where
multiple Plumed objects exist. A possible example is multithreading: loading an action
from a Plumed object used in one thread will not affect other threads.
Another example is if multiple Plumed objects are created in the C/C++ or Python interface.
If a LOAD command is used in one of these objects, the loaded action will not affect
the other objects.


*/
//+ENDPLUMEDOC

class Load :
  public virtual ActionAnyorder
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Load(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Load,"LOAD")

void Load::registerKeywords( Keywords& keys ) {
  ActionAnyorder::registerKeywords(keys);
  keys.add("compulsory","FILE","file to be loaded");
}

Load::Load(const ActionOptions&ao):
  Action(ao),
  ActionAnyorder(ao)
{
  std::string f;
  parse("FILE",f);
  checkRead();
  plumed.load(f);
}

}
}


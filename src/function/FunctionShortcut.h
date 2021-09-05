/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_function_FunctionShortcut_h
#define __PLUMED_function_FunctionShortcut_h

#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionShortcut : public ActionShortcut {
private:
/// The function that is being computed
  T myfunc;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionShortcut(const ActionOptions&);
};

template <class T>
void FunctionShortcut<T>::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ARG","the input to this function");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc; tfunc.registerKeywords( keys );
}

template <class T>
FunctionShortcut<T>::FunctionShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  std::vector<std::string> args; parseVector("ARG",args);
  plumed_assert( args.size()==1 ); plumed_assert( args[0].find_first_of('.')==std::string::npos );
  ActionWithValue* inpt = plumed.getActionSet().template selectWithLabel<ActionWithValue*>(args[0]); plumed_assert( inpt ); Value* val=inpt->copyOutput(0);
  if( val->getRank()==1 && !val->hasDerivatives() ) {
      if( actionRegister().check( getName() + "_VECTOR") ) readInputLine( getShortcutLabel() + ": " + getName() + "_VECTOR ARG=" + args[0] + " " + convertInputLineToString() );
      else plumed_merror("there is no action registered that allows you to do " + getName() + " with vectors");
  } else if( val->getRank()==2 && !val->hasDerivatives() ) {
      if( actionRegister().check( getName() + "_MATRIX") ) readInputLine( getShortcutLabel() + ": " + getName() + "_MATRIX ARG=" + args[0] + " " + convertInputLineToString() );
      else plumed_merror("there is no action registered that allows you to do " + getName() + " with matrices");
  } else plumed_error();
}

}
}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "core/Value.h"

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
  static void createAction( ActionShortcut* action, const std::vector<Value*>& vals, const std::string& allargs );
};

template <class T>
void FunctionShortcut<T>::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_ONEARG");
  keys.addActionNameSuffix("_VECTOR");
  keys.addActionNameSuffix("_MATRIX");
  keys.addActionNameSuffix("_GRID");
  T tfunc;
  tfunc.registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" || keys.getDisplayName()=="CUSTOM" || keys.getDisplayName()=="MATHEVAL" ) {
    keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix/grid","the values input to this function");
  } else {
    keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the values input to this function");
  }
  if( keys.outputComponentExists(".#!value") && keys.getDisplayName()!="DIFFERENCE" ) {
    keys.addInputKeyword("optional","MASK","vector/matrix","the label for a sparse vector/matrix that should be used to determine which elements of the vector/matrix should be computed");
  }
}

template <class T>
FunctionShortcut<T>::FunctionShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> args;
  parseVector("ARG",args);
  std::string allargs=args[0];
  for(unsigned i=1; i<args.size(); ++i) {
    allargs += "," + args[i];
  }
  std::vector<Value*> vals;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, vals );
  if( vals.size()==0 ) {
    error("found no input arguments to function");
  }
  createAction( this, vals, allargs );
}

template <class T>
void FunctionShortcut<T>::createAction( ActionShortcut* action, const std::vector<Value*>& vals, const std::string& allargs ) {
  unsigned maxrank=vals[0]->getRank();
  bool isgrid=false;
  for(unsigned i=0; i<vals.size(); ++i) {
    if( vals[i]->getRank()>0 && vals[i]->hasDerivatives() ) {
      isgrid=true;
    }
    if( vals[i]->getRank()>maxrank ) {
      maxrank=vals[i]->getRank();
    }
  }
  if( isgrid ) {
    if( actionRegister().check( action->getName() + "_GRID") ) {
      action->readInputLine( action->getShortcutLabel() + ": " + action->getName() + "_GRID ARG=" + allargs + " " + action->convertInputLineToString() );
    } else {
      plumed_merror("there is no action registered that allows you to do " + action->getName() + " with functions on a grid");
    }
  } else if( maxrank==0 ) {
    if( actionRegister().check( action->getName() + "_SCALAR") ) {
      action->readInputLine( action->getShortcutLabel() + ": " + action->getName() + "_SCALAR ARG=" + allargs + " " + action->convertInputLineToString() );
    } else {
      plumed_merror("there is no action registered that allows you to do " + action->getName() + " with scalars");
    }
  } else if( vals.size()==1 && actionRegister().check( action->getName() + "_ONEARG") ) {
    action->readInputLine( action->getShortcutLabel() + ": " + action->getName() + "_ONEARG ARG=" + allargs + " " + action->convertInputLineToString() );
  } else if( maxrank==1 ) {
    if( actionRegister().check( action->getName() + "_VECTOR") ) {
      action->readInputLine( action->getShortcutLabel() + ": " + action->getName() + "_VECTOR ARG=" + allargs + " " + action->convertInputLineToString() );
    } else {
      plumed_merror("there is no action registered that allows you to do " + action->getName() + " with vectors");
    }
  } else if( maxrank==2  ) {
    if( actionRegister().check( action->getName() + "_MATRIX") ) {
      action->readInputLine( action->getShortcutLabel() + ": " + action->getName() + "_MATRIX ARG=" + allargs + " " + action->convertInputLineToString() );
    } else {
      plumed_merror("there is no action registered that allows you to do " + action->getName() + " with matrices");
    }
  } else {
    plumed_error();
  }
}

}
}
#endif

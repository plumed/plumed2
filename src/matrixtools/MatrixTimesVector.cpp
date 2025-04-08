/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MatrixTimesVectorBase.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MCOLVAR MATRIX_VECTOR_PRODUCT
/*
Calculate the product of the matrix and the vector

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class MatrixTimesVector : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesVector(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(MatrixTimesVector,"MATRIX_VECTOR_PRODUCT")

void MatrixTimesVector::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  MatrixTimesVectorBase<Vector>::registerLocalKeywords( keys );
  keys.addActionNameSuffix("_ROWSUMS");
  keys.addActionNameSuffix("_PROPER");
}

MatrixTimesVector::MatrixTimesVector( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> args;
  parseVector("ARG",args);
  std::vector<Value*> myargs;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, myargs );

  unsigned nvectors=0, nmatrices=0;
  for(unsigned i=0; i<myargs.size(); ++i) {
    if( myargs[i]->hasDerivatives() ) {
      error("arguments should be vectors or matrices");
    }
    if( myargs[i]->getRank()<=1 ) {
      nvectors++;
    }
    if( myargs[i]->getRank()==2 ) {
      nmatrices++;
    }
  }
  std::string argstr = args[0];
  for(unsigned i=1; i<args.size(); ++i) {
    argstr += "," + args[i];
  }

  bool sumrows=false;
  if( nvectors==1 ) {
    unsigned n = myargs.size()-1;
    for(unsigned i=0; i<n; ++i) {
      if( myargs[i]->getRank()!=2 || myargs[i]->hasDerivatives() ) {
        error("all arguments other than last argument should be matrices");
      }
      if( myargs[n]->getRank()==0 ) {
        if( myargs[i]->getShape()[1]!=1 ) {
          error("number of columns in input matrix does not equal number of elements in vector");
        }
      } else if( myargs[i]->getShape()[1]!=myargs[n]->getShape()[0] ) {
        std::string str_nmat, str_nvec;
        Tools::convert( myargs[i]->getShape()[1], str_nmat);
        Tools::convert( myargs[n]->getShape()[0], str_nvec );
        error("number of columns in input matrix is " + str_nmat + " which does not equal number of elements in vector, which is " + str_nvec);
      }
    }
    if( myargs[n]->getRank()>0 ) {
      if( myargs[n]->getRank()!=1 || myargs[n]->hasDerivatives() ) {
        error("last argument to this action should be a vector");
      }
    }
    if( (myargs[n]->getPntrToAction())->getName()=="CONSTANT" ) {
      sumrows=true;
      if( myargs[n]->getRank()==0 ) {
        if( fabs( myargs[n]->get() - 1.0 )>epsilon ) {
          sumrows = false;
        }
      } else {
        for(unsigned i=0; i<myargs[n]->getShape()[0]; ++i) {
          if( fabs( myargs[n]->get(i) - 1.0 )>epsilon ) {
            sumrows=false;
            break;
          }
        }
      }
    }
    if( sumrows ) {
      readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_ROWSUMS ARG=" + argstr + " " + convertInputLineToString() );
    } else {
      readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_PROPER ARG=" + argstr + " " + convertInputLineToString() );
    }
  } else if( nmatrices==1 ) {
    if( myargs[0]->getRank()!=2 || myargs[0]->hasDerivatives() ) {
      error("first argument to this action should be a matrix");
    }
    for(unsigned i=1; i<myargs.size(); ++i) {
      if( myargs[i]->getRank()>1 || myargs[i]->hasDerivatives() ) {
        error("all arguments other than first argument should be vectors");
      }
      if( myargs[i]->getRank()==0 ) {
        if( myargs[0]->getShape()[1]!=1 ) {
          error("number of columns in input matrix does not equal number of elements in vector");
        }
      } else if( myargs[0]->getShape()[1]!=myargs[i]->getShape()[0] ) {
        error("number of columns in input matrix does not equal number of elements in vector");
      }
    }
    readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_PROPER ARG=" + argstr + " " + convertInputLineToString() );
  } else {
    error("You should either have one vector or one matrix in input");
  }
}

}
}

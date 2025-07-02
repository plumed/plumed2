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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC

#include "MatrixTimesVectorBase.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MCOLVAR MATRIX_VECTOR_PRODUCT
/*
Calculate the product of the matrix and the vector

Thiis action allows you to [multiply](https://en.wikipedia.org/wiki/Matrix_multiplication) a matrix and a vector together.
This action is primarily used to calculate coordination numbers and symmetry functions, which is what is done by the example below:

```plumed
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
PRINT ARG=cc FILE=colvar
```

Notice that you can use this action to multiply multiple matrices by a single vector as shown below:

```plumed
c1: CONTACT_MATRIX COMPONENTS GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12 D_MAX=10.0}
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1.x,c1.y,c1.z,ones
PRINT ARG=cc.x,cc.y,cc.z FILE=colvar
```

Notice that if you use this options all the input matrices must have the same sparsity pattern.  This feature
was implemented in order to making caluclating Steinhardt parameters such as [Q6](Q6.md) straightforward.

You can also multiply a single matrix by multiple vectors:

```plumed
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12 D_MAX=10.0}
ones: ONES SIZE=7
twos: CONSTANT VALUES=1,2,3,4,5,6,7
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones,twos
PRINT ARG=cc.ones,cc.twos FILE=colvar
```

This feature was implemented to make calculating local averages of the Steinhard parameters straightforward.

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
  std::string usegpustr="";
  {
    bool usegpu;
    parseFlag("USEGPU",usegpu);
    if( usegpu ) {
      usegpustr = " USEGPU";
    }
  }
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
      readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_ROWSUMS ARG="
                     + argstr + " " + convertInputLineToString() + usegpustr );
    } else {
      readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_PROPER ARG="
                     + argstr + " " + convertInputLineToString() + usegpustr  );
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
    readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT_PROPER ARG="
                   + argstr + " " + convertInputLineToString()  + usegpustr );
  } else {
    error("You should either have one vector or one matrix in input");
  }
}

}
}

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
#include "core/ActionWithVector.h"
#include "core/ActionRegister.h"

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

class MatrixTimesVector : public ActionWithVector {
private:
  bool sumrows;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesVector(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  unsigned getNumberOfDerivatives() override ;
  void prepare() override ;
  void calculate() override ;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override ;
  int checkTaskIsActive( const unsigned& itask ) const override ;
  void getNumberOfForceDerivatives( unsigned& nforces, unsigned& nderiv ) const override ;
  void gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const override ;
};

PLUMED_REGISTER_ACTION(MatrixTimesVector,"MATRIX_VECTOR_PRODUCT")

void MatrixTimesVector::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","matrix/vector/scalar","the label for the matrix and the vector/scalar that are being multiplied.  Alternatively, you can provide labels for multiple matrices and a single vector or labels for a single matrix and multiple vectors. In these cases multiple matrix vector products will be computed.");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("vector","the vector that is obtained by taking the product between the matrix and the vector that were input");
  ActionWithValue::useCustomisableComponents(keys);
}

std::string MatrixTimesVector::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( getPntrToArgument(1)->getRank()==1 ) {
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getName().find(cname)!=std::string::npos ) {
        return "the product of the matrix " + getPntrToArgument(0)->getName() + " and the vector " + getPntrToArgument(i)->getName();
      }
    }
  }
  for(unsigned i=0; i<getNumberOfArguments()-1; ++i) {
    if( getPntrToArgument(i)->getName().find(cname)!=std::string::npos ) {
      return "the product of the matrix " + getPntrToArgument(i)->getName() + " and the vector " + getPntrToArgument(getNumberOfArguments()-1)->getName();
    }
  }
  plumed_merror( "could not understand request for component " + cname );
  return "";
}

MatrixTimesVector::MatrixTimesVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  sumrows(false) {
  if( getNumberOfArguments()<2 ) {
    error("Not enough arguments specified");
  }
  unsigned nvectors=0, nmatrices=0;
  bool vectormask=false;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->hasDerivatives() ) {
      error("arguments should be vectors or matrices");
    }
    if( getPntrToArgument(i)->getRank()<=1 ) {
      nvectors++;
      ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
      if( av && av->getNumberOfMasks()>=0 ) {
        vectormask=true;
      }
    }
    if( getPntrToArgument(i)->getRank()==2 ) {
      nmatrices++;
    }
  }
  if( !vectormask ) {
    ignoreMaskArguments();
  }

  std::vector<std::size_t> shape(1);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  if( nvectors==1 ) {
    unsigned n = getNumberOfArguments()-1;
    for(unsigned i=0; i<n; ++i) {
      if( getPntrToArgument(i)->getRank()!=2 || getPntrToArgument(i)->hasDerivatives() ) {
        error("all arguments other than last argument should be matrices");
      }
      if( getPntrToArgument(n)->getRank()==0 ) {
        if( getPntrToArgument(i)->getShape()[1]!=1 ) {
          error("number of columns in input matrix does not equal number of elements in vector");
        }
      } else if( getPntrToArgument(i)->getShape()[1]!=getPntrToArgument(n)->getShape()[0] ) {
        std::string str_nmat, str_nvec;
        Tools::convert( getPntrToArgument(i)->getShape()[1], str_nmat);
        Tools::convert( getPntrToArgument(n)->getShape()[0], str_nvec );
        error("number of columns in input matrix is " + str_nmat + " which does not equal number of elements in vector, which is " + str_nvec);
      }
    }
    if( getPntrToArgument(n)->getRank()>0 ) {
      if( getPntrToArgument(n)->getRank()!=1 || getPntrToArgument(n)->hasDerivatives() ) {
        error("last argument to this action should be a vector");
      }
    }

    if( getNumberOfArguments()==2 ) {
      addValue( shape );
      setNotPeriodic();
    } else {
      for(unsigned i=0; i<getNumberOfArguments()-1; ++i) {
        std::string name = getPntrToArgument(i)->getName();
        if( name.find_first_of(".")!=std::string::npos ) {
          std::size_t dot=name.find_first_of(".");
          name = name.substr(dot+1);
        }
        addComponent( name, shape );
        componentIsNotPeriodic( name );
      }
    }
    if( (getPntrToArgument(n)->getPntrToAction())->getName()=="CONSTANT" ) {
      sumrows=true;
      if( getPntrToArgument(n)->getRank()==0 ) {
        if( fabs( getPntrToArgument(n)->get() - 1.0 )>epsilon ) {
          sumrows = false;
        }
      } else {
        for(unsigned i=0; i<getPntrToArgument(n)->getShape()[0]; ++i) {
          if( fabs( getPntrToArgument(n)->get(i) - 1.0 )>epsilon ) {
            sumrows=false;
            break;
          }
        }
      }
    }
  } else if( nmatrices==1 ) {
    if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
      error("first argument to this action should be a matrix");
    }
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()>1 || getPntrToArgument(i)->hasDerivatives() ) {
        error("all arguments other than first argument should be vectors");
      }
      if( getPntrToArgument(i)->getRank()==0 ) {
        if( getPntrToArgument(0)->getShape()[1]!=1 ) {
          error("number of columns in input matrix does not equal number of elements in vector");
        }
      } else if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(i)->getShape()[0] ) {
        error("number of columns in input matrix does not equal number of elements in vector");
      }
    }

    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      std::string name = getPntrToArgument(i)->getName();
      if( name.find_first_of(".")!=std::string::npos ) {
        std::size_t dot=name.find_first_of(".");
        name = name.substr(dot+1);
      }
      addComponent( name, shape );
      componentIsNotPeriodic( name );
    }
  } else {
    error("You should either have one vector or one matrix in input");
  }
}

unsigned MatrixTimesVector::getNumberOfDerivatives() {
  unsigned nderivatives=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    nderivatives += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nderivatives;
}

void MatrixTimesVector::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] ) {
    return;
  }
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  myval->setShape(shape);
}

void MatrixTimesVector::calculate() {
  runAllTasks();
}

int MatrixTimesVector::checkTaskIsActive( const unsigned& itask ) const {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    Value* myarg = getPntrToArgument(i);
    if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
      return 0;
    } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
      unsigned ncol = myarg->getRowLength(itask);
      unsigned base = itask*myarg->getNumberOfColumns();
      for(unsigned k=0; k<ncol; ++k) {
        if( fabs(myarg->get(base+k,false))>0 ) {
          return 1;
        }
      }
    } else {
      plumed_merror("should not be in action " + getName() );
    }
  }
  return 0;
}

void MatrixTimesVector::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  if( sumrows ) {
    unsigned n=getNumberOfArguments()-1;
    for(unsigned i=0; i<n; ++i) {
      Value* mymat = getPntrToArgument(i);
      unsigned ncol = mymat->getNumberOfColumns();
      unsigned nmat = mymat->getRowLength(task_index);
      double val=0;
      for(unsigned j=0; j<nmat; ++j) {
        val += mymat->get( task_index*ncol + j, false );
      }
      myvals.setValue( i, val );

      // And the derivatives
      if( doNotCalculateDerivatives() ) {
        continue;
      }

      unsigned dloc = task_index*ncol;
      for(unsigned j=0; j<nmat; ++j) {
        myvals.addDerivative( i, dloc + j, 1.0 );
        myvals.updateIndex( i, dloc + j );
      }
    }
  } else if( getPntrToArgument(1)->getRank()==1 ) {
    Value* mymat = getPntrToArgument(0);
    unsigned base = 0;
    unsigned ncol = mymat->getNumberOfColumns();
    unsigned nmat = mymat->getRowLength(task_index);
    unsigned dloc = task_index*ncol;
    for(unsigned i=0; i<getNumberOfArguments()-1; ++i) {
      Value* myvec = getPntrToArgument(i+1);
      base += getPntrToArgument(i)->getNumberOfStoredValues();
      double val=0;
      for(unsigned j=0; j<nmat; ++j) {
        val += mymat->get( task_index*ncol + j, false )*myvec->get( mymat->getRowIndex( task_index, j ) );
      }
      myvals.setValue( i, val );

      // And the derivatives
      if( doNotCalculateDerivatives() ) {
        continue;
      }

      for(unsigned j=0; j<nmat; ++j) {
        unsigned kind = mymat->getRowIndex( task_index, j );
        double vecval = myvec->get( kind );
        double matval = mymat->get( task_index*ncol + j, false );
        myvals.addDerivative( i, dloc + j, vecval );
        myvals.updateIndex( i, dloc + j );
        myvals.addDerivative( i, base + kind, matval );
        myvals.updateIndex( i, base + kind );
      }
    }
  } else {
    unsigned base=0, n=getNumberOfArguments()-1;
    Value* myvec = getPntrToArgument(n);
    unsigned nmat_der = 0;
    for(unsigned i=0; i<n; ++i) {
      nmat_der += getPntrToArgument(i)->getNumberOfStoredValues();
    }
    for(unsigned i=0; i<n; ++i) {
      Value* mymat = getPntrToArgument(i);
      unsigned ncol = mymat->getNumberOfColumns();
      unsigned nmat = mymat->getRowLength(task_index);
      double val=0;
      for(unsigned j=0; j<nmat; ++j) {
        val += mymat->get( task_index*ncol + j, false )*myvec->get( mymat->getRowIndex( task_index, j ) );
      }
      myvals.setValue( i, val );

      // And the derivatives
      if( doNotCalculateDerivatives() ) {
        continue;
      }

      unsigned dloc = base + task_index*ncol;
      for(unsigned j=0; j<nmat; ++j) {
        unsigned kind = mymat->getRowIndex( task_index, j );
        double vecval = myvec->get( kind );
        double matval = mymat->get( task_index*ncol + j, false );
        myvals.addDerivative( i, dloc + j, vecval );
        myvals.updateIndex( i, dloc + j );
        myvals.addDerivative( i, nmat_der + kind, matval );
        myvals.updateIndex( i, nmat_der + kind );
      }
      base += mymat->getNumberOfStoredValues();
    }
  }
}

void MatrixTimesVector::getNumberOfForceDerivatives( unsigned& nforces, unsigned& nderiv ) const {
  ActionWithVector::getNumberOfForceDerivatives( nforces, nderiv );
  if( sumrows ) {
    nderiv = getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(getNumberOfArguments()-1)->getNumberOfStoredValues();
  }
}

void MatrixTimesVector::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( !sumrows ) {
    ActionWithVector::gatherForces( itask, myvals, forces );
    return;
  }
  if( checkComponentsForForce() ) {
    unsigned base = 0;
    for(unsigned ival=0; ival<getNumberOfComponents(); ++ival) {
      const Value* myval=getConstPntrToComponent(ival);
      if( myval->forcesWereAdded() ) {
        double fforce = myval->getForce(itask);
        for(unsigned j=0; j<myvals.getNumberActive(ival); ++j) {
          unsigned jder=myvals.getActiveIndex(ival, j);
          plumed_dbg_assert( jder<forces.size() );
          forces[base+jder] += fforce*myvals.getDerivative( ival, jder );
        }
      }
      base += getPntrToArgument(ival)->getNumberOfStoredValues();
    }
  }
}

}
}

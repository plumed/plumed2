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
#include "ActionWithMatrix.h"
#include "tools/Communicator.h"

namespace PLMD {

void ActionWithMatrix::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
}

ActionWithMatrix::ActionWithMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  diagzero(false) {

  if( keywords.exists("ELEMENTS_ON_DIAGONAL_ARE_ZERO") ) {
    parseFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",diagzero);
    if( diagzero ) {
      log.printf("  setting diagonal elements equal to zero\n");
    }
  }
}


class RequiredMatrixElementsUpdater {
  RequiredMatrixElements& outmat;
public:
  RequiredMatrixElementsUpdater( RequiredMatrixElements& outmat ) : outmat(outmat) {}
  ~RequiredMatrixElementsUpdater() {
    outmat.update();
  }
};

void ActionWithMatrix::updateBookeepingArrays( RequiredMatrixElements& outmat ) {
  RequiredMatrixElementsUpdater updater(outmat);
  Value* myval = getPntrToComponent(0);
  unsigned lstart = myval->getShape()[0];
  if( getNumberOfMasks()>0 ) {
    Value* maskarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( maskarg->getNumberOfColumns() );
    }
    for(unsigned i=0; i<lstart; ++i) {
      unsigned rstart = i*(1+myval->getNumberOfColumns());
      myval->setMatrixBookeepingElement( rstart, maskarg->getRowLength(i) );
      for(unsigned j=0; j<maskarg->getRowLength(i); ++j) {
        myval->setMatrixBookeepingElement( rstart + 1 + j, maskarg->getRowIndex(i, j) );
      }
    }
  } else if ( diagzero ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( myval->getShape()[1]-1 );
    }
    for(unsigned i=0; i<lstart; ++i) {
      unsigned k=0, rstart = i*(1+myval->getNumberOfColumns());
      myval->setMatrixBookeepingElement( rstart, myval->getShape()[1]-1 );
      for(unsigned j=0; j<myval->getShape()[1]; ++j) {
        if( i!=j ) {
          myval->setMatrixBookeepingElement( rstart + 1 + k, j );
          k++;
        }
      }
    }
  } else {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( myval->getShape()[1] );
    }
  }
  Value* mycomp = getPntrToComponent(0);
  outmat.ncols = mycomp->getNumberOfColumns();
  outmat.bookeeping.resize( mycomp->matrix_bookeeping.size() );
  for(unsigned i=0; i<outmat.bookeeping.size(); ++i) {
    outmat.bookeeping[i] = mycomp->matrix_bookeeping[i];
  }
  for(unsigned i=1; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( myval );
  }
}

void ActionWithMatrix::transferStashToValues( const std::vector<double>& stash ) {
  unsigned ncols = getPntrToComponent(0)->getNumberOfColumns();
  unsigned m=0, n=0, nrows = getPntrToComponent(0)->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    for(unsigned j=0; j<ncols; ++j) {
      for(unsigned k=0; k<getNumberOfComponents(); ++k) {
        getPntrToComponent(k)->set( m, stash[n] );
        n++;
      }
      m++;
    }
  }
}

void ActionWithMatrix::transferStashToValues( const std::vector<float>& stash ) {
  unsigned ncols = getPntrToComponent(0)->getNumberOfColumns();
  unsigned m=0, n=0, nrows = getPntrToComponent(0)->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    for(unsigned j=0; j<ncols; ++j) {
      for(unsigned k=0; k<getNumberOfComponents(); ++k) {
        getPntrToComponent(k)->set( m, stash[n] );
        n++;
      }
      m++;
    }
  }
}

void ActionWithMatrix::transferForcesToStash( std::vector<double>& stash ) const {
  unsigned ncols = getConstPntrToComponent(0)->getNumberOfColumns();
  unsigned m=0, n=0, nrows = getConstPntrToComponent(0)->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    for(unsigned j=0; j<ncols; ++j) {
      for(unsigned k=0; k<getNumberOfComponents(); ++k) {
        stash[n] = getConstPntrToComponent(k)->getForce( m );
        n++;
      }
      m++;
    }
  }
}

void ActionWithMatrix::transferForcesToStash( std::vector<float>& stash ) const {
  unsigned ncols = getConstPntrToComponent(0)->getNumberOfColumns();
  unsigned m=0, n=0, nrows = getConstPntrToComponent(0)->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    for(unsigned j=0; j<ncols; ++j) {
      for(unsigned k=0; k<getNumberOfComponents(); ++k) {
        stash[n] = getConstPntrToComponent(k)->getForce( m );
        n++;
      }
      m++;
    }
  }
}

}

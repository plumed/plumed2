/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR NEIGHBORS
/*
Build a matrix with ones in for the N nearest neighbours of an atom

The following input illustrates how to use this action in tandem with [DISTANCE_MATRIX](DISTANCE_MATRIX.md) to find the six
nearest atoms to each of the first 100 atoms in the input file:

```plumed
d1: DISTANCE_MATRIX GROUP=1-100
n: NEIGHBORS ARG=d1 NLOWEST=6
```

Alternatively, if you would like to use a [CONTACT_MATRIX](CONTACT_MATRIX.md) to do something similar you would do the following:

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.5}
n: NEIGHBORS ARG=c1 NHIGHEST=6
```

This command is useful for implementing alternatives to the symmatry functions that are defined by the shortcuts
in the module symfunc.  For example, suppose that you want to calculate a variant on the [TETRAHEDRAL](TETRAHEDRAL.md) symmetry function.
In this variant on the CV the coordination sphere around each central atom is not defined using a switching function.  Instad
this coordination sphere contains only the four nearest atoms.  You can implement this CV by using the following input:

```plumed
d1: DISTANCE_MATRIX GROUP=1-100 COMPONENTS
n: NEIGHBORS ARG=d1.w NLOWEST=4
f: CUSTOM ARG=n,d1.x,d1.y,d1.z,d1.w VAR=w,x,y,z,r PERIODIC=NO FUNC=w*(((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3)
ones: ONES SIZE=100
ucv: MATRIX_VECTOR_PRODUCT ARG=f,ones
cv: CUSTOM ARG=ucv PERIODIC=NO FUNC=x/4
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class Neighbors : public ActionWithMatrix {
  bool lowest;
  unsigned number;
public:
  static void registerKeywords( Keywords& keys );
  explicit Neighbors(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  unsigned getNumberOfColumns() const override {
    return number;
  }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override {}
  void turnOnDerivatives() override ;
};

PLUMED_REGISTER_ACTION(Neighbors,"NEIGHBORS")

void Neighbors::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix","the label of an adjacency/distance matrix that will be used to find the nearest neighbors");
  keys.add("compulsory","NLOWEST","0","in each row of the output matrix set the elements that correspond to the n lowest elements in each row of the input matrix equal to one");
  keys.add("compulsory","NHIGHEST","0","in each row of the output matrix set the elements that correspond to the n highest elements in each row of the input matrix equal to one");
  keys.setValueDescription("matrix","a matrix in which the ij element is one if the ij-element of the input matrix is one of the NLOWEST/NHIGHEST elements on that row of the input matrix and zero otherwise");
}

Neighbors::Neighbors(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("found wrong number of arguments in input");
  }
  if( getPntrToArgument(0)->getRank()!=2 ) {
    error("input argument should be a matrix");
  }

  unsigned nlow;
  parse("NLOWEST",nlow);
  unsigned nhigh;
  parse("NHIGHEST",nhigh);
  if( nlow==0 && nhigh==0 ) {
    error("missing NLOWEST or NHIGHEST keyword one of these two keywords must be set in input");
  }
  if( nlow>0 && nhigh>0 ) {
    error("should only be one of NLOWEST or NHIGHEST set in input");
  }
  if( nlow>0 ) {
    number=nlow;
    lowest=true;
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d lowest elements in each row of the input matrix\n",number);
  }
  if( nhigh>0 ) {
    number=nhigh;
    lowest=false;
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d highest elements in each row of the input matrix\n",number);
  }

  // And get the shape
  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
  addValue( shape );
  setNotPeriodic();
  checkRead();
}

void Neighbors::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  warning("think about whether your symmetry functions are continuous. If the symmetry function can be calculated from distances only then you can use NEIGHBORS. If you calculate angles between vectors or use the vectors directly then the symmetry function computed using NEIGHBORS is not continuous.  It does not make sense to use such CVs when biasing");
}

unsigned Neighbors::getNumberOfDerivatives() {
  return 0;
}

void Neighbors::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  const Value* wval = getPntrToArgument(0);
  unsigned nbonds = wval->getRowLength( task_index ), ncols = wval->getShape()[1];
  if( indices.size()!=1+number ) {
    indices.resize( 1 + number );
  }
  myvals.setSplitIndex(1+number);

  unsigned nind=0;
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }
    nind++;
  }
  if( number>nind ) {
    plumed_merror("not enough matrix elements were stored");
  }

  // Now build vectors for doing sorting
  std::vector<std::pair<double,unsigned> > rows( nind );
  unsigned n=0;
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned iind = wval->getRowIndex( task_index, i );
    unsigned ipos = ncols*task_index + iind;
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }
    rows[n].first=weighti;
    rows[n].second=iind;
    n++;
  }

  // Now do the sort and clear all the stored values ready for recompute
  std::sort( rows.begin(), rows.end() );
  // This is to make this action consistent with what in other matrix actions
  unsigned start_n = getPntrToArgument(0)->getShape()[0];
  // And setup the lowest indices, which are the ones we need to calculate
  for(unsigned i=0; i<number; ++i) {
    indices[i+1] = start_n + rows[nind-1-i].second;
    if( lowest ) {
      indices[i+1] = start_n + rows[i].second;
    }
  }
}

void Neighbors::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  myvals.addValue( 0, 1.0 );
}

}
}

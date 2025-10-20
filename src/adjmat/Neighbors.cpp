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
#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "core/ActionShortcut.h"
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

You can even use these ideas with the functionality that is in the [volumes module](module_volumes.md) as shown below:

```plumed
# The atoms that are of interest
ow: GROUP ATOMS=1-16500
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=ow CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# The distance matrix
dmap: DISTANCE_MATRIX COMPONENTS GROUP=ow CUTOFF=1.0 MASK=sphere
# Find the four nearest neighbors
acv_neigh: NEIGHBORS ARG=dmap.w NLOWEST=4 MASK=sphere
# Compute a function for the atoms that are in the first coordination sphere
acv_g8: GSYMFUNC_THREEBODY ...
  WEIGHT=acv_neigh ARG=dmap.x,dmap.y,dmap.z
  FUNCTION1={FUNC=(cos(ajik)+1/3)^2 LABEL=g8}
  MASK=sphere
...
# Now compute the value of the function above for those atoms that are in the
# sphere of interest
acv: CUSTOM ARG=acv_g8.g8,sphere FUNC=y*(1-(3*x/8)) PERIODIC=NO
# And now compute the final average
acv_sum: SUM ARG=acv PERIODIC=NO
acv_norm: SUM ARG=sphere PERIODIC=NO
mean: CUSTOM ARG=acv_sum,acv_norm FUNC=x/y PERIODIC=NO
PRINT ARG=mean FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class NeighborsShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit NeighborsShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(NeighborsShortcut,"NEIGHBORS")

void NeighborsShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("optional","MASK","vector","a vector that is used to used to determine which rows of the neighbors matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","the label of an adjacency/distance matrix that will be used to find the nearest neighbors");
  keys.add("compulsory","NLOWEST","0","in each row of the output matrix set the elements that correspond to the n lowest elements in each row of the input matrix equal to one");
  keys.add("compulsory","NHIGHEST","0","in each row of the output matrix set the elements that correspond to the n highest elements in each row of the input matrix equal to one");
  keys.setValueDescription("matrix","a matrix in which the ij element is one if the ij-element of the input matrix is one of the NLOWEST/NHIGHEST elements on that row of the input matrix and zero otherwise");
  keys.addActionNameSuffix("_1LOW");
  keys.addActionNameSuffix("_NLOW");
  keys.addActionNameSuffix("_1HIGH");
  keys.addActionNameSuffix("_NHIGH");
}

NeighborsShortcut::NeighborsShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
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
  if( nlow==1 ) {
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d lowest elements in each row of the input matrix\n",nlow);
    readInputLine( getShortcutLabel() + ": NEIGHBORS_1LOW N=1 " + convertInputLineToString() );
  } else if( nlow>1 ) {
    std::string str_n;
    Tools::convert( nlow, str_n );
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d lowest elements in each row of the input matrix\n",nlow);
    readInputLine( getShortcutLabel() + ": NEIGHBORS_NLOW N=" + str_n + " " + convertInputLineToString() );
  } else if( nhigh==1 ) {
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d highest elements in each row of the input matrix\n",nhigh);
    readInputLine( getShortcutLabel() + ": NEIGHBORS_1HIGH N=1 " + convertInputLineToString() );
  } else if( nhigh>1 ) {
    std::string str_n;
    Tools::convert( nlow, str_n );
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d highest elements in each row of the input matrix\n",nhigh);
    readInputLine( getShortcutLabel() + ": NEIGHBORS_NHIGH N=" + str_n + " " + convertInputLineToString() );
  } else {
    error("do not know what I am supposed to do");
  }
}

class NeighborCalcInput {
public:
  unsigned number;
  unsigned nind;
  unsigned ncols;
  unsigned nbonds;
  View<const double> matrow;
  View<const std::size_t> bookrow;
  NeighborCalcInput( unsigned task_index,
                     unsigned n,
                     const ArgumentBookeepingHolder& arg,
                     double* d ):
    number(n),
    nind(0),
    ncols(arg.ncols),
    nbonds(arg.bookeeping[task_index*(1+ncols)]),
    matrow( d + task_index*ncols, nbonds ),
    bookrow( arg.bookeeping.data() + task_index*(1+ncols)+1, nbonds ) {
    for(unsigned i=0; i<nbonds; ++i) {
      if( matrow[i]<epsilon ) {
        continue ;
      }
      nind++;
    }
    if( number>nind ) {
      plumed_merror("not enough matrix elements were stored");
    }
  }
  void getSortedData( std::vector<std::pair<double,unsigned> >& rows ) const ;
};

void NeighborCalcInput::getSortedData( std::vector<std::pair<double,unsigned> >& rows ) const {
  unsigned n=0;
  for(unsigned i=0; i<nbonds; ++i) {
    if( matrow[i]<epsilon ) {
      continue ;
    }
    rows[n].first=matrow[i];
    rows[n].second=bookrow[i];
    n++;
  }
  std::sort( rows.begin(), rows.end() );
}

class OneLowInput {
public:
  unsigned number;
  static void calculate( const NeighborCalcInput& input, View<double>& output );
};

void OneLowInput::calculate( const NeighborCalcInput& input, View<double>& output ) {
  unsigned nv = input.bookrow[0];
  double min = input.matrow[0];
  for(unsigned i=1; i<input.nbonds; ++i) {
    if( input.matrow[i]<min ) {
      min = input.matrow[i];
      nv = input.bookrow[i];
    }
  }
  output[0] = nv;
}

class NLowInput {
public:
  unsigned number;
  static void calculate( const NeighborCalcInput& input, View<double>& output );
};

void NLowInput::calculate( const NeighborCalcInput& input, View<double>& output ) {
  std::vector<std::pair<double,unsigned> > rows( input.nind );
  input.getSortedData( rows );

  for(unsigned i=0; i<input.number; ++i) {
    output[i] = rows[i].second;
  }
}

class OneHighInput {
public:
  unsigned number;
  static void calculate( const NeighborCalcInput& input, View<double>& output );
};

void OneHighInput::calculate( const NeighborCalcInput& input, View<double>& output ) {
  unsigned nv = input.bookrow[0];
  double max = input.matrow[0];
  for(unsigned i=1; i<input.nbonds; ++i) {
    if( input.matrow[i]>max ) {
      max = input.matrow[i];
      nv = input.bookrow[i];
    }
  }
  output[0] = nv;
}

class NHighInput {
public:
  unsigned number;
  static void calculate( const NeighborCalcInput& input, View<double>& output );
};

void NHighInput::calculate( const NeighborCalcInput& input, View<double>& output ) {
  std::vector<std::pair<double,unsigned> > rows( input.nind );
  input.getSortedData( rows );

  for(unsigned i=0; i<input.number; ++i) {
    output[i] = rows[input.nind-1-i].second;
  }
}

template <class T>
class Neighbors : public ActionWithVector {
public:
  using input_type = T;
  using PTM = ParallelTaskManager<Neighbors<T>>;
private:
/// The parallel task manager
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit Neighbors(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  void turnOnDerivatives() override ;
  void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) override ;
  void prepare() override ;
  void calculate() override ;
  static void performTask( std::size_t task_index,
                           const T& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
};

typedef Neighbors<OneLowInput> olinp;
PLUMED_REGISTER_ACTION(olinp,"NEIGHBORS_1LOW")
typedef Neighbors<NLowInput> nlinp;
PLUMED_REGISTER_ACTION(nlinp,"NEIGHBORS_NLOW")
typedef Neighbors<OneHighInput> ohinp;
PLUMED_REGISTER_ACTION(ohinp,"NEIGHBORS_1HIGH")
typedef Neighbors<NHighInput> nhinp;
PLUMED_REGISTER_ACTION(nhinp,"NEIGHBORS_NHIGH")

template <class T>
void Neighbors<T>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.setDisplayName("NEIGHBORS");
  keys.addInputKeyword("optional","MASK","vector","a vector that is used to used to determine which rows of the neighbors matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","the label of an adjacency/distance matrix that will be used to find the nearest neighbors");
  keys.add("compulsory","N","the number of non-zero elements in each row of the output matrix");
  keys.setValueDescription("matrix","a matrix in which the ij element is one if the ij-element of the input matrix is one of the NLOWEST/NHIGHEST elements on that row of the input matrix and zero otherwise");
  PTM::registerKeywords( keys );
}

template <class T>
Neighbors<T>::Neighbors(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  unsigned nargs = getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    nargs = nargs - getNumberOfMasks();
  }
  if( nargs!=1 ) {
    error("found wrong number of arguments in input");
  }
  if( getPntrToArgument(0)->getRank()!=2 ) {
    error("input argument should be a matrix");
  }

  T myinp;
  parse("N",myinp.number);
  // And get the shape
  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->reshapeMatrixStore( myinp.number );
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  // Setup the parallel task manager
  taskmanager.setupParallelTaskManager( 0, 0 );
  taskmanager.setActionInput( myinp );
}

template <class T>
void Neighbors<T>::turnOnDerivatives() {
  error("If the symmetry function can be calculated from distances only then the derivatives of a function that is computed with NEIGHBORS will be continuous. If you calculate angles between vectors or use the vectors directly then the symmetry function computed using NEIGHBORS is not continuous. Out of an abundance of caution we thus forbid forces for this action.  If you are interested on applying forces with this action email gareth.tribello@gmail.com");
}

template <class T>
unsigned Neighbors<T>::getNumberOfDerivatives() {
  return 0;
}

template <class T>
void Neighbors<T>::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0]
      && myval->getShape()[1]==getPntrToArgument(0)->getShape()[1] ) {
    return;
  }
  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
  myval->setShape(shape);
  myval->reshapeMatrixStore( taskmanager.getActionInput().number );
}

template <class T>
void Neighbors<T>::calculate() {
  taskmanager.runAllTasks();
}

template <class T>
void Neighbors<T>::performTask( std::size_t task_index,
                                const T& actiondata,
                                ParallelActionsInput& input,
                                ParallelActionsOutput& output ) {
  T::calculate( NeighborCalcInput( task_index,
                                   actiondata.number,
                                   ArgumentBookeepingHolder::create( 0, input ),
                                   input.inputdata ),
                output.values );
}

template <class T>
void Neighbors<T>::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) {
  Value* myval = getPntrToComponent(0);
  // All values are set equal to one
  for(unsigned i=0; i<myval->getNumberOfStoredValues(); ++i) {
    myval->set( i, 1 );
  }

  // And we set the bookeeping data from the data in the stash
  unsigned k=0;
  std::vector<std::size_t> rowindices( taskmanager.getActionInput().number );
  for(unsigned i=0; i<myval->getShape()[0]; ++i) {
    for(unsigned j=0; j<rowindices.size(); ++j) {
      rowindices[j] = static_cast<std::size_t>( stash[k] );
      ++k;
    }
    myval->setRowIndices( i, rowindices );
  }
}

}
}

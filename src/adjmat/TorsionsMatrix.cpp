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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"
#include "core/ParallelTaskManager.h"
#include "tools/Torsion.h"

//+PLUMEDOC MCOLVAR TORSIONS_MATRIX
/*
Calculate the matrix of torsions between two vectors of molecules

This action was implemented to ensure that we can calculate the [SMAC](SMAC.md) collective variable that is discussed in
[this paper](https://www.sciencedirect.com/science/article/abs/pii/S0009250914004503?via%3Dihub). This particular action
tracks the relative orientations for all the pairs of molecules in a set much like the variables described in the crystdistrib module.

The orientations of molecules can be specified using either [PLANE](PLANE.md) or [DISTANCE](DISTANCE.md).  The example below shows how you can use
internal vectors connecting two atoms in the molecules to define the orientation of that molecule.  Three of these internal
vectors are calculated using a DISTANCE command in the input below.  The matrix of torsional angles between these various
vectors is then computed:

```plumed
d1: DISTANCE ATOMS1=1,5 ATOMS2=11,15 ATOMS3=21,25 COMPONENTS
s: VSTACK ARG=d1.x,d1.y,d1.z
sT: TRANSPOSE ARG=s
m: TORSIONS_MATRIX ARG=s,sT POSITIONS1=1,11,21 POSITIONS2=1,11,21
PRINT ARG=m FILE=matrix
```

In this example, the torsional angle in element $(1,2)$ of the matrix with label `m` is the angle between the plane containing atoms 1,5 and 10 and the plane
connecting atoms 1,10 and 15.  In other words, the elements in this matrix are the torsional angles between the vectors in the input matrices
around the vector connecting the corresponding atomic positions that are specified using the `POSTIONS` keyword.

You can also calculate a matrix of torsional angles between two different groups of molecules by using an input like the one below:

```plumed
pA: PLANE ATOMS1=1,2,3 ATOMS2=11,12,13
sA: VSTACK ARG=pA.x,pA.y,pA.z
pB: PLANE ATOMS1=21,22,23 ATOMS2=31,32,33 ATOMS3=41,42,43
sB: VSTACK ARG=pB.x,pB.y,pB.z
sBT: TRANSPOSE ARG=sB
m: TORSIONS_MATRIX ARG=sA,sBT POSITIONS1=1,11 POSITIONS2=21,31,41
PRINT ARG=m FILE=matrix
```

In this example, the orientations of the molecules are specified using the [PLANE](PLANE.md) action and is given by a normal to the plane containing the three atoms from the molecule
that was specified.  The final output is $2 \times 3$ matrix that contains all the torsional angles between the molecules defined by the two PLANE actions.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class TorsionsMatrixInput {
public:
  RequiredMatrixElements outmat;
};

class TorsionsMatrix : public ActionWithMatrix {
public:
  using input_type = TorsionsMatrixInput;
  using PTM = ParallelTaskManager<TorsionsMatrix>;
private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit TorsionsMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void prepare() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  static void performTask( std::size_t task_index, const TorsionsMatrixInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index, const TorsionsMatrixInput& actiondata, const ParallelActionsInput& input, const ForceInput& fdata, ForceOutput forces );
};

PLUMED_REGISTER_ACTION(TorsionsMatrix,"TORSIONS_MATRIX")

void TorsionsMatrix::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of this matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","an Nx3 and a 3xN matrix that contain the bond vectors that you would like to determine the torsion angles between");
  keys.add("atoms","POSITIONS1","the positions to use for the molecules specified using the first argument");
  keys.add("atoms","POSITIONS2","the positions to use for the molecules specified using the second argument");
  PTM::registerKeywords( keys );
  keys.setValueDescription("matrix","the matrix of torsions between the two vectors of input directors");
}

TorsionsMatrix::TorsionsMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  taskmanager(this) {
  int nm=getNumberOfMasks();
  if( nm<0 ) {
    nm = 0;
  }
  if( getNumberOfArguments()-nm!=2 ) {
    error("should be two arguments to this action, a matrix and a vector");
  }
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
    error("first argument to this action should be a matrix");
  }
  if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) {
    error("second argument to this action should be a matrix");
  }
  if( getPntrToArgument(0)->getShape()[1]!=3 || getPntrToArgument(1)->getShape()[0]!=3 ) {
    error("number of columns in first matrix and number of rows in second matrix should equal 3");
  }

  std::vector<AtomNumber> atoms_a;
  parseAtomList("POSITIONS1", atoms_a );
  if( atoms_a.size()!=getPntrToArgument(0)->getShape()[0] ) {
    error("mismatch between number of atoms specified using POSITIONS1 and number of arguments in vector input");
  }
  log.printf("  using positions of these atoms for vectors in first matrix \n");
  for(unsigned int i=0; i<atoms_a.size(); ++i) {
    if ( (i+1) % 25 == 0 ) {
      log.printf("  \n");
    }
    log.printf("  %d", atoms_a[i].serial());
  }
  log.printf("\n");
  std::vector<AtomNumber> atoms_b;
  parseAtomList("POSITIONS2", atoms_b );
  if( atoms_b.size()!=getPntrToArgument(1)->getShape()[1] ) {
    error("mismatch between number of atoms specified using POSITIONS2 and number of arguments in vector input");
  }
  log.printf("  using positions of these atoms for vectors in second matrix \n");
  for(unsigned i=0; i<atoms_b.size(); ++i) {
    if ( (i+1) % 25 == 0 ) {
      log.printf("  \n");
    }
    log.printf("  %d", atoms_b[i].serial());
    atoms_a.push_back( atoms_b[i] );
  }
  log.printf("\n");
  requestAtoms( atoms_a, false );

  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape );
  setPeriodic("-pi","pi");

  if( nm>0 ) {
    unsigned iarg = getNumberOfArguments()-1;
    if( getPntrToArgument(iarg)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
      error("argument passed to MASK keyword should be a matrix");
    }
    if( getPntrToArgument(iarg)->getShape()[0]!=shape[0] || getPntrToArgument(iarg)->getShape()[1]!=shape[1] ) {
      error("argument passed to MASK keyword has the wrong shape");
    }
  }
  taskmanager.setupParallelTaskManager( 1, 2, getPntrToArgument(0)->getNumberOfStoredValues() );
  taskmanager.setActionInput( TorsionsMatrixInput() );
}

unsigned TorsionsMatrix::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms() + 9 + getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues();
}

void TorsionsMatrix::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[1] ) {
    return;
  }
  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  myval->setShape(shape);
  myval->reshapeMatrixStore( shape[1] );
}

void TorsionsMatrix::calculate() {
  updateBookeepingArrays( taskmanager.getActionInput().outmat );
  unsigned nvals = getPntrToComponent(0)->getNumberOfColumns();
  taskmanager.setupParallelTaskManager( 1, 21*nvals*getPntrToArgument(0)->getNumberOfColumns(), getPntrToArgument(1)->getNumberOfStoredValues()+3*getNumberOfAtoms()+9, nvals );
  taskmanager.runAllTasks();
}

void TorsionsMatrix::getInputData( std::vector<double>& inputdata ) const {
  std::size_t total_data = getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues() + 3*getNumberOfAtoms();

  if( inputdata.size()!=total_data ) {
    inputdata.resize( total_data );
  }

  total_data = 0;
  Value* myarg = getPntrToArgument(0);
  for(unsigned j=0; j<myarg->getNumberOfStoredValues(); ++j) {
    inputdata[total_data] = myarg->get(j,false);
    total_data++;
  }
  myarg = getPntrToArgument(1);
  for(unsigned j=0; j<myarg->getNumberOfStoredValues(); ++j) {
    inputdata[total_data] = myarg->get(j,false);
    total_data++;
  }
  for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
    Vector pos( getPosition(j) );
    inputdata[total_data+0] = pos[0];
    inputdata[total_data+1] = pos[1];
    inputdata[total_data+2] = pos[2];
    total_data += 3;
  }

}

void TorsionsMatrix::performTask( std::size_t task_index,
                                  const TorsionsMatrixInput& actiondata,
                                  ParallelActionsInput& input,
                                  ParallelActionsOutput& output ) {
  std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
  std::size_t nelements = actiondata.outmat.bookeeping[fstart];
  ArgumentBookeepingHolder arg0( 0, input ), arg1( 1, input );
  std::size_t nargdata = arg1.start + arg1.shape[0]*arg1.ncols;

  // Get the position and orientation for the first molecule
  std::size_t atbase = nargdata + 3*task_index;
  Vector atom1( input.inputdata[atbase], input.inputdata[atbase+1], input.inputdata[atbase+2] );
  std::size_t agbase = arg0.ncols*task_index;
  Vector v1( input.inputdata[agbase], input.inputdata[agbase+1], input.inputdata[agbase+2] );

  // Get the distances to all the molecules in the coordination sphere
  std::vector<Vector> atom2( nelements );
  for(unsigned i=0; i<nelements; ++i) {
    std::size_t at2base = nargdata + 3*arg0.shape[0] + 3*actiondata.outmat.bookeeping[fstart+1+i];
    atom2[i] = Vector( input.inputdata[at2base], input.inputdata[at2base+1], input.inputdata[at2base+2] ) - atom1;
  }
  input.pbc->apply( atom2, nelements );

  // Now compute all the torsions
  Torsion t;
  Vector dv1, dconn, dv2 ;
  for(unsigned i=0; i<nelements; ++i) {
    if( atom2[i].modulo2()<epsilon ) {
      if( !input.noderiv ) {
        View<double, 21> der( output.derivatives.data() + 21*i );
        for(unsigned j=0; j<21; ++j) {
          der[j] = 0;
        }
      }
      continue ;
    }

    std::size_t ag2base = arg1.start + actiondata.outmat.bookeeping[fstart+1+i];
    Vector v2(input.inputdata[ag2base], input.inputdata[ag2base+arg1.ncols], input.inputdata[ag2base+2*arg1.ncols] );
    output.values[i] = t.compute( v1, atom2[i], v2, dv1, dconn, dv2 );

    if( input.noderiv ) {
      continue ;
    }

    std::size_t base = + 21*i;
    View<double, 3> deriv1( output.derivatives.data() + base );
    deriv1 = dv1;
    View<double, 3> deriv2( output.derivatives.data() + base + 3 );
    deriv2 = dv2;
    View<double, 3> deriv3( output.derivatives.data() + base + 6 );
    deriv3 = -dconn;
    View<double, 3> deriv4( output.derivatives.data() + base + 9 );
    deriv4 = dconn;
    Tensor vir( -extProduct( atom2[i], dconn ) );
    View2D<double,3,3> virial( output.derivatives.data() + base + 12 );
    for(unsigned j=0; j<3; ++j) {
      for(unsigned k=0; k<3; ++k) {
        virial[j][k] = vir[j][k];
      }
    }
  }

}

void TorsionsMatrix::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

void TorsionsMatrix::gatherForces( std::size_t task_index,
                                   const TorsionsMatrixInput& actiondata,
                                   const ParallelActionsInput& input,
                                   const ForceInput& fdata,
                                   ForceOutput forces ) {
  std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
  std::size_t nelements = actiondata.outmat.bookeeping[fstart];
  ArgumentBookeepingHolder arg0( 0, input ), arg1( 1, input );
  std::size_t arg1start = task_index*arg0.ncols;
  std::size_t atom1start = arg1.shape[0]*arg1.ncols + 3*task_index;
  std::size_t b2astart = arg1.shape[0]*arg1.ncols + 3*arg0.shape[0];
  for(unsigned i=0; i<nelements; ++i) {
    std::size_t base = 21*i;
    double force = fdata.force[i];

    View<double, 3> deriv1( fdata.deriv.data() + base );
    forces.thread_unsafe[ arg1start ] += force*deriv1[0];
    forces.thread_unsafe[ arg1start + 1 ] += force*deriv1[1];
    forces.thread_unsafe[ arg1start + 2 ] += force*deriv1[2];

    View<double, 3> deriv2( fdata.deriv.data() + base + 3 );
    forces.thread_safe[ actiondata.outmat.bookeeping[fstart+1+i] ] += force*deriv2[0];
    forces.thread_safe[ actiondata.outmat.bookeeping[fstart+1+i] + arg1.ncols ] += force*deriv2[1];
    forces.thread_safe[ actiondata.outmat.bookeeping[fstart+1+i] + 2*arg1.ncols ] += force*deriv2[2];

    View<double, 3> deriv3( fdata.deriv.data() + base + 6 );
    forces.thread_safe[ atom1start ] += force*deriv3[0];
    forces.thread_safe[ atom1start + 1 ] += force*deriv3[1];
    forces.thread_safe[ atom1start + 2 ] += force*deriv3[2];

    View<double, 3> deriv4( fdata.deriv.data() + base + 9 );
    std::size_t atom2start = b2astart + 3*actiondata.outmat.bookeeping[fstart+1+i];
    forces.thread_safe[ atom2start ] += force*deriv4[0];
    forces.thread_safe[ atom2start + 1 ] += force*deriv4[1];
    forces.thread_safe[ atom2start + 2 ] += force*deriv4[2];

    std::size_t n=0;
    View<double, 9> vir( fdata.deriv.data() + base + 12 );
    for(unsigned j=forces.thread_safe.size()-9; j<forces.thread_safe.size(); ++j) {
      forces.thread_safe[j] += force*vir[n];
      ++n;
    }
  }
}

}
}

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
#include "core/ActionRegister.h"
#include "tools/LeptonCall.h"
#include "tools/Angle.h"

namespace PLMD {
namespace symfunc {

//+PLUMEDOC COLVAR GSYMFUNC_THREEBODY
/*
Calculate functions of the coordinates of the coordinates of all pairs of bonds in the first coordination sphere of an atom

This shortcut can be used to calculate [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html) that
are like those defined by Behler in the paper that is cited in the bibliography below. The particular symmetry functions that are computed
by this shortcut are the angular ones that are functions of the set pairs of atoms in the coordination sphere of the central atom.  One of
the angular symmetry functions that Behler introduces is:

$$
G^5_i = 2^{1-\zeta} \sum_{j,k\ne i} (1 + \lambda\cos\theta_{ijk})^\zeta e^{-\nu(R_{ij}^2 + R_{ik}^2)} f_c(R_{ij}) f_c(R_{ik})
$$

In this expression $\zeta$, $\nu$ and $\lambda$ are all parameters.  $f_c$ is a switching function which acts upon $R_{ij}$, the distance between atom $i$ and atom $j$, and
$R_{ik}$, the distance between atom $i$ and atom $k$.  $\theta_{ijk}$ is then the angle between the vector that points from atom $i$ to atom $j$ and the vector that points from
atom $i$ to atom $k$.  THe input below can be used to get PLUMED to calculate the 100 values for this symmetry function for the 100 atoms in a system.

```plumed
# Calculate the contact matrix and the x,y and z components of the bond vectors
# This action calculates 4 100x100 matrices
cmat: CONTACT_MATRIX GROUP=1-100 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS

# Compute the symmetry function for the 100 atoms from the 4 100x100 matrices output
# by cmat.  The output from this action is a vector with 100 elements
beh3: GSYMFUNC_THREEBODY ...
    WEIGHT=cmat.w ARG=cmat.x,cmat.y,cmat.z
    FUNCTION1={FUNC=0.25*exp(-0.1*(rij+rik))*(1+3*cos(ajik))^3 LABEL=g5}
...

# Print the 100 symmetry function values to a file
PRINT ARG=beh3.g5 FILE=colvar
```

The GSYMFUNC_THREEBODY action sums over all the distinct triples of atoms that are identified in the contact matrix.  This action uses the same functionality as [CUSTOM](CUSTOM.md) and can thus compute any
function of the following four quantities:

* `rij` - the distance between atom $i$ and atom $j$
* `rik` - the distance between atom $i$ and atom $k$
* `rjk` - the distance between atom $j$ and atom $k$
* `ajik` - the angle between the vector connecting atom $i$ to atom $j$ and the vector connecting atom $i$ to atom $k$.

Furthermore we can calculate more than one function of these four quantities at a time as illustrated by the input below:

```plumed
# Calculate the contact matrix and the x,y and z components of the bond vectors
# This action calculates 4 100x100 matrices
cmat: CONTACT_MATRIX GROUP=1-100 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS

# Compute the 4 symmetry function below for the 100 atoms from the 4 100x100 matrices output
# by cmat.  The output from this action is a vector with 100 elements
beh3: GSYMFUNC_THREEBODY ...
    WEIGHT=cmat.w ARG=cmat.x,cmat.y,cmat.z
    FUNCTION1={FUNC=0.25*(cos(pi*sqrt(rjk)/4.5)+1)*exp(-0.1*(rij+rik+rjk))*(1+2*cos(ajik))^2 LABEL=g4}
    FUNCTION2={FUNC=0.25*exp(-0.1*(rij+rik))*(1+3.5*cos(ajik))^3 LABEL=g5}
    FUNCTION3={FUNC=0.125*(1+6.6*cos(ajik))^4 LABEL=g6}
    FUNCTION4={FUNC=sin(3.0*(ajik-1)) LABEL=g7}
...

# Print the 4 sets of 100 symmetry function values to a file
PRINT ARG=beh3.g4,beh3.g5,beh3.g6,beh3.g7 FILE=colvar
```

You can even use this action in tandem with the features that are in the [volumes module](module_volumes.md) as shown below:

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

You can read more about how to calculate more Behler-type symmetry functions [here](https://www.plumed-tutorials.org/lessons/23/001/data/Behler.html).

*/
//+ENDPLUMEDOC

class ThreeBodyGFunctionsInput {
public:
  bool multi_action_input;
  std::vector<std::string> funcstr;
  std::vector<LeptonCall> functions;
  ThreeBodyGFunctionsInput& operator=( const ThreeBodyGFunctionsInput& m ) {
    multi_action_input = m.multi_action_input;
    for(unsigned i=0; i<m.funcstr.size(); ++i) {
      addFunction( m.funcstr[i], NULL );
    }
    return *this;
  }
  void addFunction( std::string myfunc, ActionWithVector* action ) {
    funcstr.push_back( myfunc );
    functions.push_back( LeptonCall() );
    std::vector<std::string> argnames(1);
    argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) {
      argnames.push_back("rij");
    }
    if( myfunc.find("rik")!=std::string::npos ) {
      if( action && argnames.size()<2 ) {
        action->error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rik");
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
      if( action && argnames.size()<2 ) {
        action->error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rjk");
    }
    functions[functions.size()-1].set( myfunc, argnames, action, true );
  }
};

class ThreeBodyGFunctions : public ActionWithVector {
public:
  using input_type = ThreeBodyGFunctionsInput;
  using PTM = ParallelTaskManager<ThreeBodyGFunctions>;
private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  unsigned getNumberOfDerivatives() override;
  static std::size_t getIndex( std::size_t irow, std::size_t jcol, const ArgumentBookeepingHolder& mat );
  static void performTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const ThreeBodyGFunctionsInput& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addInputKeyword("optional","MASK","vector","a vector that is used to used to determine which symmetry functions should be calculated");
  keys.addInputKeyword("compulsory","ARG","matrix","three matrices containing the bond vectors of interest");
  keys.addInputKeyword("compulsory","WEIGHT","matrix","the matrix that contains the weights that should be used for each connection");
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  PTM::registerKeywords( keys );
  ActionWithValue::useCustomisableComponents( keys );
  keys.addDOI("10.1063/1.3553717");
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  unsigned nargs = getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    nargs = nargs - getNumberOfMasks();
  }
  if( nargs!=3 ) {
    error("found wrong number of arguments in input");
  }
  std::vector<Value*> wval;
  parseArgumentList("WEIGHT",wval);
  if( wval.size()!=1 ) {
    error("keyword WEIGHT should be provided with the label of a single action");
  }

  for(unsigned i=0; i<3; ++i) {
    if( getPntrToArgument(i)->getRank()!=2 ) {
      error("input argument should be a matrix");
    }
    if( wval[0]->getShape()[0]!=getPntrToArgument(i)->getShape()[0] || wval[0]->getShape()[1]!=getPntrToArgument(i)->getShape()[1] ) {
      error("mismatched shapes of matrices in input");
    }
  }
  log.printf("  using bond weights from matrix labelled %s \n",wval[0]->getName().c_str() );
  // Rerequest the arguments
  std::vector<Value*> myargs( getArguments() );
  myargs.push_back( wval[0] );
  requestArguments( myargs );
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];

  // And now read the functions to compute
  ThreeBodyGFunctionsInput input;
  for(int i=1;; i++) {
    std::string myfunc, mystr, lab, num;
    Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) {
      break;
    }
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) {
      error("found no LABEL in FUNCTION" + num + " specification");
    }
    addComponent( lab, shape );
    componentIsNotPeriodic( lab );
    if( !Tools::parse(data,"FUNC",myfunc) ) {
      error("found no FUNC in FUNCTION" + num + " specification");
    }
    log.printf("  component labelled %s is computed using %s \n",lab.c_str(), myfunc.c_str() );
    input.addFunction( myfunc, this );
  }
  checkRead();
  input.multi_action_input = getPntrToArgument(3)->getPntrToAction()!=getPntrToArgument(0)->getPntrToAction();
  taskmanager.setupParallelTaskManager( 4*wval[0]->getShape()[1], 0 );
  taskmanager.setActionInput( input );
}

std::string ThreeBodyGFunctions::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getName().find(cname)!=std::string::npos ) {
      std::string num;
      Tools::convert( i+1, num );
      return "the function defined by the FUNCTION" + num + " keyword";
    }
  }
  plumed_error();
  return "";
}

unsigned ThreeBodyGFunctions::getNumberOfDerivatives() {
  return 0;
}

void ThreeBodyGFunctions::calculate() {
  taskmanager.runAllTasks();
}

std::size_t ThreeBodyGFunctions::getIndex( std::size_t irow, std::size_t jcol, const ArgumentBookeepingHolder& mat ) {
  if( mat.shape[1]==mat.ncols ) {
    return irow*mat.ncols + jcol;
  }
  for(unsigned i=0; i<mat.bookeeping[(1+mat.ncols)*irow]; ++i) {
    if( mat.bookeeping[(1+mat.ncols)*irow+1+i]==jcol ) {
      return irow*mat.ncols+i;
    }
  }
  plumed_merror("could not find index");
  return 0;
}

void ThreeBodyGFunctions::performTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  const double* xpntr=NULL;
  const double* ypntr=NULL;
  const double* zpntr=NULL;
  /// This function uses lepton so it is unlikely that we will GPU it.  that is why I have allowed myself to create vectors here
  std::vector<double> xvals, yvals, zvals;
  auto arg0 = ArgumentBookeepingHolder::create( 0, input );
  auto arg1 = ArgumentBookeepingHolder::create( 1, input );
  auto arg2 = ArgumentBookeepingHolder::create( 2, input);
  auto arg3 = ArgumentBookeepingHolder::create( 3, input );

  std::size_t rowlen = arg3.bookeeping[(1+arg3.ncols)*task_index];
  View<const double> wval( input.inputdata + arg3.start + arg3.ncols*task_index, rowlen );
  if( actiondata.multi_action_input ) {
    xvals.resize( rowlen );
    yvals.resize( rowlen );
    zvals.resize( rowlen );
    const auto wbooks = arg3.bookeeping.subview((1+arg3.ncols)*task_index+1, rowlen);
    for(unsigned i=0; i<rowlen; ++i) {
      xvals[i] = input.inputdata[arg0.start + getIndex( task_index, wbooks[i], arg0) ];
      yvals[i] = input.inputdata[arg1.start + getIndex( task_index, wbooks[i], arg1) ];
      zvals[i] = input.inputdata[arg2.start + getIndex( task_index, wbooks[i], arg2) ];
    }
    xpntr = xvals.data();
    ypntr = yvals.data();
    zpntr = zvals.data();
  } else {
    xpntr = input.inputdata + arg0.start + arg0.ncols*task_index;
    ypntr = input.inputdata + arg1.start + arg1.ncols*task_index;
    zpntr = input.inputdata + arg2.start + arg2.ncols*task_index;
  }
  View<const double> xval( xpntr, rowlen );
  View<const double> yval( ypntr, rowlen );
  View<const double> zval( zpntr, rowlen );
  for(unsigned i=0; i<output.derivatives.size(); ++i) {
    output.derivatives[i] = 0;
  }
  View2D<double> derivatives( output.derivatives.data(), actiondata.functions.size(), 4*arg3.shape[1] );

  Angle angle;
  Vector disti, distj;
  std::vector<double> values(4);
  std::vector<Vector> der_i(4), der_j(4);
  for(unsigned i=0; i<rowlen; ++i) {
    if( wval[i]<epsilon ) {
      continue;
    }
    disti[0] = xval[i];
    disti[1] = yval[i];
    disti[2] = zval[i];
    values[1] = disti.modulo2();
    der_i[1]=2*disti;
    der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      if( wval[j]<epsilon) {
        continue;
      }
      distj[0] = xval[j];
      distj[1] = yval[j];
      distj[2] = zval[j];
      values[2] = distj.modulo2();
      der_j[1].zero();
      der_j[2]=2*distj;
      der_i[3] = ( disti - distj );
      values[3] = der_i[3].modulo2();
      der_i[3] = 2*der_i[3];
      der_j[3] = -der_i[3];
      // Compute angle between bonds
      values[0] = angle.compute( disti, distj, der_i[0], der_j[0] );
      // Compute product of weights
      double weightij = wval[i]*wval[j];
      for(unsigned n=0; n<actiondata.functions.size(); ++n) {
        double nonweight = actiondata.functions[n].evaluate( values );
        output.values[n] += nonweight*weightij;

        if( input.noderiv ) {
          continue;
        }

        for(unsigned m=0; m<actiondata.functions[n].getNumberOfArguments(); ++m) {
          double der = weightij*actiondata.functions[n].evaluateDeriv( m, values );
          derivatives[n][i] += der*der_i[m][0];
          derivatives[n][rowlen+i] += der*der_i[m][1];
          derivatives[n][2*rowlen+i] += der*der_i[m][2];
          derivatives[n][j] += der*der_j[m][0];
          derivatives[n][rowlen+j] += der*der_j[m][1];
          derivatives[n][2*rowlen+j] += der*der_j[m][2];
        }
        derivatives[n][3*rowlen+i] += nonweight*wval[j];
        derivatives[n][3*rowlen+j] += nonweight*wval[i];
      }
    }
  }
}

void ThreeBodyGFunctions::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

int ThreeBodyGFunctions::getNumberOfValuesPerTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata ) {
  return 1;
}

void ThreeBodyGFunctions::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const ThreeBodyGFunctionsInput& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  auto arg0 = ArgumentBookeepingHolder::create( 0, input );
  auto arg1 = ArgumentBookeepingHolder::create( 1, input );
  auto arg2 = ArgumentBookeepingHolder::create( 2, input);
  auto arg3 = ArgumentBookeepingHolder::create( 3, input );
  std::size_t rowlen = arg3.bookeeping[(1+arg3.ncols)*task_index];
  if( actiondata.multi_action_input ) {
    View<const std::size_t> wbooks( arg3.bookeeping.data()+(1+arg3.ncols)*task_index+1, rowlen);
    for(unsigned j=0; j<rowlen; ++j) {
      std::size_t matpos = task_index*arg3.ncols + j;
      std::size_t xpos = getIndex( task_index, wbooks[j], arg0 );
      std::size_t ypos = getIndex( task_index, wbooks[j], arg1 );
      std::size_t zpos = getIndex( task_index, wbooks[j], arg2 );
      for(unsigned i=0; i<input.ncomponents; ++i) {
        force_indices.indices[i][j] = arg0.start + xpos;
        force_indices.indices[i][rowlen+j] = arg1.start + ypos;
        force_indices.indices[i][2*rowlen+j] = arg2.start + zpos;
        force_indices.indices[i][3*rowlen+j] = arg3.start + matpos;
      }
    }
  } else {
    for(unsigned j=0; j<rowlen; ++j) {
      std::size_t matpos = task_index*arg3.ncols + j;
      for(unsigned i=0; i<input.ncomponents; ++i) {
        force_indices.indices[i][j] = arg0.start + matpos;
        force_indices.indices[i][rowlen+j] = arg1.start + matpos;
        force_indices.indices[i][2*rowlen+j] = arg2.start + matpos;
        force_indices.indices[i][3*rowlen+j] = arg3.start + matpos;
      }
    }
  }
  for(unsigned i=0; i<input.ncomponents; ++i) {
    force_indices.threadsafe_derivatives_end[i] = force_indices.tot_indices[i] = 4*rowlen;
  }
}

}
}

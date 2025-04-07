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

You can read more about how to calculate more Behler-type symmetry functions [here](https://www.plumed-tutorials.org/lessons/23/001/data/Behler.html).

*/
//+ENDPLUMEDOC

class ThreeBodyGFunctions : public ActionWithVector {
private:
  std::vector<LeptonCall> functions;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override ;
  unsigned getNumberOfDerivatives() override;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix","three matrices containing the bond vectors of interest");
  keys.addInputKeyword("compulsory","WEIGHT","matrix","the matrix that contains the weights that should be used for each connection");
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  ActionWithValue::useCustomisableComponents( keys );
  keys.addDOI("10.1063/1.3553717");
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao) {
  if( getNumberOfArguments()!=3 ) {
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
  for(unsigned i=0; i<myargs.size(); ++i) {
    myargs[i]->buildDataStore();
  }
  std::vector<unsigned> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];

  // And now read the functions to compute
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
    functions.push_back( LeptonCall() );
    std::vector<std::string> argnames(1);
    argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) {
      argnames.push_back("rij");
    }
    if( myfunc.find("rik")!=std::string::npos ) {
      if( argnames.size()<2 ) {
        error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rik");
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
      if( argnames.size()<2 ) {
        error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rjk");
    }
    functions[i-1].set( myfunc, argnames, this, true );
  }
  checkRead();
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
  runAllTasks();
}

void ThreeBodyGFunctions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  const Value* wval = getPntrToArgument(3);
  const Value* xval = getPntrToArgument(0);
  const Value* yval = getPntrToArgument(1);
  const Value* zval = getPntrToArgument(2);
  Angle angle;
  Vector disti, distj;
  unsigned matsize = wval->getNumberOfValues();
  std::vector<double> values(4);
  std::vector<Vector> der_i(4), der_j(4);
  unsigned nbonds = wval->getRowLength( task_index ), ncols = wval->getShape()[1];
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }
    disti[0] = xval->get( ipos );
    disti[1] = yval->get( ipos );
    disti[2] = zval->get( ipos );
    values[1] = disti.modulo2();
    der_i[1]=2*disti;
    der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      unsigned jpos = ncols*task_index + wval->getRowIndex( task_index, j );
      double weightj = wval->get( jpos );
      if( weightj<epsilon ) {
        continue ;
      }
      distj[0] = xval->get( jpos );
      distj[1] = yval->get( jpos );
      distj[2] = zval->get( jpos );
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
      double weightij = weighti*weightj;
      // Now compute all symmetry functions
      for(unsigned n=0; n<functions.size(); ++n) {
        unsigned ostrn = getConstPntrToComponent(n)->getPositionInStream();
        double nonweight = functions[n].evaluate( values );
        myvals.addValue( ostrn, nonweight*weightij );
        if( doNotCalculateDerivatives() ) {
          continue;
        }

        for(unsigned m=0; m<functions[n].getNumberOfArguments(); ++m) {
          double der = weightij*functions[n].evaluateDeriv( m, values );
          myvals.addDerivative( ostrn, ipos, der*der_i[m][0] );
          myvals.addDerivative( ostrn, matsize+ipos, der*der_i[m][1] );
          myvals.addDerivative( ostrn, 2*matsize+ipos, der*der_i[m][2] );
          myvals.addDerivative( ostrn, jpos, der*der_j[m][0] );
          myvals.addDerivative( ostrn, matsize+jpos, der*der_j[m][1] );
          myvals.addDerivative( ostrn, 2*matsize+jpos, der*der_j[m][2] );
        }
        myvals.addDerivative( ostrn, 3*matsize+ipos, nonweight*weightj );
        myvals.addDerivative( ostrn, 3*matsize+jpos, nonweight*weighti );
      }
    }
  }
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  // And update the elements that have derivatives
  // Needs a separate loop here as there may be forces from j
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }

    for(unsigned n=0; n<functions.size(); ++n) {
      unsigned ostrn = getConstPntrToComponent(n)->getPositionInStream();
      myvals.updateIndex( ostrn, ipos );
      myvals.updateIndex( ostrn, matsize+ipos );
      myvals.updateIndex( ostrn, 2*matsize+ipos );
      myvals.updateIndex( ostrn, 3*matsize+ipos );
    }
  }
}

}
}

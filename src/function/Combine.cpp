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
#include "Combine.h"
#include "FunctionTemplateBase.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION COMBINE
/*
Calculate a polynomial combination of a set of other variables.

The functional form of this function is
\f[
C=\sum_{i=1}^{N_{arg}} c_i (x_i-a_i)^{p_i}
\f]

The coefficients c, the parameters a and the powers p are provided as vectors.

Notice that COMBINE is not able to predict which will be periodic domain
of the computed value automatically. The user is thus forced to specify it
explicitly. Use PERIODIC=NO if the resulting variable is not periodic,
and PERIODIC=A,B where A and B are the two boundaries if the resulting variable
is periodic.



\par Examples

The following input tells plumed to print the distance between atoms 3 and 5
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\plumedfile
DISTANCE LABEL=dist      ATOMS=3,5 COMPONENTS
COMBINE  LABEL=distance2 ARG=dist.x,dist.y,dist.z POWERS=2,2,2 PERIODIC=NO
COMBINE  LABEL=distance  ARG=distance2 POWERS=0.5 PERIODIC=NO
PRINT ARG=distance,distance2
\endplumedfile
(See also \ref PRINT and \ref DISTANCE).

The following input tells plumed to add a restraint on the
cube of a dihedral angle. Notice that since the angle has a
periodic domain
-pi,pi its cube has a domain -pi**3,pi**3.
\plumedfile
t: TORSION ATOMS=1,3,5,7
c: COMBINE ARG=t POWERS=3 PERIODIC=-31.0062766802998,31.0062766802998
RESTRAINT ARG=c KAPPA=10 AT=0
\endplumedfile



*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION COMBINE_SCALAR
/*
Calculate a polynomial combination of a set of other variables.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION COMBINE_VECTOR
/*
Add together the elements of a set of vectors elementwise

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR COMBINE_MATRIX
/*
Calculate the sum of a number of matrices

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Combine> CombineShortcut;
PLUMED_REGISTER_ACTION(CombineShortcut,"COMBINE")
typedef FunctionOfScalar<Combine> ScalarCombine;
PLUMED_REGISTER_ACTION(ScalarCombine,"COMBINE_SCALAR")
typedef FunctionOfVector<Combine> VectorCombine;
PLUMED_REGISTER_ACTION(VectorCombine,"COMBINE_VECTOR")
typedef FunctionOfMatrix<Combine> MatrixCombine;
PLUMED_REGISTER_ACTION(MatrixCombine,"COMBINE_MATRIX")

void Combine::registerKeywords(Keywords& keys) {
  keys.use("PERIODIC");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","PARAMETERS","0.0","the parameters of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficients so that in total they are equal to one");
  keys.setValueDescription("a linear compbination");
}

void Combine::read( ActionWithArguments* action ) {
  coefficients.resize( action->getNumberOfArguments() );
  parameters.resize( action->getNumberOfArguments() );
  powers.resize( action->getNumberOfArguments() );
  parseVector(action,"COEFFICIENTS",coefficients);
  if(coefficients.size()!=static_cast<unsigned>(action->getNumberOfArguments())) {
    action->error("Size of COEFFICIENTS array should be the same as number for arguments");
  }
  parseVector(action,"PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(action->getNumberOfArguments())) {
    action->error("Size of PARAMETERS array should be the same as number for arguments");
  }
  parseVector(action,"POWERS",powers);
  if(powers.size()!=static_cast<unsigned>(action->getNumberOfArguments())) {
    action->error("Size of POWERS array should be the same as number for arguments");
  }

  parseFlag(action,"NORMALIZE",normalize);
  if(normalize) {
    double n=0.0;
    for(unsigned i=0; i<coefficients.size(); i++) {
      n+=coefficients[i];
    }
    for(unsigned i=0; i<coefficients.size(); i++) {
      coefficients[i]*=(1.0/n);
    }
  }

  action->log.printf("  with coefficients:");
  for(unsigned i=0; i<coefficients.size(); i++) {
    action->log.printf(" %f",coefficients[i]);
  }
  action->log.printf("\n  with parameters:");
  for(unsigned i=0; i<parameters.size(); i++) {
    action->log.printf(" %f",parameters[i]);
  }
  action->log.printf("\n  and powers:");
  for(unsigned i=0; i<powers.size(); i++) {
    action->log.printf(" %f",powers[i]);
  }
  action->log.printf("\n");
}

void Combine::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  vals[0]=0.0;
  for(unsigned i=0; i<coefficients.size(); ++i) {
    double cv = action->difference( i, parameters[i], args[i] );
    vals[0] += coefficients[i]*pow( cv, powers[i] );
    derivatives(0,i) = coefficients[i]*powers[i]*pow(cv,powers[i]-1.0);
  }
}

}
}



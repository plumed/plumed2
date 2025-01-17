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
#include "Custom.h"
#include "core/ActionRegister.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "tools/OpenMP.h"
#include "tools/LeptonCall.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION CUSTOM
/*
Calculate a combination of variables using a custom expression.

This action computes an  arbitrary function of one or more
collective variables. Arguments are chosen with the ARG keyword,
and the function is provided with the FUNC string. Notice that this
string should contain no space. Within FUNC, one can refer to the
arguments as x,y,z, and t (up to four variables provided as ARG).
This names can be customized using the VAR keyword (see examples below).

This function is implemented using the Lepton library, that allows to evaluate
algebraic expressions and to automatically differentiate them.

If you want a function that depends not only on collective variables
but also on time you can use the \subpage TIME action.

\par Examples

The following input tells plumed to perform a metadynamics
using as a CV the difference between two distances.
\plumedfile
dAB: DISTANCE ATOMS=10,12
dAC: DISTANCE ATOMS=10,15
diff: CUSTOM ARG=dAB,dAC FUNC=y-x PERIODIC=NO
# notice: the previous line could be replaced with the following
# diff: COMBINE ARG=dAB,dAC COEFFICIENTS=-1,1
METAD ARG=diff SIGMA=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100
\endplumedfile
(see also \ref DISTANCE, \ref COMBINE, and \ref METAD).
Notice that forces applied to diff will be correctly propagated
to atoms 10, 12, and 15.
Also notice that since CUSTOM is used without the VAR option
the two arguments should be referred to as x and y in the expression FUNC.
For simple functions
such as this one it is possible to use \ref COMBINE.

The following input tells plumed to print the angle between vectors
identified by atoms 1,2 and atoms 2,3
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\plumedfile
DISTANCE LABEL=d1 ATOMS=1,2 COMPONENTS
DISTANCE LABEL=d2 ATOMS=2,3 COMPONENTS
CUSTOM ...
  LABEL=theta
  ARG=d1.x,d1.y,d1.z,d2.x,d2.y,d2.z
  VAR=ax,ay,az,bx,by,bz
  FUNC=acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz)))
  PERIODIC=NO
... CUSTOM
PRINT ARG=theta
\endplumedfile
(See also \ref PRINT and \ref DISTANCE).

Notice that this action implements a large number of functions (trigonometric, exp, log, etc).
Among the useful functions, have a look at the step function (that is the Heaviside function).
`step(x)` is defined as 1 when `x` is positive and `0` when x is negative. This allows for
a straightforward implementation of if clauses.

For example, imagine that you want to implement a restraint that only acts when a
distance is larger than 0.5. You can do it with
\plumedfile
d: DISTANCE ATOMS=10,15
m: CUSTOM ARG=d FUNC=0.5*step(0.5-x)+x*step(x-0.5) PERIODIC=NO
# check the function you are applying:
PRINT ARG=d,m FILE=checkme
RESTRAINT ARG=d AT=0.5 KAPPA=10.0
\endplumedfile
(see also \ref DISTANCE, \ref PRINT, and \ref RESTRAINT)

The meaning of the function `0.5*step(0.5-x)+x*step(x-0.5)` is:
- If x<0.5 (step(0.5-x)!=0) use 0.5
- If x>0.5 (step(x-0.5)!=0) use x
Notice that the same could have been obtained using an \ref UPPER_WALLS
However, with CUSTOM you can create way more complex definitions.

\warning If you apply forces on the variable (as in the previous example) you should
make sure that the variable is continuous!
Conversely, if you are just analyzing a trajectory you can safely use
discontinuous variables.

A possible continuity check with gnuplot is
\verbatim
# this allow to step function to be used in gnuplot:
gnuplot> step(x)=0.5*(erf(x*10000000)+1)
# here you can test your function
gnuplot> p 0.5*step(0.5-x)+x*step(x-0.5)
\endverbatim

Also notice that you can easily make logical operations on the conditions that you
create. The equivalent of the AND operator is the product: `step(1.0-x)*step(x-0.5)` is
only equal to 1 when x is between 0.5 and 1.0. By combining negation and AND you can obtain an OR. That is,
`1-step(1.0-x)*step(x-0.5)` is only equal to 1 when x is outside the 0.5-1.0 interval.

CUSTOM can be used in combination with \ref DISTANCE to implement variants of the
DISTANCE keyword that were present in PLUMED 1.3 and that allowed to compute
the distance of a point from a line defined by two other points, or the progression
along that line.
\plumedfile
# take center of atoms 1 to 10 as reference point 1
p1: CENTER ATOMS=1-10
# take center of atoms 11 to 20 as reference point 2
p2: CENTER ATOMS=11-20
# take center of atoms 21 to 30 as reference point 3
p3: CENTER ATOMS=21-30

# compute distances
d12: DISTANCE ATOMS=p1,p2
d13: DISTANCE ATOMS=p1,p3
d23: DISTANCE ATOMS=p2,p3

# compute progress variable of the projection of point p3
# along the vector joining p1 and p2
# notice that progress is measured from the middle point
onaxis: CUSTOM ARG=d13,d23,d12 FUNC=(0.5*(y^2-x^2)/z) PERIODIC=NO

# compute between point p3 and the vector joining p1 and p2
fromaxis: CUSTOM ARG=d13,d23,d12,onaxis VAR=x,y,z,o FUNC=(0.5*(y^2+x^2)-o^2-0.25*z^2) PERIODIC=NO

PRINT ARG=onaxis,fromaxis

\endplumedfile

Notice that these equations have been used to combine \ref RMSD
from different snapshots of a protein so as to define
progression (S) and distance (Z) variables \cite perez2015atp.


*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MATHEVAL_SCALAR
/*
Calculate a function of a set of input scalars

See \ref MATHEVAL

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION CUSTOM_SCALAR
/*
Calculate a function of a set of input scalars

See \ref CUSTOM

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MATHEVAL_VECTOR
/*
Calculate a function of a set of input vectors elementwise

See \ref MATHEVAL

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION CUSTOM_VECTOR
/*
Calculate a function of a set of input vectors elementwise

See \ref CUSTOM

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR CUSTOM_MATRIX
/*
Calculate an arbitrary function piecewise for one or multiple input matrices.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR MATHEVAL_MATRIX
/*
Calculate an arbitrary function piecewise for one or multiple input matrices.

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Custom> CustomShortcut;
PLUMED_REGISTER_ACTION(CustomShortcut,"CUSTOM")
PLUMED_REGISTER_ACTION(CustomShortcut,"MATHEVAL")
typedef FunctionOfScalar<Custom> ScalarCustom;
PLUMED_REGISTER_ACTION(ScalarCustom,"CUSTOM_SCALAR")
PLUMED_REGISTER_ACTION(ScalarCustom,"MATHEVAL_SCALAR")
typedef FunctionOfVector<Custom> VectorCustom;
PLUMED_REGISTER_ACTION(VectorCustom,"CUSTOM_VECTOR")
PLUMED_REGISTER_ACTION(VectorCustom,"MATHEVAL_VECTOR")
typedef FunctionOfMatrix<Custom> MatrixCustom;
PLUMED_REGISTER_ACTION(MatrixCustom,"CUSTOM_MATRIX")
PLUMED_REGISTER_ACTION(MatrixCustom,"MATHEVAL_MATRIX")

//+PLUMEDOC FUNCTION MATHEVAL
/*
An alias to the CUSTOM function that can also be used to calaculate combinations of variables using a custom expression.

Documentation for this action is identical to that for \ref CUSTOM

This alias is kept in order to maintain compatibility with previous PLUMED versions.
However, notice that as of PLUMED 2.5 the libmatheval library is not linked anymore,
and the \ref MATHEVAL function is implemented using the Lepton library.

\par Examples

Just replace \ref CUSTOM with \ref MATHEVAL.

\plumedfile
d: DISTANCE ATOMS=10,15
m: MATHEVAL ARG=d FUNC=0.5*step(0.5-x)+x*step(x-0.5) PERIODIC=NO
# check the function you are applying:
PRINT ARG=d,m FILE=checkme
RESTRAINT ARG=d AT=0.5 KAPPA=10.0
\endplumedfile
(see also \ref DISTANCE, \ref PRINT, and \ref RESTRAINT)

*/
//+ENDPLUMEDOC

void Custom::registerKeywords(Keywords& keys) {
  keys.use("PERIODIC");
  keys.add("compulsory","FUNC","the function you wish to evaluate");
  keys.add("optional","VAR","the names to give each of the arguments in the function.  If you have up to three arguments in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
  keys.setValueDescription("scalar/vector/matrix/grid","an arbitrary function");
}

void Custom::read( ActionWithArguments* action ) {
  // Read in the variables
  std::vector<std::string> var;
  parseVector(action,"VAR",var);
  parse(action,"FUNC",func);
  if(var.size()==0) {
    var.resize(action->getNumberOfArguments());
    if(var.size()>3) {
      action->error("Using more than 3 arguments you should explicitly write their names with VAR");
    }
    if(var.size()>0) {
      var[0]="x";
    }
    if(var.size()>1) {
      var[1]="y";
    }
    if(var.size()>2) {
      var[2]="z";
    }
  }
  if(var.size()!=action->getNumberOfArguments()) {
    action->error("Size of VAR array should be the same as number of arguments");
  }
  // Check for operations that are not multiplication (this can probably be done much more cleverly)
  bool onlymultiplication = func.find("*")!=std::string::npos;
  // Find first bracket in expression
  if( func.find("(")!=std::string::npos ) {
    std::size_t br = func.find_first_of("(");
    std::string subexpr=func.substr(0,br);
    onlymultiplication = func.find("*")!=std::string::npos;
    if( subexpr.find("/")!=std::string::npos ) {
      std::size_t sl = func.find_first_of("/");
      std::string aa = subexpr.substr(0,sl);
      subexpr=aa;
    }
    if( subexpr.find("+")!=std::string::npos || subexpr.find("-")!=std::string::npos ) {
      onlymultiplication=false;
    }
    // Now work out which vars are in multiplication
    if( onlymultiplication ) {
      for(unsigned i=0; i<var.size(); ++i) {
        if( subexpr.find(var[i])!=std::string::npos &&
            action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
          check_multiplication_vars.push_back(i);
        }
      }
    }
  } else if( func.find("/")!=std::string::npos ) {
    onlymultiplication=true;
    if( func.find("+")!=std::string::npos || func.find("-")!=std::string::npos ) {
      onlymultiplication=false;
    }
    if( onlymultiplication ) {
      std::size_t br = func.find_first_of("/");
      std::string subexpr=func.substr(0,br);
      for(unsigned i=0; i<var.size(); ++i) {
        if( subexpr.find(var[i])!=std::string::npos &&
            action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
          check_multiplication_vars.push_back(i);
        }
      }
    }
  } else if( func.find("+")!=std::string::npos || func.find("-")!=std::string::npos ) {
    onlymultiplication=false;
  } else {
    for(unsigned i=0; i<var.size(); ++i) {
      if( action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
        check_multiplication_vars.push_back(i);
      }
    }
  }
  if( check_multiplication_vars.size()>0 ) {
    action->log.printf("  optimizing implementation as function only involves multiplication \n");
  }

  action->log.printf("  with function : %s\n",func.c_str());
  action->log.printf("  with variables :");
  for(unsigned i=0; i<var.size(); i++) {
    action->log.printf(" %s",var[i].c_str());
  }
  action->log.printf("\n");
  function.set( func, var, action );
  std::vector<double> zeros( action->getNumberOfArguments(), 0 );
  double fval = abs(function.evaluate(zeros));
  zerowhenallzero=(fval<epsilon );
  if( zerowhenallzero ) {
    action->log.printf("  not calculating when all arguments are zero \n");
  }
}

std::string Custom::getGraphInfo( const std::string& name ) const {
  return FunctionTemplateBase::getGraphInfo( name ) + + "\n" + "FUNC=" + func;
}

bool Custom::getDerivativeZeroIfValueIsZero() const {
  return check_multiplication_vars.size()>0;
}

std::vector<Value*> Custom::getArgumentsToCheck( const std::vector<Value*>& args ) {
  std::vector<Value*> fargs( check_multiplication_vars.size() );
  for(unsigned i=0; i<check_multiplication_vars.size(); ++i) {
    fargs[i] = args[check_multiplication_vars[i]];
  }
  return fargs;
}

void Custom::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  if( args.size()>1 ) {
    bool allzero=false;
    if( check_multiplication_vars.size()>0 ) {
      for(unsigned i=0; i<check_multiplication_vars.size(); ++i) {
        if( fabs(args[check_multiplication_vars[i]])<epsilon ) {
          allzero=true;
          break;
        }
      }
    } else if( zerowhenallzero ) {
      allzero=(fabs(args[0])<epsilon);
      for(unsigned i=1; i<args.size(); ++i) {
        if( fabs(args[i])>epsilon ) {
          allzero=false;
          break;
        }
      }
    }
    if( allzero ) {
      vals[0]=0;
      for(unsigned i=0; i<args.size(); i++) {
        derivatives(0,i) = 0.0;
      }
      return;
    }
  }
  vals[0] = function.evaluate( args );
  if( !noderiv ) {
    for(unsigned i=0; i<args.size(); i++) {
      derivatives(0,i) = function.evaluateDeriv( i, args );
    }
  }
}

}
}



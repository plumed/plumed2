/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION COMBINE
/*
Calculate a polynomial combination of a set of other variables.

The functional form of this function is
\f[
C=\sum_{i=1}^{N_{arg}} c_i x_i^{p_i}
\f]

The coefficients c and powers p are provided as vectors.



\par Examples
The following input tells plumed to print the distance between atoms 3 and 5
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\verbatim
DISTANCE LABEL=dist      ATOMS=3,5 COMPONENTS
COMBINE  LABEL=distance2 ARG=dist.x,dist.y,dist.z POWERS=2,2,2 PERIODIC=NO
COMBINE  LABEL=distance  ARG=distance2 POWERS=0.5 PERIODIC=NO
PRINT ARG=distance,distance2
\endverbatim
(See also \ref PRINT and \ref DISTANCE).


*/
//+ENDPLUMEDOC


class Combine :
  public Function
{
  bool normalize;
  std::vector<double> coefficients;
  std::vector<double> powers;
public:
  explicit Combine(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Combine,"COMBINE")

void Combine::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficents so that in total they are equal to one");
}

Combine::Combine(const ActionOptions&ao):
Action(ao),
Function(ao),
normalize(false),
coefficients(getNumberOfArguments(),1.0),
powers(getNumberOfArguments(),1.0)
{
  parseVector("COEFFICIENTS",coefficients);
  if(coefficients.size()!=static_cast<unsigned>(getNumberOfArguments()))
    error("Size of COEFFICIENTS array should be the same as number for arguments");

  parseVector("POWERS",powers);
  if(powers.size()!=static_cast<unsigned>(getNumberOfArguments()))
    error("Size of POWERS array should be the same as number for arguments");

  parseFlag("NORMALIZE",normalize);

  if(normalize){
    double n=0.0;
    for(unsigned i=0;i<coefficients.size();i++) n+=coefficients[i];
    for(unsigned i=0;i<coefficients.size();i++) coefficients[i]*=(1.0/n);
  }
 
  addValueWithDerivatives(); 
  checkRead();

  log.printf("  with coefficients:");
  for(unsigned i=0;i<coefficients.size();i++) log.printf(" %f",coefficients[i]);
  log.printf("\n");
  log.printf("  and powers:");
  for(unsigned i=0;i<powers.size();i++) log.printf(" %f",powers[i]);
  log.printf("\n");
}

void Combine::calculate(){
  double combine=0.0;
  for(unsigned i=0;i<coefficients.size();++i){
    combine+=coefficients[i]*pow(getArgument(i),powers[i]);
    setDerivative(i,coefficients[i]*powers[i]*pow(getArgument(i),powers[i]-1.0));
  };
  setValue(combine);
}

}
}



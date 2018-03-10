/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION MOMENTS
/*

\par Examples

*/
//+ENDPLUMEDOC


class Moments :
  public Function
{
  bool isperiodic;
  double min, max, pfactor;
  std::vector<int> powers;
public:
  explicit Moments(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculate();
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const { plumed_error(); }
};


PLUMED_REGISTER_ACTION(Moments,"MOMENTS")

void Moments::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the list of argument from which you would like to calculate moments.");
  keys.add("compulsory","POWERS","calculate the central moments of the distribution of collective variables. "
           "The \\f$m\\f$th central moment of a distribution is calculated using \\f$\\frac{1}{N} \\sum_{i=1}^N ( s_i - \\overline{s} )^m \\f$, where \\f$\\overline{s}\\f$ is "
           "the average for the distribution. The POWERS keyword takes a lists of integers as input or a range. Each integer is a value of \\f$m\\f$. The final "
           "calculated values can be referenced using moment-\\f$m\\f$.");
  keys.addOutputComponent("moment","default","the central moments of the distribution of values. The second central moment "
                          "would be referenced elsewhere in the input file using "
                          "<em>label</em>.moment-2, the third as <em>label</em>.moment-3, etc.");
}

Moments::Moments(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  isperiodic(getPntrToArgument(0)->isPeriodic())
{
  if( isperiodic ) {
    std::string str_min, str_max; getPntrToArgument(0)->getDomain( str_min, str_max );
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( !getPntrToArgument(i)->isPeriodic() ) error("cannot mix periodic and non periodic variables when calculating moments");
      std::string str_min2, str_max2; getPntrToArgument(0)->getDomain( str_min2, str_max2);
      if( str_min!=str_min2 || str_max!=str_max2 ) error("all input arguments should have same domain when calculating moments");
    }
    Tools::convert(str_min,min); Tools::convert(str_max,max); pfactor = 2*pi / ( max-min );
  } else {
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->isPeriodic() ) error("cannot mix periodic and non periodic variables when calculating moments");
    }
  }

  parseVector("POWERS",powers);
  for(unsigned i=0; i<powers.size(); ++i) {
    if( powers[i]<2 ) error("first central moment is zero do you really need to calculate that");
    log.printf("  computing %dth central moment of distribution of input cvs \n", powers[i]);
    std::string pwstr; Tools::convert( powers[i], pwstr );
    addComponentWithDerivatives( "moment-" + pwstr );
    componentIsNotPeriodic( "moment-" + pwstr );
  }
  checkRead();
}

void Moments::calculate() {
  double mean=0; double inorm = 1.0 / static_cast<double>( getNumberOfScalarArguments() );
  if( isperiodic ) {
    double sinsum=0, cossum=0, val;
    for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) { val=pfactor*( getArgumentScalar(i) - min ); sinsum+=sin(val); cossum+=cos(val); }
    mean = 0.5 + atan2( inorm*sinsum, inorm*cossum ) / (2*pi);
    mean = min + (max-min)*mean;
  } else {
    for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) mean+=getArgumentScalar(i);
    mean *= inorm;
  }

  Value* arg0 = getPntrToArgument(0);
  for(unsigned npow=0; npow<powers.size(); ++npow) {
    double dev1=0;
    if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) dev1+=pow( arg0->difference( mean, getArgumentScalar(i) ), powers[npow] - 1 );
      dev1*=inorm;
    }
    double moment=0; Value* myval = getPntrToComponent(npow); double prefactor = powers[npow]*inorm;
    for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) {
      double tmp=arg0->difference( mean, getArgumentScalar(i) );
      moment+=pow( tmp, powers[npow] );
      if( !doNotCalculateDerivatives() ) {
        myval->addDerivative( i, prefactor*(pow( tmp, powers[npow] - 1 ) - dev1) );
      }
    }
    myval->set( inorm*moment );
  }
}

}
}



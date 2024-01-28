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
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"
#include "FunctionTemplateBase.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION MOMENTS
/*
Calculate the moments of the distribution of input quantities

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MOMENTS_SCALAR
/*
Calculate the moments of the distribution of input quantities

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MOMENTS_VECTOR
/*
Calculate the moments of the distribution of input vectors

\par Examples

*/
//+ENDPLUMEDOC

class Moments : public FunctionTemplateBase {
  bool isperiodic, scalar_out;
  double min, max, pfactor;
  std::vector<int> powers;
public:
  void registerKeywords(Keywords& keys) override;
  void read( ActionWithArguments* action ) override;
  bool zeroRank() const override { return scalar_out; }
  bool doWithTasks() const override { return !scalar_out; }
  std::vector<std::string> getComponentsPerLabel() const override ;
  void setPeriodicityForOutputs( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionShortcut<Moments> MomentsShortcut;
PLUMED_REGISTER_ACTION(MomentsShortcut,"MOMENTS")
typedef FunctionOfScalar<Moments> ScalarMoments;
PLUMED_REGISTER_ACTION(ScalarMoments,"MOMENTS_SCALAR")
typedef FunctionOfVector<Moments> VectorMoments;
PLUMED_REGISTER_ACTION(VectorMoments,"MOMENTS_VECTOR")

void Moments::registerKeywords(Keywords& keys) {
  keys.add("compulsory","POWERS","calculate the central moments of the distribution of collective variables. "
           "The \\f$m\\f$th central moment of a distribution is calculated using \\f$\\frac{1}{N} \\sum_{i=1}^N ( s_i - \\overline{s} )^m \\f$, where \\f$\\overline{s}\\f$ is "
           "the average for the distribution. The POWERS keyword takes a lists of integers as input or a range. Each integer is a value of \\f$m\\f$. The final "
           "calculated values can be referenced using moment-\\f$m\\f$.");
  keys.addOutputComponent("moment","default","the central moments of the distribution of values. The second central moment "
                          "would be referenced elsewhere in the input file using "
                          "<em>label</em>.moment-2, the third as <em>label</em>.moment-3, etc.");
}

void Moments::read( ActionWithArguments* action ) {
  scalar_out = action->getNumberOfArguments()==1;
  if( scalar_out && action->getPntrToArgument(0)->getRank()==0 ) action->error("cannot calculate moments if only given one variable");

  isperiodic = (action->getPntrToArgument(0))->isPeriodic();
  if( isperiodic ) {
    std::string str_min, str_max; (action->getPntrToArgument(0))->getDomain( str_min, str_max );
    for(unsigned i=1; i<action->getNumberOfArguments(); ++i) {
      if( !(action->getPntrToArgument(i))->isPeriodic() ) action->error("cannot mix periodic and non periodic variables when calculating moments");
      std::string str_min2, str_max2; (action->getPntrToArgument(i))->getDomain( str_min2, str_max2);
      if( str_min!=str_min2 || str_max!=str_max2 ) action->error("all input arguments should have same domain when calculating moments");
    }
    Tools::convert(str_min,min); Tools::convert(str_max,max); pfactor = 2*pi / ( max-min );
  } else {
    for(unsigned i=1; i<action->getNumberOfArguments(); ++i) {
      if( (action->getPntrToArgument(i))->isPeriodic() ) action->error("cannot mix periodic and non periodic variables when calculating moments");
    }
  }

  parseVector(action,"POWERS",powers);
  for(unsigned i=0; i<powers.size(); ++i) {
    if( powers[i]<2 ) action->error("first central moment is zero do you really need to calculate that");
    action->log.printf("  computing %dth central moment of distribution of input cvs \n", powers[i]);
  }
}

std::vector<std::string> Moments::getComponentsPerLabel() const {
  std::vector<std::string> comp; std::string num;
  for(unsigned i=0; i<powers.size(); ++i) {
    Tools::convert(powers[i],num); comp.push_back( "-" + num );
  }
  return comp;
}

void Moments::setPeriodicityForOutputs( ActionWithValue* action ) {
  for(unsigned i=0; i<powers.size(); ++i) { std::string num; Tools::convert(powers[i],num); action->componentIsNotPeriodic("moment-" + num); }
}

void Moments::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  double mean=0; double inorm = 1.0 / static_cast<double>( args.size() );
  if( isperiodic ) {
    double sinsum=0, cossum=0, val;
    for(unsigned i=0; i<args.size(); ++i) { val=pfactor*( args[i] - min ); sinsum+=sin(val); cossum+=cos(val); }
    mean = 0.5 + atan2( inorm*sinsum, inorm*cossum ) / (2*pi);
    mean = min + (max-min)*mean;
  } else {
    for(unsigned i=0; i<args.size(); ++i) mean+=args[i];
    mean *= inorm;
  }

  Value* arg0 = action->getPntrToArgument(0);
  for(unsigned npow=0; npow<powers.size(); ++npow) {
    double dev1=0;
    for(unsigned i=0; i<args.size(); ++i) dev1+=pow( arg0->difference( mean, args[i] ), powers[npow] - 1 );
    dev1*=inorm; vals[npow] = 0; double prefactor = powers[npow]*inorm;
    for(unsigned i=0; i<args.size(); ++i) {
      double tmp=arg0->difference( mean, args[i] ); vals[npow] += inorm*pow( tmp, powers[npow] );
      derivatives(npow,i) = prefactor*(pow( tmp, powers[npow] - 1 ) - dev1);
    }
  }
}

}
}



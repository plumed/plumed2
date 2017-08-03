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
#include "tools/HistogramBead.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BETWEEN
/*
Use a switching function to determine how many of the input variables are within a certain range.

\par Examples

*/
//+ENDPLUMEDOC


class Between :
  public Function
{
  HistogramBead hist;
public:
  explicit Between(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Between,"BETWEEN")

void Between::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG"); 
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER, SMEAR and KERNEL keywords.");
}

Between::Between(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  std::string str_min, str_max, tstr_min, tstr_max;
  bool isPeriodic = getPntrToArgument(0)->isPeriodic();
  if( isPeriodic ) getPntrToArgument(0)->getDomain( str_min, str_max );
  for(unsigned i=1;i<getNumberOfArguments();++i){
      if( isPeriodic ){
          if( !getPntrToArgument(i)->isPeriodic() ) error("cannot mix periodic and non periodic arguments");
          getPntrToArgument(i)->getDomain( tstr_min, tstr_max ); 
          if( tstr_min!=str_min || tstr_max!=str_max ) error("cannot mix periodic arguments with different domains");
      }
  }
  std::string hinput; parse("SWITCH",hinput);
  if(hinput.length()==0) {
     std::string low, up, sme;
     parse("LOWER",low); parse("UPPER",up); parse("SMEAR",sme);
     hinput = "GAUSSIAN LOWER=" + low + " UPPER=" + up + " SMEAR=" + sme;
  }
  std::string errors; hist.set( hinput, errors );
  if( errors.size()!=0 ) error( errors );
  log.printf("  %s \n", hist.description().c_str() );

  if( !isPeriodic ) hist.isNotPeriodic();
  else {
    double min; Tools::convert( str_min, min );
    double max; Tools::convert( str_max, max );
    hist.isPeriodic( min, max );
  }

  addValueWithDerivatives();
  checkRead();
}

void Between::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  plumed_dbg_assert( args.size()==1 ); double dv, f = hist.calculate(args[0], dv);  
  addValue( 0, f, myvals ); addDerivative( 0, 0, dv, myvals );
}

}
}



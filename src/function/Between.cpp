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
  keys.add("compulsory","KERNEL","GAUSSIAN","the type of kernel function to use");
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
}

Between::Between(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  std::string str_min, str_max, tstr_min, tstr_max;
  bool isPeriodic = getPntrToComponent(0)->isPeriodic();
  if( isPeriodic ) getPntrToComponent(0)->getDomain( str_min, str_max );
  for(unsigned i=1;i<getNumberOfComponents();++i){
      if( isPeriodic ){
          if( !getPntrToComponent(i)->isPeriodic() ) error("cannot mix periodic and non periodic arguments");
          getPntrToComponent(i)->getDomain( tstr_min, tstr_max ); 
          if( tstr_min!=str_min || tstr_max!=str_max ) error("cannot mix periodic arguments with different domains");
      }
  }
  std::string ktype, low, up, sme; parse("KERNEL",ktype);
  parse("LOWER",low); parse("UPPER",up); parse("SMEAR",sme);
  std::string hinput = ktype + " LOWER=" + low + " UPPER=" + up + " SMEAR" + sme;
  std::string errors; hist.set( hinput, errors );
  if( errors.size()!=0 ) error( errors );

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
  setValue( 0, f, myvals ); addDerivative( 0, 0, dv, myvals );
}

}
}



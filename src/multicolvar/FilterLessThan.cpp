/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2016 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "MultiColvarFilter.h"

//+PLUMEDOC MFILTERS MFILTER_LESS
/*
This action can be used to filter the distribution of colvar values in a multicolvar 
so that one can compute the mean and so on for only those multicolvars less than a tolerance.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class FilterLess : public MultiColvarFilter {
private:
  SwitchingFunction sf;
public:
  static void registerKeywords( Keywords& keys );
  explicit FilterLess(const ActionOptions& ao);
  double applyFilter( const double& val, double& df ) const ;
}; 

PLUMED_REGISTER_ACTION(FilterLess,"MFILTER_LESS")

void FilterLess::registerKeywords( Keywords& keys ){
  MultiColvarFilter::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
}

FilterLess::FilterLess(const ActionOptions& ao):
Action(ao),
MultiColvarFilter(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     sf.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     sf.set(nn,mm,r_0,d_0);
  }
  log.printf("  filtering colvar values and focussing only on those less than %s\n",( sf.description() ).c_str() );

  checkRead();  
}

double FilterLess::applyFilter( const double& val, double& df ) const {
  double f = sf.calculate( val, df ); df*=val;
  return f;
}

}
}

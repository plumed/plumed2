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
#include "tools/HistogramBead.h"
#include "MultiColvarFilter.h"

//+PLUMEDOC MFILTERS MFILTER_BETWEEN
/*
This action can be used to filter the distribution of colvar values in a multicolvar 
so that one can compute the mean and so on for only those multicolvars within a certain range.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class FilterBetween : public MultiColvarFilter {
private:
  HistogramBead hb;
public:
  static void registerKeywords( Keywords& keys );
  explicit FilterBetween(const ActionOptions& ao);
  double applyFilter( const double& val, double& df ) const ;
}; 

PLUMED_REGISTER_ACTION(FilterBetween,"MFILTER_BETWEEN")

void FilterBetween::registerKeywords( Keywords& keys ){
  MultiColvarFilter::registerKeywords( keys );
  keys.add("compulsory","LOWER","the lower boundary for the range of interest");
  keys.add("compulsory","UPPER","the upper boundary for the range of interest");
  keys.add("compulsory","SMEAR","0.5","the ammount by which to smear the value for kernel density estimation");
  keys.add("optional","BEAD","This keywords is used if you want to employ an alternative to the function defeind above. "
                             "The following provides information on the \\ref histogrambead that are available. "
                             "When this keyword is present you no longer need the LOWER, UPPER and SMEAR keywords.");   
}

FilterBetween::FilterBetween(const ActionOptions& ao):
Action(ao),
MultiColvarFilter(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("BEAD",sw);
  if( getPntrToMultiColvar()->isPeriodic() ){
     std::string min, max; getPntrToMultiColvar()->retrieveDomain( min, max );
     double mlow, mhigh; Tools::convert( min,mlow ); Tools::convert( max, mhigh);
     hb.isPeriodic( mlow, mhigh );
  } else {
     hb.isNotPeriodic();
  }

  if(sw.length()>0){
     hb.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading BEAD keyword : " + errors );
  } else {
     double l, u, s; std::string ll, uu, ss;
     parse("LOWER",l); parse("UPPER",u); parse("SMEAR",s);
     Tools::convert(l,ll); Tools::convert(u,uu); Tools::convert(s,ss);
     sw="GAUSSIAN LOWER=" + ll + " UPPER=" + uu + " SMEAR=" + ss;
     hb.set(sw,errors); plumed_massert(errors.length()==0,"problems with bead" + errors);
  }
  log.printf("  filtering colvar values and focussing only on those values in range %s\n",( hb.description() ).c_str() );

  checkRead();  
}

double FilterBetween::applyFilter( const double& val, double& df ) const {
  double f = hb.calculate( val, df ); 
  return f;
}

}
}

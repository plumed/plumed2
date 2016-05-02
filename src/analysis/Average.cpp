/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
#include "vesselbase/ActionWithAveraging.h"
#include "core/ActionRegister.h"
#include "AverageVessel.h"

//+PLUMEDOC GRIDCALC AVERAGE 
/* 
Calculate the ensemble average of a collective variable

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class Average : public vesselbase::ActionWithAveraging {
private:
  AverageVessel* myaverage;
public:
  static void registerKeywords( Keywords& keys );
  explicit Average( const ActionOptions& );
  void performOperations( const bool& from_update );
  void finishAveraging();
  bool isPeriodic(){ return false; } 
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(Average,"AVERAGE")

void Average::registerKeywords( Keywords& keys ){
  vesselbase::ActionWithAveraging::registerKeywords( keys ); keys.use("ARG");
}

Average::Average( const ActionOptions& ao ):
Action(ao),
ActionWithAveraging(ao)
{
   addValue(); // Create a value so that we can output the average
   if( getNumberOfArguments()!=1 ) error("only one quantity can be averaged at a time");
   std::string instring; 
   if( getPntrToArgument(0)->isPeriodic() ){
      std::string min, max; getPntrToArgument(0)->getDomain(min,max); 
      instring = "PERIODIC=" + min + "," + max; setPeriodic( min, max );
   } else {
      setNotPeriodic();
   }
   // Create a vessel to hold the average
   vesselbase::VesselOptions da("myaverage","",-1,instring,this);
   Keywords keys; AverageVessel::registerKeywords( keys );
   vesselbase::VesselOptions dar( da, keys );
   myaverage = new AverageVessel(dar); setAveragingAction( myaverage, false );
}

void Average::performOperations( const bool& from_update ){
   myaverage->accumulate( cweight, getArgument(0) );
}

void Average::finishAveraging(){
   setValue( myaverage->getAverage() );
}

}
}

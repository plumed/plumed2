/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "FunctionVessel.h"
#include "MultiColvar.h"

namespace PLMD {

class mean : public NormedSumVessel {
private:
  MultiColvar* mycolv;
public:
  static void reserveKeyword( Keywords& keys );
  mean( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
};

PLUMED_REGISTER_VESSEL(mean,"AVERAGE")

void mean::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("AVERAGE",false,"take the average value of these variables and store it in value called average.");
}

mean::mean( const VesselOptions& da ) :
NormedSumVessel(da)
{
  if( getAction()->isPeriodic() ) error("MEAN cannot be used with periodic variables");

  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "average is used to take the average values of multi colvars");

  useNorm();
  addOutput("average");
  log.printf("  value %s.average contains the average value\n",(getAction()->getLabel()).c_str());
}

void mean::compute( const unsigned& i, const unsigned& j, Value& theval ){
  plumed_assert( j==0 );
  theval=mycolv->retreiveLastCalculatedValue(); 
}

void mean::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveColvarWeight( i, weight );
}

}

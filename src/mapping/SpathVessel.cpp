/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class SpathVessel : public vesselbase::FunctionVessel {
private:
  bool foundoneclose;
  unsigned mycoordnumber;
  Mapping* mymap;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  SpathVessel( const vesselbase::VesselOptions& da );
  std::string function_description();
  void prepare();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(SpathVessel,"SPATH")

void SpathVessel::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords(keys);
}

void SpathVessel::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("SPATH",false,"docs should not appear",true);
  keys.addOutputComponent("spath","SPATH","the position on the path");
}

SpathVessel::SpathVessel( const vesselbase::VesselOptions& da ):
FunctionVessel(da)
{
  mymap=dynamic_cast<Mapping*>( getAction() );
  plumed_massert( mymap, "SpathVessel can only be used with mappings");
  // Retrieve the index of the property in the underlying mapping
  mycoordnumber=mymap->getPropertyIndex( getLabel() );
}

std::string SpathVessel::function_description(){
  return "the position on the path";
}

void SpathVessel::prepare(){
  foundoneclose=false;
}

bool SpathVessel::calculate(){
  double weight=getAction()->getElementValue(0);
  bool addval=addValueUsingTolerance( 1, weight );
  if( addval ){
      foundoneclose=true;
      double pp=mymap->getPropertyValue( mycoordnumber );
      addValueIgnoringTolerance( 0, weight*pp );
      getAction()->chainRuleForElementDerivatives( 0, 0, pp, this );
      getAction()->chainRuleForElementDerivatives( 1, 0, 1.0, this );
  }
  return ( weight>getNLTolerance() );
}

void SpathVessel::finish(){
  if( !foundoneclose ) error("all weights in path are less than tolerance. Increase tolerance or reconsider path");
  double numerator=getFinalValue(0), denom=getFinalValue(1);
  setOutputValue( numerator / denom );
  std::vector<double> df(2);
  df[0] = 1.0 / denom;
  df[1] = - numerator / (denom*denom);
  mergeFinalDerivatives( df );
}

}
}

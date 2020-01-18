/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class SpathVessel : public vesselbase::FunctionVessel {
private:
  bool foundoneclose;
  Mapping* mymap;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit SpathVessel( const vesselbase::VesselOptions& da );
  std::string value_descriptor() override;
  void prepare() override;
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const override;
};

PLUMED_REGISTER_VESSEL(SpathVessel,"SPATH")

void SpathVessel::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords(keys);
}

void SpathVessel::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","SPATH","docs should not appear");
  keys.addOutputComponent("spath","SPATH","the position on the path");
}

SpathVessel::SpathVessel( const vesselbase::VesselOptions& da ):
  FunctionVessel(da),
  foundoneclose(false)
{
  mymap=dynamic_cast<Mapping*>( getAction() );
  plumed_massert( mymap, "SpathVessel can only be used with mappings");
  // Retrieve the index of the property in the underlying mapping
  usetol=true; norm=true;

  for(unsigned i=0; i<mymap->getFullNumberOfTasks(); ++i) {
    if( mymap->getTaskCode(i)!=mymap->getPositionInFullTaskList(i) ) error("mismatched tasks and codes");
  }
}

std::string SpathVessel::value_descriptor() {
  return "the position on the path";
}

void SpathVessel::prepare() {
  foundoneclose=false;
}

void SpathVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const {
  double pp=mymap->getPropertyValue( current, getLabel() ), weight=myvals.get(0);
  if( weight<getTolerance() ) return;
  unsigned nderivatives=getFinalValue()->getNumberOfDerivatives();
  buffer[bufstart] += weight*pp; buffer[bufstart+1+nderivatives] += weight;
  if( getAction()->derivativesAreRequired() ) {
    myvals.chainRule( 0, 0, 1, 0, pp, bufstart, buffer );
    myvals.chainRule( 0, 1, 1, 0, 1.0, bufstart, buffer );
  }
}

}
}

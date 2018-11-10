/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "SetupReferenceBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace setup {

class ReadVector : public SetupReferenceBase {
public: 
  static void registerKeywords( Keywords& keys );
  explicit ReadVector(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(ReadVector,"READ_VECTOR")

void ReadVector::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","REFERENCE","the values of the components of the vector");
}

ReadVector::ReadVector(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
  std::vector<double> reference; parseVector("REFERENCE",reference); 
  std::vector<unsigned> shape(1); shape[0] = reference.size();
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
  for(unsigned i=0;i<reference.size();++i) getPntrToComponent(0)->set( i, reference[i] );
}

}
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace multicolvar {

class SelectMassCharge : public MultiColvarBase {
private:
  bool constant;
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectMassCharge(const ActionOptions&);
// active methods:
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const override;
  unsigned getNumberOfDerivatives() const override;
};

PLUMED_REGISTER_ACTION(SelectMassCharge,"MASSES")
PLUMED_REGISTER_ACTION(SelectMassCharge,"CHARGES")

void SelectMassCharge::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys ); keys.remove("ATOMS");
  keys.add("atoms","ATOM","the index of the atom whose mass/charge you would like to use");
  keys.add("atoms","ATOMS","the indices of the atoms that you would like to use the masses/charges of");
}

SelectMassCharge::SelectMassCharge(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  constant(true)
{
  if(getNumberOfAtomsInEachCV()!=1) error("Number of specified atoms should be 1");
  addValue(); setNotPeriodic(); getPntrToOutput(0)->buildDataStore( getLabel() );
  // Determine if any of the masses/charges are not constant
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      unsigned nn, kk; getValueIndices( getAbsoluteIndex(i), nn, kk ); 
      if( getName()=="MASSES" && !masv[nn]->isConstant() ) constant=false;
      if( getName()=="CHARGES" && !chargev[nn]->isConstant() ) constant=false;
  }
  if( constant ) getPntrToOutput(0)->setConstant();
}

unsigned SelectMassCharge::getNumberOfDerivatives() const {
  return getNumberOfAtoms();
}

// calculator
void SelectMassCharge::compute( const std::vector<Vector>& pos, MultiValue& myvals ) const {
  
  if( getName()=="MASSES" ) setValue( 0, getMass( myvals.getTaskIndex() ), myvals );
  else if( chargesWereSet ) setValue( 0, getCharge( myvals.getTaskIndex() ) , myvals );

  if( constant || doNotCalculateDerivatives() ) return;
  unsigned itask=myvals.getTaskIndex(), jval=getPntrToOutput(0)->getPositionInStream();
  myvals.addDerivative( jval, ablocks[0][itask], 1.0 ); myvals.updateIndex( jval, ablocks[0][itask] );
}

}
}




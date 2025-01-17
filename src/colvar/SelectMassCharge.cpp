/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "Colvar.h"
#include "ColvarShortcut.h"
#include "core/ActionRegister.h"
#include "MultiColvarTemplate.h"
#include "tools/Pbc.h"

//+PLUMEDOC MCOLVAR CHARGE
/*
Get the charges of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR CHARGE_SCALAR
/*
Get the charges of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR CHARGE_VECTOR
/*
Get the charges of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR MASS
/*
Get the mass of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR MASS_SCALAR
/*
Get the mass of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR MASS_VECTOR
/*
Get the mass of one or multiple atoms

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

class SelectMassCharge : public Colvar {
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectMassCharge(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<SelectMassCharge> MQShortcut;
PLUMED_REGISTER_ACTION(MQShortcut,"MASS")
PLUMED_REGISTER_ACTION(MQShortcut,"CHARGE")
PLUMED_REGISTER_ACTION(SelectMassCharge,"MASS_SCALAR")
PLUMED_REGISTER_ACTION(SelectMassCharge,"CHARGE_SCALAR")
typedef MultiColvarTemplate<SelectMassCharge> MQMulti;
PLUMED_REGISTER_ACTION(MQMulti,"MASS_VECTOR")
PLUMED_REGISTER_ACTION(MQMulti,"CHARGE_VECTOR")

void SelectMassCharge::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom number");
  keys.add("atoms","ATOMS","the atom numbers that you would like to store the masses and charges of");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  std::string acname = keys.getDisplayName();
  std::size_t und = acname.find("_SCALAR");
  if( und==std::string::npos ) {
    und = acname.find("_VECTOR");
  }
  keys.setDisplayName( acname.substr(0,und) );
  keys.setValueDescription("the " + keys.getDisplayName() + " of the atom");
}

SelectMassCharge::SelectMassCharge(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  unsigned mode=getModeAndSetupValues(this);
  requestAtoms(atoms);
}

void SelectMassCharge::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOM",num,t);
  if( t.size()==1 ) {
    aa->log.printf("  for atom %d\n",t[0].serial());
  } else if( num<0 || t.size()!=0 ) {
    aa->error("Number of specified atoms should be 1");
  }
}

unsigned SelectMassCharge::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  bool constant=true;
  ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( av );
  plumed_assert( aa );
  for(unsigned i=0; i<aa->getNumberOfAtoms(); ++i) {
    std::pair<std::size_t,std::size_t> p = aa->getValueIndices( aa->getAbsoluteIndex(i) );
    if( av->getName().find("MASS")!=std::string::npos && !aa->masv[p.first]->isConstant() ) {
      constant=false;
    }
    if( av->getName().find("CHARGE")!=std::string::npos && !aa->chargev[p.first]->isConstant() ) {
      constant=false;
    }
  }
  if( !constant ) {
    av->error("cannot deal with non-constant " + av->getName() + " values");
  }
  (av->copyOutput(0))->setConstant();
  return 0;
}

// calculator
void SelectMassCharge::calculate() {
  std::vector<double> masses(1), charges(1), vals(1);
  std::vector<Vector> pos;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
  calculateCV( 0, masses, charges, pos, vals, derivs, virial, this );
  setValue( vals[0] );
}

void SelectMassCharge::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                                    const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                                    std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  if( aa->getName().find("MASSES")!=std::string::npos ) {
    vals[0]=masses[0];
  } else if( aa->chargesWereSet ) {
    vals[0]=charges[0];
  }
}

}
}




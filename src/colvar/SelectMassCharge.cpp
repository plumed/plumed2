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

namespace
{
enum class MC {Mass,Charge};
} // namespace unnamed

template <MC mq>
class SelectMassCharge : public Colvar {
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectMassCharge(const ActionOptions&);
// active methods:
  void calculate() override;
  MULTICOLVAR_DEFAULT(multiColvars::emptyMode);
private:
  AtomNumber theAtom;
};

typedef SelectMassCharge<MC::Mass> SelectMass;
typedef SelectMassCharge<MC::Charge> SelectCharge;
PLUMED_REGISTER_ACTION(SelectMass,"MASS_SCALAR")
PLUMED_REGISTER_ACTION(SelectCharge,"CHARGE_SCALAR")

typedef ColvarShortcut<SelectMass> MassShortcut;
typedef ColvarShortcut<SelectCharge> ChargeShortcut;
PLUMED_REGISTER_ACTION(MassShortcut,"MASS")
PLUMED_REGISTER_ACTION(ChargeShortcut,"CHARGE")

typedef MultiColvarTemplate<SelectMass> MassMulti;
typedef MultiColvarTemplate<SelectCharge> ChargeMulti;
PLUMED_REGISTER_ACTION(MassMulti,"MASS_VECTOR")
PLUMED_REGISTER_ACTION(ChargeMulti,"CHARGE_VECTOR")

template <MC mq>
void SelectMassCharge<mq>::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom number");
  keys.add("atoms","ATOMS","the atom numbers that you would like to store the masses and charges of");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  std::string acname = keys.getDisplayName(); std::size_t und = acname.find("_SCALAR");
  if( und==std::string::npos ) und = acname.find("_VECTOR");
  keys.setDisplayName( acname.substr(0,und) ); keys.setValueDescription("the " + keys.getDisplayName() + " of the atom");
}

template <MC mq>
SelectMassCharge<mq>::SelectMassCharge(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  theAtom=atoms[0];
  /*Modetype mode=*/getModeAndSetupValues(this);
  requestAtoms(atoms);
}

template <MC mq>
void SelectMassCharge<mq>::parseAtomList(  int const num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOM",num,t);
  if( t.size()==1 ) {
    aa->log.printf("  for atom %d\n",t[0].serial());
  } else if( num<0 || t.size()!=0 ) {
    aa->error("Number of specified atoms should be 1");
  }
}

template <MC mq>
typename SelectMassCharge<mq>::Modetype SelectMassCharge<mq>::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  bool constant=true;
  ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( av );
  plumed_assert( aa );
  for(unsigned i=0; i<aa->getNumberOfAtoms(); ++i) {
    std::pair<std::size_t,std::size_t> p = aa->getValueIndices( aa->getAbsoluteIndex(i) );
    if constexpr( mq == MC::Mass ) {
      constant = aa->isMassConstant(p.first);
    } else {
      constant = aa->isChargeConstant(p.first);
    }
  }
  if( !constant ) {
    av->error("cannot deal with non-constant " + av->getName() + " values");
  }
  (av->copyOutput(0))->setConstant();
  return {};
}

// calculator

template <MC mq>
void SelectMassCharge<mq>::calculate() {
  std::vector<Vector> posdummy;
  std::vector<std::vector<Vector> > derivsdummy;
  std::vector<Tensor> virialdummy;

  std::vector<double> massesOrCharges(1);
  if constexpr( mq == MC::Mass ) {
    massesOrCharges[0]=getMass(theAtom.index());
  } else {
    massesOrCharges[0]=getCharge(theAtom.index());
  }
  std::vector<double> vals(1);

  calculateCV( {}, massesOrCharges, massesOrCharges, posdummy, multiColvars::Ouput(vals, derivsdummy, virialdummy), this );

  setValue( vals[0] );
  // does the code above give the same result as doing:
  // if constexpr( mq == MC::Mass ) {
  //   setValue( getMass(theAtom.index) );
  // } else {
  //   setValue( getCharge(theAtom.index) );
  // }
  // ???
  //calculateCV copies the first element from masses or charges into vals[0]
}

template <MC mq>
void SelectMassCharge<mq>::calculateCV( Modetype /*mode*/, const std::vector<double>& masses, const std::vector<double>& charges,
                                        const std::vector<Vector>& pos,multiColvars::Ouput out, const ActionAtomistic* aa ) {
  auto & vals=out.vals();
  if constexpr(mq == MC::Mass)  {
    vals[0]=masses[0];
    // } else if( aa->chargesWereSet ) { this is done in the getModeAndSetupValues by isChargeConstant
  } else {
    vals[0]=charges[0];
  }
}

} // namespace colvar
} // namespace PLMD

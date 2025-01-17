/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionToPutData.h"
#include "core/DomainDecomposition.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR ENERGY
/*
Calculate the total potential energy of the simulation box.

The potential energy can be biased e.g. with umbrella sampling \cite bart-karp98jpcb or with well-tempered metadynamics \cite Bonomi:2009p17935.

Notice that this CV could be unavailable with some MD code. When
it is available, and when also replica exchange is available,
metadynamics applied to ENERGY can be used to decrease the
number of required replicas.

\bug This \ref ENERGY does not include long tail corrections.
Thus when using e.g. LAMMPS `"pair_modify tail yes"` or GROMACS `"DispCorr Ener"` (or `"DispCorr EnerPres"`),
the potential energy from \ref ENERGY will be slightly different form the one of the MD code.
You should still be able to use \ref ENERGY and then reweight your simulation with the correct MD energy value.

\bug Acceptance for replica exchange when \ref ENERGY is biased
is computed correctly only if all the replicas have the same
potential energy function. This is for instance not true when
using GROMACS with lambda replica exchange or with plumed-hrex branch.

\par Examples

The following input instructs plumed to print the energy of the system
\plumedfile
ene: ENERGY
PRINT ARG=ene
\endplumedfile

*/
//+ENDPLUMEDOC


class Energy : public ActionToPutData {
private:
/// This is used to sum the data
  DomainDecomposition* interface;
/// This is the list of forces that must be scaled
  std::vector<ActionToPutData*> forces_to_scale;
public:
  explicit Energy(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  void wait() override;
  void apply() override;
};


PLUMED_REGISTER_ACTION(Energy,"ENERGY")

Energy::Energy(const ActionOptions&ao):
  Action(ao),
  ActionToPutData(ao),
  interface(NULL) {
  plumed.setEnergyValue( getLabel() );
  std::vector<unsigned> shape;
  addValue( shape );
  setNotPeriodic();
  setUnit( "energy", "default" );
  ActionToPutData* px=plumed.getActionSet().selectWithLabel< ActionToPutData*>("posx");
  plumed_assert(px);
  forces_to_scale.push_back(px);
  addDependency( px );
  ActionToPutData* py=plumed.getActionSet().selectWithLabel< ActionToPutData*>("posy");
  plumed_assert(py);
  forces_to_scale.push_back(py);
  addDependency( py );
  ActionToPutData* pz=plumed.getActionSet().selectWithLabel< ActionToPutData*>("posz");
  plumed_assert(pz);
  forces_to_scale.push_back(pz);
  addDependency( pz );
  ActionToPutData* bx=plumed.getActionSet().selectWithLabel< ActionToPutData*>("Box");
  plumed_assert(bx);
  forces_to_scale.push_back(bx);
  addDependency( bx );
  log<<"  Bibliography ";
  log<<plumed.cite("Bartels and Karplus, J. Phys. Chem. B 102, 865 (1998)");
  log<<plumed.cite("Bonomi and Parrinello, J. Comp. Chem. 30, 1615 (2009)");
  log<<"\n";
}

void Energy::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
}

void Energy::wait() {
  if( !interface ) {
    std::vector<DomainDecomposition*> allput=plumed.getActionSet().select<DomainDecomposition*>();
    if( allput.size()>1 ) {
      warning("found more than one interface so don't know how to sum energy");
    }
    interface = allput[0];
  }
  ActionToPutData::wait();
  if( interface ) {
    interface->sumOverDomains( copyOutput(0) );
  }
}

void Energy::apply() {
  if( getPntrToValue()->forcesWereAdded() ) {
    for(unsigned i=0; i<forces_to_scale.size(); ++i) {
      forces_to_scale[i]->rescaleForces( 1.- getPntrToValue()->getForce(0));
    }
  }
}

}
}




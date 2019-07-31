/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR FRET
/*
Calculates the FRET efficiency between a pair of atoms.
The efficiency is calculated using the Forster relation:

\f[
E=\frac{1}{1+(R/R_0)^6}
\f]

where \f$R\f$ is the distance and \f$R_0\f$ is the Forster radius.

By default the distance is computed taking into account periodic
boundary conditions. This behavior can be changed with the NOPBC flag.


\par Examples

The following input tells plumed to print the FRET efficiencies
calculated as a function of the distance between atoms 3 and 5 and
the distance between atoms 2 and 4.
\plumedfile
fe1:  FRET ATOMS=3,5 R0=5.5
fe2:  FRET ATOMS=2,4 R0=5.5
PRINT ARG=fe1,fe2
\endplumedfile

The following input computes the FRET efficiency calculated on the
terminal atoms of a polymer
of 100 atoms and keeps it at a value around 0.5.
\plumedfile
WHOLEMOLECULES ENTITY0=1-100
fe: FRET ATOMS=1,100 R0=5.5 NOPBC
RESTRAINT ARG=fe KAPPA=100 AT=0.5
\endplumedfile

Notice that NOPBC is used
to be sure that if the distance is larger than half the simulation
box the distance is compute properly. Also notice that, since many MD
codes break molecules across cell boundary, it might be necessary to
use the \ref WHOLEMOLECULES keyword (also notice that it should be
_before_ FRET).
Just be sure that the ordered list provide to WHOLEMOLECULES has the following
properties:
- Consecutive atoms should be closer than half-cell throughout the entire simulation.
- Atoms required later for the distance (e.g. 1 and 100) should be included in the list

*/
//+ENDPLUMEDOC

class FretEfficiency : public Colvar {
  bool pbc;
  double R0_;

public:
  static void registerKeywords( Keywords& keys );
  explicit FretEfficiency(const ActionOptions&);
// active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(FretEfficiency,"FRET")

void FretEfficiency::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.add("compulsory","R0","The value of the Forster radius.");
}

FretEfficiency::FretEfficiency(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=2)
    error("Number of specified atoms should be 2");
  parse("R0",R0_);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
  log.printf("  with Forster radius set to %lf\n",R0_);

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  log << " Bibliography" << plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)") << "\n";

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atoms);
}


// calculator
void FretEfficiency::calculate() {

  if(pbc) makeWhole();

  Vector distance=delta(getPosition(0),getPosition(1));
  const double dist_mod=distance.modulo();
  const double inv_dist_mod=1.0/dist_mod;

  const double ratiosix=pow(dist_mod/R0_,6);
  const double fret_eff = 1.0/(1.0+ratiosix);

  const double der = -6.0*fret_eff*fret_eff*ratiosix*inv_dist_mod;

  setAtomsDerivatives(0,-inv_dist_mod*der*distance);
  setAtomsDerivatives(1, inv_dist_mod*der*distance);
  setBoxDerivativesNoPbc();
  setValue(fret_eff);

}

}
}




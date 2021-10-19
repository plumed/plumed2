/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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

//+PLUMEDOC COLVAR DIPOLE 
/*
Calculate the dipole moments for groups of atoms.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecule with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\plumedfile
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
\endplumedfile
  
\attention
If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
where N is the number of atoms. This implies that the dipole (which for a charged system depends
on the position) is computed on the geometric center of the group.


*/
//+ENDPLUMEDOC

class Dipole : public MultiColvarBase {
private:
  bool components;
public:
  static void registerKeywords( Keywords& keys );
  explicit Dipole(const ActionOptions&);
// active methods:
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(Dipole,"DIPOLE")

void Dipole::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addOutputComponent("x","COMPONENTS","the x-components of the dipoles");
  keys.addOutputComponent("y","COMPONENTS","the y-components of the dipoles");
  keys.addOutputComponent("z","COMPONENTS","the z-components of the dipoles");
}

Dipole::Dipole(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  components(false)
{
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components) {
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
    log<<"  WARNING: components will not have the proper periodicity - see manual\n";
  } else {
    addValueWithDerivatives(); setNotPeriodic();
  }
}


// calculator
void Dipole::compute( const std::vector<Vector>& pos, MultiValue& myvals ) const {
  double ctot=0.; unsigned N=pos.size();
  std::vector<double> charges(N);
  Vector dipje;
  
  for(unsigned i=0; i<N; ++i) {
      charges[i]=getCharge(i); ctot+=charges[i];
  }
  ctot/=(double)N;

  for(unsigned i=0; i<N; ++i) {
      charges[i]-=ctot; dipje += charges[i]*pos[i];
  }

  if(!components) {
    double dipole = dipje.modulo(), idip = 1./dipole;
    for(unsigned i=0; i<N; i++) addAtomsDerivatives(0, i, charges[i]*idip*dipje, myvals );
    setBoxDerivativesNoPbc( 0, pos, myvals );
    setValue(0, dipole, myvals );
  } else {
    for(unsigned i=0; i<N; i++) {
        addAtomsDerivatives(0, i, charges[i]*Vector(1.0,0.0,0.0), myvals );
        addAtomsDerivatives(1, i, charges[i]*Vector(0.0,1.0,0.0), myvals );
        addAtomsDerivatives(2, i, charges[i]*Vector(0.0,0.0,1.0), myvals );
    }
    for(unsigned i=0;i<3;++i) {
        setValue(i, dipje[i], myvals); setBoxDerivativesNoPbc( i, pos, myvals );
    }
  }
}

}
}




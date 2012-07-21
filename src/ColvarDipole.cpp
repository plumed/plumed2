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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DIPOLE 
/*
Calcualte the dipole moment for a group of atoms.

\par Examples
The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10.
\verbatim
DIPOLE GROUP=1-10
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarDipole : public Colvar {
  vector<AtomNumber> ga_lista;
public:
  ColvarDipole(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarDipole,"DIPOLE")

void ColvarDipole::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.remove("NOPBC");
}

ColvarDipole::ColvarDipole(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  parseAtomList("GROUP",ga_lista);
  checkRead();
  addValueWithDerivatives(); setNotPeriodic();

  log.printf("  of %d atoms\n",ga_lista.size());
  for(unsigned int i=0;i<ga_lista.size();++i){
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  requestAtoms(ga_lista);
}

// calculator
void ColvarDipole::calculate()
{
 double dipole=0.;
 Tensor virial;
 vector<Vector> deriv(getNumberOfAtoms());
 Vector dipje;

// deriv.resize(getPositions().size());
// deriv.resize(getNumberOfAtoms());
 for(unsigned int i=0;i<ga_lista.size();i++) {
   dipje += (getCharge(i))*getPosition(i);
 }
 dipole = dipje.modulo();

 for(unsigned int i=0;i<ga_lista.size();i++) {
   double dfunc=getCharge(i)/dipole;
   deriv[i] = deriv[i] + (dfunc)*dipje;
   virial=virial-Tensor(getPosition(i),deriv[i]);
 }

// for(unsigned i=0;i<getPositions().size();++i) setAtomsDerivatives(i,deriv[i]);
 for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (dipole);
 setBoxDerivatives  (virial);
}

}

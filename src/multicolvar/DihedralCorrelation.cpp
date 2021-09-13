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
#include "MultiColvarBase.h"
#include "tools/Torsion.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR DIHCOR
/*
Measures the degree of similarity between dihedral angles.

This colvar calculates the following quantity.

\f[
s = \frac{1}{2} \sum_i \left[ 1 + \cos( \phi_i - \psi_i ) \right]
\f]

where the \f$\phi_i\f$ and \f$\psi\f$ values and the instantaneous values for the \ref TORSION angles of interest.

\par Examples

The following provides an example input for the DIHCOR action

\plumedfile
DIHCOR ...
  ATOMS1=1,2,3,4,5,6,7,8
  ATOMS2=5,6,7,8,9,10,11,12
  LABEL=dih
... DIHCOR
PRINT ARG=dih FILE=colvar STRIDE=10
\endplumedfile

In the above input we are calculating the correlation between the torsion angle involving atoms 1, 2, 3 and 4 and the torsion angle
involving atoms 5, 6, 7 and 8.	This is then added to the correlation between the torsion angle involving atoms 5, 6, 7 and 8 and the
correlation angle involving atoms 9, 10, 11 and 12.

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the \ref MOLINFO command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
dih: DIHCOR ...
ATOMS1=@phi-3,@psi-3
ATOMS2=@psi-3,@phi-4
ATOMS3=@phi-4,@psi-4
...
PRINT ARG=dih FILE=colvar STRIDE=10
\endplumedfile

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the fourth residue of the protein.

*/
//+ENDPLUMEDOC

class DihedralCorrelation : public MultiColvarBase {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit DihedralCorrelation(const ActionOptions&);
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(DihedralCorrelation,"DIHEDRAL_CORRELATION")

void DihedralCorrelation::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
}

DihedralCorrelation::DihedralCorrelation(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  addValueWithDerivatives(); setNotPeriodic(); checkRead();
}

void DihedralCorrelation::compute( const std::vector<Vector>& pos, MultiValue& myvals ) const {
  const Vector d10=getSeparation(pos[1],pos[0]);
  const Vector d11=getSeparation(pos[2],pos[1]);
  const Vector d12=getSeparation(pos[3],pos[2]);

  Vector dd10,dd11,dd12;
  PLMD::Torsion t1;
  const double phi1  = t1.compute(d10,d11,d12,dd10,dd11,dd12);

  const Vector d20=getSeparation(pos[5],pos[4]);
  const Vector d21=getSeparation(pos[6],pos[5]);
  const Vector d22=getSeparation(pos[7],pos[6]);

  Vector dd20,dd21,dd22;
  PLMD::Torsion t2;
  const double phi2 = t2.compute( d20, d21, d22, dd20, dd21, dd22 );

  // Calculate value
  const double diff = phi2 - phi1;
  const double value = 0.5*(1.+std::cos(diff));
  // Derivatives wrt phi1
  const double dval = 0.5*std::sin(diff);
  dd10 *= dval;
  dd11 *= dval;
  dd12 *= dval;
  // And add
  addAtomsDerivatives(0, 0, dd10, myvals );
  addAtomsDerivatives(0, 1, dd11-dd10, myvals );
  addAtomsDerivatives(0, 2, dd12-dd11, myvals );
  addAtomsDerivatives(0, 3, -dd12, myvals );
  addBoxDerivatives  (0, -(extProduct(d10,dd10)+extProduct(d11,dd11)+extProduct(d12,dd12)), myvals);
  // Derivative wrt phi2
  dd20 *= -dval;
  dd21 *= -dval;
  dd22 *= -dval;
  // And add
  addAtomsDerivatives(0, 4, dd20, myvals );
  addAtomsDerivatives(0, 5, dd21-dd20, myvals);
  addAtomsDerivatives(0, 6, dd22-dd21, myvals);
  addAtomsDerivatives(0, 7, -dd22, myvals);
  addBoxDerivatives(0, -(extProduct(d20,dd20)+extProduct(d21,dd21)+extProduct(d22,dd22)), myvals);
  setValue( 0, value, myvals );
}

// We have a little helper class here to ensure that we actually do what is required by this action
class DihedralCorrelationShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  DihedralCorrelationShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(DihedralCorrelationShortcut,"DIHCOR")

void DihedralCorrelationShortcut::registerKeywords( Keywords& keys ) {
  DihedralCorrelation::registerKeywords( keys );
}

DihedralCorrelationShortcut::DihedralCorrelationShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  readInputLine( getShortcutLabel() +"_data: DIHEDRAL_CORRELATION " + convertInputLineToString() ); 
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_data PERIODIC=NO");
} 

}
}

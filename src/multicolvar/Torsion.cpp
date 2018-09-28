/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2018 The plumed team
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
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR TORSION
/*
Calculate a torsional angle.

This command can be used to compute the torsion between four atoms or alternatively
to calculate the angle between two vectors projected on the plane
orthogonal to an axis.

\par Examples

This input tells plumed to print the torsional angle between atoms 1, 2, 3 and 4
on file COLVAR.
\plumedfile
t: TORSION ATOMS=1,2,3,4
# this is an alternative, equivalent, definition:
# t: TORSION VECTOR1=2,1 AXIS=2,3 VECTOR2=3,4
PRINT ARG=t FILE=COLVAR
\endplumedfile

If you are working with a protein you can specify the special named torsion angles \f$\phi\f$, \f$\psi\f$, \f$\omega\f$ and \f$\chi_1\f$
by using TORSION in combination with the \ref MOLINFO command.  This can be done by using the following
syntax.

\plumedfile
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
\endplumedfile

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.
*/
//+ENDPLUMEDOC

class Torsion : public MultiColvarBase {
  bool do_cosine;
public:
  static void registerKeywords( Keywords& keys );
  explicit Torsion(const ActionOptions&);
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(Torsion,"TORSION")

void Torsion::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("numbered","AXIS","two atoms that define an axis.  You can use this to find the angle in the plane perpendicular to the axis between the vectors specified using the VECTOR1 and VECTOR2 keywords.");
  keys.add("numbered","VECTORA","two atoms that define a vector.  You can use this in combination with VECTORB and AXIS");
  keys.add("numbered","VECTORB","two atoms that define a vector.  You can use this in combination with VECTORA and AXIS");
  keys.addFlag("COSINE",false,"calculate cosine instead of dihedral");
  keys.reset_style("AXIS","atoms-2"); keys.reset_style("VECTORA","atoms-2"); keys.reset_style("VECTORB","atoms-2"); 
}

Torsion::Torsion(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  do_cosine(false)
{
  parseFlag("COSINE",do_cosine);
  if(do_cosine) log.printf("  calculating cosine instead of torsion\n");

  addValueWithDerivatives(); setPeriodic( "-pi", "pi" ); checkRead();
}

void Torsion::compute( const std::vector<Vector>& pos, MultiValue& myvals ) const {
  plumed_assert( pos.size()==6 );
  const Vector d0=delta(pos[1],pos[0]);
  const Vector d1=delta(pos[3],pos[2]);
  const Vector d2=delta(pos[5],pos[4]);

  Vector dd0,dd1,dd2; PLMD::Torsion t;
  double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);
  if(do_cosine) {
    dd0 *= -sin(value);
    dd1 *= -sin(value);
    dd2 *= -sin(value);
    value = cos(value);
  }
  addAtomsDerivatives(0, 0,  dd0, myvals);
  addAtomsDerivatives(0, 1, -dd0, myvals );
  addAtomsDerivatives(0, 2,  dd1, myvals);
  addAtomsDerivatives(0, 3, -dd1, myvals);
  addAtomsDerivatives(0, 4,  dd2, myvals);
  addAtomsDerivatives(0, 5, -dd2, myvals);

  addBoxDerivatives (0, -(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)), myvals);
  setValue( 0, value, myvals );
}

}
}

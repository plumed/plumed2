/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "AtomValuePack.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

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
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(DihedralCorrelation,"DIHCOR")

void DihedralCorrelation::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the dihedral correlation values you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one dihedral correlation will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of 8 atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
}

DihedralCorrelation::DihedralCorrelation(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword( "ATOMS", 8, all_atoms );
  setupMultiColvarBase( all_atoms );
  // Stuff for central atoms
  std::vector<bool> catom_ind(8, false);
  catom_ind[1]=catom_ind[2]=catom_ind[5]=catom_ind[6]=true;
  setAtomsForCentralAtom( catom_ind );

  // And setup the ActionWithVessel
  if( getNumberOfVessels()==0 ) {
    std::string fake_input;
    addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
    readVesselKeywords();  // This makes sure resizing is done
  }

  // And check everything has been read in correctly
  checkRead();
}

double DihedralCorrelation::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  const Vector d10=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  const Vector d11=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  const Vector d12=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd10,dd11,dd12;
  PLMD::Torsion t1;
  const double phi1  = t1.compute(d10,d11,d12,dd10,dd11,dd12);

  const Vector d20=getSeparation(myatoms.getPosition(5),myatoms.getPosition(4));
  const Vector d21=getSeparation(myatoms.getPosition(6),myatoms.getPosition(5));
  const Vector d22=getSeparation(myatoms.getPosition(7),myatoms.getPosition(6));

  Vector dd20,dd21,dd22;
  PLMD::Torsion t2;
  const double phi2 = t2.compute( d20, d21, d22, dd20, dd21, dd22 );

  // Calculate value
  const double diff = phi2 - phi1;
  const double value = 0.5*(1.+cos(diff));
  // Derivatives wrt phi1
  const double dval = 0.5*sin(diff);
  dd10 *= dval;
  dd11 *= dval;
  dd12 *= dval;
  // And add
  addAtomDerivatives(1, 0, dd10, myatoms );
  addAtomDerivatives(1, 1, dd11-dd10, myatoms );
  addAtomDerivatives(1, 2, dd12-dd11, myatoms );
  addAtomDerivatives(1, 3, -dd12, myatoms );
  myatoms.addBoxDerivatives  (1, -(extProduct(d10,dd10)+extProduct(d11,dd11)+extProduct(d12,dd12)));
  // Derivative wrt phi2
  dd20 *= -dval;
  dd21 *= -dval;
  dd22 *= -dval;
  // And add
  addAtomDerivatives(1, 4, dd20, myatoms );
  addAtomDerivatives(1, 5, dd21-dd20, myatoms );
  addAtomDerivatives(1, 6, dd22-dd21, myatoms );
  addAtomDerivatives(1, 7, -dd22, myatoms );
  myatoms.addBoxDerivatives(1, -(extProduct(d20,dd20)+extProduct(d21,dd21)+extProduct(d22,dd22)));

  return value;
}

}
}

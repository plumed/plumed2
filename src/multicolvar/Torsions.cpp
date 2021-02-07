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
#include "AtomValuePack.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR TORSIONS
/*
Calculate whether or not a set of torsional angles are within a particular range.

\par Examples

The following provides an example of the input for the TORSIONS command

\plumedfile
TORSIONS ...
ATOMS1=168,170,172,188
ATOMS2=170,172,188,190
ATOMS3=188,190,192,230
BETWEEN={GAUSSIAN LOWER=0 UPPER=pi SMEAR=0.1}
LABEL=ab
... TORSIONS
PRINT ARG=ab.* FILE=colvar STRIDE=10
\endplumedfile

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the \ref MOLINFO command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
TORSIONS ...
ATOMS1=@phi-3
ATOMS2=@psi-3
ATOMS3=@phi-4
BETWEEN={GAUSSIAN LOWER=0 UPPER=pi SMEAR=0.1}
LABEL=ab
... TORSIONS
PRINT ARG=ab.* FILE=colvar STRIDE=10
\endplumedfile

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the fourth residue of the protein.


*/
//+ENDPLUMEDOC

class Torsions : public MultiColvarBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit Torsions(const ActionOptions&);
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
  bool isPeriodic() override { return true; }
  void retrieveDomain( std::string& min, std::string& max ) override { min="-pi"; max="pi"; }
};

PLUMED_REGISTER_ACTION(Torsions,"TORSIONS")

void Torsions::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the torsion angles you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one torsion will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "provide the indices of four atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.use("BETWEEN"); keys.use("HISTOGRAM");
}

Torsions::Torsions(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  int natoms=4; std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword( "ATOMS", natoms, all_atoms );
  setupMultiColvarBase( all_atoms );
  std::vector<bool> catom_ind(4, false);
  catom_ind[1]=catom_ind[2]=true;
  setAtomsForCentralAtom( catom_ind );
  // Read in the vessels
  readVesselKeywords();
  // And check everything has been read in correctly
  checkRead();
}

double Torsions::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector d0,d1,d2;
  d0=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  d1=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  d2=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd0,dd1,dd2; PLMD::Torsion t;
  double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);

  addAtomDerivatives(1,0,dd0,myatoms);
  addAtomDerivatives(1,1,dd1-dd0,myatoms);
  addAtomDerivatives(1,2,dd2-dd1,myatoms);
  addAtomDerivatives(1,3,-dd2,myatoms);

  myatoms.addBoxDerivatives  (1, -(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));

  return value;
}

}
}

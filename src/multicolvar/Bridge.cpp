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
#include "AtomValuePack.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR BRIDGE
/*
Calculate the number of atoms that bridge two parts of a structure

This quantity calculates:

\f[
 f(x) = \sum_{ijk} s_A(r_{ij})s_B(r_{ik})
\f]

where the sum over \f$i\f$ is over all the ``bridging atoms" and
\f$s_A\f$ and \f$s_B\f$ are \ref switchingfunction.

\par Examples

The following example instructs plumed to calculate the number of water molecules
that are bridging between atoms 1-10 and atoms 11-20 and to print the value
to a file

\plumedfile
w1: BRIDGE BRIDGING_ATOMS=100-200 GROUPA=1-10 GROUPB=11-20 SWITCH={RATIONAL R_0=0.2}
PRINT ARG=w1 FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class Bridge : public MultiColvarBase {
private:
  Vector dij, dik;
  SwitchingFunction sf1;
  SwitchingFunction sf2;
public:
  static void registerKeywords( Keywords& keys );
  explicit Bridge(const ActionOptions&);
// active methods:
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(Bridge,"BRIDGE")

void Bridge::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("atoms-2","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("atoms-2","GROUPA","The list of atoms that are in the first interesting part of the structure");
  keys.add("atoms-2","GROUPB","The list of atoms that are in the second interesting part of the structure");
  keys.add("optional","SWITCH","The parameters of the two \\ref switchingfunction in the above formula");
  keys.add("optional","SWITCHA","The \\ref switchingfunction on the distance between bridging atoms and the atoms in "
           "group A");
  keys.add("optional","SWITCHB","The \\ref switchingfunction on the distance between the bridging atoms and the atoms in "
           "group B");
}

Bridge::Bridge(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readThreeGroups("GROUPA","GROUPB","BRIDGING_ATOMS",false,true,all_atoms);
  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms );
  // Setup Central atom atoms
  std::vector<bool> catom_ind(3, false); catom_ind[0]=true;
  setAtomsForCentralAtom( catom_ind );

  std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    sf1.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    sf2.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    parse("SWITCHA",sfinput);
    if(sfinput.length()>0) {
      weightHasDerivatives=true;
      sf1.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCHA keyword : " + errors );
      sfinput.clear(); parse("SWITCHB",sfinput);
      if(sfinput.length()==0) error("found SWITCHA keyword without SWITCHB");
      sf2.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCHB keyword : " + errors );
    } else {
      error("missing definition of switching functions");
    }
  }
  log.printf("  distance between bridging atoms and atoms in GROUPA must be less than %s\n",sf1.description().c_str());
  log.printf("  distance between bridging atoms and atoms in GROUPB must be less than %s\n",sf2.description().c_str());

  // Setup link cells
  setLinkCellCutoff( sf1.get_dmax() + sf2.get_dmax() );

  // And setup the ActionWithVessel
  if( getNumberOfVessels()!=0 ) error("should not have vessels for this action");
  std::string fake_input;
  addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
  readVesselKeywords();
  // And check everything has been read in correctly
  checkRead();
}

double Bridge::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  double tot=0;
  for(unsigned i=2; i<myatoms.getNumberOfAtoms(); ++i) {
    Vector dij=getSeparation( myatoms.getPosition(i), myatoms.getPosition(0) );
    double dw1, w1=sf1.calculateSqr( dij.modulo2(), dw1 );
    Vector dik=getSeparation( myatoms.getPosition(i), myatoms.getPosition(1) );
    double dw2, w2=sf2.calculateSqr( dik.modulo2(), dw2 );

    tot += w1*w2;
    // And finish the calculation
    addAtomDerivatives( 1, 0,  w2*dw1*dij, myatoms );
    addAtomDerivatives( 1, 1,  w1*dw2*dik, myatoms );
    addAtomDerivatives( 1, i, -w1*dw2*dik-w2*dw1*dij, myatoms );
    myatoms.addBoxDerivatives( 1, w1*(-dw2)*Tensor(dik,dik)+w2*(-dw1)*Tensor(dij,dij) );
  }
  return tot;
}

}
}

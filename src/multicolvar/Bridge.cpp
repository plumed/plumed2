/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "MultiColvar.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

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
that are bridging betweeen atoms 1-10 and atoms 11-20 and to print the value 
to a file

\verbatim
BRIDGE BRIDGING_ATOMS=100-200 GROUPA=1-10 GROUPB=11-20 LABEL=w1
PRINT ARG=a1.mean FILE=colvar 
\endverbatim

*/
//+ENDPLUMEDOC

class Bridge : public MultiColvar {
private:
  double rcut2;
  Vector dij, dik;
  SwitchingFunction sf1;
  SwitchingFunction sf2;
public:
  static void registerKeywords( Keywords& keys );
  explicit Bridge(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
  void calculateWeight( AtomValuePack& myatoms ) const ;
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Bridge,"BRIDGE")

void Bridge::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS"); 
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
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  weightHasDerivatives=true; std::vector<AtomNumber> all_atoms;
  readThreeGroups("BRIDGING_ATOMS","GROUPA","GROUPB",false, all_atoms);
  if( all_atoms.size()>0 ) ActionAtomistic::requestAtoms( all_atoms );
  // Setup the multicolvar base
  setupMultiColvarBase();
  // Setup Central atom atoms
  std::vector<bool> catom_ind(3, false); catom_ind[0]=true;
  setAtomsForCentralAtom( catom_ind ); 

  std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ){
      sf1.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
      sf2.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );   
  } else {
      parse("SWITCHA",sfinput); 
      if(sfinput.length()>0){
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
  setLinkCellCutoff( sf1.get_dmax() );
  rcut2 = sf1.get_dmax()*sf1.get_dmax();

  // And setup the ActionWithVessel
  if( getNumberOfVessels()!=0 ) error("should not have vessels for this action");
  std::string fake_input;
  addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
  readVesselKeywords();
  // And check everything has been read in correctly
  checkRead();
}

void Bridge::calculateWeight( AtomValuePack& myatoms ) const {
  Vector dij=getSeparation( myatoms.getPosition(0), myatoms.getPosition(2) );
  double ldij = dij.modulo2();
  if( ldij>rcut2 ) { myatoms.setValue(0,0); return; }
  double dw, w=sf2.calculateSqr( ldij, dw );
  myatoms.setValue( 0, w );

  myatoms.addAtomsDerivatives( 0, 0, -dw*dij );
  myatoms.addAtomsDerivatives( 0, 2, dw*dij );
  myatoms.addBoxDerivatives( 0, (-dw)*Tensor(dij,dij) );
}

double Bridge::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector dik=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dw, w=sf1.calculateSqr( dik.modulo2(), dw );

  // And finish the calculation
  myatoms.addAtomsDerivatives( 1, 0, -dw*dik );
  myatoms.addAtomsDerivatives( 1, 1,  dw*dik );
  myatoms.addBoxDerivatives( 1, (-dw)*Tensor(dik,dik) );
  return w;
}

}
}

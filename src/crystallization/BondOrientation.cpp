/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"
#include "multicolvar/AtomValuePack.h"
#include "VectorMultiColvar.h"

//+PLUMEDOC MCOLVAR BOND_DIRECTIONS
/*
Calculate the vectors connecting atoms that are within cutoff defined using a switching function.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class BondOrientation : public VectorMultiColvar {
private:
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit BondOrientation( const ActionOptions& ao );
  double calculateWeight( const unsigned& current, const double& weight, multicolvar::AtomValuePack& myatoms ) const override;
  void calculateVector( multicolvar::AtomValuePack& myatoms ) const override;
};

PLUMED_REGISTER_ACTION(BondOrientation,"BOND_DIRECTIONS")

void BondOrientation::registerKeywords( Keywords& keys ) {
  VectorMultiColvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the vectors you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one vector will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of two atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
           "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
           "in GROUPB. This must be used in conjunction with GROUPA.");
  keys.add("compulsory","NN","12","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.use("VMEAN"); keys.use("VSUM");
}

BondOrientation::BondOrientation( const ActionOptions& ao ):
  Action(ao),
  VectorMultiColvar(ao)
{
  // Read atoms
  weightHasDerivatives=true;
  std::vector<AtomNumber> all_atoms;
  readTwoGroups( "GROUP", "GROUPA", "GROUPB", all_atoms );
  if( atom_lab.size()==0 ) readAtomsLikeKeyword( "ATOMS", 2, all_atoms );
  setupMultiColvarBase( all_atoms );
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  orientation of those bonds with lengths are less than %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  double rcut = switchingFunction.get_dmax(); rcut2 = rcut*rcut;
  // Set the dimensionality of the vectors
  setVectorDimensionality(3);
}

double BondOrientation::calculateWeight( const unsigned& current, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) ); double distm=distance.modulo2();
  if( distm>rcut2 ) return 0.0;
  double df, ww=switchingFunction.calculateSqr( distm, df );
  // Derivatives of weights
  addAtomDerivatives( 0, 0, -df*weight*distance, myatoms );
  addAtomDerivatives( 0, 1, df*weight*distance, myatoms );
  myatoms.addBoxDerivatives( 0, (-df)*weight*Tensor(distance,distance) );
  return ww;
}

void BondOrientation::calculateVector( multicolvar::AtomValuePack& myatoms ) const {
  Vector distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );

  addAtomDerivatives( 2, 0, Vector(-1.0,0,0), myatoms );
  addAtomDerivatives( 2, 1, Vector(+1.0,0,0), myatoms );
  myatoms.addBoxDerivatives( 2, Tensor(distance,Vector(-1.0,0,0)) );
  myatoms.addValue( 2, distance[0] );

  addAtomDerivatives( 3, 0, Vector(0,-1.0,0), myatoms );
  addAtomDerivatives( 3, 1, Vector(0,+1.0,0), myatoms );
  myatoms.addBoxDerivatives( 3, Tensor(distance,Vector(0,-1.0,0)) );
  myatoms.addValue( 3, distance[1] );

  addAtomDerivatives( 4, 0, Vector(0,0,-1.0), myatoms );
  addAtomDerivatives( 4, 1, Vector(0,0,+1.0), myatoms );
  myatoms.addBoxDerivatives( 4, Tensor(distance,Vector(0,0,-1.0)) );
  myatoms.addValue( 4, distance[2] );
}

}
}

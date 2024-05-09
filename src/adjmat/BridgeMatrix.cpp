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
#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace adjmat {

//+PLUMEDOC MCOLVAR BRIDGE_MATRIX
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

class BridgeMatrix : public AdjacencyMatrixBase {
private:
  Vector dij, dik;
  SwitchingFunction sf1;
  SwitchingFunction sf2;
public:
  static void registerKeywords( Keywords& keys );
  explicit BridgeMatrix(const ActionOptions&);
// active methods:
  double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(BridgeMatrix,"BRIDGE_MATRIX")

void BridgeMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("optional","SWITCH","The parameters of the two switchingfunction in the above formula");
  keys.add("optional","SWITCHA","The switchingfunction on the distance between bridging atoms and the atoms in "
           "group A");
  keys.add("optional","SWITCHB","The switchingfunction on the distance between the bridging atoms and the atoms in "
           "group B");
}

BridgeMatrix::BridgeMatrix(const ActionOptions&ao):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  bool oneswitch; std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    sf1.set(sfinput,errors); oneswitch=true;
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    sf2.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    parse("SWITCHA",sfinput);
    if(sfinput.length()>0) {
      sf1.set(sfinput,errors); oneswitch=false;
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
  setLinkCellCutoff( oneswitch, sf1.get_dmax() + sf2.get_dmax() );

  // And check everything has been read in correctly
  checkRead();
}

double BridgeMatrix::calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const {
  double tot=0; if( pos2.modulo2()<epsilon ) return 0.0;
  for(unsigned i=0; i<natoms; ++i) {
    Vector dij= getPosition(i,myvals); double dijm = dij.modulo2();
    double dw1, w1=sf1.calculateSqr( dijm, dw1 ); if( dijm<epsilon ) { w1=0.0; dw1=0.0; }
    Vector dik=pbcDistance( getPosition(i,myvals), pos2 ); double dikm=dik.modulo2();
    double dw2, w2=sf2.calculateSqr( dikm, dw2 ); if( dikm<epsilon ) { w2=0.0; dw2=0.0; }

    tot += w1*w2;
    // And finish the calculation
    addAtomDerivatives( 0,  -w2*dw1*dij, myvals );
    addAtomDerivatives( 1,  w1*dw2*dik, myvals );
    addThirdAtomDerivatives( i, -w1*dw2*dik+w2*dw1*dij, myvals );
    addBoxDerivatives( w1*(-dw2)*Tensor(dik,dik)+w2*(-dw1)*Tensor(dij,dij), myvals );
  }
  return tot;
}

}
}

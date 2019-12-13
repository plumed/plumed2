/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "core/ActionRegister.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization {

//+PLUMEDOC MCOLVAR PLANES
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule.

\par Examples


*/
//+ENDPLUMEDOC


class MoleculePlane : public VectorMultiColvar {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit MoleculePlane( const ActionOptions& ao );
  AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const override;
  void calculateVector( multicolvar::AtomValuePack& myatoms ) const override;
};

PLUMED_REGISTER_ACTION(MoleculePlane,"PLANES")

void MoleculePlane::registerKeywords( Keywords& keys ) {
  VectorMultiColvar::registerKeywords( keys ); keys.use("VMEAN");
  keys.add("numbered","MOL","The numerical indices of the atoms in the molecule. If three atoms are specified the orientation "
           "of the molecule is taken as the normal to the plane containing the vector connecting the first and "
           "second atoms and the vector connecting the second and third atoms.  If four atoms are specified the "
           "orientation of the molecule is taken as the normal to the plane containing the vector connecting the "
           "first and second atoms and the vector connecting the third and fourth atoms. The molecule is always "
           "assumed to lie at the geometric center for the three/four atoms.");
  keys.reset_style("MOL","atoms");
}

MoleculePlane::MoleculePlane( const ActionOptions& ao ):
  Action(ao),
  VectorMultiColvar(ao)
{
  std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword("MOL",-1,all_atoms);
  if( ablocks.size()!=3 && ablocks.size()!=4 ) error("number of atoms in molecule specification is wrong.  Should be three or four.");

  if( all_atoms.size()==0 ) error("No atoms were specified");
  setVectorDimensionality( 3 ); setupMultiColvarBase( all_atoms );
}

AtomNumber MoleculePlane::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_lab.size() );
  plumed_assert( atom_lab[iatom].first==0 );
  return ActionAtomistic::getAbsoluteIndex( ablocks[0][atom_lab[iatom].second] );
}

void MoleculePlane::calculateVector( multicolvar::AtomValuePack& myatoms ) const {
  Vector d1, d2, cp;
  if( myatoms.getNumberOfAtoms()==3 ) {
    d1=getSeparation( myatoms.getPosition(1), myatoms.getPosition(0) );
    d2=getSeparation( myatoms.getPosition(1), myatoms.getPosition(2) );
  } else {
    d1=getSeparation( myatoms.getPosition(1), myatoms.getPosition(0) );
    d2=getSeparation( myatoms.getPosition(2), myatoms.getPosition(3) );
  }
  cp = crossProduct( d1, d2 );

  addAtomDerivatives( 2, 0, crossProduct( Vector(-1.0,0,0), d2 ), myatoms );
  if( myatoms.getNumberOfAtoms()==3 ) {
    addAtomDerivatives( 2, 1, crossProduct( Vector(+1.0,0,0), d2 ) + crossProduct( Vector(-1.0,0,0), d1 ), myatoms );
    addAtomDerivatives( 2, 2, crossProduct( Vector(+1.0,0,0), d1 ), myatoms );
  } else {
    addAtomDerivatives( 2, 1, crossProduct( Vector(+1.0,0,0), d2 ), myatoms );
    addAtomDerivatives( 2, 2, crossProduct( Vector(-1.0,0,0), d1 ), myatoms );
    addAtomDerivatives( 2, 3, crossProduct( Vector(+1.0,0,0), d1 ), myatoms );
  }
  myatoms.addBoxDerivatives( 2, Tensor(d1,crossProduct(Vector(+1.0,0,0), d2)) + Tensor( d2, crossProduct(Vector(-1.0,0,0), d1)) );
  myatoms.addValue( 2, cp[0] );

  addAtomDerivatives( 3, 0, crossProduct( Vector(0,-1.0,0), d2 ), myatoms );
  if( myatoms.getNumberOfAtoms()==3 ) {
    addAtomDerivatives( 3, 1, crossProduct( Vector(0,+1.0,0), d2 ) + crossProduct( Vector(0,-1.0,0), d1 ), myatoms );
    addAtomDerivatives( 3, 2, crossProduct( Vector(0,+1.0,0), d1 ), myatoms );
  } else {
    addAtomDerivatives( 3, 1, crossProduct( Vector(0,+1.0,0), d2 ), myatoms );
    addAtomDerivatives( 3, 2, crossProduct( Vector(0,-1.0,0), d1 ), myatoms );
    addAtomDerivatives( 3, 3, crossProduct( Vector(0,+1.0,0), d1 ), myatoms );
  }
  myatoms.addBoxDerivatives( 3, Tensor(d1,crossProduct(Vector(0,+1.0,0), d2)) + Tensor( d2, crossProduct(Vector(0,-1.0,0), d1)) );
  myatoms.addValue( 3, cp[1] );

  addAtomDerivatives( 4, 0, crossProduct( Vector(0,0,-1.0), d2 ), myatoms );
  if( myatoms.getNumberOfAtoms()==3 ) {
    addAtomDerivatives( 4, 1, crossProduct( Vector(0,0,+1.0), d2 ) + crossProduct( Vector(0,0,-1.0), d1 ), myatoms );
    addAtomDerivatives( 4, 2, crossProduct( Vector(0,0,+1.0), d1 ), myatoms);
  } else {
    addAtomDerivatives( 4, 1, crossProduct( Vector(0,0,-1.0), d2 ), myatoms);
    addAtomDerivatives( 4, 2, crossProduct( Vector(0,0,-1.0), d1 ), myatoms);
    addAtomDerivatives( 4, 3, crossProduct( Vector(0,0,+1.0), d1 ), myatoms);
  }
  myatoms.addBoxDerivatives( 4, Tensor(d1,crossProduct(Vector(0,0,+1.0), d2)) + Tensor( d2, crossProduct(Vector(0,0,-1.0), d1)) );
  myatoms.addValue( 4, cp[2] );
}

// Vector MoleculePlane::getCentralAtom(){
//   Vector com; com.zero();
//   if( getNAtoms()==3 ){
//       com+=(1.0/3.0)*getPosition(0);
//       com+=(1.0/3.0)*getPosition(1);
//       com+=(1.0/3.0)*getPosition(2);
//       addCentralAtomDerivatives( 0, (1.0/3.0)*Tensor::identity() );
//       addCentralAtomDerivatives( 1, (1.0/3.0)*Tensor::identity() );
//       addCentralAtomDerivatives( 2, (1.0/3.0)*Tensor::identity() );
//       return com;
//   }
//   com+=0.25*getPosition(0);
//   com+=0.25*getPosition(1);
//   com+=0.25*getPosition(2);
//   com+=0.25*getPosition(3);
//   addCentralAtomDerivatives( 0, 0.25*Tensor::identity() );
//   addCentralAtomDerivatives( 1, 0.25*Tensor::identity() );
//   addCentralAtomDerivatives( 2, 0.25*Tensor::identity() );
//   addCentralAtomDerivatives( 3, 0.25*Tensor::identity() );
//   return com;
// }

}
}

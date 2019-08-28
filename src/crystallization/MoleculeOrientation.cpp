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
#include "core/ActionRegister.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization {

//+PLUMEDOC MCOLVAR MOLECULES
/*
Calculate the vectors connecting a pair of atoms in order to represent the orientation of a molecule.

At its simplest this command can be used to calculate the average length of an internal vector in a
collection of different molecules.  When used in conjunction with MutiColvarFunctions in can be used
to do a variety of more complex tasks.

\par Examples

The following input tells plumed to calculate the distances between two of the atoms in a molecule.
This is done for the same set of atoms four different molecules and the average separation is then
calculated.

\plumedfile
MOLECULES MOL1=1,2 MOL2=3,4 MOL3=5,6 MOL4=7,8 MEAN LABEL=mm
PRINT ARG=mm.mean FILE=colvar
\endplumedfile


*/
//+ENDPLUMEDOC


class MoleculeOrientation : public VectorMultiColvar {
private:
  unsigned nvectors;
public:
  static void registerKeywords( Keywords& keys );
  explicit MoleculeOrientation( const ActionOptions& ao );
  AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const override;
  void calculateVector( multicolvar::AtomValuePack& myatoms ) const override;
  void normalizeVector( std::vector<double>& vals ) const override;
  void normalizeVectorDerivatives( MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(MoleculeOrientation,"MOLECULES")

void MoleculeOrientation::registerKeywords( Keywords& keys ) {
  VectorMultiColvar::registerKeywords( keys ); keys.use("MEAN"); keys.use("VMEAN");
  keys.add("numbered","MOL","The numerical indices of the atoms in the molecule. The orientation of the molecule is equal to "
           "the vector connecting the first two atoms specified.  If a third atom is specified its position "
           "is used to specify where the molecule is.  If a third atom is not present the molecule is assumed "
           "to be at the center of the vector connecting the first two atoms.");
  keys.reset_style("MOL","atoms");
}

MoleculeOrientation::MoleculeOrientation( const ActionOptions& ao ):
  Action(ao),
  VectorMultiColvar(ao)
{
  std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword("MOL",-1,all_atoms);
  nvectors = std::floor( ablocks.size() / 2 );
  if( ablocks.size()%2!=0 && 2*nvectors+1!=ablocks.size() ) error("number of atoms in molecule specification is wrong.  Should be two or three.");

  if( all_atoms.size()==0 ) error("No atoms were specified");
  setVectorDimensionality( 3*nvectors ); setupMultiColvarBase( all_atoms );

  if( ablocks.size()==2*nvectors+1  ) {
    std::vector<bool> catom_ind(ablocks.size(), false); catom_ind[ablocks.size()-1]=true;
    setAtomsForCentralAtom( catom_ind );
  }
}

AtomNumber MoleculeOrientation::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_lab.size() );
  plumed_assert( atom_lab[iatom].first==0 );
  return ActionAtomistic::getAbsoluteIndex( ablocks[0][atom_lab[iatom].second] );
}

void MoleculeOrientation::calculateVector( multicolvar::AtomValuePack& myatoms ) const {
  for(unsigned i=0; i<nvectors; ++i) {
    Vector distance; distance=getSeparation( myatoms.getPosition(2*i+0), myatoms.getPosition(2*i+1) );

    addAtomDerivatives( 2+3*i+0, 2*i+0, Vector(-1.0,0,0), myatoms );
    addAtomDerivatives( 2+3*i+0, 2*i+1, Vector(+1.0,0,0), myatoms );
    myatoms.addBoxDerivatives( 2+3*i+0, Tensor(distance,Vector(-1.0,0,0)) );
    myatoms.addValue( 2+3*i+0, distance[0] );

    addAtomDerivatives( 2+3*i+1, 2*i+0, Vector(0,-1.0,0), myatoms );
    addAtomDerivatives( 2+3*i+1, 2*i+1, Vector(0,+1.0,0), myatoms );
    myatoms.addBoxDerivatives( 2+3*i+1, Tensor(distance,Vector(0,-1.0,0)) );
    myatoms.addValue( 2+3*i+1, distance[1] );

    addAtomDerivatives( 2+3*i+2, 2*i+0, Vector(0,0,-1.0), myatoms );
    addAtomDerivatives( 2+3*i+2, 2*i+1, Vector(0,0,+1.0), myatoms );
    myatoms.addBoxDerivatives( 2+3*i+2, Tensor(distance,Vector(0,0,-1.0)) );
    myatoms.addValue( 2+3*i+2, distance[2] );
  }
}

void MoleculeOrientation::normalizeVector( std::vector<double>& vals ) const {
  for(unsigned i=0; i<nvectors; ++i) {
    double norm=0;
    for(unsigned j=0; j<3; ++j) norm += vals[2+3*i+j]*vals[2+3*i+j];
    norm = sqrt(norm);

    double inorm = 1.0; if( norm>epsilon ) inorm = 1.0 / norm;
    for(unsigned j=0; j<3; ++j) vals[2+3*i+j] = inorm*vals[2+3*i+j];
  }
}

void MoleculeOrientation::normalizeVectorDerivatives( MultiValue& myvals ) const {
  std::vector<double> weight( nvectors ), wdf( nvectors );
  for(unsigned ivec=0; ivec<nvectors; ++ivec) {
    double v=0; for(unsigned jcomp=0; jcomp<3; ++jcomp) v += myvals.get( 2+3*ivec+jcomp )*myvals.get( 2+3*ivec+jcomp );
    v=sqrt(v); weight[ivec]=1.0; wdf[ivec]=1.0;
    if( v>epsilon ) { weight[ivec] = 1.0 / v; wdf[ivec] = 1.0 / ( v*v*v ); }
  }

  for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
    unsigned jder=myvals.getActiveIndex(j);
    for(unsigned ivec=0; ivec<nvectors; ++ivec) {
      double comp2=0.0; for(unsigned jcomp=0; jcomp<3; ++jcomp) comp2 += myvals.get(2+3*ivec+jcomp)*myvals.getDerivative( 2+3*ivec+jcomp, jder );
      for(unsigned jcomp=0; jcomp<3; ++jcomp) {
        myvals.setDerivative( 2+3*ivec+jcomp, jder, weight[ivec]*myvals.getDerivative( 2+3*ivec+jcomp, jder ) - wdf[ivec]*comp2*myvals.get(2+3*ivec+jcomp) );
      }
    }
  }
}

}
}

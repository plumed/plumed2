/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR POS2VECTOR
/*

*/
//+ENDPLUMEDOC

class Pos2Vector :
  public ActionAtomistic,
  public ActionWithValue
{
private:
  bool nopbc;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  explicit Pos2Vector(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  void calculate() override;
  void apply() override;
};

PLUMED_REGISTER_ACTION(Pos2Vector,"POS2VECTOR")

void Pos2Vector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys); ActionAtomistic::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms whose positions are being used to construct the vector");
  keys.addFlag("NOPBC",false,"do not make the group of atoms whole when storing in the vector");
}

Pos2Vector::Pos2Vector(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao)
{
  std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms); parseFlag("NOPBC",nopbc); 
  log.printf("  making vector from positions of atoms : "); 
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ", atoms[i].serial() );
  log.printf("\n"); requestAtoms( atoms ); forcesToApply.resize( 3*atoms.size() + 9 );
  // Create vector to hold the atoms
  std::vector<unsigned> shape(1); shape[0] = 3*atoms.size(); addValue( shape );
}


unsigned Pos2Vector::getNumberOfDerivatives() const {
  return 0;
}

void Pos2Vector::calculate() {
  if(!nopbc) makeWhole();

  for(unsigned i=0; i<getNumberOfAtoms();++i) {
      Vector pos = getPosition(i); for(unsigned j=0; j<3; ++j) getPntrToOutput(0)->set( 3*i+j, pos[j] );
  }
} 

void Pos2Vector::apply() {
  if( doNotCalculateDerivatives() || !getPntrToOutput(0)->forcesWereAdded() ) return;

  Tensor virial; virial.zero();
  for(unsigned i=0;i<getNumberOfAtoms();++i) {
      Vector f, pos = getPosition(i);
      for(unsigned j=0;j<3;++j) { f[j] = getPntrToOutput(0)->getForce(3*i+j); forcesToApply[3*i+j] = f[j]; }
      virial -= Tensor( pos, f );
  }
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) forcesToApply[3*getNumberOfAtoms() + 3*i + j] = virial(i,j);
  unsigned mm=0; setForcesOnAtoms( forcesToApply, mm );
}

}
}

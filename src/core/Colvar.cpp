/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "Colvar.h"
#include <vector>
#include <string>

using namespace std;
namespace PLMD{

Colvar::Colvar(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
isEnergy(false)
{
}

void Colvar::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}  

void Colvar::requestAtoms(const vector<AtomNumber> & a){
  plumed_massert(!isEnergy,"request atoms should not be called if this is energy");
// Tell actionAtomistic what atoms we are getting
  ActionAtomistic::requestAtoms(a);
// Resize the derivatives of all atoms
  for(int i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
// Set the size of the forces array
  forces.resize(3*getNumberOfAtoms()+9);
}

void Colvar::apply(){
  vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());
  const unsigned    nat=getNumberOfAtoms();
  const unsigned    ncp=getNumberOfComponents();
  const unsigned    fsz=f.size();

  for(unsigned i=0;i<fsz;i++) f[i].zero();
  v.zero();

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>comm.Get_size()) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  if(!isEnergy){
    for(unsigned i=rank;i<ncp;i+=stride){
      if(getPntrToComponent(i)->applyForce(forces)){
        for(unsigned j=0;j<nat;++j){
          f[j][0]+=forces[3*j+0];
          f[j][1]+=forces[3*j+1];
          f[j][2]+=forces[3*j+2];
        }
        v(0,0)+=forces[3*nat+0];
        v(0,1)+=forces[3*nat+1];
        v(0,2)+=forces[3*nat+2];
        v(1,0)+=forces[3*nat+3];
        v(1,1)+=forces[3*nat+4];
        v(1,2)+=forces[3*nat+5];
        v(2,0)+=forces[3*nat+6];
        v(2,1)+=forces[3*nat+7];
        v(2,2)+=forces[3*nat+8];
      }
    }
    if(ncp>comm.Get_size()) {
      if(fsz>0) comm.Sum(&f[0][0],3*fsz);
      comm.Sum(&v[0][0],9);
    }
  } else if( isEnergy ){
    forces.resize(1);
    if(getPntrToComponent(0)->applyForce(forces)) modifyForceOnEnergy()+=forces[0];
  }
}

void Colvar::setBoxDerivativesNoPbc(Value* v){
  Tensor virial;
  unsigned nat=getNumberOfAtoms();
  for(unsigned i=0;i<nat;i++) virial-=Tensor(getPosition(i),
    Vector(v->getDerivative(3*i+0),
           v->getDerivative(3*i+1),
           v->getDerivative(3*i+2)));
  setBoxDerivatives(v,virial);
}
}

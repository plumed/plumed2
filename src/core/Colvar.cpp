/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "tools/OpenMP.h"
#include <vector>
#include <string>

using namespace std;
namespace PLMD {

Colvar::Colvar(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  isEnergy(false),
  isExtraCV(false)
{
}

void Colvar::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

void Colvar::requestAtoms(const vector<AtomNumber> & a) {
  plumed_massert(!isEnergy,"request atoms should not be called if this is energy");
// Tell actionAtomistic what atoms we are getting
  ActionAtomistic::requestAtoms(a);
// Resize the derivatives of all atoms
  for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
}

void Colvar::apply() {
  vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());
  const unsigned    nat=getNumberOfAtoms();
  const unsigned    ncp=getNumberOfComponents();
  const unsigned    fsz=f.size();

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>4*comm.Get_size()) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();
  if(nt>ncp/(4*stride)) nt=1;

  if(!isEnergy && !isExtraCV) {
    #pragma omp parallel num_threads(nt)
    {
      vector<Vector> omp_f(fsz);
      Tensor         omp_v;
      vector<double> forces(3*nat+9);
      #pragma omp for
      for(unsigned i=rank; i<ncp; i+=stride) {
        if(getPntrToComponent(i)->applyForce(forces)) {
          for(unsigned j=0; j<nat; ++j) {
            omp_f[j][0]+=forces[3*j+0];
            omp_f[j][1]+=forces[3*j+1];
            omp_f[j][2]+=forces[3*j+2];
          }
          omp_v(0,0)+=forces[3*nat+0];
          omp_v(0,1)+=forces[3*nat+1];
          omp_v(0,2)+=forces[3*nat+2];
          omp_v(1,0)+=forces[3*nat+3];
          omp_v(1,1)+=forces[3*nat+4];
          omp_v(1,2)+=forces[3*nat+5];
          omp_v(2,0)+=forces[3*nat+6];
          omp_v(2,1)+=forces[3*nat+7];
          omp_v(2,2)+=forces[3*nat+8];
        }
      }
      #pragma omp critical
      {
        for(unsigned j=0; j<nat; ++j) f[j]+=omp_f[j];
        v+=omp_v;
      }
    }

    if(ncp>4*comm.Get_size()) {
      if(fsz>0) comm.Sum(&f[0][0],3*fsz);
      comm.Sum(&v[0][0],9);
    }

  } else if( isEnergy ) {
    vector<double> forces(1);
    if(getPntrToComponent(0)->applyForce(forces)) modifyForceOnEnergy()+=forces[0];
  } else if( isExtraCV ) {
    vector<double> forces(1);
    if(getPntrToComponent(0)->applyForce(forces)) modifyForceOnExtraCV()+=forces[0];
  }
}

void Colvar::setBoxDerivativesNoPbc(Value* v) {
  Tensor virial;
  unsigned nat=getNumberOfAtoms();
  for(unsigned i=0; i<nat; i++) virial-=Tensor(getPosition(i),
                                          Vector(v->getDerivative(3*i+0),
                                              v->getDerivative(3*i+1),
                                              v->getDerivative(3*i+2)));
  setBoxDerivatives(v,virial);
}
}

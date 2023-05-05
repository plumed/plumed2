/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

void Colvar::requestAtoms(const std::vector<AtomNumber> & a) {
  plumed_massert(!isEnergy,"request atoms should not be called if this is energy");
// Tell actionAtomistic what atoms we are getting
  ActionAtomistic::requestAtoms(a);
// Resize the derivatives of all atoms
  for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
}

void Colvar::apply() {
  const unsigned    nat=getNumberOfAtoms();
  const unsigned    ncp=getNumberOfComponents();

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>4*comm.Get_size()) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  std::vector<double> f(3*nat+9,0);
  unsigned nt=OpenMP::getNumThreads();
  if(nt>ncp/(4*stride)) nt=1;

  if(!isEnergy && !isExtraCV) {
    #pragma omp parallel num_threads(nt)
    {
      std::vector<double> omp_f(3*nat+9,0);
      std::vector<double> forces(3*nat+9);
      #pragma omp for
      for(unsigned i=rank; i<ncp; i+=stride) {
        if(getPntrToComponent(i)->applyForce(forces)) {
          for(unsigned j=0; j<forces.size(); ++j) omp_f[j]+=forces[j];
        }
      }
      #pragma omp critical
      {
        for(unsigned j=0; j<f.size(); ++j) f[j]+=omp_f[j];
      }
    }

    if(ncp>4*comm.Get_size()) comm.Sum(&f[0],3*nat+9);
    unsigned ind=0; setForcesOnAtoms( f, ind );
  } else if( isEnergy ) {
    std::vector<double> forces(1);
    if(getPntrToComponent(0)->applyForce(forces)) modifyForceOnEnergy()+=forces[0];
  } 
  // else if( isExtraCV ) {
  //   std::vector<double> forces(1);
  //   if(getPntrToComponent(0)->applyForce(forces)) modifyForceOnExtraCV()+=forces[0];
  // }
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

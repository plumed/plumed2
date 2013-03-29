/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_manyrestraints_ManyRestraintsBase_h
#define __PLUMED_manyrestraints_ManyRestraintsBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "vesselbase/ActionWithVessel.h"

namespace PLMD {
namespace manyrestraints {

class ManyRestraintsBase :
 public ActionAtomistic,
 public ActionWithValue,
 public ActionPilot,
 public vesselbase::ActionWithVessel
{
private:
  std::vector<double> forcesToApply;
protected:
  void createRestraints( const unsigned& nrestraints );
/// Add some derivatives to a particular atom
  void addAtomsDerivatives(const int&,const Vector&);
/// Add some derivatives to the virial
  void addBoxDerivatives(const Tensor&);
public:
  static void registerKeywords( Keywords& keys );
  ManyRestraintsBase(const ActionOptions&);
  unsigned getNumberOfDerivatives();
/// Deactivate task now does nothing
  void deactivate_task(){}
  void apply();
};

inline
unsigned ManyRestraintsBase::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms() + 9;
}

inline
void ManyRestraintsBase::addAtomsDerivatives(const int& iatom, const Vector& der){
  plumed_dbg_assert( iatom<getNumberOfAtoms() );
  addElementDerivative( 3*iatom+0, der[0] );
  addElementDerivative( 3*iatom+1, der[1] );
  addElementDerivative( 3*iatom+2, der[2] );
}
  
inline
void ManyRestraintsBase::addBoxDerivatives(const Tensor& vir){
  int natoms=getNumberOfAtoms();
  addElementDerivative( 3*natoms+0, vir(0,0) );
  addElementDerivative( 3*natoms+1, vir(0,1) );
  addElementDerivative( 3*natoms+2, vir(0,2) );
  addElementDerivative( 3*natoms+3, vir(1,0) );
  addElementDerivative( 3*natoms+4, vir(1,1) );
  addElementDerivative( 3*natoms+5, vir(1,2) );
  addElementDerivative( 3*natoms+6, vir(2,0) );
  addElementDerivative( 3*natoms+7, vir(2,1) );
  addElementDerivative( 3*natoms+8, vir(2,2) );
}

}
}

#endif

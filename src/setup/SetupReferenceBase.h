/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#ifndef __PLUMED_setup_SetupReferenceBase_h
#define __PLUMED_setup_SetupReferenceBase_h

#include "core/ActionSetup.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace setup {

class SetupReferenceBase : 
public ActionSetup, 
public ActionAtomistic,
public ActionWithValue
{
friend class DRMSD;
protected:
  bool hasatoms;
  unsigned nblocks;
  std::vector<unsigned> blocks;
  std::vector<AtomNumber> myindices, mygroup;
public: 
  static void registerKeywords( Keywords& keys );
  explicit SetupReferenceBase(const ActionOptions&ao);
  ~SetupReferenceBase();
  void activate() {} 
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() const { return 0; }
/// The number of virtual atoms that are calculated by this action
  unsigned getNumberOfVirtualAtoms() const ;
  void getNatomsAndNargs( unsigned& natoms, unsigned& nargs ) const ;
  void transferDataToPlumed( const unsigned& npos, std::vector<double>& masses, std::vector<double>& charges, 
                             std::vector<Vector>& positions, const std::string& argname, PlumedMain& plmd ) const ;
};

inline
unsigned SetupReferenceBase::getNumberOfVirtualAtoms() const {
  return myindices.size();
}

}
}
#endif

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
#ifndef __PLUMED_setup_ReadReferenceConfiguration_h
#define __PLUMED_setup_ReadReferenceConfiguration_h

#include "core/ActionSetup.h"
#include "core/ActionAtomistic.h"

namespace PLMD {
namespace setup {

class ReadReferenceConfiguration : 
public ActionSetup, 
public ActionAtomistic
{
friend class DRMSD;
friend class CalculateReferenceValue;
private:
  unsigned nblocks;
  std::vector<unsigned> blocks;
  std::vector<AtomNumber> myindices;
public: 
  static void registerKeywords( Keywords& keys );
  explicit ReadReferenceConfiguration(const ActionOptions&ao);
  ~ReadReferenceConfiguration();
  void activate() {} 
  /// The number of virtual atoms that are calculated by this action
  unsigned getNumberOfVirtualAtoms() const ;
};

inline
unsigned ReadReferenceConfiguration::getNumberOfVirtualAtoms() const {
  return myindices.size();
}

}
}
#endif

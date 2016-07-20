/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#ifndef __PLUMED_multicolvar_MultiColvarFunction_h
#define __PLUMED_multicolvar_MultiColvarFunction_h

#include "tools/Matrix.h"
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "CatomPack.h"
#include "vesselbase/StoreDataVessel.h"

namespace PLMD {
namespace multicolvar {

class MultiColvarFunction : public MultiColvarBase {
private:
/// A tempory vector that is used for retrieving vectors
  std::vector<double> tvals;
protected:
/// Get the derivatives for the central atom with index ind
  CatomPack getCentralAtomPackFromInput( const unsigned& ind ) const ;
/// Build sets by taking one multicolvar from each base
  void buildSets();
public:
  explicit MultiColvarFunction(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
};

inline
CatomPack MultiColvarFunction::getCentralAtomPackFromInput( const unsigned& ind ) const {
  plumed_dbg_assert( atom_lab[ind].first>0 ); unsigned mmc=atom_lab[ind].first-1;
  unsigned basen=0;
  for(unsigned i=0;i<mmc;++i) basen+=mybasemulticolvars[i]->getNumberOfAtoms();
  return mybasemulticolvars[mmc]->getCentralAtomPack( basen, atom_lab[ind].second );
}

}
}
#endif

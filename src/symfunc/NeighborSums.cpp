/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "NeighborSums.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace symfunc {

/// The actions of the core (the plugin derives its own from NeighborSums)
template <class PF>
class NeighborSumsCPU : public NeighborSums<PF> {
public:
  static void registerKeywords( Keywords& keys ) {
    NeighborSums<PF>::registerKeywords( keys );
  }
  explicit NeighborSumsCPU( const ActionOptions& ao ):
    Action(ao),
    NeighborSums<PF>(ao) {
    this->checkRead();
  }
};

// The names are those of the shortcuts that create these actions, with a suffix
typedef NeighborSumsCPU<ContactPairs> ContactSums;
PLUMED_REGISTER_ACTION(ContactSums,"COORDINATIONNUMBER_SUMS")
PLUMED_REGISTER_ACTION(ContactSums,"Q1_DENOM_SUMS")
PLUMED_REGISTER_ACTION(ContactSums,"Q3_DENOM_SUMS")
PLUMED_REGISTER_ACTION(ContactSums,"Q4_DENOM_SUMS")
PLUMED_REGISTER_ACTION(ContactSums,"Q6_DENOM_SUMS")
typedef NeighborSumsCPU<SphericalHarmonicPairs> SphericalHarmonicSums;
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q1_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q3_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q4_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q6_SUMS")

}
}

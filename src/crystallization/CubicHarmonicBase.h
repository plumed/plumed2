/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#ifndef __PLUMED_crystallization_CubicHarmonicBase_h
#define __PLUMED_crystallization_CubicHarmonicBase_h

#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace crystallization {

class CubicHarmonicBase : public multicolvar::MultiColvarBase {
private:
//  double nl_cut;
  double rcut2;
  double rotationmatrix[3][3];
  bool unormalized;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit CubicHarmonicBase(const ActionOptions&);
// active methods:
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
  virtual double calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const = 0;
/// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
};

}
}
#endif

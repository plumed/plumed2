/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_gridtools_ActionWithIntegral_h
#define __PLUMED_gridtools_ActionWithIntegral_h

#include "ActionWithInputGrid.h"

namespace PLMD {
namespace gridtools {

class ActionWithIntegral : public ActionWithInputGrid {
private:
  double volume;
  std::vector<double> forcesToApply;
protected:
/// Get the volume of a grid box
  double getVolume() const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithIntegral(const ActionOptions&ao);
  unsigned getNumberOfDerivatives();
  void turnOnDerivatives();
/// Unless I am mistaken an integral should never be a periodic function
  bool isPeriodic() { return false; }
  void apply();
};

inline
double ActionWithIntegral::getVolume() const {
  return volume;
}

inline
unsigned ActionWithIntegral::getNumberOfDerivatives() {
  return ingrid->getNumberOfPoints();
}

}
}
#endif


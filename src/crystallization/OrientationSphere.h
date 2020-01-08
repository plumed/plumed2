/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_crystallization_OrientationSphere_h
#define __PLUMED_crystallization_OrientationSphere_h

#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace crystallization {

class OrientationSphere : public multicolvar::MultiColvarBase {
private:
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit OrientationSphere(const ActionOptions&);
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
  virtual double computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                        Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const = 0;
  virtual double calculateCoordinationPrefactor( const double& coord, double& df ) const ;
  bool isPeriodic() { return false; }
};

inline
double OrientationSphere::calculateCoordinationPrefactor( const double& coord, double& df ) const {
  df=0.0; return 1.0;
}

}
}
#endif

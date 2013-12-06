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
#ifndef __PLUMED_crystallization_OrientationSphere_h
#define __PLUMED_crystallization_OrientationSphere_h

#include "multicolvar/MultiColvarFunction.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath> 

namespace PLMD {
namespace crystallization {

class OrientationSphere : public multicolvar::MultiColvarFunction {
private:
  std::vector<double> catom_orient, catom_der, this_orient;
  std::vector<double> catom_iorient, catom_ider, this_iorient;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  OrientationSphere(const ActionOptions&);
  double compute();
  void calculateWeight();
  virtual double transformDotProduct( const double& dot, double& df );
  Vector getCentralAtom();  
  bool isPeriodic(){ return false; }
};

inline
double OrientationSphere::transformDotProduct( const double& dot, double& df ){
  df=1.0; return dot;
}

}
}
#endif

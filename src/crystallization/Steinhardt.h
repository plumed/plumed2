/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_crystallization_Steinhardt_h
#define __PLUMED_crystallization_Steinhardt_h

#include <complex>
#include "tools/SwitchingFunction.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization {

class Steinhardt : public VectorMultiColvar {
private:
  unsigned tmom;
  double rcut,rcut2;
  SwitchingFunction switchingFunction;
protected:
  std::vector<double> coeff_poly;
  std::vector<double> normaliz;
  void setAngularMomentum( const unsigned& ang );
public:
  static void registerKeywords( Keywords& keys );
  explicit Steinhardt( const ActionOptions& ao );
  void calculateVector( multicolvar::AtomValuePack& myatoms ) const ;
  double deriv_poly( const unsigned&, const double&, double& ) const ;
};

}
}
#endif

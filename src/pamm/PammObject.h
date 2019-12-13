/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#ifndef __PLUMED_pamm_PammObject_h
#define __PLUMED_pamm_PammObject_h

#include <vector>
#include "core/Value.h"
#include "tools/KernelFunctions.h"

namespace PLMD {
namespace pamm {

class PammObject {
private:
/// Regularisation parameter to use
  double regulariser;
/// Is the domain periodic
  std::vector<bool> pbc;
/// The domain of the function
  std::vector<std::string> min, max;
/// List of kernel functions involved
  std::vector<std::unique_ptr<KernelFunctions>> kernels;
public:
// Explicit definitions for constructor, copy constructor and destructor
  PammObject();
  PammObject( const PammObject& );
/// GB: I fixed this (should return PammObject&, it was returning PammObject
// However I am not sure the implementation makes sense.
  PammObject& operator=(const PammObject& po) { plumed_error(); regulariser=po.regulariser; return *this; }
/// Setup the Pamm object
  void setup( const std::string& filename, const double& reg, const std::vector<std::string>& valnames,
              const std::vector<bool>& pbcin, const std::vector<std::string>& imin, const std::vector<std::string>& imax,
              std::string& errorstr );
///
  void evaluate( const std::vector<double>& invar, std::vector<double>& outvals, std::vector<std::vector<double> >& der ) const ;
///
  unsigned getNumberOfKernels() const ;
///
  std::vector<double> getKernelCenter( const unsigned& kno ) const ;
///
  std::vector<double> getKernelSupport( const unsigned& kno ) const ;
};

inline
unsigned PammObject::getNumberOfKernels() const {
  return kernels.size();
}

inline
std::vector<double> PammObject::getKernelCenter( const unsigned& kno ) const {
  return kernels[kno]->getCenter();
}

inline
std::vector<double> PammObject::getKernelSupport( const unsigned& kno ) const {
  return kernels[kno]->getContinuousSupport();
}

}
}

#endif


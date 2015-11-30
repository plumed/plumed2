/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_multicolvar_HBPammObject_h
#define __PLUMED_multicolvar_HBPammObject_h

#include <vector>
#include "tools/Vector.h"
#include "core/Value.h"
#include "AtomValuePack.h"
#include "tools/KernelFunctions.h"

namespace PLMD {
namespace multicolvar {

class HBPammObject {
private:
/// Pointer to base class in multicolvar
  MultiColvarBase* mymulti;
/// Regularisation parameter to use
  double regulariser;
/// List of kernel functions involved
  std::vector<KernelFunctions*> kernels;
public:
// Explicit definitions for constructor, copy constructor and destructor
  HBPammObject();
  HBPammObject( const HBPammObject& );
  ~HBPammObject();
/// Setup the HBPamm object
  void setup( const std::string& filename, const double& reg, MultiColvarBase* mybase, std::string& errorstr );
/// Get the cutoff to use throughout
  double get_cutoff() const ;
/// Evaluate the HBPamm Object
  double evaluate( const unsigned& dno, const unsigned& ano, const unsigned& hno, 
                   const Vector& d_da, const double& md_da, AtomValuePack& myatoms ) const ;
};

}
}

#endif

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
#ifndef __PLUMED_reference_FakeFrame_h
#define __PLUMED_reference_FakeFrame_h

#include "ReferenceConfiguration.h"

namespace PLMD {

class FakeFrame :
  public PLMD::ReferenceConfiguration
{
public:
  explicit FakeFrame( const ReferenceConfigurationOptions& ro ) : ReferenceConfiguration(ro) {}
  void read( const PDB& ) override { plumed_merror("should not be called"); }
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const override {
    plumed_merror("should not be called"); return 1.0;
  }
};

}
#endif




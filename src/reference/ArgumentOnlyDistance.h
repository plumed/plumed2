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
#ifndef __PLUMED_reference_ArgumentOnlyDistance_h
#define __PLUMED_reference_ArgumentOnlyDistance_h

#include <vector>
#include <string>
#include "ReferenceArguments.h"

namespace PLMD {

class PDB;
class Pbc;

class ArgumentOnlyDistance : public ReferenceArguments {
public:
  explicit ArgumentOnlyDistance( const ReferenceConfigurationOptions& ro );
  void read( const PDB& pdb ) override;
  bool pcaIsEnabledForThisReference() override { return true; }
  void setupPCAStorage( ReferenceValuePack& mypack ) override { mypack.switchOnPCAOption(); }
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const override;
  double calculate( const std::vector<Value*>& vals, ReferenceValuePack& myder, const bool& squared ) const ;
};

}
#endif

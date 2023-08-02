/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "MetricRegister.h"
#include "ArgumentOnlyDistance.h"
#include "core/Value.h"
#include <limits>

namespace PLMD {

class DotProductDistance : public ArgumentOnlyDistance {
public:
  explicit DotProductDistance( const ReferenceConfigurationOptions& ro );
  void read( const PDB& ) override;
  double calculateArgumentDistance( const std::vector<Value*> & vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const override;
};

PLUMED_REGISTER_METRIC(DotProductDistance,"DOTPRODUCT")

DotProductDistance::DotProductDistance( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration(ro),
  ArgumentOnlyDistance(ro)
{
}

void DotProductDistance::read( const PDB& pdb ) {
  readArgumentsFromPDB( pdb );
}

double DotProductDistance::calculateArgumentDistance( const std::vector<Value*> & vals, const std::vector<double>& arg,
    ReferenceValuePack& myder, const bool& squared ) const {
  double dot=0.0;
  for (std::size_t i=0; i<vals.size(); ++i) dot+=getReferenceArgument(i)*arg[i];
  for (std::size_t i=0; i<vals.size(); ++i) myder.setArgumentDerivatives( i, -getReferenceArgument(i)/dot );
  if(dot==0.0) dot=std::numeric_limits<double>::min();
  return -log(dot);
}


}

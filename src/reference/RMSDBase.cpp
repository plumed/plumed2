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
#include "RMSDBase.h"

namespace PLMD {

RMSDBase::RMSDBase( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration(ro),
  SingleDomainRMSD(ro)
{
}

double RMSDBase::calculate( const std::vector<Vector>& pos, ReferenceValuePack& myder, const bool& squared ) const {
//  clearDerivatives();
  return calc( pos, myder, squared );
}

double RMSDBase::calc( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const {
  plumed_dbg_assert( pos.size()==getNumberOfAtoms() );
  return calc( pos, myder, squared );
}

}

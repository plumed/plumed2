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
#include "ArgumentOnlyDistance.h"
#include "core/Value.h"

namespace PLMD {

ArgumentOnlyDistance::ArgumentOnlyDistance( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration(ro),
  ReferenceArguments(ro)
{
}

void ArgumentOnlyDistance::read( const PDB& pdb ) {
  readArgumentsFromPDB( pdb );
}

double ArgumentOnlyDistance::calculate( const std::vector<Value*>& vals, ReferenceValuePack& myder, const bool& squared ) const {
  std::vector<double> tmparg( vals.size() );
  for(unsigned i=0; i<vals.size(); ++i) tmparg[i]=vals[i]->get();
  double d=calculateArgumentDistance( vals, tmparg, myder, squared );
  if( !myder.updateComplete() ) myder.updateDynamicLists();
  return d;
}

double ArgumentOnlyDistance::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg,
                                   ReferenceValuePack& myder, const bool& squared ) const {
  plumed_dbg_assert( pos.size()==0 );
  double d=calculateArgumentDistance( vals, arg, myder, squared );
  if( !myder.updateComplete() ) myder.updateDynamicLists();
  return d;
}

}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "MultiValue.h"

namespace PLMD{
namespace vesselbase{

MultiValue::MultiValue( const unsigned& nvals, const unsigned& nder ):
values(nvals),
derivatives(nvals,nder)
{
  hasDerivatives.clear();
  for(unsigned i=0;i<nder;++i) hasDerivatives.addIndexToList( i );
  hasDerivatives.deactivateAll();
}

void MultiValue::clearAll(){
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();
  for(unsigned i=0;i<values.size();++i) clear(i);
  hasDerivatives.deactivateAll();
}

void MultiValue::clear( const unsigned& ival ){
  values[ival]=0;
  for(unsigned i=0;i<hasDerivatives.getNumberActive();++i){
      unsigned jder=hasDerivatives[i]; derivatives(ival,jder)=0.;   
  }
}

void MultiValue::chainRule( const unsigned& ival, const unsigned& iout, const unsigned& stride, const unsigned& off, 
                            const double& df, const unsigned& bufstart, std::vector<double>& buffer ){
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  plumed_dbg_assert( off<stride );
  unsigned start=bufstart+stride*(derivatives.ncols()+1)*iout + stride + off; 
  for(unsigned i=0;i<hasDerivatives.getNumberActive();++i){
      unsigned jder=hasDerivatives[i];
      buffer[start+jder*stride] += df*derivatives(ival,jder);
  }
}

}
}

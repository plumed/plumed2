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

MultiValue::MultiValue( const unsigned& nvals, const unsigned& nder ):
values(nvals),
derivatives(nvals,nder),
atLeastOneSet(false)
{
  std::vector<unsigned> myind( nder );
  for(unsigned i=0;i<nder;++i) myind[i]=i;
  hasDerivatives.createIndexListFromVector( myind ); 
}

void MultiValue::resize( const unsigned& nvals, const unsigned& nder ){
  values.resize(nvals); derivatives.resize( nvals, nder );
  hasDerivatives.clear(); std::vector<unsigned> myind( nder ); 
  for(unsigned i=0;i<nder;++i) myind[i]=i;
  hasDerivatives.createIndexListFromVector( myind );
  atLeastOneSet=false;
}

void MultiValue::clearAll(){
  if( atLeastOneSet && !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();
  for(unsigned i=0;i<values.size();++i) clear(i);
  hasDerivatives.deactivateAll(); atLeastOneSet=false;
}

void MultiValue::clear( const unsigned& ival ){
  values[ival]=0;
  if( atLeastOneSet ){
      for(unsigned i=0;i<hasDerivatives.getNumberActive();++i){
          unsigned jder=hasDerivatives[i]; derivatives(ival,jder)=0.;   
      }
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

void MultiValue::copyValues( MultiValue& outvals ) const { 
  plumed_dbg_assert( values.size()<=outvals.getNumberOfValues() );
  for(unsigned i=0;i<values.size();++i) outvals.setValue( i, values[i] ); 
  
}

void MultiValue::copyDerivatives( MultiValue& outvals ){
  plumed_dbg_assert( values.size()<=outvals.getNumberOfValues() && derivatives.ncols()<=outvals.getNumberOfDerivatives() );
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  for(unsigned i=0;i<values.size();++i){
     for(unsigned j=0;j<hasDerivatives.getNumberActive();++j){
        unsigned jder=hasDerivatives[j];
        outvals.addDerivative( i, jder, derivatives(i,jder) );
     }
  }
}

void MultiValue::quotientRule( const unsigned& nder, const unsigned& dder, const unsigned& oder ){
  plumed_dbg_assert( nder<values.size() && dder<values.size() && oder<values.size() );
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  double weight = values[dder], pref = values[nder] / (weight*weight);
  for(unsigned j=0;j<hasDerivatives.getNumberActive();++j){
      unsigned jder=hasDerivatives[j];
      derivatives(oder,jder) = derivatives(nder,jder) / weight - pref*derivatives(dder,jder);
  }
  values[oder] = values[nder] / values[dder];
}

}

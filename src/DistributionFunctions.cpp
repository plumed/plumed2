/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ActionWithDistribution.h"
#include "DistributionFunctions.h"
#include "PlumedException.h"
#include "Value.h"

namespace PLMD {

DistributionFunction::DistributionFunction( const std::string& parameters ):
fine(true)
{
}

DistributionFunction::~DistributionFunction(){
  for(unsigned i=0;i<accumulators.size();++i){ delete thesevalues[i]; delete accumulators[i]; }
}

void DistributionFunction::addAccumulator( const bool wderiv ){
  accumulators.push_back(new Value());
  thesevalues.push_back(new Value());
  hasDeriv.push_back(wderiv);
  plumed_assert( accumulators.size()==thesevalues.size() );
  plumed_assert( accumulators.size()==hasDeriv.size() );
}

void DistributionFunction::setNumberOfDerivatives( const unsigned nder ){
  for(unsigned i=0;i<accumulators.size();++i){
      if(hasDeriv[i]) accumulators[i]->resizeDerivatives(nder); 
      else accumulators[i]->resizeDerivatives(0);
  }
}

void DistributionFunction::reset(){
  for(unsigned i=0;i<accumulators.size();++i){
      accumulators[i]->set(0); accumulators[i]->clearDerivatives();
  }
}

void DistributionFunction::clear(){
  for(unsigned i=0;i<thesevalues.size();++i){
      thesevalues[i]->set(0); thesevalues[i]->clearDerivatives();
  }
}

void DistributionFunction::copyValue( const unsigned nn, Value* value_in  ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to copy derivatives to an accumulator with no derivatives");

  copy( value_in, thesevalues[nn] );
//  unsigned nder=value_in->getNumberOfDerivatives(); 
//  if( nder!=thesevalues[nn]->getNumberOfDerivatives() ){ thesevalues[nn]->resizeDerivatives(nder); }
//  thesevalues[nn]->clearDerivatives();
//  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i) thesevalues[nn]->addDerivative( i, value_in->getDerivative(i) );
//  thesevalues[nn]->set( value_in->get() );
}

void DistributionFunction::extractDerivatives( const unsigned nn, Value* value_out  ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to copy derivatives from an accumulator with no derivatives");
  plumed_assert( accumulators[nn]->getNumberOfDerivatives()==value_out->getNumberOfDerivatives() );
  copy( accumulators[nn], value_out );
//  value_out->clearDerivatives();
//  for(unsigned i=0;i<value_out->getNumberOfDerivatives();++i) value_out->addDerivative( i, accumulators[nn]->getDerivative(i) );
}

void DistributionFunction::setValue( const unsigned nn, const double f ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  thesevalues[nn]->set(f);
} 

void DistributionFunction::chainRule( const unsigned nn, const double f ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to do chain rule on an accumulator with no derivatives");
  thesevalues[nn]->chainRule(f);
} 

void DistributionFunction::multiplyValue( const unsigned nn, Value* val ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to do product rule on an accumulator with no derivatives");
  Value* tmpvalue = new Value(); 
  tmpvalue->resizeDerivatives( val->getNumberOfDerivatives() );
  product( thesevalues[nn], val, tmpvalue ); 
  copyValue( nn, tmpvalue );  
  delete tmpvalue;
}

void DistributionFunction::mergeDerivatives( const unsigned kk, ActionWithDistribution& action ){
  for(unsigned i=0;i<accumulators.size();++i){
     if( thesevalues[i]->valueHasBeenSet() ){
         accumulators[i]->add( thesevalues[i]->get() );
         if(hasDeriv[i]){ action.mergeDerivatives( kk, thesevalues[i], accumulators[i] ); }
     }
  }
}

unsigned DistributionFunction::requiredBufferSpace() const {
  unsigned nbuf=0;
  for(unsigned i=0;i<accumulators.size();++i){
      nbuf+=1;
      if( hasDeriv[i] ) nbuf+=accumulators[i]->getNumberOfDerivatives();
  }
  return nbuf;
}

void DistributionFunction::copyDataToBuffers( unsigned& bufsize, std::vector<double>& buffers ) const {
  plumed_assert( ( bufsize+requiredBufferSpace() )<=buffers.size() );
  for(unsigned i=0;i<accumulators.size();++i){
      buffers[bufsize]=accumulators[i]->get(); bufsize++;
      if( hasDeriv[i] ){
          for(unsigned j=0;j<accumulators[i]->getNumberOfDerivatives();++j){
              buffers[bufsize]=accumulators[i]->getDerivative(j); bufsize++;
          }
      }
  }
}

void DistributionFunction::retrieveDataFromBuffers( unsigned& bufsize, const std::vector<double>& buffers ){
  plumed_assert( ( bufsize+requiredBufferSpace() )<=buffers.size() );
  for(unsigned i=0;i<accumulators.size();++i){
      accumulators[i]->set( buffers[bufsize] ); bufsize++;
      if( hasDeriv[i] ){
          accumulators[i]->clearDerivatives();
          for(unsigned j=0;j<accumulators[i]->getNumberOfDerivatives();++j){
              accumulators[i]->addDerivative( j, buffers[bufsize] ); bufsize++;
          }
      }
  }
}



}

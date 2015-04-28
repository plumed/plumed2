/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "vesselbase/ActionWithVessel.h"
#include "multicolvar/ActionVolume.h"
#include "VectorMultiColvar.h"
#include "Gradient.h"

namespace PLMD {
namespace crystallisation {

class GradientVessel : public vesselbase::FunctionVessel {
private:
  bool isdens;
  unsigned nweights, ncomponents;
  std::vector<unsigned> starts;
  std::vector<double> val_interm;
  Matrix<double> der_interm;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  GradientVessel( const vesselbase::VesselOptions& da );
  std::string function_description();
  void resize();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(GradientVessel,"GRADIENT")

void GradientVessel::registerKeywords( Keywords& keys ){
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void GradientVessel::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("GRADIENT",false,"calculate the gradient",true);
  keys.addOutputComponent("gradient","GRADIENT","the gradient");
}

GradientVessel::GradientVessel( const vesselbase::VesselOptions& da ) :
FunctionVessel(da)
{
   Gradient* vg=dynamic_cast<Gradient*>( getAction() );
   plumed_assert( vg ); isdens=(vg->getPntrToMultiColvar())->isDensity();
   nweights = vg->nbins[0] + vg->nbins[1] + vg->nbins[2];
   ncomponents = vg->vend;

   starts.push_back(0);
   if( vg->nbins[0]>0 ){
       starts.push_back( vg->nbins[0] );
       if( vg->nbins[1]>0 ){
           starts.push_back( vg->nbins[0] + vg->nbins[1] );    
           if( vg->nbins[2]>0 ) starts.push_back( nweights );  
       } else if( vg->nbins[2]>0 ){
           starts.push_back( nweights );
       }
   } else if( vg->nbins[1]>0 ){
       starts.push_back( vg->nbins[1] );
       if( vg->nbins[2]>0 ) starts.push_back( nweights );
   } else if( vg->nbins[2]>0 ){
       starts.push_back( nweights );
   }
}

std::string GradientVessel::function_description(){
  return "the gradient";
}

void GradientVessel::resize(){
  if( getAction()->derivativesAreRequired() ){
     unsigned nder=getAction()->getNumberOfDerivatives();
     resizeBuffer( (1+nder)*(ncomponents+1)*nweights );
     setNumberOfDerivatives( nder );
     val_interm.resize( ncomponents*nweights );
     der_interm.resize( ncomponents*nweights, nder );
  } else {
     setNumberOfDerivatives(0); 
     resizeBuffer( (ncomponents+1)*nweights );
     val_interm.resize( ncomponents*nweights );
  }
}

bool GradientVessel::calculate(){
  for(unsigned iw=0;iw<nweights;++iw){
      unsigned xx = (ncomponents+1)*iw;
      double weight=getAction()->getElementValue(ncomponents + iw);
      addValueIgnoringTolerance( xx, weight ); 
      getAction()->chainRuleForElementDerivatives( xx , ncomponents + iw, 1.0, this );
      for(unsigned jc=0;jc<ncomponents;++jc){
          double colvar=getAction()->getElementValue( jc );
          addValueIgnoringTolerance( xx + 1 + jc, weight*colvar );
          getAction()->chainRuleForElementDerivatives( xx + 1 + jc, jc, weight, this );
          getAction()->chainRuleForElementDerivatives( xx + 1 + jc, ncomponents + iw, colvar, this );   
      }
  }

  return true;
}

void GradientVessel::finish(){
  der_interm=0;  // Clear all interim derivatives
  unsigned nder = getAction()->getNumberOfDerivatives();

  if( isdens ){
      for(unsigned iw=0;iw<nweights;++iw){
          val_interm[iw] = getFinalValue( 2*iw );
          if( getAction()->derivativesAreRequired() ){
              unsigned wstart = 2*iw*(nder+1) + 1;
              for(unsigned jder=0;jder<nder;++jder) der_interm( iw, jder ) += getBufferElement( wstart + jder );
          }
      }
  } else {
      for(unsigned iw=0;iw<nweights;++iw){
          unsigned xx = (ncomponents+1)*iw;
          double sum=0, ww=getFinalValue( xx );
          for(unsigned jc=0;jc<ncomponents;++jc) val_interm[ iw*ncomponents + jc ] = getFinalValue( xx + 1 + jc ) / ww;
          if( getAction()->derivativesAreRequired() ){
              unsigned wstart = xx*(nder+1) + 1;
              for(unsigned jc=0;jc<ncomponents;++jc){
                  unsigned bstart = ( xx + 1 + jc )*(nder+1) + 1;
                  double val = getFinalValue( xx + 1 + jc );
                  for(unsigned jder=0;jder<nder;++jder) 
                     der_interm( iw*ncomponents + jc, jder ) = (1.0/ww)*getBufferElement( bstart + jder ) - (val/(ww*ww))*getBufferElement( wstart + jder );
              }
          }
      }
  }

  double tmp, diff2=0.0; 

  if( getAction()->derivativesAreRequired() ){
     for(unsigned j=0;j<starts.size()-1;++j){
        for(unsigned bin=starts[j];bin<starts[j+1];++bin){
           for(unsigned jc=0;jc<ncomponents;++jc){
               if( bin==starts[j] ){
                  tmp=val_interm[(starts[j+1]-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
                  for(unsigned jder=0;jder<nder;++jder){
                      addDerivativeToFinalValue( jder, +2.0*tmp*der_interm( (starts[j+1]-1)*ncomponents + jc, jder) );
                      addDerivativeToFinalValue( jder, -2.0*tmp*der_interm( bin*ncomponents + jc, jder ) );
                  }
               } else {
                  tmp=val_interm[(bin-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
                  for(unsigned jder=0;jder<nder;++jder){
                      addDerivativeToFinalValue( jder, +2.0*tmp*der_interm( (bin-1)*ncomponents + jc, jder) );
                      addDerivativeToFinalValue( jder, -2.0*tmp*der_interm( bin*ncomponents + jc, jder ) );
                  }
               }
               diff2+=tmp*tmp;
           }
        }
     }
  } else {
     for(unsigned j=0;j<starts.size()-1;++j){
        for(unsigned bin=starts[j];bin<starts[j+1];++bin){
           for(unsigned jc=0;jc<ncomponents;++jc){
               if( bin==starts[j] ) tmp=val_interm[(starts[j+1]-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
               else tmp=val_interm[(bin-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
               diff2+=tmp*tmp;
           }
        }
     }
  }
  setOutputValue( diff2 );  
}

}
}

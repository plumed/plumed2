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

namespace PLMD {
namespace crystallization {

class VectorMean : public vesselbase::FunctionVessel {
private:
  unsigned ncomp, vstart, wnum;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  VectorMean( const vesselbase::VesselOptions& da );
  std::string function_description();
  void resize();
  bool calculate( std::vector<double>& buffer );
  void finish( const std::vector<double>& buffer );
};

PLUMED_REGISTER_VESSEL(VectorMean,"VMEAN")

void VectorMean::registerKeywords( Keywords& keys ){
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void VectorMean::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("VMEAN",false,"calculate the norm of the mean vector.",true);
  keys.addOutputComponent("vmean","VMEAN","the norm of the mean vector. The output component can be refererred to elsewhere in the input "
                                          "file by using the label.vmean");
}

VectorMean::VectorMean( const vesselbase::VesselOptions& da ) :
FunctionVessel(da)
{
   usetol=true;
   multicolvar::ActionVolume* vg=dynamic_cast<multicolvar::ActionVolume*>( getAction() );
   if( vg ){
       ncomp = getAction()->getNumberOfQuantities() - 2;
       vstart=1; wnum=ncomp+1;
   } else { 
       vstart=5; wnum=1; 
       ncomp = getAction()->getNumberOfQuantities() - 5;
   }
}

std::string VectorMean::function_description(){
  return "the norm of the mean vector";
}

void VectorMean::resize(){
  if( ncomp==0 ) ncomp=getAction()->getNumberOfQuantities() - 5;

  if( getAction()->derivativesAreRequired() ){
     unsigned nder=getAction()->getNumberOfDerivatives();
     resizeBuffer( (1+nder)*(ncomp+1) );
     setNumberOfDerivatives( nder );
  } else {
     setNumberOfDerivatives(0); 
     resizeBuffer(ncomp+1);
  }
}

bool VectorMean::calculate( std::vector<double>& buffer ){
  double weight=getAction()->getElementValue(wnum);
  unsigned nder=getAction()->getNumberOfDerivatives();
  plumed_dbg_assert( weight>=getTolerance() ); 
  // bool addval = addValueUsingTolerance( 0, weight );
  buffer[bufstart] += weight;
  if( diffweight ) getAction()->chainRuleForElementDerivatives( 0, wnum, 1.0, bufstart, buffer );
  for(unsigned i=0;i<ncomp;++i){
      double colvar=getAction()->getElementValue( vstart  + i );
      buffer[bufstart + (1+i)*(1+nder)] += weight*colvar;  
      // addValueIgnoringTolerance( 1 + i, weight*colvar );
      getAction()->chainRuleForElementDerivatives( 1+i, vstart+i, weight, bufstart, buffer );
      if( diffweight ) getAction()->chainRuleForElementDerivatives( 1+i, wnum, colvar, bufstart, buffer );
  }
  return true;
}

void VectorMean::finish( const std::vector<double>& buffer ){
  double sum=0, ww=buffer[bufstart]; // getFinalValue(0);
  unsigned nder=getAction()->getNumberOfDerivatives();
  for(unsigned i=0;i<ncomp;++i){ 
     double tmp = buffer[bufstart+(nder+1)*(i+1)] / ww;
     sum+=tmp*tmp; 
  }
  double tw = 1.0 / sqrt(sum);
  setOutputValue( sqrt(sum) ); 
  if( !getAction()->derivativesAreRequired() ) return;

  for(unsigned icomp=0;icomp<ncomp;++icomp){
      double tmp = buffer[(icomp+1)*(1+nder)] / ww; // getFinalValue(icomp+1) / ww;
      unsigned bstart = bufstart + (1+icomp)*(nder+1) + 1;
      for(unsigned jder=0;jder<nder;++jder) addDerivativeToFinalValue( jder, (tw*tmp/ww)*( buffer[bstart + jder] - tmp*buffer[bufstart + 1 + jder] ) );
  }
}

}
}

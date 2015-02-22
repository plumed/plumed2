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

class VectorSum : public vesselbase::FunctionVessel {
private:
  unsigned ncomp, vstart, wnum;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  VectorSum( const vesselbase::VesselOptions& da );
  std::string function_description();
  void resize();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(VectorSum,"VSUM")

void VectorSum::registerKeywords( Keywords& keys ){
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void VectorSum::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("VSUM",false,"calculate the norm of the sum of vectors.",true);
  keys.addOutputComponent("vsum","VSUM","the norm of sum of vectors. The output component can be refererred to elsewhere in the input "
                                        "file by using the label.vsum");
}

VectorSum::VectorSum( const vesselbase::VesselOptions& da ) :
FunctionVessel(da)
{
   multicolvar::ActionVolume* vg=dynamic_cast<multicolvar::ActionVolume*>( getAction() );
   if( vg ){
       ncomp = getAction()->getNumberOfQuantities() - 2;
       vstart=1; wnum=ncomp+1;
   } else { 
       vstart=5; wnum=1; 
       ncomp = getAction()->getNumberOfQuantities() - 5;
   }
}

std::string VectorSum::function_description(){
  return "the norm of the mean vector";
}

void VectorSum::resize(){
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

bool VectorSum::calculate(){
  double weight=getAction()->getElementValue(wnum);
  plumed_dbg_assert( weight>=getTolerance() );
  bool addval = addValueUsingTolerance( 0, weight );
  if( diffweight ) getAction()->chainRuleForElementDerivatives( 0, wnum, 1.0, this );
  for(unsigned i=0;i<ncomp;++i){
      double colvar=getAction()->getElementValue( vstart  + i );
      addValueIgnoringTolerance( 1 + i, weight*colvar );
      getAction()->chainRuleForElementDerivatives( 1+i, vstart+i, weight, this );
      if( diffweight ) getAction()->chainRuleForElementDerivatives( 1+i, wnum, colvar, this );
  }
  return addval;
}

void VectorSum::finish(){
  double sum=0;
  for(unsigned i=0;i<ncomp;++i){ 
     double tmp = getFinalValue(i+1);
     sum+=tmp*tmp; 
  }
  double tw = 1.0 / sqrt(sum);
  setOutputValue( sqrt(sum) ); 
  if( !getAction()->derivativesAreRequired() ) return;

  unsigned nder = getAction()->getNumberOfDerivatives(); 
  for(unsigned icomp=0;icomp<ncomp;++icomp){
      double tmp = getFinalValue(icomp+1);
      unsigned bstart = (1+icomp)*(nder+1) + 1;
      for(unsigned jder=0;jder<nder;++jder) addDerivativeToFinalValue( jder, tw*tmp*getBufferElement( bstart + jder ) );
  }
}

}
}

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
#include "NormedSumVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

NormedSumVessel::NormedSumVessel( const VesselOptions& da ):
VesselAccumulator(da),
donorm(false)
{
  if( getAction()->isPeriodic() ){
      std::string min, max;
      getAction()->retrieveDomain( min, max );
      myvalue2.setDomain( min, max );
  } else {
      myvalue2.setNotPeriodic();
  }
}

void NormedSumVessel::useNorm(){
  donorm=true; addBufferedValue();
}

bool NormedSumVessel::calculate( const unsigned& icv, const double& tolerance ){
  bool keep=false; double val, fval, df;
  if(donorm){
     bool hasNonZeroDerivatives;
     double ww=getWeight( icv, hasNonZeroDerivatives );
     if( ww>tolerance ){
         keep=true; 
         addToBufferElement( 0, ww ); 
         if( hasNonZeroDerivatives ) addDerivativesOfWeight( icv ); 
     }
     if(!keep) return false;

     unsigned jout; 
     val=getAction()->getElementValue();
     for(unsigned j=1;j<getNumberOfValues()+1;++j){
        fval=compute( icv, j-1, val, df  );
        if( fabs( fval )>tolerance ){
            keep=true; 
            jout=value_starts[j]; 
            addToBufferElement( jout, fval); jout++;
            if( !hasNonZeroDerivatives ) getAction()->chainRuleForElementDerivatives( icv, jout, df, this );
        }  
     }
  } else {
     unsigned jout;
     val=getAction()->getElementValue();
     for(unsigned j=0;j<getNumberOfValues();++j){
        fval=compute( icv, j, val, df );
        if( fabs( fval )>tolerance ){
            keep=true; 
            jout=value_starts[j]; 
            addToBufferElement( jout, fval ); jout++;
            getAction()->chainRuleForElementDerivatives( icv, jout, df, this );
        }
     }   
  }
  return keep;
}

void NormedSumVessel::finish( const double& tolerance ){
  if( donorm ){
     getValue(0, myweight2 ); 
     for(unsigned i=0;i<getNumberOfValues();++i){
         getValue( i+1, myvalue2 );       /// periodicity is missing from here
         quotient( myvalue2, myweight2, getPntrToOutput(i) );
     }
  } else {
     for(unsigned i=0;i<getNumberOfValues();++i){
         getValue( i, myvalue2 ); copy( myvalue2, getPntrToOutput(i) );
     }
  }
}

}
}

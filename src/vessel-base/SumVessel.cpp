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
#include "SumVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

SumVessel::SumVessel( const VesselOptions& da ):
VesselAccumulator(da)
{
}

bool SumVessel::calculate( const unsigned& icv, const double& tolerance ){
  bool keep=false; double f, df; unsigned jout;
  Value myval=getAction()->retreiveLastCalculatedValue();
  for(unsigned j=0;j<getNumberOfValues();++j){
      f=compute( j, myval.get(), df );
      if( fabs(f)>tolerance ){
          keep=true; 
          jout=value_starts[j]; 
          addToBufferElement( jout, f ); jout++;
          getAction()->mergeDerivatives( icv, myval, df, jout, this );        
      }  
  }
  return keep;
}

double SumVessel::final_computations( const unsigned& ival, const double& valin, double& df ){
  df=1; return valin; 
}

void SumVessel::finish( const double& tolerance ){
  double f, df;
  for(unsigned i=0;i<getNumberOfValues();++i){
      getValue( i, myvalue2 ); 
      f=final_computations( i, myvalue2.get(), df );
      myvalue2.chainRule(df); myvalue2.set(f);
      copy( myvalue2, getPntrToOutput(i) );
  }
}
}

}

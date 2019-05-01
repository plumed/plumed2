/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "vesselbase/ActionWithVessel.h"
#include "multicolvar/ActionVolume.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization {

class VectorMean : public vesselbase::FunctionVessel {
private:
  unsigned nder;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit VectorMean( const vesselbase::VesselOptions& da );
  std::string value_descriptor();
  void resize();
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  void finish( const std::vector<double>& buffer );
};

PLUMED_REGISTER_VESSEL(VectorMean,"VMEAN")

void VectorMean::registerKeywords( Keywords& keys ) {
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void VectorMean::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","VMEAN","calculate the norm of the mean vector.");
  keys.addOutputComponent("vmean","VMEAN","the norm of the mean vector. The output component can be referred to elsewhere in the input "
                          "file by using the label.vmean");
}

VectorMean::VectorMean( const vesselbase::VesselOptions& da ) :
  FunctionVessel(da),
  nder(0)
{
}

std::string VectorMean::value_descriptor() {
  return "the norm of the mean vector";
}

void VectorMean::resize() {
  unsigned ncomp=getAction()->getNumberOfQuantities() - 2;

  if( getAction()->derivativesAreRequired() ) {
    nder=getAction()->getNumberOfDerivatives();
    resizeBuffer( (1+nder)*(ncomp+1) ); getFinalValue()->resizeDerivatives( nder );
  } else {
    nder=0; resizeBuffer(ncomp+1);
  }
}

void VectorMean::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  unsigned ncomp=getAction()->getNumberOfQuantities()-2;

  double weight=myvals.get(0); plumed_dbg_assert( weight>=getTolerance() );
  buffer[bufstart] += weight;
  for(unsigned i=0; i<ncomp; ++i) buffer[bufstart + (1+i)*(1+nder)] += weight*myvals.get(2+i);
  if( !getAction()->derivativesAreRequired() ) return;

  if( diffweight ) myvals.chainRule( 0, 0, 1, 0, 1.0, bufstart, buffer );
  for(unsigned i=0; i<ncomp; ++i) {
    double colvar=myvals.get(2+i);
    myvals.chainRule( 2+i, 1+i, 1, 0, weight, bufstart, buffer );
    if( diffweight ) myvals.chainRule( 0, 1+i, 1, 0, colvar, bufstart, buffer );
  }
  return;
}

void VectorMean::finish( const std::vector<double>& buffer ) {
  unsigned ncomp=getAction()->getNumberOfQuantities()-2;
  double sum=0, ww=buffer[bufstart];
  for(unsigned i=0; i<ncomp; ++i) {
    double tmp = buffer[bufstart+(nder+1)*(i+1)] / ww;
    sum+=tmp*tmp;
  }
  setOutputValue( sqrt(sum) );
  if( !getAction()->derivativesAreRequired() ) return;

  Value* fval=getFinalValue(); double tw = 1.0 / sqrt(sum);
  for(unsigned icomp=0; icomp<ncomp; ++icomp) {
    double tmp = buffer[bufstart + (icomp+1)*(1+nder)] / ww;
    unsigned bstart = bufstart + (1+icomp)*(nder+1) + 1;
    for(unsigned jder=0; jder<nder; ++jder) fval->addDerivative( jder, (tw*tmp/ww)*( buffer[bstart + jder] - tmp*buffer[bufstart + 1 + jder] ) );
  }
}

}
}

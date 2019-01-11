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

class VectorSum : public vesselbase::FunctionVessel {
private:
  unsigned nder;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit VectorSum( const vesselbase::VesselOptions& da );
  std::string value_descriptor();
  void resize();
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  void finish( const std::vector<double>& buffer );
};

PLUMED_REGISTER_VESSEL(VectorSum,"VSUM")

void VectorSum::registerKeywords( Keywords& keys ) {
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void VectorSum::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","VSUM","calculate the norm of the sum of vectors.");
  keys.addOutputComponent("vsum","VSUM","the norm of sum of vectors. The output component can be referred to elsewhere in the input "
                          "file by using the label.vsum");
}

VectorSum::VectorSum( const vesselbase::VesselOptions& da ) :
  FunctionVessel(da),
  nder(0)
{
}

std::string VectorSum::value_descriptor() {
  return "the norm of the mean vector";
}

void VectorSum::resize() {
  unsigned ncomp=getAction()->getNumberOfQuantities() - 2;

  if( getAction()->derivativesAreRequired() ) {
    nder=getAction()->getNumberOfDerivatives();
    resizeBuffer( (1+nder)*ncomp ); getFinalValue()->resizeDerivatives( nder );
  } else {
    nder=0; resizeBuffer(ncomp);
  }
}

void VectorSum::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  unsigned ncomp=getAction()->getNumberOfQuantities()-2;

  double weight=myvals.get(0);
  plumed_dbg_assert( weight>=getTolerance() );
  for(unsigned i=0; i<ncomp; ++i) buffer[bufstart + i*(1+nder)] += weight*myvals.get(2+i);
  if( !getAction()->derivativesAreRequired() ) return;

  for(unsigned i=0; i<ncomp; ++i) {
    double colvar=myvals.get(2+i);
    myvals.chainRule( 2+i, i, 1, 0, weight, bufstart, buffer );
    if( diffweight ) myvals.chainRule( 0, i, 1, 0, colvar, bufstart, buffer );
  }
  return;
}

void VectorSum::finish( const std::vector<double>& buffer ) {
  unsigned ncomp=getAction()->getNumberOfQuantities()-2;

  double sum=0;
  for(unsigned i=0; i<ncomp; ++i) {
    double tmp = buffer[bufstart+(nder+1)*i];
    sum+=tmp*tmp;
  }
  double tw = 1.0 / sqrt(sum);
  setOutputValue( sqrt(sum) );
  if( !getAction()->derivativesAreRequired() ) return;

  Value* fval=getFinalValue();
  for(unsigned icomp=0; icomp<ncomp; ++icomp) {
    double tmp = buffer[bufstart + icomp*(1+nder)];
    unsigned bstart = bufstart + icomp*(nder+1) + 1;
    for(unsigned jder=0; jder<nder; ++jder) fval->addDerivative( jder, tw*tmp*buffer[bstart + jder] );
  }
}

}
}

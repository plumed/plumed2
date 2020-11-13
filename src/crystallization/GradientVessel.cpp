/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "Gradient.h"

namespace PLMD {
namespace crystallization {

class GradientVessel : public vesselbase::FunctionVessel {
private:
  bool isdens;
  size_t nweights, ncomponents;
  std::vector<unsigned> starts;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit GradientVessel( const vesselbase::VesselOptions& da );
  std::string value_descriptor();
  void resize();
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  void finish( const std::vector<double>& buffer );
};

PLUMED_REGISTER_VESSEL(GradientVessel,"GRADIENT")

void GradientVessel::registerKeywords( Keywords& keys ) {
  vesselbase::FunctionVessel::registerKeywords(keys);
}

void GradientVessel::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","GRADIENT","calculate the gradient");
  keys.addOutputComponent("gradient","GRADIENT","the gradient");
}

GradientVessel::GradientVessel( const vesselbase::VesselOptions& da ) :
  FunctionVessel(da)
{
  Gradient* vg=dynamic_cast<Gradient*>( getAction() );
  plumed_assert( vg ); isdens=(vg->getPntrToMultiColvar())->isDensity();
  nweights = vg->nbins[0] + vg->nbins[1] + vg->nbins[2];
  if( (vg->getPntrToMultiColvar())->getNumberOfQuantities()>2 ) {
    ncomponents = (vg->getPntrToMultiColvar())->getNumberOfQuantities() - 2;
  } else {
    ncomponents = 1;
  }

  starts.push_back(0);
  if( vg->nbins[0]>0 ) {
    starts.push_back( vg->nbins[0] );
    if( vg->nbins[1]>0 ) {
      starts.push_back( vg->nbins[0] + vg->nbins[1] );
      if( vg->nbins[2]>0 ) starts.push_back( nweights );
    } else if( vg->nbins[2]>0 ) {
      starts.push_back( nweights );
    }
  } else if( vg->nbins[1]>0 ) {
    starts.push_back( vg->nbins[1] );
    if( vg->nbins[2]>0 ) starts.push_back( nweights );
  } else if( vg->nbins[2]>0 ) {
    starts.push_back( nweights );
  }
}

std::string GradientVessel::value_descriptor() {
  return "the gradient";
}

void GradientVessel::resize() {
  if( getAction()->derivativesAreRequired() ) {
    unsigned nder=getAction()->getNumberOfDerivatives();
    resizeBuffer( (1+nder)*(ncomponents+1)*nweights );
    getFinalValue()->resizeDerivatives( nder );
  } else {
    resizeBuffer( (ncomponents+1)*nweights );
  }
}

void GradientVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  unsigned nder;
  if( getAction()->derivativesAreRequired() ) nder=getAction()->getNumberOfDerivatives();
  else nder=0;
  unsigned wstart, cstart; if( ncomponents==1 ) { cstart=1; wstart=2; } else { cstart=2; wstart=2+ncomponents; }

  for(unsigned iw=0; iw<nweights; ++iw) {
    unsigned xx = (ncomponents+1)*iw;
    double weight=myvals.get(wstart+iw);
    buffer[bufstart+xx*(nder+1)] += weight;
    myvals.chainRule( wstart + iw, xx, 1, 0, 1.0, bufstart, buffer );
    for(unsigned jc=0; jc<ncomponents; ++jc) {
      double colvar=myvals.get( cstart + jc );
      buffer[bufstart+(xx+1+jc)*(nder+1) ] += weight*colvar;
      myvals.chainRule( cstart + jc, xx + 1 + jc, 1, 0, weight, bufstart, buffer );
      myvals.chainRule( wstart + iw, xx + 1 + jc, 1, 0, colvar, bufstart, buffer );
    }
  }
}

void GradientVessel::finish( const std::vector<double>& buffer ) {
  std::vector<double> val_interm( ncomponents*nweights );
  unsigned nder;
  if( getAction()->derivativesAreRequired() ) nder=getAction()->getNumberOfDerivatives();
  else nder=0;
  Matrix<double> der_interm( ncomponents*nweights, nder ); der_interm=0;

  if( isdens ) {
    for(unsigned iw=0; iw<nweights; ++iw) {
      val_interm[iw] = buffer[bufstart + 2*iw*(1+nder)];
      if( getAction()->derivativesAreRequired() ) {
        unsigned wstart = bufstart + 2*iw*(nder+1) + 1;
        for(unsigned jder=0; jder<nder; ++jder) der_interm( iw, jder ) += buffer[ wstart + jder ];
      }
    }
  } else {
    for(unsigned iw=0; iw<nweights; ++iw) {
      unsigned xx = (ncomponents+1)*iw;
      double ww=buffer[bufstart + xx*(1+nder)];
      for(unsigned jc=0; jc<ncomponents; ++jc) val_interm[ iw*ncomponents + jc ] = buffer[bufstart + (xx+1+jc)*(1+nder)] / ww;
      if( getAction()->derivativesAreRequired() ) {
        unsigned wstart = bufstart + xx*(nder+1) + 1;
        for(unsigned jc=0; jc<ncomponents; ++jc) {
          unsigned bstart = bufstart + ( xx + 1 + jc )*(nder+1) + 1;
          double val = buffer[bufstart + (nder+1)*(xx+1+jc)];
          for(unsigned jder=0; jder<nder; ++jder)
            der_interm( iw*ncomponents + jc, jder ) = (1.0/ww)*buffer[bstart + jder] - (val/(ww*ww))*buffer[wstart + jder];
        }
      }
    }
  }

  double tmp, diff2=0.0;

  if( getAction()->derivativesAreRequired() ) {
    Value* fval=getFinalValue();
    for(unsigned j=0; j<starts.size()-1; ++j) {
      for(unsigned bin=starts[j]; bin<starts[j+1]; ++bin) {
        for(unsigned jc=0; jc<ncomponents; ++jc) {
          if( bin==starts[j] ) {
            tmp=val_interm[(starts[j+1]-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
            for(unsigned jder=0; jder<nder; ++jder) {
              fval->addDerivative( jder, +2.0*tmp*der_interm( (starts[j+1]-1)*ncomponents + jc, jder) );
              fval->addDerivative( jder, -2.0*tmp*der_interm( bin*ncomponents + jc, jder ) );
            }
          } else {
            tmp=val_interm[(bin-1)*ncomponents + jc] - val_interm[bin*ncomponents + jc];
            for(unsigned jder=0; jder<nder; ++jder) {
              fval->addDerivative( jder, +2.0*tmp*der_interm( (bin-1)*ncomponents + jc, jder) );
              fval->addDerivative( jder, -2.0*tmp*der_interm( bin*ncomponents + jc, jder ) );
            }
          }
          diff2+=tmp*tmp;
        }
      }
    }
  } else {
    for(unsigned j=0; j<starts.size()-1; ++j) {
      for(unsigned bin=starts[j]; bin<starts[j+1]; ++bin) {
        for(unsigned jc=0; jc<ncomponents; ++jc) {
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

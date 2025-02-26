/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_colvar_ColvarOutputh
#define __PLUMED_colvar_ColvarOutputh

#include <vector>

#include "Pbc.h"
#include "View.h"
#include "View2D.h"
#include "Vector.h"
#include "Tensor.h"
#include "ColvarInput.h"

namespace PLMD {

class Colvar;

namespace colvar {

class ColvarOutput {
private:
  class DerivHelper {
  private:
    std::size_t nderivPerComponent;
    double* derivatives;
  public:
    DerivHelper(double* d, std::size_t n ) : nderivPerComponent(n), derivatives(d) {}
    View2D<double, helpers::dynamic_extent, 3> operator[](std::size_t i) {
      return { derivatives + i*nderivPerComponent, nderivPerComponent };
    }

    Vector getAtomDerivatives( std::size_t valueID, std::size_t atomID ) {
      std::size_t base = valueID*nderivPerComponent + 3*atomID;
      return Vector( derivatives[base], derivatives[base+1], derivatives[base+2] );
    }

    View<double,3> getDerivativesView( std::size_t valueID, std::size_t atomID) {
      std::size_t base = valueID*nderivPerComponent + 3*atomID;
      return { derivatives +base};
    }
  };
  class VirialHelper {
  private:
    std::size_t nderivPerComponent;
    double* derivatives;
  public:
    VirialHelper(double* d, std::size_t n ) : nderivPerComponent(n), derivatives(d) {}
    Tensor operator[](std::size_t i) const {
      std::size_t n=(i+1)*nderivPerComponent;
      return Tensor( derivatives[n-9],
                     derivatives[n-8],
                     derivatives[n-7],
                     derivatives[n-6],
                     derivatives[n-5],
                     derivatives[n-4],
                     derivatives[n-3],
                     derivatives[n-2],
                     derivatives[n-1] );
    }
    View<double,9> getView(std::size_t i) const {
      std::size_t n=(i+1)*nderivPerComponent;
      return {derivatives+n};
    }
    void set( std::size_t i, const Tensor& v ) {
      std::size_t n=(i+1)*nderivPerComponent;
      derivatives[n-9]=v[0][0]; derivatives[n-8]=v[0][1]; derivatives[n-7]=v[0][2];
      derivatives[n-6]=v[1][0]; derivatives[n-5]=v[1][1]; derivatives[n-4]=v[1][2];
      derivatives[n-3]=v[2][0]; derivatives[n-2]=v[2][1]; derivatives[n-1]=v[2][2];
    }
  };
  size_t ncomponents;
public:
  View<double> values;
  DerivHelper derivs;
  VirialHelper virial;
  ColvarOutput( View<double> v, std::size_t nderivPerComponent, double *derivatives ):
    ncomponents(v.size()),
    values(v),
    derivs(derivatives,nderivPerComponent),
    virial(derivatives,nderivPerComponent)
  {}

  static ColvarOutput createColvarOutput( std::vector<double>& v,
                                          std::vector<double>& d,
                                          Colvar* action );

  Vector getAtomDerivatives( std::size_t i, std::size_t a ) {
    return derivs.getAtomDerivatives(i,a);
  }

#pragma acc routine seq
  void setBoxDerivativesNoPbc( const ColvarInput& inpt ) {
    //now with no extra allocated memory: (actually I was searching for a bug...that was not here...)
    unsigned nat=inpt.pos.size();
    for(unsigned i=0; i<ncomponents; ++i) {
      auto v = virial.getView(i);
      LoopUnroller<9>::_zero(v.data());
      for(unsigned j=0; j<nat; j++) {
        const auto deriv =  derivs.getDerivativesView(i,j);
        v[0] -= inpt.pos[j][0]*deriv[0];
        v[1] -= inpt.pos[j][0]*deriv[1];
        v[2] -= inpt.pos[j][0]*deriv[2];

        v[3] -= inpt.pos[j][1]*deriv[0];
        v[4] -= inpt.pos[j][1]*deriv[1];
        v[5] -= inpt.pos[j][1]*deriv[2];

        v[6] -= inpt.pos[j][2]*deriv[0];
        v[7] -= inpt.pos[j][2]*deriv[1];
        v[8] -= inpt.pos[j][2]*deriv[2];
      }
    }
  }
};

} // namespace colvar
} //namespace PLMD
#endif // __PLUMED_colvar_ColvarOutputh

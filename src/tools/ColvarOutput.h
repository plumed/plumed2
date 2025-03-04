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
#ifndef __PLUMED_tools_ColvarOutput_h
#define __PLUMED_tools_ColvarOutput_h

#include <vector>

#include "Pbc.h"
#include "View.h"
#include "View2D.h"
#include "Vector.h"
#include "Tensor.h"

namespace PLMD {

class Colvar;

class ColvarOutput {
private:
  class DerivHelper {
  private:
    std::size_t nderivPerComponent;
    double* derivatives;
  public:
    DerivHelper(double* d, std::size_t n ) : nderivPerComponent(n), derivatives(d) {}
    View2D<double, helpers::dynamic_extent, 3> operator[](std::size_t i) {
      //the -9 is to "exclude" the virial (even if tecnically is still accessible)
      return { derivatives + i*nderivPerComponent, nderivPerComponent-9 };
    }

    Vector getAtomDerivatives( std::size_t valueID, std::size_t atomID ) {
      std::size_t base = valueID*nderivPerComponent + 3*atomID;
      return Vector( derivatives[base], derivatives[base+1], derivatives[base+2] );
    }

    View<double,3> getView( std::size_t valueID, std::size_t atomID) {
      std::size_t base = valueID*nderivPerComponent + 3*atomID;
      return View<double,3> { derivatives +base};
    }
  };
public:
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
      std::size_t n=(i+1)*nderivPerComponent-9;
      return View<double,9> {derivatives+n};
    }
    void set( std::size_t i, const Tensor& v ) {
      std::size_t n=(i+1)*nderivPerComponent;
      derivatives[n-9]=v[0][0];
      derivatives[n-8]=v[0][1];
      derivatives[n-7]=v[0][2];
      derivatives[n-6]=v[1][0];
      derivatives[n-5]=v[1][1];
      derivatives[n-4]=v[1][2];
      derivatives[n-3]=v[2][0];
      derivatives[n-2]=v[2][1];
      derivatives[n-1]=v[2][2];
    }
  };
  std::size_t ncomponents;
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
};

} //namespace PLMD

#endif // __PLUMED_tools_ColvarOutput_h

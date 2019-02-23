/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#ifndef __PLUMED_tools_ConjugateGradient_h
#define __PLUMED_tools_ConjugateGradient_h

#include "MinimiseBase.h"

namespace PLMD {

template <class FCLASS>
class ConjugateGradient : public MinimiseBase<FCLASS> {
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
  const unsigned ITMAX;
  const double EPS;
public:
  ConjugateGradient( FCLASS* funcc ) : MinimiseBase<FCLASS>(funcc), ITMAX(200), EPS(1E-10) {}
  void minimise( const double& ftol, std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
void ConjugateGradient<FCLASS>::minimise( const double& ftol, std::vector<double>& p, engf_pointer myfunc ) {
  std::vector<double> xi( p.size() ), g( p.size() ), h( p.size() );
  double fp = this->calcDerivatives( p, xi, myfunc );
  for(unsigned j=0; j<p.size(); ++j) { g[j] = -xi[j]; xi[j]=h[j]=g[j]; }

  for(unsigned its=0; its<ITMAX; ++its) {
    double fret=this->linemin( xi, p, myfunc );
    // The exit condition
    if( 2.0*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS)) { return; }
    fp = fret; double igeng = this->calcDerivatives( p, xi, myfunc );
    double ddg=0., gg=0.;
    for(unsigned j=0; j<p.size(); ++j) { gg += g[j]*g[j]; ddg += (xi[j]+g[j])*xi[j]; }

    if( gg==0.0 ) return;

    double gam=ddg/gg;
    for(unsigned j=0; j<p.size(); ++j) { g[j] = -xi[j]; xi[j]=h[j]=g[j]+gam*h[j]; }
  }
  plumed_merror("Too many interactions in conjugate gradient");
}

}

#endif

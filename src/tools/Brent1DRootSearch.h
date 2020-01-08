/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#ifndef __PLUMED_tools_Brent1DRootSearch_h
#define __PLUMED_tools_Brent1DRootSearch_h

#include "Tools.h"

#include <vector>
#include <string>

namespace PLMD {

/// A class for doing parabolic interpolation and minimisation of
/// 1D functions using Brent's method.
template <class FCLASS>
class Brent1DRootSearch {
private:
/// Has the minimum been bracketed
  bool bracketed;
/// The tolerance for the line minimiser
  double tol;
/// Maximum number of interactions in line minimiser
  const unsigned ITMAX;
/// A small number that protects against trying to achieve fractional
/// accuracy for a minimum that happens to be exactly zero
  const double EPS;
/// The factor by which to expand the range when bracketing
  const double EXPAND;
/// This is the type specifier for the function to minimise
  typedef double(FCLASS::*eng_pointer)( const double& val );
/// Three points bracketting the minimum and the corresponding function values
  double ax,bx,fa,fb;
/// The class containing the function we are trying to minimise
  FCLASS myclass_func;
public:
  explicit Brent1DRootSearch( const FCLASS& pf,  const double& t=3.0E-8 );
/// Bracket the minium
  void bracket( const double& ax, const double& xx, eng_pointer eng );
/// Find the minimum between two brackets
  double search( eng_pointer eng );
};

template <class FCLASS>
Brent1DRootSearch<FCLASS>::Brent1DRootSearch( const FCLASS& pf, const double& t ):
  bracketed(false),
  tol(t),
  ITMAX(100),
  EPS(3.0E-8),
  EXPAND(1.6),
  ax(0), bx(0),
  fa(0), fb(0),
  myclass_func(pf)
{
}

template <class FCLASS>
void Brent1DRootSearch<FCLASS>::bracket( const double& a, const double& b, eng_pointer eng ) {
  plumed_assert( a!=b ); ax=a; bx=b; fa=(myclass_func.*eng)(a); fb=(myclass_func.*eng)(b);
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) plumed_merror("input points do not bracket root");
  bracketed=true;
}

template <class FCLASS>
double Brent1DRootSearch<FCLASS>::search( eng_pointer eng ) {
  plumed_dbg_assert( bracketed );

  double cx=bx, d, e, min1, min2, fc=fb, p, q, r, s, tol1, xm;
  for(unsigned iter=0; iter<ITMAX; iter++) {
    if ( (fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0) ) { cx=ax; fc=fa; e=d=bx-ax; }
    if( fabs(fc) < fabs(fb) ) { ax=bx; bx=cx; cx=ax; fa=fb; fb=fc; fc=fa; }
    tol1=2*EPS*fabs(bx)+0.5*tol; xm=0.5*(cx-bx);
    if( fabs(xm) <= tol1 || fb == 0.0 ) return bx;
    if( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ) {
      s=fb/fa;
      if( ax==cx ) {
        p=2.0*xm*s; q=1.0-s;
      } else {
        q=fa/fc; r=fb/fc; p=s*(2.0*xm*q*(q-r)-(bx-ax)*(r-1.0)); q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p); min1=3.0*xm*q-fabs(tol1*q); min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d; d=p/q;
      } else {
        d=xm; e=d;
      }
    } else {
      d=xm; e=d;
    }
    ax=bx; fa=fb;
    if( fabs(d) > tol1 ) bx+=d;
    else if(xm<0 ) bx -= fabs(tol1); // SIGN(tol1,xm);
    else bx += tol1;
    fb = (myclass_func.*eng)(bx);
  }

  plumed_merror("Too many interactions in zbrent");
}

}
#endif

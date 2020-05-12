/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_tools_Minimise1DBrent_h
#define __PLUMED_tools_Minimise1DBrent_h

#include "Tools.h"

#include <vector>
#include <string>

namespace PLMD {

/// A class for doing parabolic interpolation and minimisation of
/// 1D functions using Brent's method.
template <class FCLASS>
class Minimise1DBrent {
private:
/// Has the minimum been bracketed
  bool bracketed;
/// Has the function been minimised
  bool minimised;
/// The tolerance for the line minimiser
  double tol;
/// The default ration by which successive intervals are magnified
  const double GOLD;
/// The maximum magnification allowed for a parabolic fit step
  const double GLIMIT;
/// Use to prevent any possible division by zero
  const double TINY;
/// Maximum number of interactions in line minimiser
  const unsigned ITMAX;
/// The value of the golden ratio
  const double CGOLD;
/// A small number that protects against trying to achieve fractional
/// accuracy for a minimum that happens to be exactly zero
  const double ZEPS;
/// This is the type specifier for the function to minimise
  typedef double(FCLASS::*eng_pointer)( const double& val );
/// Three points bracketting the minimum and the corresponding function values
  double ax,bx,cx,fa,fb,fc,fmin;
/// The class containing the function we are trying to minimise
  FCLASS myclass_func;
public:
  explicit Minimise1DBrent( const FCLASS& pf,  const double& t=3.0E-8 );
/// Bracket the minium
  void bracket( const double& ax, const double& xx, eng_pointer eng );
/// Find the minimum between two brackets
  double minimise( eng_pointer eng );
/// Return the value of the function at the minimum
  double getMinimumValue() const ;
};

template <class FCLASS>
Minimise1DBrent<FCLASS>::Minimise1DBrent( const FCLASS& pf, const double& t ):
  bracketed(false),
  minimised(false),
  tol(t),
  GOLD(1.618034),
  GLIMIT(100.0),
  TINY(1.0E-20),
  ITMAX(100),
  CGOLD(0.3819660),
  ZEPS(epsilon*1.0E-3),
  ax(0),bx(0),cx(0),
  fa(0),fb(0),fc(0),
  fmin(0),
  myclass_func(pf)
{
}

template <class FCLASS>
void Minimise1DBrent<FCLASS>::bracket( const double& a, const double& b, eng_pointer eng ) {
  ax=a; bx=b; double fu;
  fa=(myclass_func.*eng)(ax); fb=(myclass_func.*eng)(bx);
  if( fb>fa ) {
    double tmp;
    tmp=ax; ax=bx; bx=tmp;
    tmp=fa; fa=fb; fb=tmp;
  }
  cx=bx+GOLD*(bx-ax);
  fc=(myclass_func.*eng)(cx);
  while ( fb > fc ) {
    double r=(bx-ax)*(fb-fc);
    double q=(bx-cx)*(fb-fa);
    double u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*(fabs(q-r)>TINY?fabs(q-r):TINY)*(q-r>=0?1:-1));
    double ulim=bx+GLIMIT*(cx-bx);
    if((bx-u)*(u-cx) > 0.0 ) {
      fu=(myclass_func.*eng)(u);
      if( fu < fc ) {
        ax=bx; bx=u; fa=fb; fb=fu; bracketed=true; return;
      } else if( fu > fb ) {
        cx=u; fc=fu; bracketed=true; return;
      }
      u=cx+GOLD*(cx-bx); fu=(myclass_func.*eng)(u);
    } else if((cx-u)*(u-ulim) > 0.0 ) {
      fu=(myclass_func.*eng)(u);
      if( fu<fc ) {
        bx=cx; cx=u; u+=GOLD*(u-bx);
        fb=fc; fc=fu; fu=(myclass_func.*eng)(u);
      }
    } else if( (u-ulim)*(ulim-cx) >= 0.0 ) {
      u=ulim;
      fu=(myclass_func.*eng)(u);
    } else {
      u=cx+GOLD*(cx-bx);
      fu=(myclass_func.*eng)(u);
    }
    ax=bx; bx=cx; cx=u;
    fa=fb; fb=fc; fc=fu;
  }
  bracketed=true;
}

template <class FCLASS>
double Minimise1DBrent<FCLASS>::minimise( eng_pointer eng ) {
  plumed_dbg_assert( bracketed );

  double a,b,d=0.0,etemp,fu,fv,fw,fx;
  double p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx );
  b=(ax >= cx ? ax : cx );
  x=w=v=bx;
  fw=fv=fx=(myclass_func.*eng)(x);
  for(unsigned iter=0; iter<ITMAX; ++iter) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if( fabs(x-xm) <= (tol2-0.5*(b-a))) {
      fmin=fx; minimised=true; return x;
    }
    if( fabs(e) > tol1 ) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if( q > 0.0 ) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if( fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) ) {
        d = CGOLD*(e=(x >= xm ? a-x : b-x ));
      } else {
        d=p/q; u=x+d;
        if(u-a < tol2 || b-u < tol2 ) d=(xm-x>=0?fabs(tol1):-fabs(tol1));
      }
    } else {
      d=CGOLD*(e=( x >= xm ? a-x : b-x ));
    }
    if( fabs(d)>=tol1) u=x+d; else u=x+(d>=0?fabs(tol1):-fabs(tol1));
    fu=(myclass_func.*eng)(u);
    if( fu <= fx ) {
      if( u >= x ) a=x; else b=x;
      v=w; fv=fw;
      w=x; fw=fx;
      x=u; fx=fu;
    } else {
      if( u < x ) a=u; else b=u;
      if( fu <=fw || w==x ) {
        v=w; w=u; fv=fw; fw=fu;
      } else if( fu <= fv || v==x || v==w ) {
        v=u; fv=fu;
      }
    }
  }
  plumed_merror("Too many interactions in brent");
}

template <class FCLASS>
double Minimise1DBrent<FCLASS>::getMinimumValue() const {
  plumed_dbg_assert( minimised );
  return fmin;
}

}
#endif

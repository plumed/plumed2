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
#ifndef __PLUMED_tools_MinimiseBase_h
#define __PLUMED_tools_MinimiseBase_h

#include "Minimise1DBrent.h"

namespace PLMD {

template <class FCLASS>
class F1dim {
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der ) const ;
  typedef double(FCLASS::*engfnc_pointer)( const std::vector<double>& p, std::vector<double>& der ) ;
/// Pointer to the vector containing an initial position on the vector
  const std::vector<double>& p;
/// The direction of the vector we are minimising along
  const std::vector<double>& dir;
/// Tempory vector that holds a point at which we want to calculate the energy
  std::vector<double> pt;
/// Vector that holds the derivatives at the point at which we calculate the energy (these are not used)
  std::vector<double> fake_der;
/// Class containging the function in the class
  FCLASS* func;
/// Member of class that calculates the energy we are trying to mnimise
  engf_pointer calc;
/// Member of class that calcualtes the energy we are trying to minimise
  engfnc_pointer calc2;
public:
  explicit F1dim( const std::vector<double>& pp, const std::vector<double>& dd, FCLASS* ff, engf_pointer cc, engfnc_pointer cc2 );
/// Calculate the energy at \f$\mathbf{p} + xt*\mathbf{dir}\f$
  double getEng( const double& xt );
};

template <class FCLASS>
F1dim<FCLASS>::F1dim( const std::vector<double>& pp, const std::vector<double>& dd, FCLASS* ff, engf_pointer cc, engfnc_pointer cc2 ):
  p(pp),
  dir(dd),
  pt(pp.size()),
  fake_der(pp.size()),
  func(ff),
  calc(cc),
  calc2(cc2)
{
  plumed_assert( calc || calc2 );
}

template <class FCLASS>
double F1dim<FCLASS>::getEng( const double& xt ) {
  for(unsigned j=0; j<pt.size(); ++j) pt[j] = p[j] + xt*dir[j];
  if( calc ) return (func->*calc)(pt,fake_der);
  return (func->*calc2)(pt,fake_der);
}

template <class FCLASS>
class MinimiseBase {
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
/// The class that calculates the energy given a position
  FCLASS* myclass_func;
protected:
/// This calculates the derivatives at a point
  double calcDerivatives( const std::vector<double>& p, std::vector<double>& der, engf_pointer myfunc );
public:
  explicit MinimiseBase( FCLASS* funcc ) : myclass_func(funcc) {}
/// This is the line minimiser
  double linemin( const std::vector<double>& dir, std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
double MinimiseBase<FCLASS>::linemin( const std::vector<double>& dir, std::vector<double>& p, engf_pointer myfunc ) {
  // Construct the object that turns points on a line into vectors
  F1dim<FCLASS> f1dim( p, dir, myclass_func, NULL, myfunc );

  // Construct an object that will do the line search for the minimum
  Minimise1DBrent<F1dim<FCLASS> > bb(f1dim);

  // This does the actual line minimisation
  double ax=0.0, xx=1.0;
  bb.bracket( ax, xx, &F1dim<FCLASS>::getEng );
  double xmin=bb.minimise( &F1dim<FCLASS>::getEng );
  for(unsigned i=0; i<p.size(); ++i) p[i] += xmin*dir[i];
  return bb.getMinimumValue();
}

template <class FCLASS>
double MinimiseBase<FCLASS>::calcDerivatives( const std::vector<double>& p, std::vector<double>& der, engf_pointer myfunc ) {
  return (myclass_func->*myfunc)( p, der );
}

}
#endif

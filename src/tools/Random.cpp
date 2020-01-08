/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "Random.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <iterator>
#include <functional>

namespace PLMD {

const double Random::fact=5.9604644775390625e-8;     /* 1 / 2^24  */
const double Random::EPS=3.0e-16;
const double Random::AM=1.0/IM;
const double Random::RNMX=(1.0-EPS); // 1.0-EPS;
const std::string Random::noname="noname";

Random::Random(const std::string & name):
  switchGaussian(false),
  saveGaussian(0.0),
// this is required because if a Random object is created during
// initialization than Random::noname could still be initialized.
// In practice: without this it is not possible to declare
// a static Random object without enforcing the order of the
// static constructors.
  name(&name!=&noname?name:"noname")
{
  iy=0;
  iv[0]=0;
  setSeed(0);
}

void Random::setSeed(int idum_) {
  if(idum_>0) idum_=-idum_;
  idum=idum_;
  incPrec=false;
}

double Random::RandU01 ()
{
  if (incPrec)
    return U01d();
  else
    return U01();
}


double Random::U01d ()
{
  double u;
  u = U01();
  u += U01() * fact;
  return (u < 1.0) ? u : (u - 1.0);
}

double Random::U01() {
  int j,k;
  double temp;
  if (idum <= 0 || !iy) {
    if (-idum < 1) idum=1;
    else idum = -idum;
    for (j=NTAB+7; j>=0; j--) {
      k=idum/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

void Random::WriteStateFull(std::ostream & out)const {
  out<<name<<std::endl;
  out<<idum<<" "<<iy;
  for(int i=0; i<NTAB; i++) {
    out<<" "<<iv[i];
  };
  out<<" "<<switchGaussian;
  out<<" "<<saveGaussian;
  out<<std::endl;
}

void Random::ReadStateFull (std::istream & in) {
  getline(in,name);
  in>>idum>>iy;
  for (int i = 0; i < NTAB; i++) in>>iv[i];
  in>>switchGaussian;
  in>>saveGaussian;
}

void Random::toString(std::string & str)const {
  std::ostringstream ostr;
  ostr<<idum<<"|"<<iy;
  for(int i=0; i<NTAB; i++) {
    ostr<<"|"<<iv[i];
  };
  str=ostr.str();
}

void Random::fromString(const std::string & str) {
  std::string s=str;
  for(unsigned i=0; i<s.length(); i++) if(s[i]=='|') s[i]=' ';
  std::istringstream istr(s.c_str());
  istr>>idum>>iy;
  for (int i = 0; i < NTAB; i++) istr>>iv[i];
}

// This allows to have the same stream of random numbers
// with different compilers:
#ifdef __INTEL_COMPILER
#pragma intel optimization_level 0
#endif

double Random::Gaussian() {
  double v1,v2,rsq;
  if(switchGaussian) {
    switchGaussian=false;
    return saveGaussian;
  }
  while(true) {
    v1=2.0*RandU01()-1.0;
    v2=2.0*RandU01()-1.0;
    rsq=v1*v1+v2*v2;
    if(rsq<1.0 && rsq>0.0) break;
  }
  double fac=sqrt(-2.*std::log(rsq)/rsq);
  saveGaussian=v1*fac;
  switchGaussian=true;
  return v2*fac;
}

void Random::Shuffle(std::vector<unsigned>& vec) {
  std::iterator_traits<std::vector<unsigned>::iterator >::difference_type i, n;
  n = vec.end() - vec.begin();
  for(i=n-1; i>0; --i) {
    std::swap(vec[i], vec[(int)round(RandU01() * IM) % i]);
  }
}


}

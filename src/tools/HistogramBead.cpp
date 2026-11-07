/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "HistogramBead.h"
#include <vector>
#include "Tools.h"

#pragma GCC diagnostic error "-Wswitch"

/*
IMPORTANT NOTE FOR DEVELOPERS:

If you add a new type of function in this file please add documentation for your new switching function type in function/Between.cpp

*/
namespace PLMD {

constexpr double DP2CUTOFF=6.25;
//Non periodic constructor
HistogramBead::HistogramBead(KernelType kt,const  double l,
                             const double h,
                             const double w):
  type(kt),
  periodicity(Periodicity::notperiodic) {
  set(l,h,w);
}
//periodic constructor
HistogramBead::HistogramBead(KernelType kt,
                             double mlow,
                             double mhigh,const  double l,
                             const double h,
                             const double w):
  type(kt) {
  isPeriodic(mlow,mhigh);
  set(l,h,w);
}

HistogramBead::HistogramBead(const HistogramBead&)=default;
HistogramBead::HistogramBead(HistogramBead&&)=default;
HistogramBead& HistogramBead::operator=(const HistogramBead&)=default;
HistogramBead& HistogramBead::operator=(HistogramBead&&)=default;

std::string HistogramBead::description() const {
  std::ostringstream ostr;
  switch (type) {
  case KernelType::gaussian:
    ostr<<"between "<<lowb<<" and "
        <<highb<<" width of gaussian window equals "<<width;
    break;
  case KernelType::triangular:
    ostr<<"between "<<lowb<<" and "
        <<highb<<" width of triangular window equals "<<width;
    break;
  }
  return ostr.str();
}

void HistogramBead::set( const std::string& params, std::string& errormsg ) {
  std::vector<std::string> data=Tools::getWords(params);
  if(data.size()<1) {
    errormsg="No input has been specified";
    return;
  }

  std::string name=data[0];

  if(name=="GAUSSIAN") {
    type=KernelType::gaussian;
    cutoff=std::sqrt(2.0*DP2CUTOFF);
  } else if(name=="TRIANGULAR") {
    type=KernelType::triangular;
    cutoff=1.0;
  } else {
    plumed_merror("cannot understand kernel type " + name );
  }

  double smear;
  bool found_r=Tools::parse(data,"LOWER",lowb);
  if( !found_r ) {
    errormsg="Lower bound has not been specified use LOWER";
  }
  found_r=Tools::parse(data,"UPPER",highb);
  if( !found_r ) {
    errormsg="Upper bound has not been specified use UPPER";
  }
  if( lowb>=highb ) {
    errormsg="Lower bound is higher than upper bound";
  }

  smear=0.5;
  Tools::parse(data,"SMEAR",smear);
  width=smear*(highb-lowb);
}

void HistogramBead::set(const  double l, const double h, const double w) {
  lowb=l;
  highb=h;
  width=w;
  switch (type) {
  case KernelType::gaussian : {
    cutoff=std::sqrt(2.0*DP2CUTOFF);
  }
  break;
  case KernelType::triangular : {
    cutoff=1.0;
  }
  }
}

HistogramBead::KernelType HistogramBead::getKernelType( const std::string& ktype ) {
  if(ktype=="gaussian") {
    return KernelType::gaussian;
  } else if(ktype=="triangular") {
    return KernelType::triangular;
  } else {
    plumed_merror("cannot understand kernel type " + ktype );
  }
}

double HistogramBead::calculate(const double x, double& df ) const {
  double res=0.0;
  switch (type) {
  case KernelType::gaussian : {
    const double lowB = difference( x, lowb ) / ( std::sqrt(2.0) * width );
    const double upperB = difference( x, highb ) / ( std::sqrt(2.0) * width );
    df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( std::sqrt(PLMD::twopi)*width );
    res = 0.5*( erf( upperB ) - erf( lowB ) );
  }
  break;
  case KernelType::triangular : {
    const double lowB = ( difference( x, lowb ) / width );
    const double upperB = ( difference( x, highb ) / width );
    df=0.0;
    if( std::fabs(lowB)<1. ) {
      df = (1 - std::fabs(lowB)) / width;
    }
    if( std::fabs(upperB)<1. ) {
      df -= (1 - std::fabs(upperB)) / width;
    }
    if (upperB<=-1. || lowB >=1.) {
      break;
    }
    const double ia = (lowB>-1.0) ? lowB : -1.0;
    const double ib = (upperB<1.0) ? upperB : 1.0;
    res = (ib*(2.-std::fabs(ib))-ia*(2.-std::fabs(ia)))*0.5;
  }
  }
  return res;
}

double HistogramBead::calculateWithCutoff( double x, double& df ) const {
  return calculateWithCutoff( x, lowb, highb, width, df );
}

double HistogramBead::calculateWithCutoff(double x,
    double l,
    double u,
    double s,
    double& df ) const {
  double lowB = difference( x, l ) / s;
  double upperB = difference( x, u ) / s;
  if( upperB<=-cutoff || lowB>=cutoff ) {
    df=0.0;
    return 0.0;
  }
  double f=0.0;
  switch (type) {
  case KernelType::gaussian : {
    lowB /= std::sqrt(2.0);
    upperB /= std::sqrt(2.0);
    df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( std::sqrt(2*pi)*s );
    f = 0.5*( erf( upperB ) - erf( lowB ) );
  }
  break;
  case KernelType::triangular : {
    df=0;
    if( std::fabs(lowB)<1. ) {
      df = (1.0 - std::fabs(lowB)) / s;
    }
    if( std::fabs(upperB)<1. ) {
      df -= (1.0 - std::fabs(upperB)) / s;
    }
    if (upperB<=-1. || lowB >=1.) {
      f=0.0;
    } else {
      double ia, ib;
      if( lowB>-1.0 ) {
        ia=lowB;
      } else {
        ia=-1.0;
      }
      if( upperB<1.0 ) {
        ib=upperB;
      } else {
        ib=1.0;
      }
      f = (ib*(2.-std::fabs(ib))-ia*(2.-std::fabs(ia)))*0.5;
    }
  }
  }
  return f;
}

double HistogramBead::lboundDerivative( const double x ) const {
  switch (type) {
  case KernelType::gaussian : {
    double lowB = difference( x, lowb ) / ( std::sqrt(2.0) * width );
    return exp( -lowB*lowB ) / ( std::sqrt(2*pi)*width );
  }
  break;
  case  KernelType::triangular : {
    assert(false);
    //TODO: something about this:
    // double lowB = fabs( difference( x, lowb ) / width );
    // if( lowB<1 ) {
    //   return ( 1 - (lowB) ) / 2*width;
    // } else {
    //   return 0;
    // }
  }
  }
  //warning::unreachable
  return 0.0;
}

double HistogramBead::uboundDerivative( const double x ) const {
  switch (type) {
  case KernelType::gaussian : {
    double upperB = difference( x, highb ) / ( std::sqrt(2.0) * width );
    return exp( -upperB*upperB ) / ( std::sqrt(2*pi)*width );
  }
  break;
  case KernelType::triangular : {
    assert(false);
    //TODO: something about this:
    // double upperB = fabs( difference( x, highb ) / width );
    // if( upperB<1 ) {
    //   return ( 1 - (upperB) ) / 2*width;
    // } else {
    //   return 0;
    // }
  }
  }
  //warning::unreachable
  return 0.0;
}

inline
double HistogramBead::difference( double d1, const double d2 ) const {
  switch (periodicity) {
  case Periodicity::notperiodic :
    break;
  case Periodicity::periodic: {
    // Make sure the point is in the target range
    d1*=inv_max_minus_min;
    d1=Tools::pbc(d1);
    d1*=max_minus_min;
  }
  }
  return d2-d1;
}

}

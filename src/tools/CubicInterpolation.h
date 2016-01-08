/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#ifndef __PLUMED_tools_CubicInterpolation_h
#define __PLUMED_tools_CubicInterpolation_h

#include <vector>
#include "Matrix.h"
#include "core/Value.h"

namespace PLMD {

// Abstract base class for cubic interpolation
class CInterpolation {
private:
  unsigned bold;
  Matrix<double> splinepoints;
  unsigned search1( const unsigned& kk, const double& x, const unsigned& jold ) const ;
protected:
  std::vector<double> lb, ub;
  std::vector<unsigned> np, stride;
  unsigned findBox( const std::vector<double>& pos );
  double getPointSpacing( const unsigned dir, const unsigned k ) const ;
  double getCrossTermDenominator( const unsigned i, const unsigned j ) const ;
public:
  CInterpolation( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax );
  virtual ~CInterpolation();
  unsigned getNumberOfSplinePoints() const ;
  void getNumbersOfPoints( std::vector<unsigned>& nspline ) const ;
  void getSplinePoint( const unsigned nn, std::vector<double>& pp ) const ;
  void getGridBoundaries( std::vector<double>& gmin, std::vector<double>& gmax ) const ;
  virtual void set_table( const std::vector<Value>& ff )=0;
  virtual double get_fdf( const std::vector<double>& pos )=0;
};

inline
unsigned CInterpolation::getNumberOfSplinePoints() const {
  return splinepoints.nrows();
}

inline
void CInterpolation::getSplinePoint( const unsigned nn, std::vector<double>& pp ) const {
  plumed_dbg_assert( nn<splinepoints.nrows() && pp.size()==np.size() );
  for(unsigned i=0;i<np.size();++i) pp[i]=splinepoints(nn,i); 
}

inline
double CInterpolation::getPointSpacing( const unsigned dir, const unsigned k ) const {
  unsigned i=k*stride[dir];
  return splinepoints(i+stride[dir], dir) - splinepoints(i, dir); 
}

inline
double CInterpolation::getCrossTermDenominator( const unsigned i, const unsigned j ) const {
  plumed_dbg_assert( splinepoints.ncols()==2 );
  unsigned iplus, iminus; iplus=(i+1)*stride[0]; iminus=(i-1)*stride[0];
  return ( splinepoints(iplus,0) - splinepoints(iminus,0) ) * ( splinepoints(iplus+j+1,1) - splinepoints(iplus+j-1,1) );
}

inline
void CInterpolation::getGridBoundaries( std::vector<double>& gmin, std::vector<double>& gmax ) const {
  getSplinePoint( 0, gmin ); getSplinePoint( splinepoints.nrows()-1, gmax );
}

class InterpolateCubic : public CInterpolation {
private:
  std::vector<double> clist;
public:
  InterpolateCubic( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax );
  void set_table( const std::vector<Value>& ff );
  double get_fdf( const std::vector<double>& pos );
};

class InterpolateBicubic : public CInterpolation {
private:
  Matrix<int> wt;
  std::vector<double> t1, t2;
  Matrix<double> dcross;
  std::vector<double> clist;
  void IBicCoeff( const std::vector<double>& y, const std::vector<double>& dy1, const std::vector<double>& dy2,
                  const std::vector<double>& d2y12, const double& d1, const double& d2, Matrix<double>& c );
public:
  InterpolateBicubic( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax );
  void set_table( const std::vector<Value>& ff );
  double get_fdf( const std::vector<double>& pos );
};

}


#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
  
   See http://www.plumed-code.org for more information.
  
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
#include "BicubicInterpolation.h"

namespace PLMD {
namespace vesselbase {

BicubicInterpolation::BicubicInterpolation( GridVesselBase* gg, const unsigned dstart ):
InterpolationBase(gg,dstart)
{
  plumed_dbg_assert( getDimension()==2 );
  static int wt_d[16*16]=
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
     2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
     0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
     -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
     9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
     -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
     2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
     -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
     4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};

  unsigned l=0; wt.resize(16,16); t1.resize(16); t2.resize(16);
  for (unsigned i=0;i<16;i++) for (unsigned j=0;j<16;j++){ wt(i,j)=wt_d[l++]; }
  std::vector<unsigned> nbin(2); getNumberOfGridPoints( nbin );
  dcross.resize( nbin[0], nbin[1] ); clist.resize( 16*getNumberOfSplinePoints() );
} 

void BicubicInterpolation::setInterpolationTables(){
  std::vector<unsigned> np(2); std::vector<double> der(2);

  // Get the denominator for the calculation of the mixed derivative
  getGridPointSpacing(der); double denom = (2*der[0])*(2*der[1]);

   for(unsigned i=1;i<dcross.nrows()-1;++i){
      for(unsigned j=1;j<dcross.ncols()-1;++j){
          np[0] = i+1; np[1] = j+1; double vpipj = getValue( np );
          np[0] = i-1; np[1] = j-1; double vmijm = getValue( np );
          np[0] = i+1; np[1] = j-1; double vpijm = getValue( np );
          np[0] = i-1; np[1] = j+1; double vmijp = getValue( np );
          dcross(i,j) = ( vpipj + vmijm - vpijm - vmijp ) / denom;
      }
  }

  double d1=der[0], d2=der[1]; Matrix<double> tc(4,4);
  std::vector<double> y(4), dy1(4), dy2(4), d2y12(4);

  unsigned pij=0; 
  for (unsigned i=0;i<dcross.nrows()-1;++i){
      for (unsigned j=0; j<dcross.ncols()-1;++j){
         np[0]=i; np[1]=j; y[0]=getValueAndDerivatives( np, der ); dy1[0]=der[0]; dy2[0]=der[1]; d2y12[0] = dcross( i, j );
         np[0]=i+1; np[1]=j; y[1]=getValueAndDerivatives( np, der ); dy1[1]=der[0]; dy2[1]=der[1]; d2y12[1] = dcross( i+1, j );
         np[0]=i+1; np[1]=j+1; y[2]=getValueAndDerivatives( np, der ); dy1[2]=der[0]; dy2[2]=der[1]; d2y12[2] = dcross( i+1, j+1 );
         np[0]=i; np[1]=j+1; y[3]=getValueAndDerivatives( np, der ); dy1[3]=der[0]; dy2[3]=der[1]; d2y12[3] = dcross( i, j+1 );
         IBicCoeff( y, dy1, dy2, d2y12, d1, d2, tc);

         np[0]=i; np[1]=j; pij=16*getBoxIndex( np ); 
         for(unsigned k=0; k<4; ++k){ for(unsigned n=0; n<4; ++n){ clist[pij++]=tc(k,n); } }
      }
  }
}

void BicubicInterpolation::IBicCoeff( const std::vector<double>& y, const std::vector<double>& dy1, const std::vector<double>& dy2,
                                      const std::vector<double>& d2y12, const double& d1, const double& d2, Matrix<double>& c ){
  double xx, d1d2=d1*d2;
  for(unsigned i=0;i<4;i++){ t1[i] = y[i]; t1[i+4] = dy1[i]*d1; t1[i+8] = dy2[i]*d2; t1[i+12] = d2y12[i]*d1d2; }
  for(unsigned i=0;i<16;i++){ xx=0.0; for(unsigned k=0;k<16;k++){ xx += wt(i,k)*t1[k]; } t2[i]=xx; }
  unsigned l=0; for(unsigned i=0;i<4;i++){ for(unsigned j=0;j<4;j++){ c(i,j)=t2[l++]; } }
}

double BicubicInterpolation::interpolateFunction( const unsigned& mybox, const std::vector<double>& dd ){
  plumed_dbg_assert( dd.size()==2 );

  double t=dd[0], u=dd[1];
  //faster access by pointer arithmetic (dirty dirty dirty) 
  double *cbase=&clist[(mybox+1)*16-1], *c3, *c2, *c1, *c0; double f=0.;
  for (int i=3; i>=0; i--) {    // Note to self - this has to be an int as unsigned cannot be less than zero - duh!!
      c3=cbase; c2=c3-1; c1=c2-1; c0=c1-1; cbase=c0-1;
      f = t*f + ( ( (*c3)*u + (*c2) )*u + (*c1) )*u + (*c0);
  }
  return f;
}

}
}

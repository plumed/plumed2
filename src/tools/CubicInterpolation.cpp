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
#include "CubicInterpolation.h"

namespace PLMD {

CInterpolation::CInterpolation( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax ) : 
bold(0)
{
  plumed_assert( fmin.size()==dd.size() && fmax.size()==dd.size() );
 
  np.resize( dd.size() ); stride.resize( dd.size() ); unsigned totalp=1;
  for(unsigned i=0;i<dd.size();++i){ np[i]=dd[i]; stride[dd.size()-1-i]=totalp; totalp*=np[i]; }

  unsigned ii,kk; 
  splinepoints.resize( totalp, np.size() );
  std::vector<double> delr( np.size() );
  for(unsigned j=0;j<np.size();++j) delr[j] = ( fmax[j] - fmin[j] )/static_cast<double>(np[j]-1);  

  for(unsigned i=0;i<totalp;++i){
     ii=i;
     for(unsigned j=0;j<np.size();++j){
        kk=std::floor( double(ii) / double(stride[j]) );
        ii-=kk*stride[j];
        splinepoints(i,j)=fmin[j] + kk*delr[j];
     }
     plumed_assert(ii==0);
  }
  lb.resize( np.size() ); ub.resize( np.size() );
}

CInterpolation::~CInterpolation(){
  splinepoints.resize(0,0); lb.resize(0); ub.resize(0); np.resize(0); stride.resize(0);
}

void CInterpolation::getNumbersOfPoints( std::vector<unsigned>& nspline ) const {
  nspline.resize( np.size() );
  for(unsigned i=0;i<np.size();++i) nspline[i]=np[i];
}

unsigned CInterpolation::findBox( const std::vector<double>& pos ){
  plumed_dbg_massert( pos.size()==np.size(), "position size does not match the size of the grid");

  unsigned jold, ccf_box, bnew=0;
  for(unsigned i=0;i<np.size();++i){
     jold=static_cast<int>( std::floor( double(bold)/double(stride[i]) ) );
     bold-=jold*stride[i];
     ccf_box=search1( i, pos[i], jold );
     bnew+=ccf_box; 
  }
  plumed_dbg_assert( bold==0 ); bold=bnew;
  for(unsigned i=0;i<np.size();++i){ lb[i]=splinepoints(bold,i); ub[i]=splinepoints(bold+stride[i],i); }
  return bold;
}

unsigned CInterpolation::search1( const unsigned& kk, const double& x, const unsigned& jold ) const {
    int inc=stride[kk], jl=jold*stride[kk], ju=(jold+1)*stride[kk], jm; 
    if ( x>=splinepoints(jl,kk) && x<splinepoints( ju, kk ) ) return jl;
    else {
        if( x>=splinepoints(jl, kk ) ){
            while(true){
                ju=jl+inc;
                if( x<splinepoints( ju, kk ) ) break;
                else if( ju>=(np[kk]-1)*inc ){ju=(np[kk]-1)*inc; break; } 
                jl=ju; 
            }
        }
        else{
            ju=jl;
            while(true){
                jl=jl-inc;
                if( x>=splinepoints( jl, kk ) ) break;
                else if( jl<=0 ){ jl=0; break; }
                ju=jl; 
            }
        }
    }
    while( ju-jl>inc ){
      jm = (ju+jl) / (2*inc) ;
      if ( x>splinepoints(jm*inc,kk) ) jl=jm*inc; else ju=jm*inc;
    }
    plumed_dbg_assert( jl%stride[kk]==0 && ju==jl+stride[kk] );
    return jl;
}

InterpolateCubic::InterpolateCubic( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax ) :
CInterpolation(dd,fmin,fmax)
{
  plumed_massert(np.size()==1,"should be one dimensional data");
  clist.resize( 4*np[0] );
}

void InterpolateCubic::set_table( const std::vector<Value>& ff ){
  plumed_assert( getNumberOfSplinePoints()==ff.size() ); 
  plumed_assert( ff[0].getNumberOfDerivatives()==1 );

  double d1, norm; unsigned pij;
  for(unsigned i=0;i<np[0]-1;++i){
      d1 = getPointSpacing( 0, i );   
      norm=(d1*d1)/6.0; pij=i*4;
      clist[pij]=ff[i].get(); pij++;
      clist[pij]=ff[i+1].get(); pij++;
      clist[pij]=ff[i].getDerivative(0)*norm; pij++;
      clist[pij]=ff[i+1].getDerivative(0)*norm;   
  }
}

double InterpolateCubic::get_fdf( const std::vector<double>& pos ){
  plumed_dbg_assert( pos.size()==1 );
  
  unsigned mybox=findBox( pos );
  double d1=ub[0] - lb[0]; 
  double b=( pos[0] - lb[0] ) / d1, a=( ub[0] - pos[0] ) / d1;
  
  double *cbase=&clist[(mybox*4)+3], *c3=cbase-1, *c2=c3-1, *c1=c2-1;
  double f=a*(*c1) + b*(*c2) + (a*a*a-a)*(*c3) + (b*b*b-b)*(*cbase);
  return f;  
}

InterpolateBicubic::InterpolateBicubic( const std::vector<unsigned>& dd, const std::vector<double>& fmin, const std::vector<double>& fmax ) :
CInterpolation(dd,fmin,fmax)
{
  plumed_massert(np.size()==2,"should be two dimensional data");
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

  // This is to set up the coefficient matrix
  unsigned l=0; wt.resize(16,16); t1.resize(16); t2.resize(16);
  for (unsigned i=0;i<16;i++) for (unsigned j=0;j<16;j++){ wt(i,j)=wt_d[l++]; }
  // Resize everything
  dcross.resize( np[0], np[1] ); clist.resize( np[0] * np[1] * 4 * 4 );
}

void InterpolateBicubic::set_table( const std::vector<Value>& ff ){
  plumed_assert( getNumberOfSplinePoints()==ff.size() ); 
  plumed_assert( ff[0].getNumberOfDerivatives()==2 );

  dcross=0.0; unsigned iplus, iminus;
  for(unsigned i=1;i<np[0]-1;++i){
      iplus=(i+1)*stride[0]; iminus=(i-1)*stride[0];
      for(unsigned j=1;j<np[1]-1;++j){
          dcross(i,j) = ( ff[iplus+j+1].get() + ff[iminus+j-1].get() - ff[iplus+j-1].get() - ff[iminus+j+1].get() ) /
                          getCrossTermDenominator( i, j );
      }
  } 

  double d1, d2; Matrix<double> tc(4,4);
  std::vector<double> y(4), dy1(4), dy2(4), d2y12(4);

  unsigned pij=0; unsigned ipos;
  for (unsigned i=0;i<np[0]-1;++i){
      ipos=i*stride[0]; d1 = getPointSpacing( 0, i );    
      for (unsigned j=0; j<np[1]-1;++j){
         d2 = getPointSpacing( 1, j );                   
         y[0] = ff[ipos+j].get(); y[1] = ff[ipos+stride[0]+j].get(); y[2] = ff[ipos+stride[0]+j+1].get(); y[3] = ff[ipos+j+1].get();
         dy1[0] = ff[ipos+j].getDerivative(0); dy1[1] = ff[ipos+stride[0]+j].getDerivative(0); 
         dy1[2] = ff[ipos+stride[0]+j+1].getDerivative(0); dy1[3] = ff[ipos+j+1].getDerivative(0);
         dy2[0] = ff[ipos+j].getDerivative(1); dy2[1] = ff[ipos+stride[0]+j].getDerivative(1); 
         dy2[2] = ff[ipos+stride[0]+j+1].getDerivative(1); dy2[3] = ff[ipos+j+1].getDerivative(1);
         d2y12[0] = dcross( i, j ); d2y12[1] = dcross( i+1, j ); d2y12[2] = dcross( i+1, j+1 ); d2y12[3] = dcross( i, j+1 );
         IBicCoeff( y, dy1, dy2, d2y12, d1, d2, tc);

         pij=( ipos+j )*16;
         for(unsigned k=0; k<4; ++k){ for(unsigned n=0; n<4; ++n){ clist[pij++]=tc(k,n); } }
      }
  }
}

void InterpolateBicubic::IBicCoeff( const std::vector<double>& y, const std::vector<double>& dy1, const std::vector<double>& dy2, 
                                    const std::vector<double>& d2y12, const double& d1, const double& d2, Matrix<double>& c ){
  double xx, d1d2=d1*d2;
  for(unsigned i=0;i<4;i++){ t1[i] = y[i]; t1[i+4] = dy1[i]*d1; t1[i+8] = dy2[i]*d2; t1[i+12] = d2y12[i]*d1d2; }
  for(unsigned i=0;i<16;i++){ xx=0.0; for(unsigned k=0;k<16;k++){ xx += wt(i,k)*t1[k]; } t2[i]=xx; }
  unsigned l=0; for(unsigned i=0;i<4;i++){ for(unsigned j=0;j<4;j++){ c(i,j)=t2[l++]; } }
}

double InterpolateBicubic::get_fdf( const std::vector<double>& pos ){

   plumed_dbg_assert( pos.size()==2 );
   unsigned mybox=findBox( pos );
   double d1 = ub[0] - lb[0], d2 = ub[1] - lb[1];
   double t = (pos[0] - lb[0]) / d1, u = (pos[1] - lb[1]) / d2;

   //faster access by pointer arithmetic (dirty dirty dirty)
   double *cbase=&clist[(mybox+1)*16-1], *c3, *c2, *c1, *c0;

   double f=0.;
   for (int i=3; i>=0; i--) {    // Note to self - this has to be an int as unsigned cannot be less than zero - duh!!
       c3=cbase; c2=c3-1; c1=c2-1; c0=c1-1; cbase=c0-1;
       f= t*f + ( ( (*c3)*u + (*c2) )*u + (*c1) )*u + (*c0);
   }
   delete cbase; delete c3; delete c2; delete c1; delete c0;
   return f;
}

}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "KernelFunctions.h"

namespace PLMD {

//+PLUMEDOC INTERNAL kernelfunctions
/*
Constructing histograms is something you learnt to do relatively early in life. You perform an experiment a number of times, 
count the number of times each result comes up and then draw a bar graph that describes how often each of the results came up.  
This only works when there are a finite number of possible results.  If the result a number between 0 and 1 the bar chart is 
less easy to draw as there are as many possible results as there are numbers between zero and one - an infinite number.  
To resolve this problem we replace probability, \f$P\f$ with probability density, \f$\pi\f$, and write the probability of getting 
a number between \f$a\f$ and \f$b\f$ as:

\f[
P = \int_{a}^b \textrm{d}x \pi(x)   
\f]

To calculate probability densities from a set of results we use a process called kernel density estimation.  
Histograms are accumulated by adding up kernel functions, \f$K\f$, with finite spatial extent, that integrate to one.
These functions are centered on each of the \f$n\f$-dimensional data points, \f$\mathbf{x}_i\f$. The overall effect of this
is that each result we obtain in our experiments contributes to the probability density in a finite sized region of the space.  

Expressing all this mathematically in kernel density estimation we write the probability density as:

\f[
\pi(\mathbf{x}) =  \sum_i K\left[ (\mathbf{x} - \mathbf{x}_i)^T \Sigma (\mathbf{x} - \mathbf{x}_i) \right]
\f]

where \f$\Sigma\f$ is an \f$n \times n\f$ matrix called the bandwidth that controls the spatial extent of 
the kernel. Whenever we accumulate a histogram (e.g. in \ref HISTOGRAM or in \ref METAD) we use this 
technique. 

There is thus some flexibility in the particular function we use for \f$K[\mathbf{r}]\f$ in the above.
The following variants are available.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td> TYPE </td> <td> FUNCTION </td> 
</tr> <tr> 
<td> gaussian </td> <td> \f$f(r) = \frac{1}{(2 \pi)^{n} \sqrt{|\Sigma^{-1}|}} \exp\left(-0.5 r^2 \right)\f$ </td>
</tr> <tr>
<td> triangular </td> <td> \f$f(r) = \frac{3}{V} ( 1 - | r | )H(1-|r|) \f$ </td>
</tr> <tr>
<td> uniform </td> <td> \f$f(r) = \frac{1}{V}H(1-|r|)\f$ </td>
</tr>
</table>

In the above \f$H(y)\f$ is a function that is equal to one when \f$y>0\f$ and zero when \f$y \le 0\f$. \f$n\f$ is
the dimensionality of the vector \f$\mathbf{x}\f$ and \f$V\f$ is the volume of an elipse in an \f$n\f$ dimensional
space which is given by:

\f{eqnarray*}{
V &=& | \Sigma^{-1} | \frac{ \pi^{\frac{n}{2}} }{\left( \frac{n}{2} \right)! } \qquad \textrm{for even} \quad n \\
V &=& | \Sigma^{-1} | \frac{ 2^{\frac{n+1}{2}} \pi^{\frac{n-1}{2}} }{ n!! } 
\f}

In \ref METAD the normalization constants are ignored so that the value of the function at \f$r=0\f$ is equal
to one.  In addition in \ref METAD we must be able to differentiate the bias in order to get forces.  This limits 
the kernels we can use in this method.
*/
//+ENDPLUMEDOC

KernelFunctions::KernelFunctions( const std::vector<double>& at, const std::vector<double>& sig, const std::string& type, const double& w, const bool& norm ):
center(at),
width(sig)
{
  unsigned ncv=center.size();
  if( width.size()==ncv ) diagonal=true;
  else if( width.size()==(ncv*(ncv+1))/2 ) diagonal=false;
  else plumed_massert(0,"specified sigma is neither diagonal or full covariance matrix");

  // Setup the kernel type
  if(type=="GAUSSIAN" || type=="gaussian"){
      ktype=gaussian;
  } else if(type=="UNIFORM" || type=="uniform"){
      ktype=uniform;
  } else if(type=="TRIANGULAR" || type=="triangular"){
      ktype=triangular;
  }

  if( norm ){
    double det;
    if(diagonal){
       det=1; for(unsigned i=0;i<width.size();++i) det*=width[i];
    } else {
       unsigned ncv=ndim(); 
       Matrix<double> mymatrix( getMatrix() ), myinv( ncv, ncv );
       Invert(mymatrix,myinv); double logd;
       logdet( myinv, logd );
       det=exp(logd);
    }
    double volume;
    if( ktype==gaussian ){
       for(unsigned i=0;i<width.size();++i) det*=width[i];
       volume=pow( 2*pi, 0.5*ncv ) * pow( det, 0.5 );
    } else if( ktype==uniform || ktype==triangular ){
       if( ncv%2==1 ){
          double dfact=1;
          for(unsigned i=1;i<ncv;i+=2) dfact*=static_cast<double>(i);
          volume=( pow( pi, (ncv-1)/2 ) ) * ( pow( 2., (ncv+1)/2 ) ) / dfact;
       } else {
          double fact=1.;
          for(unsigned i=1;i<ncv/2;++i) fact*=static_cast<double>(i);
          volume=pow( pi,ncv/2 ) / fact;
       }
       if(ktype==uniform) volume*=det;
       else if(ktype==triangular) volume*=det / 3.;
    } else {
       plumed_merror("not a valid kernel type");
    } 
    height=w / volume;  
  } else {
    height=w;
  }
}

double KernelFunctions::getCutoff( const double& width ) const {
  double DP2CUTOFF=6.25;
  if( ktype==gaussian ) return sqrt(2.0*DP2CUTOFF)*width;
  else if(ktype==triangular ) return width;
  else if(ktype==uniform) return width;
  else plumed_merror("No valid kernel type");
  return 0.0;
}

std::vector<unsigned> KernelFunctions::getSupport( const std::vector<double>& dx ) const {
  plumed_assert( ndim()==dx.size() );
  std::vector<unsigned> support( dx.size() );
  if(diagonal){
     for(unsigned i=0;i<dx.size();++i) support[i]=static_cast<unsigned>(ceil( getCutoff(width[i])/dx[i] ));
  } else {
     unsigned ncv=ndim(); 
     Matrix<double> mymatrix( getMatrix() ), myinv( ncv,ncv );
     Invert(mymatrix,myinv);
     Matrix<double> myautovec(ncv,ncv); std::vector<double> myautoval(ncv);  
     diagMat(myinv,myautoval,myautovec);
     for(unsigned i=0;i<dx.size();++i){
         double extent=fabs(sqrt(myautoval[0])*myautovec(i,0)); 
         support[i]=static_cast<unsigned>(ceil( getCutoff( extent )/dx[i] ));
     }
  }
  return support;
}

double KernelFunctions::evaluate( const std::vector<Value*>& pos, std::vector<double>& derivatives, bool usederiv ) const {
  plumed_assert( pos.size()==ndim() && derivatives.size()==ndim() );
  if( usederiv ) plumed_massert( ktype!=uniform, "step function can not be differentiated" ); 

  double r2=0;
  if(diagonal){ 
     for(unsigned i=0;i<ndim();++i){
         derivatives[i]=pos[i]->difference( center[i] ) / width[i]; 
         r2+=derivatives[i]*derivatives[i];
     }
  } else {
     Matrix<double> mymatrix( getMatrix() ); 
     for(unsigned i=0;i<mymatrix.nrows();++i){
        double dp_i, dp_j; derivatives[i]=0;
        dp_i=pos[i]->difference( center[i] ); 
        for(unsigned j=0;j<mymatrix.ncols();++j){
          if(i==j) dp_j=dp_i;
          else dp_j=pos[j]->difference( center[j] );

          derivatives[i]+=mymatrix(i,j)*dp_j;
          r2+=dp_i*dp_j*mymatrix(i,j);
        }
     }
  }
  double kderiv, kval;
  if(ktype==gaussian){
     kval=height*exp(-0.5*r2); kderiv=-kval;
  } else {
     double r=sqrt(r2);
     if(ktype==triangular){
        if( r<1.0 ){
            if(r==0) kderiv=0;
            kderiv=-1; kval=height*( 1. - fabs(r) );
        } else {
            kval=0.; kderiv=0.;
        }
     } else if(ktype==uniform){ 
        kderiv=0.;
        if(r<1.0) kval=height;
        else kval=0;
     } else {
         plumed_merror("Not a valid kernel type");
     }
     kderiv*=height / r ;
  }  
  for(unsigned i=0;i<ndim();++i) derivatives[i]*=kderiv;
  return kval;
}

KernelFunctions* KernelFunctions::read( PlumedIFile* ifile, const std::vector<std::string>& valnames ){
  std::string sss; ifile->scanField("multivariate",sss);
  std::vector<double> cc( valnames.size() ), sig;
  if( sss=="true" ){
     sig.resize( valnames.size() );
     for(unsigned i=0;i<valnames.size();++i){
         ifile->scanField(valnames[i],cc[i]);
         ifile->scanField("sigma_"+valnames[i],sig[i]);
     }
  } else if( sss=="false" ){
     unsigned ncv=valnames.size();
     sig.resize( (ncv*(ncv+1))/2 );
     Matrix<double> upper(ncv,ncv), lower(ncv,ncv);
     for(unsigned i=0;i<ncv;++i){
         ifile->scanField(valnames[i],cc[i]);
         for(unsigned j=0;j<ncv-i;j++){ ifile->scanField("sigma_" +valnames[j+i] + "_" + valnames[j], lower(j+i,j) ); upper(j,j+i)=lower(j+i,j); }
     }
     Matrix<double> mymult( ncv, ncv ), invmatrix(ncv,ncv);
     mult(lower,upper,mymult); Invert( mymult, invmatrix );
     unsigned k=0; 
     for(unsigned i=0;i<ncv;i++){
         for(unsigned j=i;j<ncv;j++){ sig[k]=invmatrix(i,j); k++; }
     }
  } else {
      plumed_merror("multivariate flag should equal true or false");
  } 
  double h; ifile->scanField("height",h);
  return new KernelFunctions( cc, sig, "gaussian", h, false );
}

}

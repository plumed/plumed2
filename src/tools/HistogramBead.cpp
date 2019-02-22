/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include <limits>
#include "Tools.h"
#include "Keywords.h"

namespace PLMD {

//+PLUMEDOC INTERNAL histogrambead
/*
A function that can be used to calculate whether quantities are between fixed upper and lower bounds.

If we have multiple instances of a variable we can estimate the probability density function
for that variable using a process called kernel density estimation:

\f[
P(s) = \sum_i K\left( \frac{s - s_i}{w} \right)
\f]

In this equation \f$K\f$ is a symmetric function that must integrate to one that is often
called a kernel function and \f$w\f$ is a smearing parameter.  From a probability density function calculated using
kernel density estimation we can calculate the number/fraction of values between an upper and lower
bound using:

\f[
w(s) = \int_a^b \sum_i K\left( \frac{s - s_i}{w} \right)
\f]

All the input to calculate a quantity like \f$w(s)\f$ is generally provided through a single
keyword that will have the following form:

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ SMEAR=\f$\frac{w}{b-a}\f$}

This will calculate the number of values between \f$a\f$ and \f$b\f$.  To calculate
the fraction of values you add the word NORM to the input specification.  If the
function keyword SMEAR is not present \f$w\f$ is set equal to \f$0.5(b-a)\f$. Finally,
type should specify one of the kernel types that is present in plumed. These are listed
in the table below:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td> TYPE </td> <td> FUNCTION </td>
</tr> <tr>
<td> GAUSSIAN </td> <td> \f$\frac{1}{\sqrt{2\pi}w} \exp\left( -\frac{(s-s_i)^2}{2w^2} \right)\f$ </td>
</tr> <tr>
<td> TRIANGULAR </td> <td> \f$ \frac{1}{2w} \left( 1. - \left| \frac{s-s_i}{w} \right| \right) \quad \frac{s-s_i}{w}<1 \f$ </td>
</tr>
</table>

Some keywords can also be used to calculate a discrete version of the histogram.  That
is to say the number of values between \f$a\f$ and \f$b\f$, the number of values between
\f$b\f$ and \f$c\f$ and so on.  A keyword that specifies this sort of calculation would look
something like

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ NBINS=\f$n\f$ SMEAR=\f$\frac{w}{n(b-a)}\f$}

This specification would calculate the following vector of quantities:

\f[
w_j(s) = \int_{a + \frac{j-1}{n}(b-a)}^{a + \frac{j}{n}(b-a)} \sum_i K\left( \frac{s - s_i}{w} \right)
\f]

*/
//+ENDPLUMEDOC

void HistogramBead::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
}

HistogramBead::HistogramBead():
  init(false),
  lowb(0.0),
  highb(0.0),
  width(0.0),
  cutoff(std::numeric_limits<double>::max()),
  type(gaussian),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
}

std::string HistogramBead::description() const {
  std::ostringstream ostr;
  ostr<<"betweeen "<<lowb<<" and "<<highb<<" width of gaussian window equals "<<width;
  return ostr.str();
}

void HistogramBead::generateBins( const std::string& params, std::vector<std::string>& bins ) {
  std::vector<std::string> data=Tools::getWords(params);
  plumed_massert(data.size()>=1,"There is no input for this keyword");

  std::string name=data[0];

  unsigned nbins; std::vector<double> range(2); std::string smear;
  bool found_nb=Tools::parse(data,"NBINS",nbins);
  plumed_massert(found_nb,"Number of bins in histogram not found");
  bool found_r=Tools::parse(data,"LOWER",range[0]);
  plumed_massert(found_r,"Lower bound for histogram not specified");
  found_r=Tools::parse(data,"UPPER",range[1]);
  plumed_massert(found_r,"Upper bound for histogram not specified");
  plumed_massert(range[0]<range[1],"Range specification is dubious");
  bool found_b=Tools::parse(data,"SMEAR",smear);
  if(!found_b) { Tools::convert(0.5,smear); }

  std::string lb,ub; double delr = ( range[1]-range[0] ) / static_cast<double>( nbins );
  for(unsigned i=0; i<nbins; ++i) {
    Tools::convert( range[0]+i*delr, lb );
    Tools::convert( range[0]+(i+1)*delr, ub );
    bins.push_back( name + " " +  "LOWER=" + lb + " " + "UPPER=" + ub + " " + "SMEAR=" + smear );
  }
  plumed_assert(bins.size()==nbins);
}

void HistogramBead::set( const std::string& params, std::string& errormsg ) {
  std::vector<std::string> data=Tools::getWords(params);
  if(data.size()<1) {
    errormsg="No input has been specified";
    return;
  }

  std::string name=data[0]; const double DP2CUTOFF=6.25;
  if(name=="GAUSSIAN") { type=gaussian; cutoff=sqrt(2.0*DP2CUTOFF); }
  else if(name=="TRIANGULAR") { type=triangular; cutoff=1.; }
  else plumed_merror("cannot understand kernel type " + name );

  double smear;
  bool found_r=Tools::parse(data,"LOWER",lowb);
  if( !found_r ) errormsg="Lower bound has not been specified use LOWER";
  found_r=Tools::parse(data,"UPPER",highb);
  if( !found_r ) errormsg="Upper bound has not been specified use UPPER";
  if( lowb>=highb ) errormsg="Lower bound is higher than upper bound";

  smear=0.5; Tools::parse(data,"SMEAR",smear);
  width=smear*(highb-lowb); init=true;
}

void HistogramBead::set( double l, double h, double w) {
  init=true; lowb=l; highb=h; width=w;
  const double DP2CUTOFF=6.25;
  if( type==gaussian ) cutoff=sqrt(2.0*DP2CUTOFF);
  else if( type==triangular ) cutoff=1.;
  else plumed_error();
}

void HistogramBead::setKernelType( const std::string& ktype ) {
  if(ktype=="gaussian") type=gaussian;
  else if(ktype=="triangular") type=triangular;
  else plumed_merror("cannot understand kernel type " + ktype );
}

double HistogramBead::calculate( double x, double& df ) const {
  plumed_dbg_assert(init && periodicity!=unset );
  double lowB, upperB, f;
  if( type==gaussian ) {
    lowB = difference( x, lowb ) / ( sqrt(2.0) * width );
    upperB = difference( x, highb ) / ( sqrt(2.0) * width );
    df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2*pi)*width );
    f = 0.5*( erf( upperB ) - erf( lowB ) );
  } else if( type==triangular ) {
    lowB = ( difference( x, lowb ) / width );
    upperB = ( difference( x, highb ) / width );
    df=0;
    if( fabs(lowB)<1. ) df = (1 - fabs(lowB)) / width;
    if( fabs(upperB)<1. ) df -= (1 - fabs(upperB)) / width;
    if (upperB<=-1. || lowB >=1.) {
      f=0.;
    } else {
      double ia, ib;
      if( lowB>-1.0 ) { ia=lowB; } else { ia=-1.0; }
      if( upperB<1.0 ) { ib=upperB; } else { ib=1.0; }
      f = (ib*(2.-fabs(ib))-ia*(2.-fabs(ia)))*0.5;
    }
  } else {
    plumed_merror("function type does not exist");
  }
  return f;
}

double HistogramBead::calculateWithCutoff( double x, double& df ) const {
  plumed_dbg_assert(init && periodicity!=unset );

  double lowB, upperB, f;
  lowB = difference( x, lowb ) / width ; upperB = difference( x, highb ) / width;
  if( upperB<=-cutoff || lowB>=cutoff ) { df=0; return 0; }

  if( type==gaussian ) {
    lowB /= sqrt(2.0); upperB /= sqrt(2.0);
    df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2*pi)*width );
    f = 0.5*( erf( upperB ) - erf( lowB ) );
  } else if( type==triangular ) {
    df=0;
    if( fabs(lowB)<1. ) df = (1 - fabs(lowB)) / width;
    if( fabs(upperB)<1. ) df -= (1 - fabs(upperB)) / width;
    if (upperB<=-1. || lowB >=1.) {
      f=0.;
    } else {
      double ia, ib;
      if( lowB>-1.0 ) { ia=lowB; } else { ia=-1.0; }
      if( upperB<1.0 ) { ib=upperB; } else { ib=1.0; }
      f = (ib*(2.-fabs(ib))-ia*(2.-fabs(ia)))*0.5;
    }
  } else {
    plumed_merror("function type does not exist");
  }
  return f;
}

double HistogramBead::lboundDerivative( const double& x ) const {
  if( type==gaussian ) {
    double lowB = difference( x, lowb ) / ( sqrt(2.0) * width );
    return exp( -lowB*lowB ) / ( sqrt(2*pi)*width );
  } else if ( type==triangular ) {
    plumed_error();
//      lowB = fabs( difference( x, lowb ) / width );
//      if( lowB<1 ) return ( 1 - (lowB) ) / 2*width;
//      else return 0;
  } else {
    plumed_merror("function type does not exist");
  }
  return 0;
}

double HistogramBead::uboundDerivative( const double& x ) const {
  plumed_dbg_assert(init && periodicity!=unset );
  if( type==gaussian ) {
    double upperB = difference( x, highb ) / ( sqrt(2.0) * width );
    return exp( -upperB*upperB ) / ( sqrt(2*pi)*width );
  } else if ( type==triangular ) {
    plumed_error();
//      upperB = fabs( difference( x, highb ) / width );
//      if( upperB<1 ) return ( 1 - (upperB) ) / 2*width;
//      else return 0;
  } else {
    plumed_merror("function type does not exist");
  }
  return 0;
}

}

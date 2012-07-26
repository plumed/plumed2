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
#include "FieldVessel.h"
#include "ActionWithDistribution.h"

namespace PLMD {

class DistributionField : public FieldVessel {
private:
  double sigma, fsigma2, fnorm;
  Value tmpder;
public:
  static void reserveKeyword( Keywords& keys );
  DistributionField( const VesselOptions& da );
  void calculateEnergy( const unsigned&, const std::vector<double>&, const Value& , Value& , std::vector<Value>& );
  void printKeywords();
  double calculateField( const std::vector<double>& pp );
  void calculateFieldDerivatives( const std::vector<double>& pp, std::vector<double>& tmpforce );
};

PLUMED_REGISTER_VESSEL(DistributionField,"DISTRIBUTION")

void DistributionField::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","DISTRIBUTION","create a field cv from the distribution of the collective variables. From a distribution of collective variables "
                                  "we calculate a 1D-field using "
                                  " \\f$\\psi(z)=\\frac{\\sum_i G(z-s_i,\\sigma)}{\\int \\textrm{d}z \\sum_i G(z-s_i,\\sigma)}\\f$ where \\f$ G(z-s_i,\\sigma)\\f$ "
                                  " is a normalized Gaussian function with standard deviation \\f$\\sigma\\f$ centered at the value of the \\f$i\\f$th collective variable"
                                  " \\f$s_i\\f$.  In other words we use a histogram calculated from the values of all the constituent collective variables. "
                                  + FieldVessel::documentation() );
}

void DistributionField::printKeywords(){
  Keywords field_keys;
  FieldVessel::getKeywords( field_keys );
  field_keys.add("compulsory","SIGMA","the ammount by which to smear each of the quantities in your distribution"); 
  field_keys.print(log);
}

DistributionField::DistributionField( const VesselOptions& da ):
FieldVessel(da)
{
  plumed_massert(get_Ndx()==1, "Can only create 1D distributions");

  std::vector<std::string> data=Tools::getWords(da.parameters);
  bool found_sigma=Tools::parse( data, "SIGMA", sigma );
  if(!found_sigma) error("did not find SIGMA keyword");
  fsigma2=sigma*sigma; fnorm = 1.0 / ( sqrt(2*pi)*sigma );

  setLabel("distribution");
  log.printf("  %s.distribution can be used to reference a 1 dimensional field that describes the histogram of values\n",(getAction()->getLabel()).c_str() );
  std::vector<unsigned> nspline( get_nspline() ); std::vector<double> min( get_min() ),max( get_max() );
  log.printf("  histogram is calculated at %d evenly space points between %f and %f. Each value is represented by a Gaussian of width %f\n",nspline[0],min[0],max[0],sigma);  
}

void DistributionField::calculateEnergy( const unsigned& jcv, const std::vector<double>& at, const Value& value_in, Value& stress, std::vector<Value>& der ){
  plumed_assert( at.size()==1 && stress.getNumberOfDerivatives()==1 );

  double dd, ss, du; dd=value_in.difference( at[0] ); 
  ss=fnorm*exp( -(dd*dd)/(2*fsigma2) ); du=dd/fsigma2;
  stress.set( ss ); 
  stress.clearDerivatives();
  stress.addDerivative( 0, -du*ss );

  if( !mergeBeforeInterpolation() ){
     for(unsigned icv=0;icv<get_NdX();++icv){ der[icv].set(0); der[icv].clearDerivatives(); }
     der[jcv].set(du*ss);
     der[jcv].addDerivative( 0, du*ss*(1./(fsigma2*du) - du) );
  } else {
     double d2x=1./(fsigma2*du) - du;
     getAction()->mergeDerivatives( jcv, value_in, du*ss, tmpder );
     for(unsigned i=0;i<tmpder.getNumberOfDerivatives();++i){
         der[i].set( tmpder.getDerivative(i) );
         der[i].clearDerivatives();
         der[i].addDerivative(0, d2x*tmpder.getDerivative(i) );
     }
  }
}

double DistributionField::calculateField( const std::vector<double>& pp ){
  return getFieldAt( pp );
}

void DistributionField::calculateFieldDerivatives( const std::vector<double>& pp, std::vector<double>& tmpforce ){
  getDerivativesAt( pp, tmpforce ); 
}

}

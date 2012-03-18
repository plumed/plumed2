#include "DistributionFunctions.h"

namespace PLMD {

std::string cvdens::documentation(){
  std::ostringstream ostr;
  ostr<<"To make this quantity continuous it is calculated using \\f$ a = \\frac{\\sum_{i=1}^N s_i w(x_i)w(y_i)w(z_i)}{\\sum_{i=1}^N w(x_i)w(y_i)w(z_i)}\\f$ ";
  ostr<<"where the sum is over the collective variables \\f$ s_i \\f$ which can be thought to be at \\f$ (x_i,y_i,z_i)\\f$. The three \\f$w(x)\\f$ functions ";
  ostr<<"specify the extent of the subregion in which we are calculating the average CV in each direction. To make these assignment functions continuous they ";
  ostr<<"are calculated using "<<HistogramBead::documentation(true)<<".  If no parameters for the \\f$ w(d) \\f$ function in a particular direction is specified ";
  ostr<<"then the subcell is assumed to incorporate the entirity of the box in that direction."; 
  return ostr.str();
}

cvdens::cvdens( const std::string& parameters ) :
DistributionFunction(parameters)
{
  HistogramBead xbin; xbin.set(parameters, "X"); 
  if ( xbin.hasBeenSet() ){ beads.push_back(xbin); dir.push_back(0); }
  HistogramBead ybin; ybin.set(parameters, "Y");
  if ( ybin.hasBeenSet() ){ beads.push_back(ybin); dir.push_back(1); }
  HistogramBead zbin; zbin.set(parameters, "Z");
  if ( zbin.hasBeenSet() ){ beads.push_back(zbin); dir.push_back(2); }  
  plumed_assert( beads.size()>=1 );

  isDensity=(parameters.find("density")!=std::string::npos);
  addAccumulator( true );    // This holds the numerator - value of cv times "location in box" 
  addAccumulator( true );    // This holds the denominator - number of atoms in box
}

std::string cvdens::message(){
  std::ostringstream ostr;
  ostr<<"average value of cv for ";
  for(unsigned i=0;i<dir.size();++i){
     if( dir[i]==0 ) ostr<<"x "<<beads[i].description();
     if( dir[i]==1 ) ostr<<"y "<<beads[i].description();
     if( dir[i]==2 ) ostr<<"z "<<beads[i].description();
     if( dir.size()>1 && i!=(dir.size()-1) ) ostr<<", "; 
  }
  return ostr.str();
}

void cvdens::calculate( Value* value_in, std::vector<Value>& aux ){
  plumed_assert( aux.size()==3 );

  double f, df;
  copyValue( 0, value_in );
  for(unsigned i=0;i<dir.size();++i){
     f=beads[i].calculate( aux[ dir[i] ].get(), df );
     aux[ dir[i] ].chainRule(df); aux[ dir[i] ].set(f);
     multiplyValue( 0, &aux[ dir[i] ] );
     if( i==0 ){ copyValue( 1, &aux[ dir[i] ] ); }
     else{ multiplyValue( 1, &aux[ dir[i] ] ); }
  }
}

void cvdens::finish( Value* value_out ){
  if( !isDensity ){ quotient( getPntrToAccumulator(0), getPntrToAccumulator(1), value_out ); }
  else{ extractDerivatives( 1, value_out ); value_out->set( getPntrToAccumulator(1)->get() ); }
}

}

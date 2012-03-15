#include "DistributionFunctions.h"

namespace PLMD {

void within::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par WITHIN/HISTOGRAM \n\n";
  ostr<<"Calculates the number of colvars that are within a certain range.  To make this quantity continuous it is calculated using: \n";
  ostr<<"\\f[\n";
  ostr<<"S = \\sum_i \\int_a^b G( s_i, \\sigma*(b-a) )";
  ostr<<"\\f]\n";
  ostr<<"where \\f$G( s_i, \\sigma )\\f$ is a normalized Gaussian function of width \\f$\\sigma\\f$ centered at the value of the colvar \\f$s_i\\f$.  \\f$a\\f$ and \\f$b\\f$ are\n";
  ostr<<"the lower and upper bounds of the range of interest respectively.  The values of \\f$a\\f$ and \\f$b\\f$ must be specified using the WITHIN keyword.  If this keyword\n";
  ostr<<"has three input numbers then the third is assumed to be the value of \\f$\\sigma\\f$.  You can specify that you want to investigate multiple rangles by using multiple instances\n"; 
  ostr<<"of the WITHIN keyword (WITHIN1,WITHIN2 etc).  Alternatively, if you want to calculate a discretized distribution function you can use the HISTOGRAM keyword in \n";
  ostr<<"tandem with the RANGE keyword.  HISTOGRAM takes as input the number of bins in your distribution and (optionally) the value of the smearing parameter \\f$\\sigma\\f$.\n";
  ostr<<"RANGE takes the upper and lower bound of the histogram.  The interval between the upper and lower bound specified using RANGE is then divided into \n";
  ostr<<"equal-width bins.\n";  
  docs=ostr.str();
}

within::within( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{ 
  if( parameters.size()==3 ){
     Tools::convert(parameters[0],a);
     Tools::convert(parameters[1],b); 
     Tools::convert(parameters[2],sigma);
  } else if(parameters.size()==2){
     Tools::convert(parameters[0],a);
     Tools::convert(parameters[1],b);
     sigma=0.5;
  } else {
     error("WITHIN keyword takes two or three arguments");
  }
  if(a>=b) error("For WITHIN keyword upper bound is greater than lower bound");
  hist.set(a,b,sigma);
  addAccumulator( true );
}

std::string within::message(){
  std::ostringstream ostr;
  ostr<<"number of values between "<<a<<" and "<<b<<" The value of the smearing parameter is "<<sigma;
  return ostr.str();
}

void within::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  double df, f; f=hist.calculate( value_in->get() , df );
  chainRule(0, df); setValue(0, f);
}

void within::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}

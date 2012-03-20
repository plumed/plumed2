#include "DistributionFunctions.h"

namespace PLMD {

std::string min::documentation(){
  std::ostringstream ostr;
  ostr<<"To make this quantity continuous the minimum is calculated using ";
  ostr<<"\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ ";
  ostr<<"The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)";
  return ostr.str();
}

min::min( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::vector<std::string> data=Tools::getWords(parameters);
  if( data.size()!=1 ){ error("There should only be a value for beta in the input to MIN"); return; }
  bool found_beta=Tools::parse(data,"BETA",beta);
  if (!found_beta){ error("No value for BETA specified in call to MIN"); return; } 
  addAccumulator( true );
}

std::string min::message(){
  std::ostringstream ostr;
  ostr<<"the minimum value. beta is equal to "<<beta; 
  return ostr.str();
}

void min::printKeywords( Log& log ){
  Keywords mkeys; 
  mkeys.add("compulsory","BETA","the value of beta for the equation in the manual");
  mkeys.print(log);
}

std::string min::getLabel(){
  return "min";
}

void min::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in );
  double p, df, tmp; p=value_in->get();
  tmp=exp( beta/p ); df=tmp/(p*p); 
  chainRule( 0, df ); setValue( 0, tmp ); 
}

void min::finish( Value* value_out ){
  double a=getPntrToAccumulator(0)->get();
  if( a==0 ) plumed_massert(0,"Error in calculating minimum.  There is probably an NL_CUTOFF in your input that is too small");
  extractDerivatives( 0, value_out );
  double df, dist=beta/std::log( a ); df=dist*dist/a;
  value_out->chainRule(df); value_out->set(dist);
}

}

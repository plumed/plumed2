#include "DistributionFunctions.h"

namespace PLMD {

void min::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par MIN"<<std::endl;
  ostr<<std::endl;
  ostr<<"Calculate the minimum value.  To make this quantity continuous the minimum is calculated using:"<<std::endl;
  ostr<<std::endl;
  ostr<<"\\f["<<std::endl;
  ostr<<"\\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) }"<<std::endl;
  ostr<<"\\f]"<<std::endl;
  ostr<<std::endl;
  ostr<<"The keyword MIN=\\f$\\beta\\f$ tells plumed to calculate this quantity and sets \\f$\\beta\\f$ to an appropriate value.  Once"<<std::endl;
  ostr<<"calcualted the final value is referneced using label.min"<<std::endl;  
  docs=ostr.str();
}

min::min( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  if( parameters.size()==1 ){
     Tools::convert(parameters[0],beta);
  } else {
     error("MIN keyword should have one argument - the value for beta");
  }
  addAccumulator( true );
}

std::string min::message(){
  std::ostringstream ostr;
  ostr<<"the minimum value. beta is equal to "<<beta; 
  return ostr.str();
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

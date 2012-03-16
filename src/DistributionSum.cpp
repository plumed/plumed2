#include "DistributionFunctions.h"

namespace PLMD {

void sum::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par SUM \n\n";
  ostr<<"Calculate the sum of all the colvars in the distribution.  Once calculated the final value is referenced\n";
  ostr<<"using lable.sum.\n";
  docs=ostr.str();
}

sum::sum( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],nval); 
  addAccumulator( true );
  addAccumulator( false );
}

std::string sum::message(){
  std::ostringstream ostr;
  ostr<<"the sum of all the values"; 
  return ostr.str();
}

void sum::calculate( Value* value_in, std::vector<Value>& aux ){
  copyDerivatives( 0, value_in ); 
  setValue( 0, value_in->get() );
  setValue( 1, 1.0 );
}

void sum::finish( Value* value_out ){
  if ( getPntrToAccumulator(1)->get()!=nval ) printf("WARNING: A neighbor list is causing discontinuities in a sum");
  extractDerivatives( 0, value_out ); value_out->set( getPntrToAccumulator(0)->get() );
}

}

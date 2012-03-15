#include "DistributionFunctions.h"

namespace PLMD {

void mean::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par AVERAGE \n\n";
  ostr<<"Calculate the average of all the colvars in the distribution.  Once calculated the final value is referenced\n";
  ostr<<"using lable.average.\n";
  docs=ostr.str();
}

mean::mean( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],nval); 
  addAccumulator( true );
  addAccumulator( false );
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

void mean::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  setValue( 1, 1.0 );
}

void mean::finish( Value* value_out ){
  if ( getPntrToAccumulator(1)->get()!=nval ) printf("WARNING: A neighbor list is causing discontinuities in an average");
  extractDerivatives( 0, value_out );
  value_out->chainRule(1.0/nval); value_out->set(getPntrToAccumulator(0)->get()/nval); 
}

}

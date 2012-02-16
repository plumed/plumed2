#include "DistributionFunctions.h"

namespace PLMD {

DistributionFunction::DistributionFunction( const std::vector<std::string>& parameters ):
fine(true)
{
}

double DistributionFunction::calculate( Value* value_in, Value* value_out ){
  plumed_massert(value_in->getNumberOfDerivatives()==value_out->getNumberOfDerivatives(), "mismatch for number of derivatives");

  // Copy derivative to out array
  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i) value_out->addDerivative( i, value_in->getDerivative(i) ); 
  // Compute the function and derivatives
  double f, df; f=compute( value_in->get(), df ); 
  value_out->chainRule(df); value_out->set(f);
  return f;
}

void DistributionFunction::finish( const double& total, Value* value_out ){
  double f, df; f=last_step( total, df );
  value_out->chainRule(df); value_out->set(f); 
}

}

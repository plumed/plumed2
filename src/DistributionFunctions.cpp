#include "DistributionFunctions.h"

namespace PLMD {

DistributionFunction::DistributionFunction( const std::vector<std::string>& parameters ):
fine(true)
{
}

void DistributionFunction::copyDerivatives( Value* value_in, Value* value_out ){
  plumed_massert(value_in->getNumberOfDerivatives()==value_out->getNumberOfDerivatives(), "mismatch for number of derivatives"); 
  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i) value_out->addDerivative( i, value_in->getDerivative(i) );
}

}

#include "DistributionFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER SUM
/**

Calculate the sum of all the colvars in the distribution.  Once calculated the final value is referenced
using lable.sum.

*/
//+ENDPLUMEDOC


sum::sum( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==0,"parameters should have zero size");
}

std::string sum::message(){
  std::ostringstream ostr;
  ostr<<"the sum of all the values"; 
  return ostr.str();
}

double sum::compute( const double p, double& df ){
  df=1.0; return p;
}

double sum::last_step( const double p, double& df ){
  df=1.0; return p;
}

}

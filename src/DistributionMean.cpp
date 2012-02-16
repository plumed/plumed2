#include "DistributionFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER AVERAGE
/**

Calculate the average of all the colvars in the distribution.  Once calculated the final value is referenced
using lable.average.

*/
//+ENDPLUMEDOC


mean::mean( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],nvalues);
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

double mean::compute( const double p, double& df ){
  df=1.0; return p;
}

double mean::last_step( const double p, double& df ){
  df=1.0/nvalues; return p/nvalues;
}

}

#include "ColvarWithoutModifiers.h"

namespace PLMD {

ColvarWithoutModifiers::ColvarWithoutModifiers(const ActionOptions& ao) :
Colvar(ao)
{
}

void ColvarWithoutModifiers::finishColvarSetup( const unsigned ncv, const double min, const double max ){
  readActionColvar();
  std::vector<double> domain(2); domain[0]=min; domain[1]=max;
  readActionWithExternalArguments( 3*getNumberOfAtoms()+9, domain ); 

  std::string ss;
  for(unsigned i=0;i<ncv;++i){
      Tools::convert(i,ss); addValue("value " + ss, false, true); 
  }
  checkRead();
}

}

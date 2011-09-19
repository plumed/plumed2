#include "Atoms.h"
#include "ActionWithExternalArguments.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR VOLUME
/**
Calculate the volume of the simulation box.

\par Syntax
\verbatim
VOLUME
\endverbatim

\par Example
The following input is printing the volume of the system
\verbatim
VOLUME LABEL=vol
PRINT ARG=vol
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC


class ColvarVolume : public ActionWithExternalArguments {
private:
  bool components;
  Tensor box;
public:
  ColvarVolume(const ActionOptions&);
// active methods:
  virtual void clearOutputForces(){};
  virtual void retrieveData();
  virtual void calculate();
  virtual void apply();
};

PLUMED_REGISTER_ACTION(ColvarVolume,"VOLUME")

ColvarVolume::ColvarVolume(const ActionOptions&ao):
ActionWithExternalArguments(ao),
components(false)
{

  registerKeyword( 0, "COMPONENTS", "use xx, yy, zz, alpha, beta, gamma as the colvars rather than the box volume");
  readAction();
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components){
// todo
  } else {
     addValueWithDerivatives("volume");
     getValue("volume")->setPeriodicity(false);
  }
}

void ColvarVolume::retrieveData(){
  box=plumed.getAtoms().box;
}

// calculator
void ColvarVolume::calculate(){
  if(components){
// todo
  } else {
    setValue( box.determinant() );
  }
}

void ColvarVolume::apply(){
  Tensor v; 
  if(components){
    // todo
  } else { 
      if (!getValue(0)->checkForced()) return;
      const double force=getValue(0)->getForce();
      v(0,0)-=force; v(1,1)-=force; v(2,2)-=force;
  }
  Tensor & vir(plumed.getAtoms().virial); vir+=v;
}

}




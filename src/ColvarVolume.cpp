#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR VOLUME
/*
Calculate the volume of the simulation box.

\par Examples
The following input tells plumed to print the volume of the system
\verbatim
VOLUME LABEL=vol
PRINT ARG=vol
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC


class ColvarVolume : public Colvar {
  bool components;

public:
  ColvarVolume(const ActionOptions&);
// active methods:
  virtual void calculate();
/// Register all the keywords for this action
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ColvarVolume,"VOLUME")

ColvarVolume::ColvarVolume(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false)
{
  std::vector<AtomNumber> atoms;
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components){
// todo
  }
  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
}

void ColvarVolume::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
  keys.addFlag("COMPONENTS",false,"use xx, yy, zz, alpha, beta, gamma as the colvars rather than the box volume");
}


// calculator
void ColvarVolume::calculate(){
  if(components){
// todo
  };

  setBoxDerivatives(-1.0*Tensor::identity());
  setValue         (getBox().determinant());
}

}




#include "Colvar.h"
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


class ColvarVolume : public Colvar {
  bool components;

public:
  ColvarVolume(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarVolume,"VOLUME")

ColvarVolume::ColvarVolume(const ActionOptions&ao):
Colvar(ao),
components(false)
{
  std::vector<AtomNumber> atoms;
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components){
// todo
  }
  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  requestAtoms(atoms);
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




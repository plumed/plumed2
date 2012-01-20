#include "Colvar.h"
#include "ActionRegister.h"
#include "PlumedMain.h"

#include <string>
#include <cmath>
#include <cassert>

namespace PLMD{

//+PLUMEDOC COLVAR ENERGY
/**
Calculate the total energy of the simulation box.

\par Syntax
\verbatim
ENERGY
\endverbatim

\par Example
The following input is printing the energy of the system
\verbatim
ENERGY LABEL=ene
PRINT ARG=ene
\endverbatim
(See also \ref PRINT).

\bug It does not work with variable cell simulations. Should include also long range corrections?

*/
//+ENDPLUMEDOC


class ColvarEnergy : public Colvar {
  bool components;

public:
  ColvarEnergy(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};


using namespace std;


PLUMED_REGISTER_ACTION(ColvarEnergy,"ENERGY")

ColvarEnergy::ColvarEnergy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false)
{
  assert(!checkNumericalDerivatives());
  std::vector<AtomNumber> atoms;
  requestAtoms(atoms);
  isEnergy=true;
  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);
}

void ColvarEnergy::registerKeywords( Keywords& keys ){
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
}


// calculator
void ColvarEnergy::calculate(){
  setValue(getEnergy());
}

}




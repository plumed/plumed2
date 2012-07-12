#include "Colvar.h"
#include "ActionRegister.h"
#include "PlumedMain.h"

#include <string>
#include <cmath>
#include <cassert>

namespace PLMD{

//+PLUMEDOC COLVAR ENERGY
/*
Calculate the total energy of the simulation box.

\par Examples
The following input instructs plumed to print the energy of the system
\verbatim
ENERGY LABEL=ene
PRINT ARG=ene
\endverbatim
(See also \ref PRINT).

\bug This command does not work with variable cell simulations. Should it also include long range corrections?

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
  isEnergy=true;
  addValueWithDerivatives(); setNotPeriodic();
  getPntrToValue()->resizeDerivatives(1);
}

void ColvarEnergy::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
}


// calculator
void ColvarEnergy::calculate(){
  setValue( getEnergy() );
  getPntrToComponent(0)->addDerivative(0,1.0);
}

}




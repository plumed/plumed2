#include "ActionRegister.h"
#include "ColvarEnergy.h"
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

using namespace std;

PLUMED_REGISTER_ACTION(ColvarEnergy,"ENERGY")

ColvarEnergy::ColvarEnergy(const ActionOptions&ao):
ActionWithExternalArguments(ao)
{
  forbidKeyword("NUMERICAL_DERIVATIVES"); forbidKeyword("STRIDE");
  readAction();
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 1, domain);
  if( checkNumericalDerivatives() ) error("numerical derivatives cannot be used with energy cvs");
  addValue("pe", true, false);
  checkRead();
}

void ColvarEnergy::retrieveData(){
  energy=plumed.getAtoms().getEnergy();
}

// calculator
void ColvarEnergy::calculate(){
  setValue( 0, energy, 1.0 );
}

void ColvarEnergy::apply(){
  std::vector<double> f(1); 
  if( getForces( 0, f ) ) plumed.getAtoms().forceOnEnergy+=f[0];
}

}




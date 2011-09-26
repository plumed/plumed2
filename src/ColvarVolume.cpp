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
/// Are we using all components of the box or just the volume
  bool components;
/// The box coordinates
  Tensor box;
public:
  ColvarVolume(const ActionOptions&);
/// Get rid of any output forces
  void clearOutputForces(){};
/// Get the cell from atoms
  void retrieveData();
/// Transfer the cell box to the value
  void calculate();
/// Apply any forces to the cell box
  void apply();
/// You can't calculate numerical derivatives
  void calculateNumericalDerivatives(){ assert(false); }
};

PLUMED_REGISTER_ACTION(ColvarVolume,"VOLUME")

ColvarVolume::ColvarVolume(const ActionOptions&ao):
ActionWithExternalArguments(ao),
components(false)
{
  forbidKeyword("NUMERICAL_DERIVATIVES"); forbidKeyword("STRIDE");
  registerKeyword( 0, "COMPONENTS", "use xx, yy, zz, alpha, beta, gamma as the colvars rather than the box volume");
  readAction();
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 1, domain );
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components){
// todo
  } else {
     addValue("volume", true, false );
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
    setValue( 0, box.determinant(), 1.0 );
  }
}

void ColvarVolume::apply(){
  Tensor v; v.clear(); std::vector<double> forces(1);
  if( components ){
    // todo
  } else if( getForces( 0, forces ) ){
       v(0,0)-=forces[0]; v(1,1)-=forces[0]; v(2,2)-=forces[0];
  }
  Tensor & vir(plumed.getAtoms().virial); vir+=v;
}

}




#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC VATOM COM
/**
Calculate the center of mass of a group of atoms

\par Example
The following input instructs plumed to claculate the distance between the
center of mass of atoms 1,2,3,4,5,6,7 and the center of mass of atoms 15,20:
\verbatim
COM ATOMS=1-7         LABEL=c1
COM ATOMS=15,20       LABEL=c2
DISTANCE ATOMS=c1,c2  LABEL=d1
PRINT ARG=d1
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC


class GenericCOM:
  public ActionWithVirtualAtom
{
public:
  GenericCOM(const ActionOptions&ao);
  void updateDynamicContent( const double& cutoff, std::vector<bool>& skips ){ assert(false); }
  void calculate();
};

PLUMED_REGISTER_ACTION(GenericCOM,"COM")

GenericCOM::GenericCOM(const ActionOptions&ao):
  ActionWithVirtualAtom(ao)
{
  allowKeyword("ATOMS"); allowKeyword("GROUP");
  forbidKeyword("UPDATE"); forbidKeyword("NL_CUTOFF");
  readActionWithVirtualAtom();
  checkRead();
}

void GenericCOM::calculate(){
  Vector pos;
  double mass(0.0),charge(0.0);
  unsigned natoms=getNumberOfAtoms();
  vector<Tensor> deriv(natoms);
  for(unsigned i=0;i<natoms;i++) mass+=getMasses(i);
  for(unsigned i=0;i<natoms;i++) charge+=getCharges(i);
  for(unsigned i=0;i<natoms;i++){
    pos+=(getMasses(i)/mass)*getPositions(i);
    deriv[i]=(getMasses(i)/mass)*Tensor::identity();
  }
  setPosition(pos);
  setMass(mass);
  setCharge(charge);
  setAtomsDerivatives(deriv);
}

}

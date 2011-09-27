#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC COM
/**
Calculate the center of mass of a group of atoms

\par Syntax
\verbatim
COM LABEL=label ATOMS=x,y,z,...
\endverbatim
The center of mass of atoms x,y,z,... is computed and stored in a virtual
atom which can be accessed through the label "label"

\par Example
The following input is printing the distance between the
center of mass of atoms 1,2,3,4,5,6,7 and that of atoms 15,20:
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
  void updateNeighbourList( const double& cutoff, std::vector<bool>& skips ){ assert(false); }
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

#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.  The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the the label for the COM action that creates it.

\par Examples
The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
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
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(GenericCOM,"COM")

void GenericCOM::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
}

GenericCOM::GenericCOM(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
  requestAtoms(atoms);
}

void GenericCOM::calculate(){
  Vector pos;
  double mass(0.0),charge(0.0);
  vector<Tensor> deriv(getNumberOfAtoms());
  for(unsigned i=0;i<getNumberOfAtoms();i++) mass+=getMass(i);
  for(unsigned i=0;i<getNumberOfAtoms();i++) charge+=getCharge(i);
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    pos+=(getMass(i)/mass)*getPosition(i);
    deriv[i]=(getMass(i)/mass)*Tensor::identity();
  }
  setPosition(pos);
  setMass(mass);
  setCharge(charge);
  setAtomsDerivatives(deriv);
}

}

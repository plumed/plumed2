#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR TEMPLATE
/**
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

*/
//+ENDPLUMEDOC
   
class ColvarTemplate : public Colvar {
  bool pbc;

public:
  ColvarTemplate(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarTemplate,"TEMPLATE")

void ColvarTemplate::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.addFlag("TEMPLATE_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
  keys.addFlag("TEMPLATE_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  keys.add("compulsory","TEMPLATE_COMPULSORY","all compulsory keywords should be added like this with a description here");
  keys.add("optional","TEMPLATE_OPTIONAL","all optional keywords that have input should be added like a description here");
  keys.add("atoms","TEMPLATE_INPUT","the keyword with which you specify what atoms to use should be added like this");
}

ColvarTemplate::ColvarTemplate(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  assert(atoms.size()==2);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");

  requestAtoms(atoms);
}


// calculator
void ColvarTemplate::calculate(){

  Vector distance;
  if(pbc){
    distance=pbcDistance(getPositions(0),getPositions(1));
  } else {
    distance=delta(getPositions(0),getPositions(1));
  }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  setAtomsDerivatives(0,-invvalue*distance);
  setAtomsDerivatives(1,invvalue*distance);
  setBoxDerivatives  (-invvalue*Tensor(distance,distance));
  setValue           (value);
}

}




#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DIPOLE 
/**
This is just a template variable

*/
//+ENDPLUMEDOC
   
class ColvarDipole : public Colvar {
  vector<AtomNumber> ga_lista;
public:
  ColvarDipole(const ActionOptions&);
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarDipole,"DIPOLE")

ColvarDipole::ColvarDipole(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  parseAtomList("GROUP",ga_lista);
  checkRead();
  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  log.printf("  of %d atoms\n",ga_lista.size());
  for(unsigned int i=0;i<ga_lista.size();++i){
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  requestAtoms(ga_lista);
}

// calculator
void ColvarDipole::calculate()
{
 double dipole=0.;
 Tensor virial;
 vector<Vector> deriv;
 Vector dipje;

 deriv.resize(getPositions().size());
 for(unsigned int i=0;i<ga_lista.size();i++) {
   dipje += (getCharges(ga_lista[i].index()))*getPositions(i);
 }
 dipole = dipje.modulo();

 for(unsigned int i=0;i<ga_lista.size();i++) {
   double dfunc=getCharges(ga_lista[i].index())/dipole;
   deriv[i] = deriv[i] + (dfunc)*dipje;
   virial=virial-Tensor(getPositions(i),deriv[i]);
 }

 for(unsigned i=0;i<getPositions().size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (dipole);
 setBoxDerivatives  (virial);
}

}

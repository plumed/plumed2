#include "Colvar.h"
#include "PlumedMain.h"
#include "ActionRegister.h"

#include "PDB.h"

#include <RMSD.h>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/**
Calculate the RMSD with respect to a reference structure

\par Syntax
\verbatim
RMSD REFERENCE=file.pdb
\endverbatim

...

*/
//+ENDPLUMEDOC
   
class ColvarRMSD : public Colvar {
  RMSD rmsd;
  vector<Vector> derivs;

public:
  ColvarRMSD(const ActionOptions&);
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  string reference;
  parse("REFERENCE",reference);
  checkRead();


  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  PDB pdb;
  pdb.read(reference,0.1/plumed.getAtoms().getInternalLengthUnits());

  rmsd.setFromPDB(pdb);

  requestAtoms(pdb.getAtomNumbers());

  derivs.resize(getNatoms());
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNatoms());
}


// calculator
void ColvarRMSD::calculate(){

  double r=rmsd.calculate(getPositions(),derivs);

  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setAtomsDerivatives(i,derivs[i]);
  Tensor virial;
  for(unsigned i=0;i<derivs.size();i++) virial=virial+(-1.0*Tensor(getPositions(i),derivs[i]));
  setBoxDerivatives(virial);
}

}




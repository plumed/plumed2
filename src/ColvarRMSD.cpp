#include "Colvar.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/**
Calculate the RMSD with respect to a reference structure

\par Syntax
for an optimal alignment following Kearsley algorithm then 
\verbatim
RMSD REFERENCE=file.pdb TYPE=OPTIMAL
\endverbatim
else, for simple case (no optimal alignment but only translation) 
\verbatim
RMSD REFERENCE=file.pdb TYPE=SIMPLE
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
PLUMED_COLVAR_INIT(ao),rmsd(log)
{
  string reference;
  parse("REFERENCE",reference);
  string type;	
  type.assign("SIMPLE");
  parse("TYPE",type);

  checkRead();


  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  PDB pdb;

  // read everything in ang and transform to nm

  pdb.read(reference,0.1/plumed.getAtoms().getUnits().length);

  rmsd.setFromPDB(pdb,type);

  requestAtoms(pdb.getAtomNumbers());

  derivs.resize(getNatoms());
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNatoms());
  log.printf("  method for alignment : %s \n",rmsd.getMethod().c_str() );

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




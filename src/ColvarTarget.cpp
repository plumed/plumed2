#include "Function.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "TargetDist.h"
#include "Atoms.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC FUNCTION TARGET
/**
This function measures the pythagorean distance from a particular structure measured in the space defined by some 
set of collective variables.

\par Examples


*/
//+ENDPLUMEDOC

class TargetFrame : public Function {
private:
  TargetDist target;
  std::vector<double> derivs;
public:
  TargetFrame(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys );
};

PLUMED_REGISTER_ACTION(TargetFrame,"TARGET")

void TargetFrame::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure. " + PDB::documentation() );
  keys.add("optional","REFERENCE_VEC","the vector of values for the CVs at the reference point (if you use this you don't need REFERENCE)");
}

TargetFrame::TargetFrame(const ActionOptions&ao):
Action(ao),
Function(ao),
target(log)
{
  std::vector<double> targ;
  parseVector("REFERENCE_VEC",targ);
  if( targ.size()!=0 ){
    target.read( targ, getArguments() );
  } else {
    string reference;
    parse("REFERENCE",reference);
    PDB pdb; 
    pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().length);
    printf("Read pdb file with %d atoms inside\n",pdb.size());
    target.read( pdb, getArguments() );
  }
  checkRead();
  derivs.resize( getNumberOfArguments() );
  addValueWithDerivatives(); setNotPeriodic();
}

void TargetFrame::calculate(){
  double r=target.calculate( derivs );
  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setDerivative(i,derivs[i]);
}

}

#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include <cstdio>

using namespace PLMD;
using namespace std;

namespace PLMD
{

//+PLUMEDOC GENERIC DUMPATOMS
/**
  Dump atom positions

\par syntax
\verbatim
DUMPATOMS [STRIDE=s] FILE=file.xyz ATOMS=list
\endverbatim

Listed atoms are written on file.xyz (currently only xyz format).
The positions are those stored in plumed exactly at the point where
the directive is put in the input file. This is relevant for directives
editing atom positions (e.g. \ref WHOLEMOLECULES).
*/
//+ENDPLUMEDOC

class GenericDumpAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  FILE*fp;
public:
  GenericDumpAtoms(const ActionOptions&);
  ~GenericDumpAtoms();
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericDumpAtoms,"DUMPATOMS")

GenericDumpAtoms::GenericDumpAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao)
{
  vector<AtomNumber> atoms;
  string file;
  parse("FILE",file);
  parseAtomList("ATOMS",atoms);
  checkRead();
  assert(file.length()>0);
  fp=fopen(file.c_str(),"w");
  requestAtoms(atoms);
}

void GenericDumpAtoms::calculate(){
  fprintf(fp,"%d\n",getNatoms());
  fprintf(fp,"here we should write the box\n");
  for(unsigned i=0;i<getNatoms();++i){
    fprintf(fp,"X %f %f %f\n",getPositions(i)(0),getPositions(i)(1),getPositions(i)(2));
  }
}

GenericDumpAtoms::~GenericDumpAtoms(){
  fclose(fp);
}
  

}

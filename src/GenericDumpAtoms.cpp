#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "Pbc.h"
#include <cstdio>
#include <cassert>

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
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){};
  void update();
};

PLUMED_REGISTER_ACTION(GenericDumpAtoms,"DUMPATOMS")

void GenericDumpAtoms::registerKeywords( Keywords& keys ){
  ActionAtomistic::registerKeywords( keys );
  keys.add("input", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates");
}

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

void GenericDumpAtoms::update(){
  fprintf(fp,"%d\n",getNatoms());
  const Tensor & t(getPbc().getBox());
  if(getPbc().isOrthorombic()){
    fprintf(fp," %f %f %f\n",t(0,0),t(1,1),t(2,2));
  }else{
    fprintf(fp," %f %f %f %f %f %f %f %f %f\n",
                 t(0,0),t(0,1),t(0,2),
                 t(1,0),t(1,1),t(1,2),
                 t(2,0),t(2,1),t(2,2)
           );
  }
  for(unsigned i=0;i<getNatoms();++i){
    fprintf(fp,"X %f %f %f\n",getPositions(i)(0),getPositions(i)(1),getPositions(i)(2));
  }
}

GenericDumpAtoms::~GenericDumpAtoms(){
  fclose(fp);
}
  

}

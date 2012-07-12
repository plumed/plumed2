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
/*

This command can be used to output the positions of a particular set of atoms.
The atoms required are ouput in a xyz formatted file.  Importantly, if your
input file contains actions that edit the atoms position (e.g. \ref WHOLEMOLECULES)
and the DUMPDERIVATIVES command appears after this instruction, then the eddited
atom positions are output.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples

The following input instructs plumed to print out the positions of atoms
1-10 together with the position of the center of mass of atoms 11-20 every
10 steps to a file called file.xyz.
\verbatim
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
\endverbatim

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
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
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
  log.printf("  printing the following atoms :");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
}

void GenericDumpAtoms::update(){
  fprintf(fp,"%d\n",getNumberOfAtoms());
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
  for(unsigned i=0;i<getNumberOfAtoms();++i){
    fprintf(fp,"X %f %f %f\n",getPosition(i)(0),getPosition(i)(1),getPosition(i)(2));
  }
}

GenericDumpAtoms::~GenericDumpAtoms(){
  fclose(fp);
}
  

}

#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC DUMPFORCES
/**

Dump the force acting on a value on a file. E.g., for a CV it dumps
the force on the CV itself. Notice that to have the forces on atoms
one should multiply this times the output of DUMPDERIVATIVES.

\par Syntax
\verbatim
DUMPFORCES ARG=what [STRIDE=s] [FILE=file]
\endverbatim
It can print at the same time forces on more than one object.

*/
//+ENDPLUMEDOC

class GenericDumpForces :
public ActionPilot,
public ActionWithArguments
{
  string file;
  FILE* fp;
public:
  void calculate(){};
  GenericDumpForces(const ActionOptions&);
  void apply(){};
  void update();
  ~GenericDumpForces();
};

PLUMED_REGISTER_ACTION(GenericDumpForces,"DUMPFORCES")

GenericDumpForces::GenericDumpForces(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fp(NULL)
{
  parse("FILE",file);
  assert(file.length()>0);
  if(comm.Get_rank()==0){
    fp=fopen(file.c_str(),"wa");
    log.printf("  on file %s\n",file.c_str());
    fprintf(fp,"%s","#! FIELDS time parameter");
    const std::vector<Value*>& arguments(getArguments());
    assert(arguments.size()>0);
    for(unsigned i=0;i<arguments.size();i++){
      fprintf(fp," %s",arguments[i]->getFullName().c_str());
    };
    fprintf(fp,"%s","\n");
  }
  checkRead();
}


void GenericDumpForces::update(){
  if(comm.Get_rank()!=0)return;
  const std::vector<Value*>& arguments(getArguments());
  fprintf(fp," %f",getTime());
  for(unsigned i=0;i<getNumberOfArguments();i++){
    fprintf(fp," %15.10f",arguments[i]->getForce());
  };
  fprintf(fp,"\n");
}

GenericDumpForces::~GenericDumpForces(){
  if(fp) fclose(fp);
}

}



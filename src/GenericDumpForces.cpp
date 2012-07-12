#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC ANALYSIS DUMPFORCES
/*
Dump the force acting on one of a values in a file.  

For a CV this command will dump
the force on the CV itself. Be aware that in order to have the forces on the atoms
you should multiply the output from this argument by the output from DUMPDERIVATIVES.
Furthermore, also note that you can output the forces on multiple quantities simultaneously
by specifying more than one argument. You can control the buffering of output using the \ref FLUSH keyword.


\par Examples
The following input instructs plumed to write a file called forces that contains
the force acting on the distance between atoms 1 and 2. 
\verbatim
DISTANCE ATOM=1,2 LABEL=distance
DUMPFORCES ARG=distance STRIDE=1 FILE=forces
\endverbatim

(See also \ref DISTANCE)

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
  static void registerKeywords(Keywords& keys);
  void apply(){};
  void update();
  ~GenericDumpForces();
};

PLUMED_REGISTER_ACTION(GenericDumpForces,"DUMPFORCES")

void GenericDumpForces::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the forces should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the forces");
}

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
    if( getNumberOfArguments()==0 ) error("no arguments have been specified");
    fprintf(fp,"%s","#! FIELDS time parameter");
    for(unsigned i=0;i<getNumberOfArguments();i++){
      fprintf(fp," %s",getPntrToArgument(i)->getName().c_str());
    };
    fprintf(fp,"%s","\n");
  }
  checkRead();
}


void GenericDumpForces::update(){
  if(comm.Get_rank()!=0)return;
  fprintf(fp," %f",getTime());
  for(unsigned i=0;i<getNumberOfArguments();i++){
    fprintf(fp," %15.10f",getPntrToArgument(i)->getForce());
  };
  fprintf(fp,"\n");
}

GenericDumpForces::~GenericDumpForces(){
  if(fp) fclose(fp);
}

}



#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC DUMPDERIVATIVES
/**
Dump the derivatives of one or more actions with respect to their input parameters
on a file. For collective variables, these are the derivatives of the collective
variable (or of one of its components) with respect to atom positions and to cell
vectors (virial-like form). For functions, there are the the derivatives with respect
the function arguments.

\verbatim
DISTANCE ATOM=1,2 LABEL=distance
# this is writing on file deriv the derivatives of the
# distance wrt atomic and cell coordinates
DUMPDERIVATIVES ARG=distance STRIDE=100 FILE=deriv
\endverbatim

It can print at the same time derivatives of more than one object, but they should
have the same number of parameters. Typically, this can be used to test numerical
derivatives againts analytical ones

\verbatim
DISTANCE ATOM=1,2 LABEL=distance
DISTANCE ATOM=1,2 LABEL=distanceN NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=distance,distanceN STRIDE=1 FILE=deriv
\endverbatim


*/
//+ENDPLUMEDOC

class DumpDerivatives :
public ActionPilot,
public ActionWithArguments
{
  string file;
  FILE* fp;
public:
  void calculate();
  DumpDerivatives(const ActionOptions&);
  void apply(){};
  ~DumpDerivatives();
};

PLUMED_REGISTER_ACTION(DumpDerivatives,"DUMPDERIVATIVES")

DumpDerivatives::DumpDerivatives(const ActionOptions&ao):
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
    fprintf(fp,"#! FIELDS time parameter derivative");
    const std::vector<Value*>& arguments(getArguments());
    assert(arguments.size()>0);
    unsigned npar=arguments[0]->getDerivatives().size();
    assert(npar>0);
    for(unsigned i=1;i<arguments.size();i++){
      assert(npar==arguments[i]->getDerivatives().size());
    }
    for(unsigned i=0;i<arguments.size();i++){
      fprintf(fp," %s",arguments[i]->getFullName().c_str());
    };
    fprintf(fp,"\n");
  }
  checkRead();
}


void DumpDerivatives::calculate(){
  if(comm.Get_rank()!=0)return;
  const std::vector<Value*>& arguments(getArguments());
  unsigned npar=arguments[0]->getDerivatives().size();
  for(unsigned ipar=0;ipar<npar;ipar++){
    fprintf(fp," %f",getTime());
    fprintf(fp," %u",ipar);
    for(unsigned i=0;i<getNumberOfArguments();i++){
      fprintf(fp," %f",arguments[i]->getDerivatives()[ipar]);
    };
    fprintf(fp,"\n");
  }
}

DumpDerivatives::~DumpDerivatives(){
  if(fp) fclose(fp);
}

}



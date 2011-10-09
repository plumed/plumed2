#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC ANALYSIS DUMPDERIVATIVES
/**

Dump the derivatives of one or more objects with respect to their input parameters
on a file. For collective variables, these are the derivatives of the collective
variable (or of one of its components) with respect to atom positions and to cell
vectors (virial-like form). For functions, there are the the derivatives with respect
the function arguments.

\par Syntax
\verbatim
DUMPDERIVATIVES ARG=what [STRIDE=s] [FILE=file]
\endverbatim
It can print at the same time derivatives of more than one object, but they should
have the same number of parameters. Typically, this can be used to test numerical
derivatives againts analytical ones

\par Example
The following input is writing on file deriv distanceB both
the analytical and numerical derivatives of distance between atoms 1 and 2.
\verbatim
DISTANCE ATOM=1,2 LABEL=distance
DISTANCE ATOM=1,2 LABEL=distanceN NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=distance,distanceN STRIDE=1 FILE=deriv
\endverbatim

(See also \ref DISTANCE)

*/
//+ENDPLUMEDOC

class GenericDumpDerivatives :
public ActionWithArguments
{
  string file;
  FILE* fp;
public:
  void calculate();
  GenericDumpDerivatives(const ActionOptions&);
  void apply(){};
  ~GenericDumpDerivatives();
};

PLUMED_REGISTER_ACTION(GenericDumpDerivatives,"DUMPDERIVATIVES")

GenericDumpDerivatives::GenericDumpDerivatives(const ActionOptions&ao):
ActionWithArguments(ao),
fp(NULL)
{
  registerKeyword(1,"FILE","the name of the file on which to write out the derivative information");
  readAction();
  std::vector<double> domain(2,0.0);
  readActionWithArguments( domain );

  parse("FILE",file);
  assert(file.length()>0);
  if(comm.Get_rank()==0){
    fp=fopen(file.c_str(),"wa");
    log.printf("  on file %s\n",file.c_str());
    unsigned nargs=getNumberOfArguments();
    if( nargs==0 ) error("no arguments specified");
    unsigned npar=getNumberOfDerivatives(0);
    for(unsigned i=1;i<nargs;i++){
        if( npar!=getNumberOfDerivatives(i) ) error("the number of derivatives must be the same in all values being dumped");
    } 
    fprintf(fp,"%s","#! FIELDS time parameter"); 
    printArgumentNames(fp);
    fprintf(fp,"%s","\n");
  }
  checkRead();
}


void GenericDumpDerivatives::calculate(){
  if(comm.Get_rank()!=0)return;
  unsigned npar=getNumberOfDerivatives(0);
  for(unsigned ipar=0;ipar<npar;ipar++){
    fprintf(fp," %f",getTime());
    fprintf(fp," %u",ipar);
    for(unsigned i=0;i<getNumberOfArguments();i++){
      fprintf(fp," %15.10f", getArgumentDerivative(i,ipar) );   
    };
    fprintf(fp,"\n");
  }
}

GenericDumpDerivatives::~GenericDumpDerivatives(){
  if(fp) fclose(fp);
}

}



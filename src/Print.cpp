#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC PRINT
/**
Print quantities on file

This action is used to periodically print quantities on a file.
Similarly to Actions of type \ref Bias, it accepts keywords
ARG and STRIDE to specify which quantities should be printed and
how frequently. It also accepts a keyword FMT specifying (in printf() style)
the format for the written numbers and a keyword FILE specifying the
name of the output file (if omitted, plumed log will be used).
This directive can be used multiple times to write multiple files,
perhaps at different stride.

Example
\verbatim
DISTANCE ATOMS=0,10 LABEL=distance
ENERGY   LABEL=energy

# this is writing distance on file COLVAR every 10 steps
PRINT ARG=distance STRIDE=10  FILE=COLVAR
# this is writing distance and energy on file COLVAR_ALL every 1000 steps
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endverbatim

*/
//+ENDPLUMEDOC

class Print :
public ActionPilot,
public ActionWithArguments
{
  string file;
  FILE* fp;
  string fmt;
public:
  void calculate();
  Print(const ActionOptions&);
  void apply(){};
  ~Print();
};

PLUMED_REGISTER_ACTION(Print,"PRINT")

Print::Print(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fp(NULL),
fmt("%f")
{
  parse("FILE",file);
  if(file.length()>0){
    if(comm.Get_rank()==0){
      fp=fopen(file.c_str(),"wa");
      log.printf("  on file %s\n",file.c_str());
      fprintf(fp,"#! FIELDS time");
      const std::vector<Value*>& arguments(getArguments());
      for(unsigned i=0;i<arguments.size();i++){
        fprintf(fp," %s",arguments[i]->getFullName().c_str());
      };
      fprintf(fp,"\n");
    }
  } else {
    log.printf("  on plumed log file\n");
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  checkRead();
}


void Print::calculate(){
    if(comm.Get_rank()!=0)return;
    if(!fp){
      log.printf("PRINT:");
      for(unsigned i=0;i<getNumberOfArguments();i++){
        log.printf(fmt.c_str(),getArgument(i));
      };
      log.printf("\n");
    } else {
      fprintf(fp," %f",getTime());
      for(unsigned i=0;i<getNumberOfArguments();i++){
        fprintf(fp,fmt.c_str(),getArgument(i));
      };
      fprintf(fp,"\n");
    }
}

Print::~Print(){
  if(fp) fclose(fp);
}

}



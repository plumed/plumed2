#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC PRINT
/**
Print quantities to a file.  This directive can be used multiple times
in the input so you can print files with different strides or print different quantities 
to different files.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples
The following input instructs plumed to print the distance between atoms 3 and 5 on a file 
called COLVAR every 10 steps, and the distance and total energy on a file called COLVAR_ALL
every 1000 steps.
\verbatim
DISTANCE ATOMS=2,5 LABEL=distance
ENERGY             LABEL=energy
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endverbatim
(See also \ref DISTANCE and \ref ENERGY).

*/
//+ENDPLUMEDOC

class GenericPrint :
public ActionPilot,
public ActionWithArguments
{
  string file;
  FILE* fp;
  string fmt;
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  int rotate;
  int rotateCountdown;
  int rotateLast;
  vector<Value*> rotateArguments;
/////////////////////////////////////////
public:
  void calculate(){};
  void prepare();
  GenericPrint(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply(){};
  void update();
  ~GenericPrint();
};

PLUMED_REGISTER_ACTION(GenericPrint,"PRINT")

void GenericPrint::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.add("compulsory","FILE","the name of the file on which to output these quantities");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.add("hidden","_ROTATE","some funky thing implemented by GBussi");
}

GenericPrint::GenericPrint(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fp(NULL),
fmt("%f"),
rotate(0)
{
  parse("FILE",file);
  if(file.length()>0){
    if(comm.Get_rank()==0){
      fp=fopen(file.c_str(),"a");
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
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  parse("_ROTATE",rotate);
  if(rotate>0){
    rotateCountdown=rotate;
    rotateArguments=getArguments();
    vector<Value*> a(1,rotateArguments[0]);
    requestArguments(vector<Value*>(1,rotateArguments[0]));
    rotateLast=0;
  }
/////////////////////////////////////////
  checkRead();
}

void GenericPrint::prepare(){
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  if(rotate>0){
    rotateCountdown--;
    if(rotateCountdown==0){
      rotateCountdown=rotate;
      rotateLast++;
      rotateLast%=rotateArguments.size();
      requestArguments(vector<Value*>(1,rotateArguments[rotateLast]));
    }
  }
/////////////////////////////////////////
}

void GenericPrint::update(){
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

GenericPrint::~GenericPrint(){
  if(fp) fclose(fp);
}

}



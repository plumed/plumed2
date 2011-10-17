#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC PRINT
/**
Print quantities on a file.

\par Example
The following instructs plumed to print the distance between atoms 3 and 5 on file COLVAR
every 10 steps, and the distance and total energy on file COLVAR_ALL
every 1000 steps.
\verbatim
DISTANCE ATOMS=2,5 LABEL=distance
ENERGY             LABEL=energy
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endverbatim
(See also \ref DISTANCE and \ref ENERGY).

When using PRINT on machines that buffer the output you can use the \ref FLUSH keyword to
force all files to dump all the output that has been buffered with a particular frequency.

*/
//+ENDPLUMEDOC

class GenericPrint :
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
  void calculate();
  void prepare();
  GenericPrint(const ActionOptions&);
  void apply(){};
  ~GenericPrint();
};

PLUMED_REGISTER_ACTION(GenericPrint,"PRINT")

GenericPrint::GenericPrint(const ActionOptions&ao):
ActionWithArguments(ao),
fp(NULL),
fmt("%f"),
rotate(0)
{
  registerKeyword(1,"FILE","the name of the file on which to write out the information");
  registerKeyword(0,"FMT","the format with which to output numbers"); 
  readAction();
  std::vector<double> domain(2,0.0);
  readActionWithArguments( domain );

  parse("FILE",file);
  if(file.length()>0){
    if(comm.Get_rank()==0){
      fp=fopen(file.c_str(),"wa");
      log.printf("  on file %s\n",file.c_str());
      fprintf(fp,"#! FIELDS time");
      printArgumentNames(fp);
      //const std::vector<Value*>& arguments(getArguments());
      //for(unsigned i=0;i<arguments.size();i++){
      //  fprintf(fp," %s",arguments[i]->getFullName().c_str());
      //};
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
//  parse("_ROTATE",rotate);
//  if(rotate>0){
//    rotateCountdown=rotate;
//    rotateArguments=getArguments();
//    vector<Value*> a(1,rotateArguments[0]);
//    requestArguments(vector<Value*>(1,rotateArguments[0]));
//    rotateLast=0;
//  }
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
//      requestArguments(vector<Value*>(1,rotateArguments[rotateLast]));
    }
  }
/////////////////////////////////////////
}

void GenericPrint::calculate(){
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



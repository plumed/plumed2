#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC ANALYSIS DUMPPROJECTIONS
/**
Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).

*/
//+ENDPLUMEDOC

class GenericDumpProjections :
public ActionPilot,
public ActionWithArguments
{
  string file;
  string fmt;
  FILE* fp;
public:
  void calculate(){};
  GenericDumpProjections(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply(){};
  void update();
  bool checkNeedsGradients()const{return true;}
  ~GenericDumpProjections();
};

PLUMED_REGISTER_ACTION(GenericDumpProjections,"DUMPPROJECTIONS")

void GenericDumpProjections::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the derivatives should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the derivatives");
  keys.add("compulsory","FMT","%15.10f","the format with which the derivatives should be output");
}

GenericDumpProjections::GenericDumpProjections(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fmt("%15.10f"),
fp(NULL)
{
  parse("FILE",file);
  assert(file.length()>0);
  parse("FMT",fmt);
  fmt=" "+fmt;
  if(comm.Get_rank()==0){
    fp=fopen(file.c_str(),"wa");
    log.printf("  on file %s\n",file.c_str());
    log.printf("  with format %s\n",fmt.c_str());
    fprintf(fp,"%s","#! FIELDS time ");
    for(unsigned i=0;i<getNumberOfArguments();i++){
      for(unsigned j=0;j<getNumberOfArguments();j++){
         fprintf(fp," %s-%s",getPntrToArgument(i)->getName().c_str(),getPntrToArgument(j)->getName().c_str());
      }
    };
    fprintf(fp,"%s","\n");
  }
  checkRead();
}


void GenericDumpProjections::update(){
  if(comm.Get_rank()!=0)return;
  fprintf(fp," %f",getTime());
  for(unsigned i=0;i<getNumberOfArguments();i++){
    for(unsigned j=0;j<getNumberOfArguments();j++){
      fprintf(fp,fmt.c_str(),getProjection(i,j));
    }
  };
  fprintf(fp,"\n");
}

GenericDumpProjections::~GenericDumpProjections(){
  if(fp) fclose(fp);
}

}



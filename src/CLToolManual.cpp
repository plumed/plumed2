#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include "PlumedConfig.h"
#include "ActionRegister.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS manual
/**
manual is a tool that you can use to construct the manual page for 
a particular action
*/
//+ENDPLUMEDOC

class CLToolManual:
public CLTool
{
public:
  int main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "print out a description of the keywords for an action in html";
  }
};

PLUMED_REGISTER_CLTOOL(CLToolManual,"manual")

int CLToolManual::main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc){

// to avoid warnings:
 (void) in;

 std::string action="none";
 bool printhelp(false);

// Start parsing options
 string prefix("");
 std::string a("");
 for(int i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help"){
      printhelp=true;
      break;
    }
    if(a.find("--action=")==0){
      a.erase(0,a.find("=")+1);
      action=a; 
    } 
 } 

 if(printhelp){
   fprintf(out,"%s",
   "Usage: info [options]\n"
   "Options:\n"
   "  [--help|-h]           : prints this help\n"
   "  [--action]            : print the manual for this action\n"
);
   return 0; 
 }

 if(action=="none"){
    fprintf(stderr,"missing --action flag\n");
    return 1;
 }

 std::cerr<<actionRegister(); 
 if( !actionRegister().printManual(action) ){
    fprintf(stderr,"specified action is not registered\n");
    return 1; 
 }

 return 0;
}

} // End of namespace

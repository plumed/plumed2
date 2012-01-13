#include "PlumedException.h"
#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace PLMD;

/**
This is the entry point to the command line tools
included in the plumed library.
*/

int CLTool::globalMain(int argc, char **argv,FILE*in,FILE*out,PlumedCommunicator& pc){
  int i;
  bool printhelp=false;

// Check if PLUMED_ROOT is defined
  string root;
  {
    char* croot=getenv("PLUMED_ROOT");
    if(!croot){
      string msg=
         "ERROR: I cannot find PLUMED\n"
         "Please set PLUMED_ROOT environment variable\n";
      fprintf(stderr,"%s",msg.c_str());
      return 1;
    }
    root=croot;
  }

// Check if PLUMED_ROOT/patches/ directory exists (as a further check)
  {
    vector<string> files=Tools::ls(root);
    if(find(files.begin(),files.end(),"patches")==files.end()) {
      string msg=
         "ERROR: I cannot find $PLUMED_ROOT/patches/ directory\n"
         "Check your PLUMED_ROOT variable\n";
      fprintf(stderr,"%s",msg.c_str());
      return 1;
    }
  }

// Build list of available C++ tools:
  vector<string> availableCxx=cltoolRegister().list();
// Build list of available shell tools:
  vector<string> availableShell;
  {
    vector<string> tmp;
    tmp=Tools::ls(string(root+"/scripts"));
    for(unsigned j=0;j<tmp.size();++j){
      unsigned ff=tmp[j].find(".sh");
      if(ff==string::npos) tmp[j].erase();
      else                 tmp[j].erase(ff);
    }
    for(unsigned j=0;j<tmp.size();++j) if(tmp[j].length()>0) availableShell.push_back(tmp[j]);
  }

  
// Start parsing options
  string prefix("");
  string a("");
  for(i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="help" || a=="-h" || a=="--help"){
      printhelp=true;
      break;
    } else if(a=="--has-mpi"){
      if(PlumedCommunicator::initialized()) return 0;
      else return 1;
    } else if(a[0]=='-') {
      string msg="ERROR: Unknown option " +a;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    } else break;
  }
  if(printhelp){
    string msg=
        "Usage: plumed [options] [command] [command options]\n"
        "  plumed [command] -h : to print help for a specific command\n"
        "Options:\n"
        "  [help|-h|--help] : to print this help\n"
        "  [--has-mpi]      : fails if plumed is running with MPI\n"
        "Commands:\n";
    fprintf(out,"%s",msg.c_str());
    for(unsigned j=0;j<availableCxx.size();++j){
      CLTool *cl=cltoolRegister().create(availableCxx[j]);
      plumed_assert(cl);
      string manual=availableCxx[j]+" : "+cl->description();
      delete cl;
      fprintf(out,"  plumed %s\n", manual.c_str());
    }
    for(unsigned j=0;j<availableShell.size();++j){
      string cmd=root+"/scripts/"+availableShell[j]+".sh --description";
      FILE *fp=popen(cmd.c_str(),"r");
      string line,manual;
      while(Tools::getline(fp,line))manual+=line;
      pclose(fp);
      manual= availableShell[j]+" : "+manual;
      fprintf(out,"  plumed %s\n", manual.c_str());
    }
    return 0;
  }
  if(i==argc) {
    fprintf(out,"%s","Nothing to do. Use 'plumed help' for help\n");
    return 0;
  }

// this is the command to be executed:
  string command(argv[i]);

  if(find(availableCxx.begin(),availableCxx.end(),command)!=availableCxx.end()){
    CLTool *cl=cltoolRegister().create(command);
    plumed_assert(cl);
    int ret=cl->main(argc-i,&argv[i],in,out,pc);
    delete cl;
    return ret;
  }

  if(find(availableShell.begin(),availableShell.end(),command)!=availableShell.end()){
    plumed_massert(in==stdin,"shell tools can only work on stdin");
    plumed_massert(out==stdout,"shell tools can only work on stdin");
    string cmd=root+"/scripts/"+command+".sh";
    for(int j=i+1;j<argc;j++) cmd+=string(" ")+argv[j];
    system(cmd.c_str());
    return 0;
  }

  string msg="ERROR: unknown command " + command;
  fprintf(stderr,"%s\n",msg.c_str());
  return 1;

}

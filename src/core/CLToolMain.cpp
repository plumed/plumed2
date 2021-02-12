/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLToolMain.h"
#include "config/Config.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/DLLoader.h"
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <memory>
#include <unordered_map>

namespace PLMD {

CLToolMain::CLToolMain():
  argc(0),
  in(stdin),
  out(stdout)
{
}

CLToolMain::~CLToolMain() {
// empty destructor to delete unique_ptr
}

#define CHECK_NULL(val,word) plumed_massert(val,"NULL pointer received in cmd(\"CLTool " + word + "\")");

void CLToolMain::cmd(const std::string& word,void*val) {

// Enumerate all possible commands:
  enum {
#include "CLToolMainEnum.inc"
  };

// Static object (initialized once) containing the map of commands:
  const static std::unordered_map<std::string, int> word_map = {
#include "CLToolMainMap.inc"
  };

  std::vector<std::string> words=Tools::getWords(word);
  unsigned nw=words.size();
  if(nw==0) {
    // do nothing
  } else {
    int iword=-1;
    char**v;
    char*vv;
    const auto it=word_map.find(words[0]);
    if(it!=word_map.end()) iword=it->second;
    switch(iword) {
    case cmd_setArgc:
      CHECK_NULL(val,word);
      argc=*static_cast<int*>(val);
      break;
    case cmd_setArgv:
      CHECK_NULL(val,word);
      v=static_cast<char**>(val);
      for(int i=0; i<argc; ++i) argv.push_back(std::string(v[i]));
      break;
    case cmd_setArgvLine:
      CHECK_NULL(val,word);
      vv=static_cast<char*>(val);
      argv=Tools::getWords(vv);
      break;
    case cmd_setIn:
      CHECK_NULL(val,word);
      in=static_cast<FILE*>(val);
      break;
    case cmd_setOut:
      CHECK_NULL(val,word);
      out=static_cast<FILE*>(val);
      break;
    case cmd_setMPIComm:
      comm.Set_comm(val);
      break;
    case cmd_setMPIFComm:
      comm.Set_fcomm(val);
      break;
    case cmd_run:
      CHECK_NULL(val,word);
      argc=argv.size();
      {
        int n=0; for(int i=0; i<argc; ++i) n+=argv[i].length()+1;
        std::vector<char> args(n);
        std::vector<char*> vvv(argc);
        char* ptr=&args[0];
        for(int i=0; i<argc; ++i) {
          vvv[i]=ptr;
          for(unsigned c=0; c<argv[i].length(); ++c) {
            *ptr=argv[i][c]; ptr++;
          }
          *ptr=0; ptr++;
        }
        int ret=run(argc,&vvv[0],in,out,comm);
        *static_cast<int*>(val)=ret;
      }
      break;
    default:
      plumed_merror("cannot interpret cmd(\"CLTool " + word + "\"). check plumed developers manual to see the available commands.");
      break;
    }
  }
}

/**
This is the entry point to the command line tools
included in the plumed library.
*/

int CLToolMain::run(int argc, char **argv,FILE*in,FILE*out,Communicator& pc) {
  int i;
  bool printhelp=false;

  DLLoader dlloader;

  std::string root=config::getPlumedRoot();

  bool standalone_executable=false;

// Start parsing options
  std::string prefix("");
  std::string a("");
  for(i=1; i<argc; i++) {
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="help" || a=="-h" || a=="--help") {
      printhelp=true;
      break;
    } else if(a=="--has-mpi") {
      if(Communicator::initialized()) return 0;
      else return 1;
    } else if(a=="--has-cregex") {
      return (config::hasCregex()?0:1);
    } else if(a=="--has-dlopen") {
      return (config::hasDlopen()?0:1);
    } else if(a=="--has-molfile") {
      return (config::hasMolfile()?0:1);
    } else if(a=="--has-external-molfile") {
      return (config::hasExternalMolfile()?0:1);
    } else if(a=="--has-zlib") {
      return (config::hasZlib()?0:1);
    } else if(a=="--has-xdrfile") {
      return (config::hasXdrfile()?0:1);
    } else if(a=="--is-installed") {
      return (config::isInstalled()?0:1);
    } else if(a=="--no-mpi") {
// this is ignored, as it is parsed in main
      continue;
    } else if(a=="--mpi") {
// this is ignored, as it is parsed in main
      continue;
    } else if(a=="--standalone-executable") {
      standalone_executable=true;
    } else if(Tools::startWith(a,"--load=")) {
      a.erase(0,a.find("=")+1);
      prefix="";
      void *p=dlloader.load(a);
      if(!p) {
        fprintf(stderr,"ERROR: cannot load library %s\n",a.c_str());
        fprintf(stderr,"ERROR: %s\n",dlloader.error().c_str());
        return 1;
      }
    } else if(a=="--load") {
      prefix="--load=";
    } else if(a[0]=='-') {
      std::string msg="ERROR: Unknown option " +a;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    } else break;
  }

// Check if plumedRoot/patches/ directory exists (as a further check)
  if(!standalone_executable) {
    std::vector<std::string> files=Tools::ls(root);
    if(find(files.begin(),files.end(),"patches")==files.end()) {
      std::string msg=
        "WARNING: I cannot find "+root+"/patches/ directory. Set PLUMED_ROOT or reinstall PLUMED\n\n";
      fprintf(stderr,"%s",msg.c_str());
    }
  }

// Build list of available C++ tools:
  std::vector<std::string> availableCxx=cltoolRegister().list();
// Build list of available shell tools:
  std::vector<std::string> availableShell;
  if(!standalone_executable) {
    std::vector<std::string> tmp;
    tmp=Tools::ls(std::string(root+"/scripts"));
    for(unsigned j=0; j<tmp.size(); ++j) {
      size_t ff=tmp[j].find(".sh");
      if(ff==std::string::npos) tmp[j].erase();
      else                 tmp[j].erase(ff);
    }
    for(unsigned j=0; j<tmp.size(); ++j) if(tmp[j].length()>0) availableShell.push_back(tmp[j]);
  }

  if(printhelp) {
    std::string msg=
      "Usage: plumed [options] [command] [command options]\n"
      "  plumed [command] -h|--help: to print help for a specific command\n"
      "Options:\n"
      "  [help|-h|--help]          : to print this help\n"
      "  [--is-installed]          : fails if plumed is not installed\n"
      "  [--has-mpi]               : fails if plumed is running without MPI\n"
      "  [--has-dlopen]            : fails if plumed is compiled without dlopen\n"
      "  [--load LIB]              : loads a shared object (typically a plugin library)\n"
      "  [--standalone-executable] : tells plumed not to look for commands implemented as scripts\n"
      "Commands:\n";
    fprintf(out,"%s",msg.c_str());
    for(unsigned j=0; j<availableCxx.size(); ++j) {
      auto cl=cltoolRegister().create(CLToolOptions(availableCxx[j]));
      plumed_assert(cl);
      std::string manual=availableCxx[j]+" : "+cl->description();
      fprintf(out,"  plumed %s\n", manual.c_str());
    }
    for(unsigned j=0; j<availableShell.size(); ++j) {
      std::string cmd=config::getEnvCommand()+" \""+root+"/scripts/"+availableShell[j]+".sh\" --description";
      FILE *fp=popen(cmd.c_str(),"r");
      std::string line,manual;
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
  std::string command(argv[i]);

  if(find(availableCxx.begin(),availableCxx.end(),command)!=availableCxx.end()) {
    auto cl=cltoolRegister().create(CLToolOptions(command));
    plumed_assert(cl);
    // Read the command line options (returns false if we are just printing help)
    if( !cl->readInput( argc-i,&argv[i],in,out ) ) { return 0; }
    int ret=cl->main(in,out,pc);
    return ret;
  }

  if(find(availableShell.begin(),availableShell.end(),command)!=availableShell.end()) {
    plumed_massert(in==stdin,"shell tools can only work on stdin");
    plumed_massert(out==stdout,"shell tools can only work on stdin");
    std::string cmd=config::getEnvCommand()+" \""+root+"/scripts/"+command+".sh\"";
    for(int j=i+1; j<argc; j++) cmd+=std::string(" ")+argv[j];
    int r=system(cmd.c_str());
// this is necessary since system seems to return numbers which are multiple
// of 256. this would make the interpretation by the shell wrong
// I just return 1 in case of failure and 0 in case of success
    if(r!=0) return 1;
    else return 0;
  }

  std::string msg="ERROR: unknown command " + command + ". Use 'plumed help' for help";
  fprintf(stderr,"%s\n",msg.c_str());
  return 1;

}
}

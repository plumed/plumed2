/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS completion
/*
Dumps the body of a bash function to be used for auto completion.

Users will typically not need this command.
See more at \ref BashAutocompletion

\par Examples

\verbatim
plumed completion
\endverbatim


*/
//+ENDPLUMEDOC

class Completion:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Completion(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  string description()const override {
    return "dump a function usable for programmable completion";
  }
};

PLUMED_REGISTER_CLTOOL(Completion,"completion")

void Completion::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
}

Completion::Completion(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int Completion::main(FILE* in, FILE*out,Communicator& pc) {
  static const char completion [] = {
#include "completion.xxd"
// cppcheck-suppress syntaxError
    , 0x00
  };
  fprintf(out,"local cmds=\"help -h --help");
// Build list of available C++ tools:
  std::vector<string> availableCxx=cltoolRegister().list();
// Build list of available shell tools:
  vector<string> tmp=Tools::ls(string(config::getPlumedRoot()+"/scripts"));
  for(unsigned j=0; j<tmp.size(); ++j) {
    size_t ff=tmp[j].find(".sh");
    if(ff==string::npos) tmp[j].erase();
    else                 tmp[j].erase(ff);
  }
  for(unsigned j=0; j<availableCxx.size(); j++) fprintf(out," %s",availableCxx[j].c_str());
  for(unsigned j=0; j<tmp.size(); ++j) if(tmp[j].length()>0) fprintf(out," %s",tmp[j].c_str());
  fprintf(out,"\"\n");

  for(unsigned j=0; j<availableCxx.size(); j++) {
    std::string s=availableCxx[j];
// handle - sign (convert to underscore)
    for(;;) {
      size_t n=s.find("-");
      if(n==std::string::npos) break;
      s[n]='_';
    }
    fprintf(out,"local cmd_keys_%s=\"",s.c_str());
    std::vector<std::string> keys=cltoolRegister().getKeys(availableCxx[j]);
    for(unsigned k=0; k<keys.size(); k++) {
// handle --help/-h
      std::string s=keys[k];
      for(;;) {
        size_t n=s.find("/");
        if(n==std::string::npos) break;
        s[n]=' ';
      }
      fprintf(out," %s",s.c_str());
    }
    fprintf(out,"\"\n");
  }

////  ALTERNATIVE IMPLEMENTATION
////  checking tools on the fly
////     for(unsigned j=0; j<tmp.size(); j++) {
////       std::string s=tmp[j];
////   // handle - sign (convert to underscore)
////       for(;;) {
////         size_t n=s.find("-");
////         if(n==std::string::npos) break;
////         s[n]='_';
////       }
////       fprintf(out,"local cmd_keys_%s=\"",s.c_str());
////       std::string cmd=config::getEnvCommand()+" \""+config::getPlumedRoot()+"/scripts/"+s+".sh\" --options";
////       FILE *fp=popen(cmd.c_str(),"r");
////       std::string line,manual;
////       while(Tools::getline(fp,line))manual+=line;
////       pclose(fp);
////       std::vector<std::string> keys=Tools::getWords(manual);
////       for(unsigned k=0; k<keys.size(); k++) {
////   // handle --help/-h
////         std::string s=keys[k];
////         for(;;) {
////           size_t n=s.find("/");
////           if(n==std::string::npos) break;
////           s[n]=' ';
////         }
////         fprintf(out," %s",s.c_str());
////       }
////       fprintf(out,"\"\n");
////     }

  fprintf(out,"%s\n",completion);
  std::string name=config::getPlumedProgramName();

  fprintf(out,
          "############################################\n"
          "## ADD THESE COMMANDS TO YOUR .bashrc FILE:\n"
          "############################################\n"
          "# _%s() { eval \"$(%s --no-mpi completion 2>/dev/null)\";}\n"
          "# complete -F _%s -o default %s\n"
          "############################################\n",
          name.c_str(),name.c_str(),name.c_str(),name.c_str());

  return 0;
}

} // End of namespace
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2023 The plumed team
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
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/CLToolMain.h"
#include "core/ModuleLoader.h"
#include "tools/DLLoader.h"
#include <sstream>
#include <cstdio>
#include <string>
#include <vector>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS completion
/*
Generate a Bash completion function body or list completion candidates.

With no options this command emits the function body used by the profile snippet
described [here](module_cltools.md). The snippet queries the executable at completion
time, including explicit global module selections. `--list` emits plain candidates
for the completion function, and `--command` selects a tool for its option list.

## Examples

```plumed
plumed completion
```


*/
//+ENDPLUMEDOC

class Completion:
  public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Completion(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "dump a function usable for programmable completion";
  }
};

PLUMED_REGISTER_CLTOOL(Completion,"completion")

void Completion::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("optional","--list","list commands, global-options, options, modules or installed-modules, one per line");
  keys.add("optional","--command","tool whose options should be listed");
}

Completion::Completion(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

int Completion::main(FILE* in, FILE*out,Communicator& pc) {
  static const char completion [] = {
#include "completion.xxd"
// cppcheck-suppress syntaxError
    , 0x00
  };
  std::string list, command;
  parse("--list",list);
  parse("--command",command);
  if(!list.empty()) {
    std::vector<std::string> candidates;
    if(list=="commands") {
      candidates=cltoolRegister().list();
      candidates.push_back("help");
      candidates.push_back("modules");
      for(auto script:Tools::ls(config::getPlumedRoot()+"/scripts")) {
        if(script.size()>3 && script.compare(script.size()-3,3,".sh")==0) {
          candidates.push_back(script.substr(0,script.size()-3));
        }
      }
    } else if(list=="global-options") {
      candidates=CLToolMain::getGlobalOptions();
    } else if(list=="options") {
      if(command=="modules") {
        candidates= {"--module", "--list", "--has", "--has-feature", "-q", "--help"};
      } else {
        candidates=cltoolRegister().getKeys(command);
        static const char scripts [] = {
#include "completion-scripts.h"
          , 0x00
        };
        std::istringstream rows(scripts);
        std::string row;
        while(std::getline(rows,row)) {
          const auto tab=row.find('\t');
          if(tab!=std::string::npos && row.substr(0,tab)==command) {
            std::istringstream options(row.substr(tab+1));
            std::string option;
            while(options>>option) {
              candidates.push_back(option);
            }
          }
        }
      }
    } else if(list=="modules") {
      DLLoader loader;
      ModuleLoader::inheritContext(loader);
      candidates=ModuleLoader::modules(loader);
    } else if(list=="installed-modules") {
      candidates=ModuleLoader::installedModules();
    } else {
      std::fprintf(stderr,"Unknown completion list: %s\n",list.c_str());
      return 1;
    }
    for(auto candidate:candidates) {
      if(list!="options") {
        std::fprintf(out,"%s\n",candidate.c_str());
        continue;
      }
      std::replace(candidate.begin(),candidate.end(),'/',' ');
      std::istringstream aliases(candidate);
      std::string alias;
      while(aliases>>alias) {
        std::fprintf(out,"%s\n",alias.c_str());
      }
    }
    return 0;
  }

  std::fprintf(out,"%s\n",completion);
  std::string plumedName=config::getPlumedProgramName();

  std::fprintf(out,
               "############################################\n"
               "## ADD THESE COMMANDS TO YOUR .bashrc FILE:\n"
               "############################################\n"
               "# _%s() { eval \"$(%s --no-mpi completion 2>/dev/null)\";}\n"
               "# complete -F _%s -o default -o filenames %s\n"
               "############################################\n",
               plumedName.c_str(),
               plumedName.c_str(),
               plumedName.c_str(),
               plumedName.c_str());

  return 0;
}

} // End of namespace
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2022-2026 The plumed team
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
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "core/CLToolRegister.h"
#include "core/InputBuilderMetadata.h"
#include "tools/IFile.h"

#include <cstdio>
#include <map>
#include <sstream>
#include <string>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS gen_json
/*
gen_json constructs a json file that includes a dictionary of actions, the keywords for those actions and the components and outputs this to standard output

This command is used during the build process of plumed. It constructs a json file that contains documentation information on all actions and cltools. This json
file is then used by the python package [PlumedToHTML](https://github.com/plumed/PlumedToHTML). This package is used to construct the example input files that
appear in [the plumed nest](https://www.plumed-nest.org), [plumed tutorials](https://www.plumed-tutorials.org) and this manual.

You will likely not ever need to run this command.  When it runs during the build process we use the following command:

```plumed
plumed --no-mpi gen_json > syntax.json
```

The optional `--actions` file can override the embedded action descriptions. Each line has the form `ACTION:description`.

*/
//+ENDPLUMEDOC

class GenJson : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit GenJson(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "print out a json file that contains the pluemd syntax";
  }
};

PLUMED_REGISTER_CLTOOL(GenJson,"gen_json")

void GenJson::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("optional","--actions","a file containing one line descriptions of the various actions");
}

GenJson::GenJson(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

int GenJson::main(FILE*,FILE*out,Communicator&) {
  std::map<std::string,std::string> actionDescriptions;
  std::string actionFile;
  const bool hasActionFile=parse("--actions",actionFile);
  if(hasActionFile) {
    IFile input;
    input.open(actionFile);
    std::string line;
    while(input.getline(line)) {
      const std::size_t colon=line.find(':');
      actionDescriptions.insert({line.substr(0,colon),line.substr(colon+1)});
    }
    input.close();
  }

  std::ostringstream json;
  writeInputBuilderJSON(json,hasActionFile ? &actionDescriptions : nullptr,getDLHandles());
  const std::string result=json.str();
  if(std::fwrite(result.data(),1,result.size(),out)!=result.size()) {
    error("failed to write JSON metadata to the command output file");
  }
  return 0;
}

} // namespace cltools
} // namespace PLMD

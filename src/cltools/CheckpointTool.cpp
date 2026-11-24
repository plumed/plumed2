/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#include "tools/Checkpoint.h"
#include "tools/Tools.h"
#include <cstdio>
#include <string>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS checkpoint
/*
Inspect a PLUMED checkpoint file.

The checkpoint files written by the [CHECKPOINT](CHECKPOINT.md) directive are
binary. This tool prints their content in a readable form: the header (step,
PLUMED version, replica), the table of the output files with the sizes
recorded at the checkpoint, and the state stored by every action. Arrays are
printed with their shape; with `--full` all their elements are printed.

```plumed
plumed checkpoint --dump plumed.cpt
```

*/
//+ENDPLUMEDOC

class CheckpointTool:
  public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit CheckpointTool(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "inspect a PLUMED checkpoint file";
  }
};

PLUMED_REGISTER_CLTOOL(CheckpointTool,"checkpoint")

void CheckpointTool::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--dump","the checkpoint file to print");
  keys.addFlag("--full",false,"print all the elements of the arrays, not only the first ones");
}

CheckpointTool::CheckpointTool(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

int CheckpointTool::main(FILE* in, FILE*out,Communicator& pc) {
  std::string file;
  parse("--dump",file);
  bool full;
  parseFlag("--full",full);
  const unsigned maxElements=full ? 0 : 8;

  CheckpointReader reader(CheckpointReader::readFile(file));
  std::fprintf(out,"file: %s\n",file.c_str());
  std::fprintf(out,"format version: %u\n",reader.getFormatVersion());
  std::fprintf(out,"written by PLUMED %s\n",reader.getPlumedVersion().c_str());
  for(const auto & entryName : reader.getEntryNames()) {
    const auto & e=reader.getEntry(entryName);
    std::string shape;
    if(!e.scalar()) {
      shape="[";
      for(unsigned i=0; i<e.dims.size(); ++i) {
        shape+=(i>0?"x":"")+std::to_string(e.dims[i]);
      }
      shape+="]";
    }
    std::fprintf(out,"%s%s = ",entryName.c_str(),shape.c_str());
    switch(e.type) {
    case CheckpointReader::Type::String:
      std::fprintf(out,"\"%s\"",e.str.c_str());
      break;
    case CheckpointReader::Type::Int:
      for(unsigned i=0; i<e.ints.size(); ++i) {
        if(maxElements>0 && i>=maxElements) {
          std::fprintf(out," ...");
          break;
        }
        std::fprintf(out,"%s%lld",(i>0?" ":""),e.ints[i]);
      }
      break;
    case CheckpointReader::Type::Double:
      for(unsigned i=0; i<e.doubles.size(); ++i) {
        if(maxElements>0 && i>=maxElements) {
          std::fprintf(out," ...");
          break;
        }
        std::fprintf(out,"%s%.17g",(i>0?" ":""),e.doubles[i]);
      }
      break;
    }
    std::fprintf(out,"\n");
  }
  return 0;
}

}
}

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
#include "tools/H5Table.h"
#include "tools/IFile.h"
#include "tools/Tools.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS h5convert
/*
Convert a PLUMED file between text and HDF5.

The files that PLUMED writes field by field (e.g. with [PRINT](PRINT.md), the hills of
[METAD](METAD.md), the kernels of [OPES_METAD](OPES_METAD.md)) are written in HDF5 format
when their name ends with `.h5` (see [Files](Files.md)). This tool converts such a file
to the equivalent text file, with the numbers written with all their digits, and a text
file to HDF5, so that the tools that read text files can be used on HDF5 files and vice versa.
The direction of the conversion is given by the extensions of the two files.

```plumed
plumed h5convert --input COLVAR.h5 --output COLVAR
plumed h5convert --input HILLS --output HILLS.h5
```

When a text file is converted, the values that are numbers are stored as double precision
numbers (with the digits that are in the text file) and the others as strings.

*/
//+ENDPLUMEDOC

class H5Convert:
  public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit H5Convert(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "convert a PLUMED file between text and HDF5 (.h5)";
  }
};

PLUMED_REGISTER_CLTOOL(H5Convert,"h5convert")

void H5Convert::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--input","the file to convert");
  keys.add("compulsory","--output","the converted file (HDF5 if its extension is .h5, text otherwise)");
}

H5Convert::H5Convert(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

namespace {
/// a value read from a text file: a number if it is one (without evaluating expressions such as pi)
h5::Value textValue(const std::string& s) {
  h5::Value v;
  const char* c=s.c_str();
  char* end=nullptr;
  const double d=std::strtod(c,&end);
  if(!s.empty() && end && *end=='\0') {
    v.type=h5::Value::Double;
    v.d=d;
  } else {
    v.type=h5::Value::String;
    v.s=s;
  }
  return v;
}
}

int H5Convert::main(FILE* in, FILE*out,Communicator& pc) {
  std::string input,output;
  parse("--input",input);
  parse("--output",output);
  const bool fromH5=h5::isH5(input);
  const bool toH5=h5::isH5(output);
  if(fromH5==toH5) {
    error("one of the two files (and only one) must be an HDF5 file (extension .h5)");
  }
  // the lines of the input file, as text (an HDF5 table gives the lines of the equivalent text file)
  IFile ifile;
  ifile.open(input);
  ifile.allowNoEOL();
  std::string line;
  if(!toH5) {
    std::ofstream ofs(output);
    if(!ofs) {
      error("cannot open "+output);
    }
    while(ifile.getline(line)) {
      ofs<<line<<"\n";
    }
    return 0;
  }
  h5::Writer writer(output,false);
  std::vector<std::string> names;
  std::vector<h5::Field> constants;
  std::vector<h5::Field> row;
  while(ifile.getline(line)) {
    const std::vector<std::string> words=Tools::getWords(line);
    if(words.size()>=2 && words[0]=="#!" && words[1]=="FIELDS") {
      // a new header: its SET lines follow
      names.assign(words.begin()+2,words.end());
      constants.clear();
    } else if(words.size()>=4 && words[0]=="#!" && words[1]=="SET") {
      // the value as it is written (with its spaces)
      const std::size_t p=line.find(words[2],line.find("SET")+3)+words[2].size()+1;
      constants.push_back(h5::Field{words[2],textValue(line.substr(p))});
    } else if(!words.empty() && words[0][0]!='#') {
      if(words.size()!=names.size()) {
        error("file "+input+": the number of values does not match the number of fields in line:\n"+line);
      }
      row.resize(names.size());
      for(unsigned i=0; i<names.size(); i++) {
        row[i].name=names[i];
        row[i].value=textValue(words[i]);
      }
      writer.addRow(row,constants);
    }
  }
  writer.flush();
  return 0;
}

}
}

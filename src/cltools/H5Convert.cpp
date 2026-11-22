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
#include "tools/H5Archive.h"
#include "tools/FileTools.h"
#include "tools/IFile.h"
#include "tools/Tools.h"
#include <cstdio>
#include <filesystem>
#include <memory>
#include <limits>
#include <cstdlib>
#include <string>
#include <vector>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS h5convert
/*
Convert a PLUMED file between text, HDF5 and HDF5 batches in TAR.

The files that PLUMED writes field by field (e.g. with [PRINT](PRINT.md), the hills of
[METAD](METAD.md), the kernels of [OPES_METAD](OPES_METAD.md)) are written in HDF5 format
when their name ends with `.h5` (see [Files](Files.md)). This tool converts such a file
to the equivalent text file, with the numbers written with all their digits, and a text
file to HDF5, so that the tools that read text files can be used on HDF5 files and vice versa.
Input tables are detected by content; the output extension selects text, `.h5`
or `.h5.tar`. HDF5-to-HDF5 and TAR rebatching preserve exact stored types.

```plumed
plumed h5convert --input COLVAR.h5 --output COLVAR
plumed h5convert --input HILLS --output HILLS.h5
plumed h5convert --input HILLS --output HILLS.h5.tar --batch-size 10000
plumed h5convert --input HILLS.h5.tar --output HILLS.h5
```

Use `-` as the text filename to read from standard input or write to standard output:

```plumed
plumed h5convert --input - --output HILLS.h5 < HILLS
plumed h5convert --input COLVAR.h5 --output - > COLVAR
```

Standard input and output support text only; the output must be binary when reading text from stdin.
`--batch-size` is a positive maximum data-record count for `.h5.tar` output only;
schema changes and the memory limit can produce shorter batches.
For quick inspection and concatenation of HDF5 tables, see [h5cat](h5cat.md).

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
  Completion complete(const std::vector<std::string>& arguments) const override {
    auto result=CLTool::complete(arguments);
    if(result.option=="--input" || result.option=="--output") {
      result.candidates.push_back("-");
    }
    if(result.option=="--batch-size") {
      result.files=false;
      result.candidates= {"4096"};
    }
    return result;
  }
  std::string description()const override {
    return "convert PLUMED text, HDF5 and TAR batch tables";
  }
};

PLUMED_REGISTER_CLTOOL(H5Convert,"h5convert")

void H5Convert::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--input","the file to convert (- for text standard input)");
  keys.add("optional","--batch-size","maximum data records per TAR batch (positive integer; default 4096)");
  keys.add("compulsory","--output","the converted file (HDF5 for .h5, TAR batches for .h5.tar, text otherwise; - for text standard output)");
}

H5Convert::H5Convert(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

namespace {
/// Write the full-precision text representation, without taking ownership of out.
void writeH5Text(h5table::Reader& reader,const std::string& input,FILE* out) {
  std::string line;
  while(reader.getline(line)) {
    plumed_massert(std::fprintf(out,"%s\n",line.c_str())>=0,"cannot write converted text for "+input);
  }
  plumed_massert(std::fflush(out)==0,"cannot flush converted text for "+input);
}

/// a value read from a text file: a number if it is one (without evaluating expressions such as pi)
h5table::Value textValue(const std::string& s) {
  h5table::Value v;
  const char* c=s.c_str();
  char* end=nullptr;
  const double d=std::strtod(c,&end);
  if(!s.empty() && end && *end=='\0') {
    v.type=h5table::Value::Double;
    v.d=d;
  } else {
    v.type=h5table::Value::String;
    v.s=s;
  }
  return v;
}
}

int H5Convert::main(FILE* in, FILE*out,Communicator& pc) {
  std::string input,output;
  parse("--input",input);
  parse("--output",output);
  const bool fromH5=input!="-" && h5table::isTable(input);
  const bool toH5=h5table::isH5(output);
  std::string batchOption;
  parse("--batch-size",batchOption);
  std::size_t batchSize=4096;
  if(!batchOption.empty()) {
    unsigned long long count=0;
    if(!h5table::isArchive(output) || batchOption[0]=='-' ||
        !Tools::convertNoexcept(batchOption,count) || count==0 || count>std::numeric_limits<std::size_t>::max()) {
      error("--batch-size requires a positive integer and .h5.tar output");
    }
    batchSize=static_cast<std::size_t>(count);
  }
  if(!fromH5 && !toH5) {
    error("at least one file must be an HDF5 table or .h5.tar archive");
  }
  if(input!="-" && output!="-") {
    std::error_code ec;
    if(std::filesystem::equivalent(pathFromUTF8(input),pathFromUTF8(output),ec) || input==output) {
      error("input and output must be different files");
    }
  }
  if(!toH5) {
    h5table::Reader reader(input);
    unique_FILE outputFile;
    std::unique_ptr<h5table::BatchFile> staging;
    if(output!="-") {
      const auto parent=pathFromUTF8(output).parent_path();
      staging=std::make_unique<h5table::BatchFile>(parent.empty()?".":pathToUTF8(parent));
      outputFile.reset(fopenUTF8(staging->path,"wb"));
      if(!outputFile) {
        error("cannot open "+output);
      }
    }
    writeH5Text(reader,input,outputFile ? outputFile.get() : out);
    if(outputFile && std::fclose(outputFile.release())!=0) {
      error("cannot close "+output);
    }
    if(staging) {
      staging->publish(output);
    }
    return 0;
  }
  const auto parent=pathFromUTF8(output).parent_path();
  h5table::BatchFile staging(parent.empty()?".":pathToUTF8(parent));
  if(h5table::isArchive(output)) {
    staging.path+=".h5.tar";
  }
  auto writer=std::make_unique<h5table::Writer>(staging.path,false,batchSize);
  if(fromH5) {
    h5table::Reader reader(input);
    std::vector<h5table::Field> fields,constants;
    while(reader.readRow(fields,constants)) {
      writer->addRow(fields,constants);
    }
    writer->flush();
    writer.reset();
    staging.publish(output);
    return 0;
  }
  IFile ifile;
  if(input=="-") {
    ifile.link(in);
  } else {
    ifile.open(input);
  }
  ifile.allowNoEOL();
  std::string line;

  std::vector<std::string> names;
  std::vector<h5table::Field> constants;
  std::vector<h5table::Field> row;
  while(ifile.getline(line)) {
    const std::vector<std::string> words=Tools::getWords(line);
    if(words.size()>=2 && words[0]=="#!" && words[1]=="FIELDS") {
      // a new header: its SET lines follow
      names.assign(words.begin()+2,words.end());
      constants.clear();
    } else if(words.size()>=4 && words[0]=="#!" && words[1]=="SET") {
      // the value as it is written (with its spaces)
      const std::size_t p=line.find(words[2],line.find("SET")+3)+words[2].size()+1;
      constants.push_back(h5table::Field{words[2],textValue(line.substr(p))});
    } else if(!words.empty() && words[0][0]!='#') {
      if(words.size()!=names.size()) {
        error("file "+input+": the number of values does not match the number of fields in line:\n"+line);
      }
      row.resize(names.size());
      for(unsigned i=0; i<names.size(); i++) {
        row[i].name=names[i];
        row[i].value=textValue(words[i]);
      }
      writer->addRow(row,constants);
    }
  }
  if(std::ferror(ifile.getFILE())) {
    error("cannot read "+input);
  }
  writer->flush();
  writer.reset();
  staging.publish(output);
  return 0;
}

//+PLUMEDOC TOOLS h5cat
/*
Write PLUMED HDF5 tables and `.h5.tar` batch archives to standard output as text.

Supply one or more filenames as positional arguments. Each table is converted to
text with full numerical precision, using the same representation as
[h5convert](h5convert.md). Files are concatenated in argument order, including
each file's field headers and constants. Rows are neither sorted nor deduplicated.

```plumed
plumed h5cat COLVAR.h5
plumed h5cat COLVAR.h5.tar
plumed h5cat COLVAR.0.h5 COLVAR.1.h5 > combined.txt
```

The output can also be piped to tools such as `less` or `awk`. Use `--` before a
filename starting with `-`. Input must be named HDF5 tables written by PLUMED;
standard input and arbitrary HDF5 containers are not supported. Filenames need
not end in `.h5`. An unreadable or invalid file stops the command with a nonzero
exit status; text from preceding files may already have been written.

*/
//+ENDPLUMEDOC

class H5Cat:
  public CLTool {
  std::vector<std::string> inputs;
public:
  static void registerKeywords(Keywords& keys) {
    CLTool::registerKeywords(keys);
    keys.addFlag("--",false,"end option processing; remaining arguments are filenames");
  }
  explicit H5Cat(const CLToolOptions& co): CLTool(co) {
    inputdata=inputType::commandline;
  }
  bool readInput(int argc,char** argv,FILE* in,FILE* out) override {
    inputs.clear();
    bool options=true;
    for(int i=1; i<argc; ++i) {
      const std::string arg=argv[i];
      if(options && (arg=="--help" || arg=="-h")) {
        std::fprintf(out,"Usage: plumed h5cat [--] FILE...\n\n"
                     "Write one or more PLUMED HDF5 tables as text to standard output,\n"
                     "in argument order, preserving headers, constants and all rows.\n"
                     "Use -- before filenames starting with -. Standard input is not supported.\n");
        return false;
      }
      if(options && arg=="--") {
        options=false;
      } else if(options && !arg.empty() && arg[0]=='-') {
        error("unknown option "+arg+"; use -- before filenames starting with -");
      } else {
        inputs.push_back(arg);
      }
    }
    if(inputs.empty()) {
      error("at least one HDF5 filename is required; use h5cat --help for help");
    }
    return true;
  }
  int main(FILE* in,FILE* out,Communicator& pc) override {
    for(const auto& input : inputs) {
      h5table::Reader reader(input);
      writeH5Text(reader,input,out);
    }
    return 0;
  }
  std::string description()const override {
    return "concatenate PLUMED HDF5 tables as text";
  }
};

PLUMED_REGISTER_CLTOOL(H5Cat,"h5cat")

}
}

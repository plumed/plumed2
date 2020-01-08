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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS info
/*
This tool allows you to obtain information about your plumed version

You can specify the information you require using the following command line
arguments

\par Examples

The following command returns the root directory for your plumed distribution.
\verbatim
plumed info --root
\endverbatim

*/
//+ENDPLUMEDOC

class Info:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Info(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  string description()const override {
    return "provide informations about plumed";
  }
};

PLUMED_REGISTER_CLTOOL(Info,"info")

void Info::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.addFlag("--configuration",false,"prints the configuration file");
  keys.addFlag("--root",false,"print the location of the root directory for the plumed source");
  keys.addFlag("--user-doc",false,"print the location of user manual (html)");
  keys.addFlag("--developer-doc",false,"print the location of user manual (html)");
  keys.addFlag("--version",false,"print the version number");
  keys.addFlag("--long-version",false,"print the version number (long version)");
  keys.addFlag("--git-version",false,"print the version number (git version, if available)");
  keys.addFlag("--include-dir",false,"print the location of the include dir");
  keys.addFlag("--soext",false,"print the extension of shared libraries (so or dylib)");
}

Info::Info(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int Info::main(FILE* in, FILE*out,Communicator& pc) {

  bool printconfiguration; parseFlag("--configuration",printconfiguration);
  bool printroot; parseFlag("--root",printroot);
  bool printuserdoc; parseFlag("--user-doc",printuserdoc);
  bool printdeveloperdoc; parseFlag("--developer-doc",printdeveloperdoc);
  bool printversion; parseFlag("--version",printversion);
  bool printlongversion; parseFlag("--long-version",printlongversion);
  bool printgitversion; parseFlag("--git-version",printgitversion);
  bool printincludedir; parseFlag("--include-dir",printincludedir);
  bool printsoext; parseFlag("--soext",printsoext);
  if(printroot) fprintf(out,"%s\n",config::getPlumedRoot().c_str());
  if(printconfiguration) fprintf(out,"%s",config::getMakefile().c_str());
  if(printincludedir) fprintf(out,"%s\n",config::getPlumedIncludedir().c_str());
  if(printuserdoc) {
    std::string userdoc=config::getPlumedHtmldir()+"/user-doc/html/index.html";
    FILE *ff=std::fopen(userdoc.c_str(),"r");
    if(ff) std::fclose(ff);
    else userdoc="http://www.plumed.org/doc-v" + config::getVersion() + "/user-doc/html/index.html";
    fprintf(out,"%s\n",userdoc.c_str());
  }
  if(printdeveloperdoc) {
    std::string developerdoc=config::getPlumedHtmldir()+"/developer-doc/html/index.html";
    FILE *ff=std::fopen(developerdoc.c_str(),"r");
    if(ff) std::fclose(ff);
    else developerdoc="http://www.plumed.org/doc-v" + config::getVersion() + "/developer-doc/html/index.html";
    fprintf(out,"%s\n",developerdoc.c_str());
  }
  if(printversion) fprintf(out,"%s\n",config::getVersion().c_str());
  if(printlongversion) fprintf(out,"%s\n",config::getVersionLong().c_str());
  if(printgitversion) fprintf(out,"%s\n",config::getVersionGit().c_str());
  if(printsoext) fprintf(out,"%s\n",config::getSoExt().c_str());

  return 0;
}



}
}

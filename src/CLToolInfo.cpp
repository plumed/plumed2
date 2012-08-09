/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "Tools.h"
#include "PlumedConfig.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

namespace PLMD {

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

class CLToolInfo:
public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  CLToolInfo(const CLToolOptions& co );
  int main(FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "provide informations about plumed";
  }
};

PLUMED_REGISTER_CLTOOL(CLToolInfo,"info")

void CLToolInfo::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.addFlag("--configuration",false,"prints the configuration file");
  keys.addFlag("--root",false,"print the location of the root directory for the plumed source");
}

CLToolInfo::CLToolInfo(const CLToolOptions& co ):
CLTool(co)
{
  inputdata=commandline;
}

int CLToolInfo::main(FILE*out,PlumedCommunicator& pc){

 bool printconfiguration; parseFlag("--configuration",printconfiguration);
 bool printroot; parseFlag("--root",printroot);
 if(printroot) fprintf(out,"%s\n",plumedRoot.c_str());

 if(printconfiguration){
    static const unsigned char conf [] ={
#include "Makefile.conf.xxd"
    , 0x00 };
    fprintf(out,"%s",conf);
    return 0;
  }

  return 0;
}



}

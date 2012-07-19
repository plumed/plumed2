/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
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
/**
This is a tool that allows you to obtain information about your plumed version
*/
//+ENDPLUMEDOC

class CLToolInfo:
public CLTool
{
public:
  int main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "provide informations about plumed";
  }
};


PLUMED_REGISTER_CLTOOL(CLToolInfo,"info")

int CLToolInfo::main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc){

// to avoid warnings:
 (void) in;

 bool printconfiguration(false);
 bool printhelp(false);

// Start parsing options
  string prefix("");
  string a("");
  for(int i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help"){
      printhelp=true;
      break;
    }
    if(a=="--configuration"){
      printconfiguration=true;
      break;
    } else if(a=="--root"){
      fprintf(out,"%s\n",plumedRoot.c_str());
    } else {
      string msg="ERROR: maximum one file at a time";
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    }
  }

  if(printhelp){
    fprintf(out,"%s",
 "Usage: info [options]\n"
 "Options:\n"
 "  [--help|-h]             : prints this help\n"
 "  [--configuration]       : prints the configuration file\n"
);
    return 0;
  }

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

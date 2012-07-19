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
#include "ActionRegister.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS manual
/**
manual is a tool that you can use to construct the manual page for 
a particular action
*/
//+ENDPLUMEDOC

class CLToolManual:
public CLTool
{
public:
  int main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "print out a description of the keywords for an action in html";
  }
};

PLUMED_REGISTER_CLTOOL(CLToolManual,"manual")

int CLToolManual::main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc){

// to avoid warnings:
 (void) in;

 std::string action="none";
 bool printhelp(false);

// Start parsing options
 string prefix("");
 std::string a("");
 for(int i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help"){
      printhelp=true;
      break;
    }
    if(a.find("--action=")==0){
      a.erase(0,a.find("=")+1);
      action=a; 
      prefix="";
    } else if(a=="--action"){
      prefix="--action=";
    } else {
      string msg="ERROR: unknown option"+a;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    }
 } 

 if(printhelp){
   fprintf(out,"%s",
   "Usage: info [options]\n"
   "Options:\n"
   "  [--help|-h]           : prints this help\n"
   "  [--action]            : print the manual for this action\n"
);
   return 0; 
 }

 if(action=="none"){
    fprintf(stderr,"missing --action flag\n");
    return 1;
 }

 std::cerr<<actionRegister(); 
 if( !actionRegister().printManual(action) ){
    fprintf(stderr,"specified action is not registered\n");
    return 1; 
 }

 return 0;
}

} // End of namespace

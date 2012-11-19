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
#include "PlumedException.h"
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace PLMD;

std::string PlumedException::format(const std::string&msg,const std::string&file,unsigned line,const std::string&function){
  std::string message;
  message="\n+++ Internal PLUMED error";
  if(file.length()>0){
    char cline[1000];
    sprintf(cline,"%u",line);
    message += "\n+++ file "+file+", line "+cline;
    if(function.length()>0) message +=", function "+function;
  }
  if(msg.length()>0) message +="\n+++ message: "+msg;
  return message;
}


PlumedException::PlumedException()
{
  this->msg=format("","",0,"");
  abortIfExceptionsAreDisabled();
}

PlumedException::PlumedException(const std::string&msg)
{
  this->msg=format(msg,"",0,"");
  abortIfExceptionsAreDisabled();
}

PlumedException::PlumedException(const std::string&msg,const std::string&file,unsigned line,const std::string&function)
{
  this->msg=format(msg,file,line,function);
  abortIfExceptionsAreDisabled();
}

void PlumedException::abortIfExceptionsAreDisabled(){
#if ! defined(__PLUMED_EXCEPTIONS)
  fprintf(stderr,"%s",what());
  fprintf(stderr,"\n");

  void* callstack[128];
  int i, frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);
  for (i = 0; i < frames; ++i) {
     fprintf(stderr,"%s\n", strs[i]);
  }
  free(strs);

  std::abort();
#endif
}



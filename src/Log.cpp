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
#include <cstdarg>
#include <cstring>
#include "Log.h"
#include "PlumedCommunicator.h"
#include "PlumedException.h"

using namespace std;
using namespace PLMD;

int Log::printf(const char*fmt,...){
  if(comm.Get_rank()>0)return 0;
  int pointer=strlen(buffer);
  va_list arg;
  va_start(arg, fmt);
  int r=vsnprintf(&buffer[pointer],buflen-pointer,fmt,arg);
  va_end(arg);
  plumed_massert(r>-1 && r<buflen-pointer,"error using fmt string " + std::string(fmt));

// Line is buffered until newline, then written with a PLUMED: prefix
  char*p1=buffer;
  char*p2;
  while((p2=strchr(p1,'\n'))){
    *p2='\0';
    fprintf(fp,"PLUMED: %s\n",p1);
    p1=p2+1;
  };
  memmove(buffer,p1,strlen(p1)+1);
  return r;
}

Log::Log(PlumedCommunicator &comm):
  fp(stdout),
  toBeClosed(false),
  comm(comm){
  buffer=new char[buflen];
  for(int i=0;i<buflen;i++) buffer[i]='\0';
}

Log::~Log(){
  if(comm.Get_rank()>0)return;
  if(!fp && toBeClosed)fclose(fp);
  delete [] buffer;
}

void Log::setFile(string str){
  if(comm.Get_rank()>0)return;
  fp=fopen(str.c_str(),"w");
  toBeClosed=true;
}

void Log::set(FILE*f){
  if(comm.Get_rank()>0)return;
  fp=f;
}

void Log::flush(){
  fflush(fp);
}


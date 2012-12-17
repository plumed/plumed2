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
#include "File.h"
#include "PlumedException.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>

#include <iostream>
#include <string>

namespace PLMD{

size_t OFile::llwrite(const char*ptr,size_t s){
  size_t r;
  if(linked) return linked->llwrite(ptr,s);
  if(! (comm && comm->Get_rank()>0)){
    if(!fp) plumed_merror("writing on uninitilized File");
    r=fwrite(ptr,1,s,fp);
  }
  if(comm) comm->Bcast(&r,1,0);
  return r;
}

OFile::OFile():
  linked(NULL),
  fieldChanged(false)
{
  fmtField();
  buflen=1;
  buffer=new char[buflen];
// these are set to zero to avoid valgrind errors
  for(unsigned i=0;i<buflen;++i) buffer[i]=0;
  buffer_string=new char [1000];
// these are set to zero to avoid valgrind errors
  for(unsigned i=0;i<1000;++i) buffer_string[i]=0;
}

OFile::~OFile(){
  delete [] buffer_string;
  delete [] buffer;
}

OFile& OFile::link(OFile&l){
  fp=NULL;
  linked=&l;
  return *this;
}

OFile& OFile::setLinePrefix(const std::string&l){
  linePrefix=l;
  return *this;
}

int OFile::printf(const char*fmt,...){
  size_t pointer=strlen(buffer);
  va_list arg;
  va_start(arg, fmt);
  int r=std::vsnprintf(&buffer[pointer],buflen-pointer,fmt,arg);
  va_end(arg);
  if(r>=buflen-pointer){
    int newlen=buflen;
    while(newlen<=r+pointer) newlen*=2;
    char* newbuf=new char [newlen];
    memmove(newbuf,buffer,buflen);
    for(int k=buflen;k<newlen;k++) newbuf[k]=0;
    delete [] buffer;
    buffer=newbuf;
    buflen=newlen;
    va_list arg;
    va_start(arg, fmt);
    r=std::vsnprintf(&buffer[pointer],buflen-pointer,fmt,arg);
    va_end(arg);
  }
  plumed_massert(r>-1 && r<buflen-pointer,"error using fmt string " + std::string(fmt));

// Line is buffered until newline, then written with a PLUMED: prefix
  char*p1=buffer;
  char*p2;
  while((p2=strchr(p1,'\n'))){
    *p2='\0';
    if(linePrefix.length()>0) llwrite(linePrefix.c_str(),linePrefix.length());
    llwrite(p1,std::strlen(p1));
    llwrite("\n",1);
    p1=p2+1;
  };
  memmove(buffer,p1,strlen(p1)+1);
  return r;
}

OFile& OFile::addConstantField(const std::string&name){
  Field f;
  f.name=name;
  const_fields.push_back(f);
  return *this;
}


OFile& OFile::clearFields(){
  fields.clear();
  const_fields.clear();
  previous_fields.clear();
  return *this;
}

OFile& OFile::fmtField(const std::string&fmt){
  this->fieldFmt=fmt;
  return *this;
}

OFile& OFile::fmtField(){
  this->fieldFmt="%23.16lg";
  return *this;
}

OFile& OFile::printField(const std::string&name,double v){
  sprintf(buffer_string,fieldFmt.c_str(),v);
  printField(name,buffer_string);
  return *this;
}

OFile& OFile::printField(const std::string&name,int v){
  sprintf(buffer_string," %d",v);
  printField(name,buffer_string);
  return *this;
}

OFile& OFile::printField(const std::string&name,const std::string & v){
  unsigned i;
  for(i=0;i<const_fields.size();i++) if(const_fields[i].name==name) break;
  if(i>=const_fields.size()){
    Field field;
    field.name=name;
    field.value=v;
    fields.push_back(field);
  } else {
    if(const_fields[i].value!=v) fieldChanged=true;
    const_fields[i].value=v;
  }
  return *this;
}

OFile& OFile::setupPrintValue( Value *val ){
  if( val->isPeriodic() ){
      addConstantField("min_" + val->getName() );
      addConstantField("max_" + val->getName() );
  }
  return *this;
}

OFile& OFile::printField( Value* val, const double& v ){
  printField( val->getName(), v );
  if( val->isPeriodic() ){
      std::string min, max; val->getDomain( min, max );
      printField( "min_" + val->getName(), min );
      printField("max_" + val->getName(), max ); 
  }  
  return *this;
}

OFile& OFile::printField(){
  bool reprint=false;
  if(fieldChanged || fields.size()!=previous_fields.size()){
    reprint=true;
  } else for(unsigned i=0;i<fields.size();i++){
    if( previous_fields[i].name!=fields[i].name ||
        (fields[i].constant && fields[i].value!=previous_fields[i].value) ){
      reprint=true;
      break;
    }
  }
  if(reprint){
    printf("#! FIELDS");
    for(unsigned i=0;i<fields.size();i++) printf(" %s",fields[i].name.c_str());
    printf("\n");
    for(unsigned i=0;i<const_fields.size();i++){
        printf("#! SET %s %s",const_fields[i].name.c_str(),const_fields[i].value.c_str());
        printf("\n");
    }
  }
  for(unsigned i=0;i<fields.size();i++) printf("%s",fields[i].value.c_str());
  printf("\n");
  previous_fields=fields;
  fields.clear();
  fieldChanged=false;
  return *this;
}

OFile& OFile::open(const std::string&path){
  plumed_assert(!cloned);
  eof=false;
  err=false;
  fp=NULL;
  this->path=path;
  if(plumed){
    this->path+=plumed->getSuffix();
  }
  if(plumed && plumed->getRestart()){
     fp=std::fopen(const_cast<char*>(this->path.c_str()),"a");
  } else {
     if(!comm || comm->Get_rank()==0){
       FILE* ff=std::fopen(const_cast<char*>(this->path.c_str()),"r");
       FILE* fff=NULL;
       if(ff){
         std::string backup;
         size_t found=this->path.find_last_of("/\\");
         std::string directory=this->path.substr(0,found+1);
         std::string file=this->path.substr(found+1);
         for(int i=0;;i++){
           std::string num;
           Tools::convert(i,num);
           backup=directory+"bck."+num+"."+file;
           fff=std::fopen(backup.c_str(),"r");
           if(!fff) break;
         }
         int check=rename(this->path.c_str(),backup.c_str());
         plumed_massert(check==0,"renaming "+this->path+" into "+backup+" failed for some reason");
       }
       if(ff) std::fclose(ff);
       if(fff) std::fclose(fff);
     }
     comm->Barrier();
     fp=std::fopen(const_cast<char*>(this->path.c_str()),"w");
  }
  if(plumed) plumed->insertFile(*this);
  return *this;
}

}



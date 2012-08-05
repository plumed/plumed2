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
#include "PlumedFile.h"
#include "PlumedException.h"
#include "Action.h"
#include "PlumedMain.h"
#include "PlumedCommunicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>

using namespace PLMD;

void PlumedFileBase::test(){
  PLMD::PlumedOFile pof;
  pof.open("ciao","w");
  pof.printf("%s\n","test1");
  pof.setLinePrefix("plumed: ");
  pof.printf("%s\n","test2");
  pof.setLinePrefix("");
  pof.addField("x1");
  pof.addConstantField("x2").printField("x2",67.0);
  pof.addField("x3");
  pof.printField("x1",10.0).printField("x3",20.12345678901234567890).printField();
  pof.printField("x1",10.0).printField("x3",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",10.0).printField("x2",777.0).printField("x1",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",67.0).printField("x1",18.0).printField();
  pof.close();

  PLMD::PlumedIFile pif;
  std::string s;
  pif.open("ciao","r");
  pif.getline(s); std::printf("%s\n",s.c_str());
  pif.getline(s); std::printf("%s\n",s.c_str());
  
  int x1,x2,x3;
  while(pif.scanField("x1",x1).scanField("x3",x2).scanField("x2",x3).scanField()){
    std::cout<<"CHECK "<<x1<<" "<<x2<<" "<<x3<<"\n";
  }
  pif.close();
}

size_t PlumedOFile::llwrite(const char*ptr,size_t s){
  size_t r;
  if(fp){
    if(! (comm && comm->Get_rank()>0)){
      r=fwrite(ptr,1,s,fp);
    }
    if(comm) comm->Bcast(&r,1,0);
  } else if(linked){
    (*linked)<<ptr;
    r=strlen(ptr);
  } else plumed_merror("writing on uninitilized PlumedFile");
  return r;
}

size_t PlumedIFile::llread(char*ptr,size_t s){
  plumed_assert(fp);
  size_t r;
  int leof=0;
  if(! (comm && comm->Get_rank()>0)){
    r=fread(ptr,1,s,fp);
    if(feof(fp)) leof=1;
  }
  if(comm) comm->Bcast(&r,1,0);
  if(comm) comm->Bcast(ptr,r,0);
  if(comm) comm->Bcast(&leof,r,0);
  if(leof) eof=true;
  return r;
}

PlumedFileBase& PlumedFileBase::link(FILE*fp){
  this->fp=fp;
  cloned=true;
  return *this;
}

PlumedFileBase& PlumedFileBase::flush(){
  fflush(fp);
  return *this;
}

PlumedFileBase& PlumedFileBase::link(PlumedCommunicator&comm){
  this->comm=&comm;
  return *this;
}

PlumedFileBase& PlumedFileBase::link(PlumedMain&plumed){
  this->plumed=&plumed;
  link(plumed.comm);
  return *this;
}

PlumedFileBase& PlumedFileBase::link(Action&action){
  this->action=&action;
  link(action.plumed);
  return *this;
}

PlumedFileBase& PlumedFileBase::open(const std::string& path,const std::string& mode){
  plumed_assert(!cloned);
  eof=false;
  err=false;
  fp=NULL;
  if(plumed){
    const std::string pathsuf=path+plumed->getSuffix();
    fp=std::fopen(const_cast<char*>(pathsuf.c_str()),const_cast<char*>(mode.c_str()));
  }
  if(!fp) fp=std::fopen(const_cast<char*>(path.c_str()),const_cast<char*>(mode.c_str()));
  plumed_massert(fp,"file " + path + "cannot be found");
  return *this;
}

void        PlumedFileBase::close(){
  plumed_assert(!cloned);
  eof=false;
  err=false;
  std::fclose(fp);
  fp=NULL;
}

PlumedFileBase::PlumedFileBase():
  fp(NULL),
  comm(NULL),
  plumed(NULL),
  action(NULL),
  cloned(false),
  eof(false),
  err(false)
{
}

PlumedFileBase::~PlumedFileBase()
{
  if(!cloned) fclose(fp);
}

PlumedFileBase::operator bool()const{
  return !eof;
}


PlumedOFile::PlumedOFile():
  linked(NULL),
  fieldChanged(false)
{
  fmtField();
  buffer=new char[buflen];
  buffer_string=new char [1000];
}

PlumedOFile::~PlumedOFile(){
  delete [] buffer_string;
  delete [] buffer;
}

PlumedOFile& PlumedOFile::link(PlumedOFile&l){
  fp=NULL;
  linked=&l;
  return *this;
}

PlumedOFile& PlumedOFile::setLinePrefix(const std::string&l){
  linePrefix=l;
  return *this;
}

int PlumedOFile::printf(const char*fmt,...){
  int pointer=strlen(buffer);
  va_list arg;
  va_start(arg, fmt);
  int r=std::vsnprintf(&buffer[pointer],buflen-pointer,fmt,arg);
  va_end(arg);
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

PlumedOFile& PlumedOFile::addField(const std::string&name){
  Field f;
  f.name=name;
  fields.push_back(f);
  fieldChanged=true;
  return *this;
}

PlumedOFile& PlumedOFile::addConstantField(const std::string&name){
  Field f;
  f.name=name;
  f.constant=true;
  fields.push_back(f);
  fieldChanged=true;
  return *this;
}


PlumedOFile& PlumedOFile::clearFields(){
  fields.clear();
  fieldChanged=true;
  return *this;
}

PlumedOFile& PlumedOFile::fmtField(const std::string&fmt){
  this->fieldFmt=fmt;
  return *this;
}

PlumedOFile& PlumedOFile::fmtField(){
  this->fieldFmt="%23.16lg";
  return *this;
}

unsigned PlumedOFile::findField(const std::string&name)const{
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  if(i>=fields.size()) plumed_merror(name);
  return i;
}

PlumedOFile& PlumedOFile::printField(const std::string&name,double v){
  sprintf(buffer_string,fieldFmt.c_str(),v);
  printField(name,buffer_string);
  return *this;
}

PlumedOFile& PlumedOFile::printField(const std::string&name,int v){
  sprintf(buffer_string," %d",v);
  printField(name,buffer_string);
  return *this;
}

PlumedOFile& PlumedOFile::printField(const std::string&name,const std::string & v){
  unsigned i=findField(name);
  if(fields[i].constant) fieldChanged=true;
  fields[i].value=v;
  fields[i].set=true;
  return *this;
}

PlumedOFile& PlumedOFile::printField(){
  if(fieldChanged){
    printf("#! FIELDS");
    for(unsigned i=0;i<fields.size();i++){
      if(!fields[i].constant)
        printf(" %s",fields[i].name.c_str());
    }
    printf("\n");
    for(unsigned i=0;i<fields.size();i++)
      if(fields[i].constant){
        printf("#! SET %s %s",fields[i].name.c_str(),fields[i].value.c_str());
        printf("\n");
    }
  }
  fieldChanged=false;
  for(unsigned i=0;i<fields.size();i++){
    plumed_assert(fields[i].set);
    if(!fields[i].constant){
      printf("%s",fields[i].value.c_str());
      fields[i].set=false;
    }
  }
  printf("\n");
  return *this;
}

PlumedIFile& PlumedIFile::advanceField(){
  plumed_assert(!inMiddleOfField);
  std::string line;
  bool done=false;
  while(!done){
    getline(line);
    if(!*this) return *this;
    std::vector<std::string> words=Tools::getWords(line);
    if(words.size()>=2 && words[0]=="#!" && words[1]=="FIELDS"){
      fields.clear();
      for(unsigned i=2;i<words.size();i++){
        Field field;
        field.name=words[i];
        fields.push_back(field);
      }
    } else if(words.size()==4 && words[0]=="#!" && words[1]=="SET"){
      Field field;
      field.name=words[2];
      field.value=words[3];
      field.constant=true;
      fields.push_back(field);
    } else {
      unsigned nf=0;
      for(unsigned i=0;i<fields.size();i++) if(!fields[i].constant) nf++;
      Tools::trimComments(line);
      words=Tools::getWords(line);
      plumed_assert(nf==words.size());
      unsigned j=0;
      for(unsigned i=0;i<fields.size();i++){
        if(fields[i].constant) continue;
        fields[i].value=words[j];
        fields[i].read=false;
        j++;
      }
      done=true;
    }
  }
  inMiddleOfField=true;
  return *this;
}

PlumedIFile& PlumedIFile::scanFieldList(std::vector<std::string>&s){
  if(!inMiddleOfField) advanceField();
  if(!*this) return *this;
  s.clear();
  for(unsigned i=0;i<fields.size();i++)
    s.push_back(fields[i].name);
  return *this;
}

PlumedIFile& PlumedIFile::scanField(const std::string&name,std::string&str){
  if(!inMiddleOfField) advanceField();
  if(!*this) return *this;
  unsigned i=findField(name);
  str=fields[i].value;
  fields[i].read=true;
  return *this;
}

PlumedIFile& PlumedIFile::scanField(const std::string&name,double &x){
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

PlumedIFile& PlumedIFile::scanField(const std::string&name,int &x){
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}


PlumedIFile& PlumedIFile::scanField(){
  for(unsigned i=0;i<fields.size();i++){
    plumed_assert(fields[i].read);
  }
  inMiddleOfField=false;
  return *this;
}

PlumedIFile::PlumedIFile():
  inMiddleOfField(false)
{
}

PlumedIFile::~PlumedIFile(){
  plumed_assert(!inMiddleOfField);
}

PlumedIFile& PlumedIFile::getline(std::string &str){
  char tmp;
  str="";
  while(llread(&tmp,1)==1 && tmp && tmp!='\n' && !eof){
    str+=tmp;
  }
  return *this;
}

unsigned PlumedIFile::findField(const std::string&name)const{
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  if(i>=fields.size()) plumed_merror(name);
  return i;
}




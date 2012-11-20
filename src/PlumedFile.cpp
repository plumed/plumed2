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
#include "Value.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>

#include <iostream>
#include <string>

using namespace PLMD;

void PlumedFileBase::test(){
  PLMD::PlumedOFile pof;
  pof.open("ciao");
  pof.printf("%s\n","test1");
  pof.setLinePrefix("plumed: ");
  pof.printf("%s\n","test2");
  pof.setLinePrefix("");
  pof.addConstantField("x2").printField("x2",67.0);
  pof.printField("x1",10.0).printField("x3",20.12345678901234567890).printField();
  pof.printField("x1",10.0).printField("x3",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",10.0).printField("x2",777.0).printField("x1",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",67.0).printField("x1",18.0).printField();
  pof.close();

  PLMD::PlumedIFile pif;
  std::string s;
  pif.open("ciao");
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
  if(linked) return linked->llwrite(ptr,s);
  if(! (comm && comm->Get_rank()>0)){
    if(!fp) plumed_merror("writing on uninitilized PlumedFile");
    r=fwrite(ptr,1,s,fp);
  }
  if(comm) comm->Bcast(&r,1,0);
  return r;
}

size_t PlumedIFile::llread(char*ptr,size_t s){
  plumed_assert(fp);
  size_t r;
  r=fread(ptr,1,s,fp);
  if(feof(fp))   eof=true;
  if(ferror(fp)) err=true;
  return r;
}

PlumedFileBase& PlumedFileBase::link(FILE*fp){
  this->fp=fp;
  cloned=true;
  return *this;
}

PlumedFileBase& PlumedFileBase::flush(){
  fflush(fp);
  if(heavyFlush){
    fclose(fp);
    fp=std::fopen(const_cast<char*>(path.c_str()),"a");
  }
  return *this;
}

PlumedFileBase& PlumedFileBase::link(PlumedCommunicator&comm){
  plumed_massert(!fp,"cannot link an already open file");
  this->comm=&comm;
  return *this;
}

PlumedFileBase& PlumedFileBase::link(PlumedMain&plumed){
  plumed_massert(!fp,"cannot link an already open file");
  this->plumed=&plumed;
  link(plumed.comm);
  return *this;
}

PlumedFileBase& PlumedFileBase::link(Action&action){
  plumed_massert(!fp,"cannot link an already open file");
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
    this->path=path+plumed->getSuffix();
    fp=std::fopen(const_cast<char*>(this->path.c_str()),const_cast<char*>(mode.c_str()));
  }
  if(!fp){
    this->path=path;
    fp=std::fopen(const_cast<char*>(this->path.c_str()),const_cast<char*>(mode.c_str()));
  }
  if(plumed) plumed->insertFile(*this);
  plumed_massert(fp,"file " + path + "cannot be found");
  return *this;
}


bool PlumedFileBase::FileExist(const std::string& path){
  FILE *ff=NULL;
  bool do_exist=false;
  if(plumed){
    this->path=path+plumed->getSuffix();
    ff=std::fopen(const_cast<char*>(this->path.c_str()),"r");
  }
  if(!ff){
    this->path=path;
    ff=std::fopen(const_cast<char*>(this->path.c_str()),"r");
  }
  if(ff) {do_exist=true; fclose(ff);}
  return do_exist; 
}

bool PlumedFileBase::isOpen(){
  bool isopen=false;
  if(fp) isopen=true;
  return isopen; 
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
  err(false),
  heavyFlush(false)
{
}

PlumedFileBase::~PlumedFileBase()
{
  if(plumed) plumed->eraseFile(*this);
  if(!cloned && fp) fclose(fp);
}

PlumedFileBase::operator bool()const{
  return !eof;
}


PlumedOFile::PlumedOFile():
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

PlumedOFile& PlumedOFile::addConstantField(const std::string&name){
  Field f;
  f.name=name;
  const_fields.push_back(f);
  return *this;
}


PlumedOFile& PlumedOFile::clearFields(){
  fields.clear();
  const_fields.clear();
  previous_fields.clear();
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

PlumedOFile& PlumedOFile::setupPrintValue( Value *val ){
  if( val->isPeriodic() ){
      addConstantField("min_" + val->getName() );
      addConstantField("max_" + val->getName() );
  }
  return *this;
}

PlumedOFile& PlumedOFile::printField( Value* val, const double& v ){
  printField( val->getName(), v );
  if( val->isPeriodic() ){
      std::string min, max; val->getDomain( min, max );
      printField( "min_" + val->getName(), min );
      printField("max_" + val->getName(), max ); 
  }  
  return *this;
}

PlumedOFile& PlumedOFile::printField(){
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

PlumedOFile& PlumedOFile::open(const std::string&path){
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


PlumedIFile& PlumedIFile::advanceField(){
  plumed_assert(!inMiddleOfField);
  std::string line;
  bool done=false;
  while(!done){
    getline(line);
    if(!*this){return *this;}
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

PlumedIFile& PlumedIFile::open(const std::string&name){
  PlumedFileBase::open(name,"r");
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

bool PlumedIFile::FieldExist(const std::string& s){
     std::vector<std::string> slist;
     scanFieldList(slist);
     int mycount = (int) std::count(slist.begin(), slist.end(), s);
     if(mycount>0) return true;
     else return false;
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

PlumedIFile& PlumedIFile::scanField(Value* val){
  double ff; scanField(  val->getName(), ff );
  val->set( ff );
  if( FieldExist("min_" + val->getName() ) ){ 
      std::string min, max;
      scanField("min_" + val->getName(), min );
      scanField("max_" + val->getName(), max );
      val->setDomain( min, max ); 
  } else {
      val->setNotPeriodic();
  }
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
  fpos_t pos;
  fgetpos(fp,&pos);
  while(llread(&tmp,1)==1 && tmp && tmp!='\n' && !eof && !err){
    str+=tmp;
  }
  if(tmp!='\n' || err){
    eof = true;
    str="";
    fsetpos(fp,&pos);
  }
  return *this;
}

unsigned PlumedIFile::findField(const std::string&name)const{
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  if(i>=fields.size()) plumed_merror(name);
  return i;
}

void PlumedIFile::reset(bool reset){
 eof = reset;
 err = reset;
 if(!reset) clearerr(fp);
 return;
} 

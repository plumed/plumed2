#include "PlumedFile.h"
#include "PlumedException.h"
#include "Action.h"
#include "PlumedMain.h"
#include "PlumedCommunicator.h"
#include <cstdarg>
#include <cstring>

using namespace PLMD;

void PlumedFileBase::test(){
  PLMD::PlumedOFile pof;
  pof.open("ciao","w");
  pof.printf("%s\n","test1");
  pof.setLinePrefix("plumed: ");
  pof.printf("%s\n","test2");
  pof.addField("x1");
  pof.addField("x2",67.0);
  pof.addField("x3");
  pof.printField("x1",10.0).printField("x3",20.12345678901234567890).printField();
  pof.printField("x1",10.0).printField("x3",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",10.0).printField("x2",777).printField("x1",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",67.0).printField("x1",18.0).printField();
  pof.close();

  PLMD::PlumedIFile pif;
  std::string s;
  pif.open("ciao","r");
  pif.getline(s); std::printf("%s\n",s.c_str());
  pif.getline(s); std::printf("%s\n",s.c_str());
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
  if(! (comm && comm->Get_rank()>0)){
    r=fread(ptr,1,s,fp);
  }
  if(comm) comm->Bcast(&r,1,0);
  if(comm) comm->Bcast(ptr,r,0);
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
  std::fclose(fp);
  fp=NULL;
}

PlumedFileBase::PlumedFileBase():
  fp(NULL),
  comm(NULL),
  plumed(NULL),
  action(NULL),
  cloned(false)
{
}

PlumedFileBase::~PlumedFileBase()
{
  if(!cloned) fclose(fp);
}

PlumedOFile::PlumedOFile():
  linked(NULL),
  fieldChanged(false),
  fieldFmt("%22.16lg")
{
  buffer=new char[buflen];
}

PlumedOFile::~PlumedOFile(){
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

PlumedOFile& PlumedOFile::addField(const std::string&name,double v){
  Field f;
  f.name=name;
  f.value=v;
  f.constant=true;
  f.set=true;
  fields.push_back(f);
  fieldChanged=true;
  return *this;
}

PlumedOFile& PlumedOFile::clearFields(){
  fields.clear();
  fieldChanged=true;
  return *this;
}

PlumedOFile& PlumedOFile::fmtFields(const std::string&fmt){
  fieldFmt=fmt;
  return *this;
}

PlumedOFile& PlumedOFile::fmtField(const std::string&name,const std::string&fmt){
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  plumed_assert(i<fields.size());
  fields[i].fmt=fmt;
  return *this;
}

PlumedOFile& PlumedOFile::printField(const std::string&name,double v){
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  plumed_assert(i<fields.size());
  if(fields[i].constant) fieldChanged=true;
  fields[i].value=v;
  fields[i].set=true;
  return *this;
}

PlumedOFile& PlumedOFile::printField(){
  if(fieldChanged){
    printf("#! FIELDS");
    for(unsigned i=0;i<fields.size();i++){
      printf(" %s",fields[i].name.c_str());
    }
    printf("\n");
    for(unsigned i=0;i<fields.size();i++)
      if(fields[i].constant){
        std::string fmt;
        if(fields[i].fmt.length()>0) fmt=fields[i].fmt;
        else fmt=fieldFmt;
        printf("#! SET %s ",fields[i].name.c_str());
        printf(fmt.c_str(),fields[i].value);
        printf("\n");
    }
  }
  fieldChanged=false;
  for(unsigned i=0;i<fields.size();i++){
    plumed_assert(fields[i].set);
    if(!fields[i].constant){
      std::string fmt;
      if(fields[i].fmt.length()>0) fmt=fields[i].fmt;
      else fmt=fieldFmt;
//      printf(" ");
      printf(fmt.c_str(),fields[i].value);
      fields[i].set=false;
    }
  }
  printf("\n");
  return *this;
}

PlumedIFile& PlumedIFile::scanFieldList(std::vector<std::string>&s){
  s=fields;
  return *this;
}

PlumedIFile& PlumedIFile::scanField(const std::string&,double&){
  plumed_error();
  return *this;
}

PlumedIFile& PlumedIFile::scanField(){
  plumed_error();
  return *this;
}

PlumedIFile::PlumedIFile(){
}

PlumedIFile::~PlumedIFile(){
}

PlumedIFile& PlumedIFile::getline(std::string &str){
  char tmp;
  str="";
  while(llread(&tmp,1)==1 && tmp && tmp!='\n'){
    str+=tmp;
  }
  return *this;
}




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
#include "IFile.h"
#include "Exception.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>
#include <cmath>

#include <iostream>
#include <string>
#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif

namespace PLMD {

size_t IFile::llread(char*ptr,size_t s) {
  plumed_assert(fp);
  size_t r;
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    int rr=gzread(gzFile(gzfp),ptr,s);
    if(rr==0)   eof=true;
    if(rr<0)    err=true;
    r=rr;
#else
    plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
  } else {
    r=fread(ptr,1,s,fp);
    if(feof(fp))   eof=true;
    if(ferror(fp)) err=true;
  }
  return r;
}

IFile& IFile::advanceField() {
  plumed_assert(!inMiddleOfField);
  std::string line;
  bool done=false;
  while(!done) {
    getline(line);
// using explicit conversion not to confuse cppcheck 1.86
    if(!bool(*this)) {return *this;}
    std::vector<std::string> words=Tools::getWords(line);
    if(words.size()>=2 && words[0]=="#!" && words[1]=="FIELDS") {
      fields.clear();
      for(unsigned i=2; i<words.size(); i++) {
        Field field;
        field.name=words[i];
        fields.push_back(field);
      }
    } else if(words.size()==4 && words[0]=="#!" && words[1]=="SET") {
      Field field;
      field.name=words[2];
      field.value=words[3];
      field.constant=true;
      fields.push_back(field);
    } else {
      unsigned nf=0;
      for(unsigned i=0; i<fields.size(); i++) if(!fields[i].constant) nf++;
      Tools::trimComments(line);
      words=Tools::getWords(line);
      if( words.size()==nf ) {
        unsigned j=0;
        for(unsigned i=0; i<fields.size(); i++) {
          if(fields[i].constant) continue;
          fields[i].value=words[j];
          fields[i].read=false;
          j++;
        }
        done=true;
      } else if( !words.empty() ) {
        plumed_merror("file " + getPath() + ": mismatch between number of fields in file and expected number");
      }
    }
  }
  inMiddleOfField=true;
  return *this;
}

IFile& IFile::open(const std::string&path) {
  plumed_massert(!cloned,"file "+path+" appears to be cloned");
  eof=false;
  err=false;
  fp=NULL;
  gzfp=NULL;
  bool do_exist=FileExist(path);
  plumed_massert(do_exist,"file " + path + " cannot be found");
  fp=std::fopen(const_cast<char*>(this->path.c_str()),"r");
  if(Tools::extension(this->path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
    gzfp=(void*)gzopen(const_cast<char*>(this->path.c_str()),"r");
#else
    plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
  }
  if(plumed) plumed->insertFile(*this);
  return *this;
}

IFile& IFile::scanFieldList(std::vector<std::string>&s) {
  if(!inMiddleOfField) advanceField();
// using explicit conversion not to confuse cppcheck 1.86
  if(!bool(*this)) return *this;
  s.clear();
  for(unsigned i=0; i<fields.size(); i++)
    s.push_back(fields[i].name);
  return *this;
}

bool IFile::FieldExist(const std::string& s) {
  std::vector<std::string> slist;
  scanFieldList(slist);
  int mycount = (int) std::count(slist.begin(), slist.end(), s);
  if(mycount>0) return true;
  else return false;
}

IFile& IFile::scanField(const std::string&name,std::string&str) {
  if(!inMiddleOfField) advanceField();
// using explicit conversion not to confuse cppcheck 1.86
  if(!bool(*this)) return *this;
  unsigned i=findField(name);
  str=fields[i].value;
  fields[i].read=true;
  return *this;
}

IFile& IFile::scanField(const std::string&name,double &x) {
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(const std::string&name,int &x) {
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(const std::string&name,long int &x) {
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(const std::string&name,unsigned &x) {
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(const std::string&name,long unsigned &x) {
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(Value* val) {
  double ff=std::numeric_limits<double>::quiet_NaN(); // this is to be sure a NaN value is replaced upon failure
  scanField(  val->getName(), ff );
  val->set( ff );
  if( FieldExist("min_" + val->getName() ) ) {
    std::string min, max;
    scanField("min_" + val->getName(), min );
    scanField("max_" + val->getName(), max );
    val->setDomain( min, max );
  } else {
    val->setNotPeriodic();
  }
  return *this;
}

IFile& IFile::scanField() {
  if(!ignoreFields) {
    for(unsigned i=0; i<fields.size(); i++) {
      plumed_massert(fields[i].read,"field "+fields[i].name+" was not read: all the fields need to be read otherwise you could miss important infos" );
    }
  }
  inMiddleOfField=false;
  return *this;
}

IFile::IFile():
  inMiddleOfField(false),
  ignoreFields(false),
  noEOL(false)
{
}

IFile::~IFile() {
  if(inMiddleOfField) std::cerr<<"WARNING: IFile closed in the middle of reading. seems strange!\n";
}

IFile& IFile::getline(std::string &str) {
  char tmp=0;
  str="";
  fpos_t pos;
  fgetpos(fp,&pos);
  while(llread(&tmp,1)==1 && tmp && tmp!='\n' && tmp!='\r' && !eof && !err) {
    str+=tmp;
  }
  if(tmp=='\r') {
    llread(&tmp,1);
    plumed_massert(tmp=='\n',"plumed only accepts \\n (unix) or \\r\\n (dos) new lines");
  }
  if(eof && noEOL) {
    if(str.length()>0) eof=false;
  } else if(eof || err || tmp!='\n') {
    eof = true;
    str="";
    if(!err) fsetpos(fp,&pos);
// there was a fsetpos here that apparently is not necessary
//  fsetpos(fp,&pos);
// I think it was necessary to have rewind working correctly
// after end of file. Since rewind is not used now anywhere,
// it should be ok not to reset position.
// This is necessary so that eof works properly for emacs files
// with no endline at end of file.
  }
  return *this;
}

unsigned IFile::findField(const std::string&name)const {
  unsigned i;
  for(i=0; i<fields.size(); i++) if(fields[i].name==name) break;
  if(i>=fields.size()) {
    plumed_merror("file " + getPath() + ": field " + name + " cannot be found");
  }
  return i;
}

void IFile::reset(bool reset) {
  eof = reset;
  err = reset;
  if(!reset && fp) clearerr(fp);
#ifdef __PLUMED_HAS_ZLIB
  if(!reset && gzfp) gzclearerr(gzFile(gzfp));
#endif
  return;
}

void IFile::allowIgnoredFields() {
  ignoreFields=true;
}

void IFile::allowNoEOL() {
  noEOL=true;
}

}

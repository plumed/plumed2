/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "OFile.h"
#include "Exception.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>

#include <iostream>
#include <string>
#include <cstdlib>
#include <cerrno>

#include <memory>
#include <utility>

#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif

namespace PLMD {

size_t OFile::llwrite(const char*ptr,size_t s) {
  size_t r;
  if(linked) return linked->llwrite(ptr,s);
  if(! (comm && comm->Get_rank()>0)) {
    if(!fp) plumed_merror("writing on uninitilized File");
    if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
      r=gzwrite(gzFile(gzfp),ptr,s);
#else
      plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
    } else {
      r=fwrite(ptr,1,s,fp);
    }
  }
//  This barrier is apparently useless since it comes
//  just before a Bcast.
//
//  Anyway, it looks like it is solving an issue that appeared on
//  TRAVIS (at least on my laptop) so I add it here.
//  GB
  if(comm) comm->Barrier();


  if(comm) comm->Bcast(r,0);
  return r;
}

OFile::OFile():
  linked(NULL),
  fieldChanged(false),
  backstring("bck"),
  enforceRestart_(false),
  enforceBackup_(false)
{
  fmtField();
  buflen=1;
  actual_buffer_length=0;
  buffer.reset(new char [buflen]);
// these are set to zero to avoid valgrind errors
  for(int i=0; i<buflen; ++i) buffer[i]=0;
  buffer_string.reset(new char [1000]);
// these are set to zero to avoid valgrind errors
  for(unsigned i=0; i<1000; ++i) buffer_string[i]=0;
}

OFile& OFile::link(OFile&l) {
  fp=NULL;
  gzfp=NULL;
  linked=&l;
  return *this;
}

OFile& OFile::setLinePrefix(const std::string&l) {
  linePrefix=l;
  return *this;
}

int OFile::printf(const char*fmt,...) {
  va_list arg;
  va_start(arg, fmt);
  int r=std::vsnprintf(&buffer[actual_buffer_length],buflen-actual_buffer_length,fmt,arg);
  va_end(arg);
  if(r>=buflen-actual_buffer_length) {
    int newlen=buflen;
    while(newlen<=r+actual_buffer_length) newlen*=2;
    std::unique_ptr<char[]> newbuf{new char [newlen]};
    memmove(newbuf.get(),buffer.get(),buflen);
    for(int k=buflen; k<newlen; k++) newbuf[k]=0;
    buffer=std::move(newbuf);
    buflen=newlen;
    va_list arg;
    va_start(arg, fmt);
    r=std::vsnprintf(&buffer[actual_buffer_length],buflen-actual_buffer_length,fmt,arg);
    va_end(arg);
  }
  plumed_massert(r>-1 && r<buflen-actual_buffer_length,"error using fmt string " + std::string(fmt));

// Line is buffered until newline, then written with a PLUMED: prefix
  char*p1=buffer.get();
  char*p2;
// newline is only searched in the just added portion:
  char*psearch=p1+actual_buffer_length;
  actual_buffer_length+=r;
  while((p2=strchr(psearch,'\n'))) {
    if(linePrefix.length()>0) llwrite(linePrefix.c_str(),linePrefix.length());
    llwrite(p1,p2-p1+1);
    actual_buffer_length-=(p2-p1)+1;
    p1=p2+1;
    psearch=p1;
  };
  if(buffer.get()!=p1) memmove(buffer.get(),p1,actual_buffer_length);
  return r;
}

OFile& OFile::addConstantField(const std::string&name) {
  Field f;
  f.name=name;
  const_fields.push_back(f);
  return *this;
}


OFile& OFile::clearFields() {
  fields.clear();
  const_fields.clear();
  previous_fields.clear();
  return *this;
}

OFile& OFile::fmtField(const std::string&fmt) {
  this->fieldFmt=fmt;
  return *this;
}

OFile& OFile::fmtField() {
  this->fieldFmt="%23.16lg";
  return *this;
}

OFile& OFile::printField(const std::string&name,double v) {
// When one tries to print -nan we print nan instead.
// The distinction between +nan and -nan is not well defined
// Always printing nan simplifies some regtest (special functions computed our of range).
  if(std::isnan(v)) v=std::numeric_limits<double>::quiet_NaN();
  sprintf(buffer_string.get(),fieldFmt.c_str(),v);
  printField(name,buffer_string.get());
  return *this;
}

OFile& OFile::printField(const std::string&name,int v) {
  sprintf(buffer_string.get()," %d",v);
  printField(name,buffer_string.get());
  return *this;
}

OFile& OFile::printField(const std::string&name,const std::string & v) {
  unsigned i;
  for(i=0; i<const_fields.size(); i++) if(const_fields[i].name==name) break;
  if(i>=const_fields.size()) {
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

OFile& OFile::setupPrintValue( Value *val ) {
  if( val->isPeriodic() ) {
    addConstantField("min_" + val->getName() );
    addConstantField("max_" + val->getName() );
  }
  return *this;
}

OFile& OFile::printField( Value* val, const double& v ) {
  printField( val->getName(), v );
  if( val->isPeriodic() ) {
    std::string min, max; val->getDomain( min, max );
    printField( "min_" + val->getName(), min );
    printField("max_" + val->getName(), max );
  }
  return *this;
}

OFile& OFile::printField() {
  bool reprint=false;
  if(fieldChanged || fields.size()!=previous_fields.size()) {
    reprint=true;
  } else for(unsigned i=0; i<fields.size(); i++) {
      if( previous_fields[i].name!=fields[i].name ||
          (fields[i].constant && fields[i].value!=previous_fields[i].value) ) {
        reprint=true;
        break;
      }
    }
  if(reprint) {
    printf("#! FIELDS");
    for(unsigned i=0; i<fields.size(); i++) printf(" %s",fields[i].name.c_str());
    printf("\n");
    for(unsigned i=0; i<const_fields.size(); i++) {
      printf("#! SET %s %s",const_fields[i].name.c_str(),const_fields[i].value.c_str());
      printf("\n");
    }
  }
  for(unsigned i=0; i<fields.size(); i++) printf("%s",fields[i].value.c_str());
  printf("\n");
  previous_fields=fields;
  fields.clear();
  fieldChanged=false;
  return *this;
}

void OFile::setBackupString( const std::string& str ) {
  backstring=str;
}

void OFile::backupAllFiles( const std::string& str ) {
  if(str=="/dev/null") return;
  plumed_assert( backstring!="bck" && !checkRestart());
  size_t found=str.find_last_of("/\\");
  std::string filename = appendSuffix(str,getSuffix());
  std::string directory=filename.substr(0,found+1);
  std::string file=filename.substr(found+1);
  if( FileExist(filename) ) backupFile("bck", filename);
  for(int i=0;; i++) {
    std::string num; Tools::convert(i,num);
    std::string filestr = directory + backstring + "." + num + "." + file;
    if( !FileExist(filestr) ) break;
    backupFile( "bck", filestr);
  }
}

void OFile::backupFile( const std::string& bstring, const std::string& fname ) {
  if(fname=="/dev/null") return;
  int maxbackup=100;
  if(std::getenv("PLUMED_MAXBACKUP")) Tools::convert(std::getenv("PLUMED_MAXBACKUP"),maxbackup);
  if(maxbackup>0 && (!comm || comm->Get_rank()==0)) {
    FILE* ff=std::fopen(const_cast<char*>(fname.c_str()),"r");
    if(ff) {
      std::fclose(ff);
      std::string backup;
      size_t found=fname.find_last_of("/\\");
      std::string directory=fname.substr(0,found+1);
      std::string file=fname.substr(found+1);
      for(int i=0;; i++) {
        std::string num;
        Tools::convert(i,num);
        if(i>maxbackup) plumed_merror("cannot backup file "+file+" maximum number of backup is "+num+"\n");
        backup=directory+bstring +"."+num+"."+file;
        FILE* fff=std::fopen(backup.c_str(),"r");
        if(!fff) break;
        else std::fclose(fff);
      }
      int check=rename(fname.c_str(),backup.c_str());
      plumed_massert(check==0,"renaming "+fname+" into "+backup+" failed for reason: "+strerror(errno));
    }
  }
}

OFile& OFile::open(const std::string&path) {
  plumed_assert(!cloned);
  eof=false;
  err=false;
  fp=NULL;
  gzfp=NULL;
  this->path=path;
  this->path=appendSuffix(path,getSuffix());
  if(checkRestart()) {
    fp=std::fopen(const_cast<char*>(this->path.c_str()),"a");
    mode="a";
    if(Tools::extension(this->path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
      gzfp=(void*)gzopen(const_cast<char*>(this->path.c_str()),"a9");
#else
      plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
    }
  } else {
    backupFile( backstring, this->path );
    if(comm)comm->Barrier();
    fp=std::fopen(const_cast<char*>(this->path.c_str()),"w");
    mode="w";
    if(Tools::extension(this->path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
      gzfp=(void*)gzopen(const_cast<char*>(this->path.c_str()),"w9");
#else
      plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
    }
  }
  if(plumed) plumed->insertFile(*this);
  return *this;
}

OFile& OFile::rewind() {
// we use here "hard" rewind, which means close/reopen
// the reason is that normal rewind does not work when in append mode
// moreover, we can take a backup of the file
  plumed_assert(fp);
  clearFields();
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    gzclose((gzFile)gzfp);
#endif
  } else fclose(fp);
  if(!comm || comm->Get_rank()==0) {
    std::string fname=this->path;
    size_t found=fname.find_last_of("/\\");
    std::string directory=fname.substr(0,found+1);
    std::string file=fname.substr(found+1);
    std::string backup=directory+backstring +".last."+file;
    int check=rename(fname.c_str(),backup.c_str());
    plumed_massert(check==0,"renaming "+fname+" into "+backup+" failed for reason: "+strerror(errno));
  }
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    gzfp=(void*)gzopen(const_cast<char*>(this->path.c_str()),"w9");
#endif
  } else fp=std::fopen(const_cast<char*>(path.c_str()),"w");
  return *this;
}

FileBase& OFile::flush() {
  if(heavyFlush) {
    if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
      gzclose(gzFile(gzfp));
      gzfp=(void*)gzopen(const_cast<char*>(path.c_str()),"a");
#endif
    } else {
      fclose(fp);
      fp=std::fopen(const_cast<char*>(path.c_str()),"a");
    }
  } else {
    FileBase::flush();
    // if(gzfp) gzflush(gzFile(gzfp),Z_FINISH);
    // for some reason flushing with Z_FINISH has problems on linux
    // I thus use this (incomplete) flush
#ifdef __PLUMED_HAS_ZLIB
    if(gzfp) gzflush(gzFile(gzfp),Z_FULL_FLUSH);
#endif
  }
  return *this;
}

bool OFile::checkRestart()const {
  if(enforceRestart_) return true;
  else if(enforceBackup_) return false;
  else if(action) return action->getRestart();
  else if(plumed) return plumed->getRestart();
  else return false;
}

OFile& OFile::enforceRestart() {
  enforceRestart_=true;
  enforceBackup_=false;
  return *this;
}

OFile& OFile::enforceBackup() {
  enforceBackup_=true;
  enforceRestart_=false;
  return *this;
}


}

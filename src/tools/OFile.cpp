/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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

#include <string>
#include <cstdlib>
#include <cerrno>
#include <utility>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif

namespace PLMD {

size_t OFile::llwrite(const char*ptr,size_t s) {
  if(linked) {
    return linked->llwrite(ptr,s);
  }
  // Only the first process writes. Until v2.10 this was followed by a Barrier and by a
  // Bcast of the number of written bytes, i.e. two collective operations per line,
  // which slowed down runs with many processes and hung if the processes did not
  // write the same number of lines. The return value is not used.
  if(comm && comm->Get_rank()>0) {
    return s;
  }
  if(!fp) {
    plumed_merror("writing on uninitialized File");
  }
  size_t r=0;
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    r=gzwrite(gzFile(gzfp),ptr,s);
#else
    plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
  } else {
    r=std::fwrite(ptr,1,s,fp);
  }
  return r;
}

OFile::OFile():
  linked(NULL),
  fieldChanged(false),
  backstring("bck"),
  enforceRestart_(false),
  enforceBackup_(false) {
  fmtField();
  buflen=1;
  actual_buffer_length=0;
  buffer.resize(buflen);
// these are set to zero to avoid valgrind errors
  for(int i=0; i<buflen; ++i) {
    buffer[i]=0;
  }
// these are set to zero to avoid valgrind errors
  buffer_string.resize(1000,0);
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
    while(newlen<=r+actual_buffer_length) {
      newlen*=2;
    }
    std::vector<char> newbuf(newlen);
    std::memmove(newbuf.data(),buffer.data(),buflen);
    for(int k=buflen; k<newlen; k++) {
      newbuf[k]=0;
    }
    std::swap(buffer,newbuf);
    buflen=newlen;
    va_list args;
    va_start(args, fmt);
    r=std::vsnprintf(&buffer[actual_buffer_length],buflen-actual_buffer_length,fmt,args);
    va_end(args);
  }
  plumed_massert(r>-1 && r<buflen-actual_buffer_length,"error using fmt string " + std::string(fmt));

// Line is buffered until newline, then written with a PLUMED: prefix
  char*p1=buffer.data();
  char*p2;
// newline is only searched in the just added portion:
  char*psearch=p1+actual_buffer_length;
  actual_buffer_length+=r;
  while((p2=std::strchr(psearch,'\n'))) {
    if(linePrefix.length()>0) {
      llwrite(linePrefix.c_str(),linePrefix.length());
    }
    llwrite(p1,p2-p1+1);
    actual_buffer_length-=(p2-p1)+1;
    p1=p2+1;
    psearch=p1;
  };
  if(buffer.data()!=p1) {
    std::memmove(buffer.data(),p1,actual_buffer_length);
  }
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
  if(std::isnan(v)) {
    v=std::numeric_limits<double>::quiet_NaN();
  }
  std::snprintf(buffer_string.data(),buffer_string.size(),fieldFmt.c_str(),v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %d",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,long int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %ld",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,long long int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %lld",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %u",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,long unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %lu",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,long long unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %llu",v);
  printField(name,buffer_string.data());
  return *this;
}

OFile& OFile::printField(const std::string&name,const std::string & v) {
  unsigned i;
  for(i=0; i<const_fields.size(); i++)
    if(const_fields[i].name==name) {
      break;
    }
  if(i>=const_fields.size()) {
    Field field;
    field.name=name;
    field.value=v;
    fields.push_back(field);
  } else {
    if(const_fields[i].value!=v) {
      fieldChanged=true;
    }
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
    std::string min, max;
    val->getDomain( min, max );
    printField( "min_" + val->getName(), min );
    printField("max_" + val->getName(), max );
  }
  return *this;
}

OFile& OFile::printField() {
  bool reprint=false;
  if(fieldChanged || fields.size()!=previous_fields.size()) {
    reprint=true;
  } else
    for(unsigned i=0; i<fields.size(); i++) {
      if( previous_fields[i].name!=fields[i].name ||
          (fields[i].constant && fields[i].value!=previous_fields[i].value) ) {
        reprint=true;
        break;
      }
    }
  if(reprint) {
    printf("#! FIELDS");
    for(unsigned i=0; i<fields.size(); i++) {
      printf(" %s",fields[i].name.c_str());
    }
    printf("\n");
    for(unsigned i=0; i<const_fields.size(); i++) {
      printf("#! SET %s %s",const_fields[i].name.c_str(),const_fields[i].value.c_str());
      printf("\n");
    }
  }
  for(unsigned i=0; i<fields.size(); i++) {
    printf("%s",fields[i].value.c_str());
  }
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
  if(str=="/dev/null") {
    return;
  }
  plumed_assert( backstring!="bck" && !checkRestart());
  size_t found=str.find_last_of("/\\");
  std::string filename = appendSuffix(str,getSuffix());
  std::string directory=filename.substr(0,found+1);
  std::string file=filename.substr(found+1);
  if( FileExist(filename) ) {
    backupFile("bck", filename);
  }
  for(int i=0;; i++) {
    std::string num;
    Tools::convert(i,num);
    std::string filestr = directory + backstring + "." + num + "." + file;
    if( !FileExist(filestr) ) {
      break;
    }
    backupFile( "bck", filestr);
  }
}

void OFile::backupFile( const std::string& bstring, const std::string& fname ) {
  if(fname=="/dev/null") {
    return;
  }
  int maxbackup=100;
  if(std::getenv("PLUMED_MAXBACKUP")) {
    Tools::convert(std::getenv("PLUMED_MAXBACKUP"),maxbackup);
  }
  if(maxbackup>0 && (!comm || comm->Get_rank()==0)) {
    FILE* ff=std::fopen(const_cast<char*>(fname.c_str()),"r");
    if(ff) {
      // no exception here
      std::fclose(ff);
      std::string backup;
      size_t found=fname.find_last_of("/\\");
      std::string directory=fname.substr(0,found+1);
      std::string file=fname.substr(found+1);
      for(int i=0;; i++) {
        std::string num;
        Tools::convert(i,num);
        if(i>maxbackup) {
          plumed_merror("cannot backup file "+file+" maximum number of backup is "+num+"\n");
        }
        backup=directory+bstring +"."+num+"."+file;
        FILE* fff=std::fopen(backup.c_str(),"r");
        // no exception here
        if(!fff) {
          break;
        } else {
          std::fclose(fff);
        }
      }
      int check=rename(fname.c_str(),backup.c_str());
      plumed_massert(check==0,"renaming "+fname+" into "+backup+" failed for reason: "+std::strerror(errno));
    }
  }
}

OFile& OFile::open(const std::string&setpath) {
  plumed_assert(!cloned);
  eof=false;
  err=false;
  fp=NULL;
  gzfp=NULL;
  path=setpath;
  path=appendSuffix(setpath,getSuffix());
  // this must be done before backing up the file, otherwise
  // the backup would rename the file that is being written by another OFile
  checkNotOpenElsewhere();
  if(checkRestart()) {
    fp=std::fopen(const_cast<char*>(path.c_str()),"a");
    mode="a";
    if(Tools::extension(path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
      gzfp=(void*)gzopen(const_cast<char*>(path.c_str()),"a9");
#else
      plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
    }
  } else {
    backupFile( backstring, path );
    if(comm) {
      comm->Barrier();
    }
    fp=std::fopen(const_cast<char*>(path.c_str()),"w");
    mode="w";
    if(Tools::extension(path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
      gzfp=(void*)gzopen(const_cast<char*>(path.c_str()),"w9");
#else
      plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
    }
    // All processes open (and truncate) the file: the first process must not
    // write before the others have opened it, or what it wrote could be truncated
    if(comm) {
      comm->Barrier();
    }
  }
  if(plumed) {
    plumed->insertFile(*this);
  }
  return *this;
}

void OFile::checkNotOpenElsewhere()const {
  if(!plumed) {
    return;
  }
  if(path=="/dev/null") {
    return;
  }
  struct stat st;
  const bool exists=(stat(path.c_str(),&st)==0);
// sharing devices, pipes, etc is legitimate
// (S_ISREG is not defined by MSVC, so the file type is tested explicitly)
  if(exists && (st.st_mode & S_IFMT)!=S_IFREG) {
    return;
  }
  for(const auto & f : plumed->getFiles()) {
    const auto other=dynamic_cast<const OFile*>(f);
    // input files are ignored: reading a file while writing it is allowed (e.g. METAD restart)
    if(!other || other==this || !other->fp) {
      continue;
    }
    bool same=(other->path==path);
#ifndef _WIN32
    // On Windows st_ino is always zero, so two different files on the same drive
    // would look the same: there only the paths are compared
    if(!same && exists) {
      // catch the same file accessed with a different path (e.g. ./FILE, absolute paths, or links)
      struct stat ost;
      if(fstat(fileno(other->fp),&ost)==0) {
        same=(ost.st_dev==st.st_dev && ost.st_ino==st.st_ino);
      }
    }
#endif
    if(same) {
      std::string msg="file " + path + " is already open for writing";
      if(other->action) {
        msg+=" by action " + other->action->getLabel() + " (" + other->action->getName() + ")";
      }
      if(action) {
        msg+=" and cannot be opened again by action " + action->getLabel() + " (" + action->getName() + ")";
      }
      msg+=". Writing the same file twice would result in a corrupted file. Please use a different file name";
      plumed_merror(msg);
    }
  }
}

OFile& OFile::rewind() {
// we use here "hard" rewind, which means close/reopen
// the reason is that normal rewind does not work when in append mode
// moreover, we can take a backup of the file
  plumed_assert(fp);
  clearFields();

  if(!comm || comm->Get_rank()==0) {
    std::string fname=this->path;
    size_t found=fname.find_last_of("/\\");
    std::string directory=fname.substr(0,found+1);
    std::string file=fname.substr(found+1);
    std::string backup=directory+backstring +".last."+file;
    int check=rename(fname.c_str(),backup.c_str());
    plumed_massert(check==0,"renaming "+fname+" into "+backup+" failed for reason: "+std::strerror(errno));
  }

  if(comm) {
    comm->Barrier();
  }

  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    gzclose((gzFile)gzfp);
    // no exception here
    gzfp=(void*)gzopen(const_cast<char*>(this->path.c_str()),"w9");
#endif
  } else {
    std::fclose(fp);
    // no exception here
    fp=std::fopen(const_cast<char*>(path.c_str()),"w");
  }
  // see OFile::open
  if(comm) {
    comm->Barrier();
  }
  return *this;
}

FileBase& OFile::flush() {
  if(heavyFlush) {
    if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
      gzclose(gzFile(gzfp));
      // no exception here
      gzfp=(void*)gzopen(const_cast<char*>(path.c_str()),"a");
#endif
    } else {
      std::fclose(fp);
      // no exception here
      fp=std::fopen(const_cast<char*>(path.c_str()),"a");
    }
  } else {
    FileBase::flush();
    // if(gzfp) gzflush(gzFile(gzfp),Z_FINISH);
    // for some reason flushing with Z_FINISH has problems on linux
    // I thus use this (incomplete) flush
#ifdef __PLUMED_HAS_ZLIB
    if(gzfp) {
      gzflush(gzFile(gzfp),Z_FULL_FLUSH);
    }
#endif
  }
  return *this;
}

bool OFile::checkRestart()const {
  if(enforceRestart_) {
    return true;
  } else if(enforceBackup_) {
    return false;
  } else if(action) {
    return action->getRestart();
  } else if(plumed) {
    return plumed->getRestart();
  } else {
    return false;
  }
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

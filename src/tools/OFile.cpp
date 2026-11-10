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
#include "FileTools.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Checkpointing.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>
#include <cctype>
#include <filesystem>
#include <sstream>

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
  if(trimPending_) {
    applyPendingTrim();
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
  if(h5mode) {
    // only blank lines (e.g. those between the rows of a grid) can be printed in a table: they are ignored
    va_list arg;
    va_start(arg, fmt);
    std::vector<char> text(std::vsnprintf(nullptr,0,fmt,arg)+1);
    va_end(arg);
    va_start(arg, fmt);
    std::vsnprintf(text.data(),text.size(),fmt,arg);
    va_end(arg);
    for(const char* c=text.data(); *c; c++) {
      if(!std::isspace(static_cast<unsigned char>(*c))) {
        plumed_merror("file "+getPath()+": HDF5 tables (.h5) can only be written field by field (as PRINT does): this file is written as free text, use another extension");
      }
    }
    return text.size()-1;
  }
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
  // after a checkpoint restart, the value it had in the last header (so that the header is not printed again)
  for(const auto & r : restored_const_fields) {
    if(r.name==name) {
      f.value=r.value;
    }
  }
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
  h5table::Value t;
  t.type=h5table::Value::Double;
  t.d=v;
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %d",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,long int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %ld",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,long long int v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %lld",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %u",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,long unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %lu",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,long long unsigned v) {
  std::snprintf(buffer_string.data(),buffer_string.size()," %llu",v);
  h5table::Value t;
  t.type=h5table::Value::Int;
  t.i=static_cast<long long>(v);
  setField(name,buffer_string.data(),t);
  return *this;
}

OFile& OFile::printField(const std::string&name,const std::string & v) {
  h5table::Value t;
  t.type=h5table::Value::String;
  t.s=v;
  setField(name,v,t);
  return *this;
}

void OFile::setField(const std::string&name,const std::string& v,const h5table::Value& typed) {
  unsigned i;
  for(i=0; i<const_fields.size(); i++)
    if(const_fields[i].name==name) {
      break;
    }
  if(i>=const_fields.size()) {
    Field field;
    field.name=name;
    field.value=v;
    field.typed=typed;
    fields.push_back(field);
  } else {
    if(const_fields[i].value!=v) {
      fieldChanged=true;
    }
    const_fields[i].value=v;
    const_fields[i].typed=typed;
  }
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

std::string OFile::getLastFields() const {
  std::string str;
  // name, constant flag, value (used only for constant fields), one field per line
  for(const auto & f : previous_fields) {
    str+=f.name+"\t0\t\n";
  }
  for(const auto & f : const_fields) {
    str+=f.name+"\t1\t"+f.value+"\n";
  }
  return str;
}

void OFile::setLastFields(const std::string& str) {
  previous_fields.clear();
  restored_const_fields.clear();
  std::istringstream is(str);
  std::string line;
  while(std::getline(is,line)) {
    const size_t t1=line.find('\t');
    const size_t t2=line.find('\t',t1+1);
    plumed_massert(t1!=std::string::npos && t2!=std::string::npos,"file "+getPath()+": cannot parse the field header stored in the checkpoint");
    Field f;
    f.name=line.substr(0,t1);
    f.value=line.substr(t2+1);
    if(line.substr(t1+1,t2-t1-1)=="1") {
      restored_const_fields.push_back(f);
    } else {
      previous_fields.push_back(f);
    }
  }
  // constant fields already added (when the header is restored at the first write, after
  // a trim at restart): they take the restored value if they have none yet, and the header
  // is printed again if a value differs or a field was not in the restored header
  bool changed=false;
  for(auto & c : const_fields) {
    bool found=false;
    for(const auto & r : restored_const_fields) {
      if(r.name==c.name) {
        found=true;
        if(c.value.empty()) {
          c.value=r.value;
        } else if(c.value!=r.value) {
          changed=true;
        }
      }
    }
    if(!found) {
      changed=true;
    }
  }
  // (with no constant field yet, at a checkpoint restart, they are added later by addConstantField)
  if(!const_fields.empty() && const_fields.size()!=restored_const_fields.size()) {
    changed=true;
  }
  fieldChanged=changed;
}

OFile& OFile::printField() {
  if(h5mode) {
    // HDF5 table: the writer starts a new block when the fields change
    if(isWriter()) {
      h5open();
      std::vector<h5table::Field> row,constants;
      for(const auto & f : fields) {
        row.push_back(h5table::Field{f.name,f.typed});
      }
      for(const auto & f : const_fields) {
        constants.push_back(h5table::Field{f.name,f.typed});
      }
      h5writer->addRow(row,constants);
    }
    previous_fields=fields;
    fields.clear();
    fieldChanged=false;
    return *this;
  }
  // before deciding whether the header is printed (see applyPendingTrim)
  if(trimPending_) {
    applyPendingTrim();
  }
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
    FILE* ff=fopenUTF8(fname,"r");
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
        FILE* fff=fopenUTF8(backup,"r");
        // no exception here
        if(!fff) {
          break;
        } else {
          std::fclose(fff);
        }
      }
      int check=renameUTF8(fname,backup);
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
  // the null device must not be backed up (renamed) nor get a suffix
  path=isWindowsNullDevice(setpath) ? std::string("/dev/null") : setpath;
  path=appendSuffix(path,getSuffix());
  // this must be done before backing up the file, otherwise
  // the backup would rename the file that is being written by another OFile
  checkNotOpenElsewhere();
  // Only the writer opens the file (the other processes open /dev/null):
  // with many processes on a parallel file system, opening (and truncating)
  // the same file from all of them is slow.
  if(h5table::isH5(path)) {
    // HDF5 table: fp is the null device on all the processes, the table is written by h5writer
    // (opened with "r+": the null device exists, nothing is created)
    h5mode=true;
    const bool restart=checkRestart();
    h5append=restart;
    if(isWriter()) {
      if(restart) {
        // trimmed at the first row (see h5open)
        trimPending_=plumed && !enforceNoTrim_ && plumed->checkpointing.trimsAtRestart() && !(action && action->getNoTrim());
      } else {
        backupFile( backstring, path );
      }
    }
    fp=std::fopen(nullDevice,"r+");
    mode=restart ? "a" : "w";
    if(comm) {
      comm->Barrier();
    }
  } else if(checkRestart()) {
    openOnWriter("a","a9");
    mode="a";
    // trimmed at the first write (see trimPending_)
    trimPending_=plumed && isWriter() && !enforceNoTrim_ && plumed->checkpointing.trimsAtRestart() && !(action && action->getNoTrim());
  } else {
    backupFile( backstring, path );
    if(comm) {
      comm->Barrier();
    }
    openOnWriter("w","w9");
    mode="w";
  }
  if(!h5mode) {
    syncOpen();
  }
  if(plumed) {
    plumed->insertFile(*this);
    plumed->checkpointing.registerOutputFile(*this);
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
  // The files are compared on the writer only: the other processes have the null device
  // as file pointer (see openOnWriter), and the files could be different if paths
  // are local to nodes. The verdict is then broadcast, so that all processes
  // stop with the same error (if only the writer stopped, the others would hang).
  const auto findDuplicate=[&]() {
    std::string msg;
    struct stat st;
    const bool exists=(stat(path.c_str(),&st)==0);
    // sharing devices, pipes, etc is legitimate
    // (S_ISREG is not defined by MSVC, so the file type is tested explicitly)
    if(exists && (st.st_mode & S_IFMT)!=S_IFREG) {
      return std::string();
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
        msg="file " + path + " is already open for writing";
        if(other->action) {
          msg+=" by action " + other->action->getLabel() + " (" + other->action->getName() + ")";
        }
        if(action) {
          msg+=" and cannot be opened again by action " + action->getLabel() + " (" + action->getName() + ")";
        }
        msg+=". Writing the same file twice would result in a corrupted file. Please use a different file name";
        break;
      }
    }
    return msg;
  };
  std::string msg;
  if(isWriter()) {
    msg=findDuplicate();
  }
  if(comm) {
    unsigned len=msg.length();
    comm->Bcast(len,0);
    std::vector<char> chars(msg.begin(),msg.end());
    chars.resize(len);
    if(len>0) {
      comm->Bcast(chars,0);
    }
    msg.assign(chars.begin(),chars.end());
  }
  if(!msg.empty()) {
    plumed_merror(msg);
  }
}

void OFile::applyPendingTrim() {
  if(!trimPending_) {
    return;
  }
  trimPending_=false;
  // the file is open in append mode and nothing has been written yet:
  // the next write goes to the new end of the file
  plumed->checkpointing.trimAtRestart(path);
  // the header of the last line kept is not printed again if the fields are the same
  plumed->checkpointing.restoreLastFields(*this);
}

OFile::~OFile() {
  h5writer.reset();
}

void OFile::closeBackend() {
  h5writer.reset();
  h5mode=false;
}

void OFile::h5open() {
  if(h5writer) {
    return;
  }
  if(trimPending_) {
    applyPendingTrim();
  }
  h5writer=std::make_unique<h5table::Writer>(path,h5append);
  h5append=true;
}

std::vector<unsigned long long> OFile::getH5Rows() {
  if(h5writer) {
    return h5writer->rows();
  }
  // nothing written yet by this run: the table as it is
  std::error_code ec;
  if(isWriter() && h5append && std::filesystem::exists(pathFromUTF8(path),ec)) {
    return h5table::rows(path);
  }
  return {};
}

bool OFile::isWriter()const {
  return !comm || comm->Get_rank()==0;
}

void OFile::openOnWriter(const char* fmode,const char* gzmode,bool reopenFp) {
  // "/dev/null" is a valid name on all platforms (on Windows fopen would fail,
  // and a NULL fp crashed at the first write or flush: WALKERS_MPI hills files)
  const std::string p=(isWriter() && path!="/dev/null") ? path : std::string(nullDevice);
  if(reopenFp) {
#ifdef _WIN32
    // binary mode: write LF line endings as on the other platforms
    // (text mode would turn every \n into \r\n)
    fp=fopenUTF8(p,(std::string(fmode)+"b").c_str());
#else
    fp=fopenUTF8(p,fmode);
#endif
  }
  if(Tools::extension(path)=="gz") {
#ifdef __PLUMED_HAS_ZLIB
    gzfp=gzopenUTF8(p,gzmode);
#else
    (void) gzmode;
    plumed_merror("file " + getPath() + ": trying to use a gz file without zlib being linked");
#endif
  }
}

void OFile::syncOpen() {
  if(!comm) {
    return;
  }
  int ok=(fp!=NULL);
  comm->Bcast(ok,0);
  if(!ok && fp) {
    // the writer could not open the file: behave as if no process could
    std::fclose(fp);
    fp=NULL;
#ifdef __PLUMED_HAS_ZLIB
    if(gzfp) {
      gzclose(gzFile(gzfp));
    }
#endif
    gzfp=NULL;
  }
}

OFile& OFile::rewind() {
// we use here "hard" rewind, which means close/reopen
// the reason is that normal rewind does not work when in append mode
// moreover, we can take a backup of the file
  plumed_assert(fp);
  clearFields();
  // a snapshot: it is rewritten, not trimmed
  trimPending_=false;
  if(plumed) {
    plumed->checkpointing.markSnapshotFile(path);
  }
  if(h5mode) {
    // HDF5 table: the previous one is kept as bck.last.name, a new one is created at the next row
    if(isWriter()) {
      h5writer.reset();
      std::string fname=this->path;
      size_t found=fname.find_last_of("/\\");
      std::string backup=fname.substr(0,found+1)+backstring+".last."+fname.substr(found+1);
      std::error_code ec;
      if(std::filesystem::exists(pathFromUTF8(fname),ec)) {
        std::filesystem::rename(pathFromUTF8(fname),pathFromUTF8(backup),ec);
        plumed_massert(!ec,"renaming "+fname+" into "+backup+" failed: "+ec.message());
      }
      h5append=false;
    }
    if(comm) {
      comm->Barrier();
    }
    return *this;
  }

  // close before renaming: on Windows an open file cannot be renamed
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    gzclose((gzFile)gzfp);
#endif
  } else {
    std::fclose(fp);
  }

  if(!comm || comm->Get_rank()==0) {
    std::string fname=this->path;
    size_t found=fname.find_last_of("/\\");
    std::string directory=fname.substr(0,found+1);
    std::string file=fname.substr(found+1);
    std::string backup=directory+backstring +".last."+file;
    // std::filesystem::rename replaces an existing backup also on Windows,
    // where std::rename fails if the target exists
    std::error_code ec;
    std::filesystem::rename(pathFromUTF8(fname),pathFromUTF8(backup),ec);
    plumed_massert(!ec,"renaming "+fname+" into "+backup+" failed: "+ec.message());
  }

  if(comm) {
    comm->Barrier();
  }

  // no exception here
  if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
    openOnWriter("w","w9",false);
#endif
  } else {
    openOnWriter("w","w9");
  }
  // see OFile::open
  syncOpen();
  return *this;
}

FileBase& OFile::flush() {
  if(h5mode) {
    if(h5writer) {
      h5writer->flush();
    }
    return *this;
  }
  if(heavyFlush) {
    if(gzfp) {
#ifdef __PLUMED_HAS_ZLIB
      gzclose(gzFile(gzfp));
      // no exception here
      openOnWriter("a","a",false);
#endif
    } else {
      std::fclose(fp);
      // no exception here
      openOnWriter("a","a");
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

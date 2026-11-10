/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#include "Checkpoint.h"
#include "Exception.h"
#include "FileTools.h"
#include <cerrno>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#ifndef _WIN32
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#endif

/*
File layout (all integers little-endian, independent of the host):

  magic            8 bytes   "PLUMEDCP"
  formatVersion    uint32
  byteOrderMark    uint32    0x01020304
  plumedVersion    string    (uint32 length + bytes)
  entries          ...
  endMarker        uint8     0
  checksum         uint32    CRC32 of everything before it

Every entry starts with a type byte and a name:
  1  begin section   name
  2  end section     (no name)
  3  integer array   name, uint8 ndims, uint64 dims[ndims], int64 data[prod(dims)]
  4  double array    name, uint8 ndims, uint64 dims[ndims], float64 data[prod(dims)]
  5  string          name, uint64 length, bytes
Scalars are arrays with ndims=0 (one element). Entry names are stored relative
to their section; the reader keys them by the full path "section/subsection/name".
*/

namespace PLMD {

namespace {

const char magic[9]="PLUMEDCP";
const uint32_t formatVersionCurrent=1;
const uint32_t byteOrderMark=0x01020304;

enum : uint8_t { tEnd=0, tBegin=1, tEndSection=2, tInt=3, tDouble=4, tString=5 };

void put8(std::string& b,uint8_t v) {
  b.push_back(static_cast<char>(v));
}

void put32(std::string& b,uint32_t v) {
  for(unsigned i=0; i<4; i++) {
    b.push_back(static_cast<char>((v>>(8*i))&0xff));
  }
}

void put64(std::string& b,uint64_t v) {
  for(unsigned i=0; i<8; i++) {
    b.push_back(static_cast<char>((v>>(8*i))&0xff));
  }
}

void putDouble(std::string& b,double v) {
  static_assert(sizeof(double)==8,"double must be 8 bytes");
  uint64_t u;
  std::memcpy(&u,&v,8);
  put64(b,u);
}

void putString(std::string& b,const std::string& s) {
  put32(b,static_cast<uint32_t>(s.size()));
  b.append(s);
}

uint32_t crc32(const std::string& data) {
  static uint32_t table[256];
  static bool init=false;
  if(!init) {
    for(uint32_t i=0; i<256; i++) {
      uint32_t c=i;
      for(unsigned k=0; k<8; k++) {
        c=(c&1) ? (0xEDB88320u^(c>>1)) : (c>>1);
      }
      table[i]=c;
    }
    init=true;
  }
  uint32_t c=0xFFFFFFFFu;
  for(unsigned char ch : data) {
    c=table[(c^ch)&0xff]^(c>>8);
  }
  return c^0xFFFFFFFFu;
}

// Sequential parser with bounds checking
class Parser {
  const std::string& b;
  size_t pos=0;
  void need(size_t n) const {
    if(pos+n>b.size()) {
      plumed_merror("checkpoint: unexpected end of data (truncated file?)");
    }
  }
public:
  explicit Parser(const std::string& bb): b(bb) {}
  size_t position() const {
    return pos;
  }
  uint8_t get8() {
    need(1);
    return static_cast<uint8_t>(b[pos++]);
  }
  uint32_t get32() {
    need(4);
    uint32_t v=0;
    for(unsigned i=0; i<4; i++) {
      v|=static_cast<uint32_t>(static_cast<uint8_t>(b[pos++]))<<(8*i);
    }
    return v;
  }
  uint64_t get64() {
    need(8);
    uint64_t v=0;
    for(unsigned i=0; i<8; i++) {
      v|=static_cast<uint64_t>(static_cast<uint8_t>(b[pos++]))<<(8*i);
    }
    return v;
  }
  double getDouble() {
    uint64_t u=get64();
    double v;
    std::memcpy(&v,&u,8);
    return v;
  }
  std::string getString32() {
    uint32_t n=get32();
    need(n);
    std::string s=b.substr(pos,n);
    pos+=n;
    return s;
  }
  std::string getString64() {
    uint64_t n=get64();
    need(n);
    std::string s=b.substr(pos,n);
    pos+=n;
    return s;
  }
};

}

// CheckpointWriter

CheckpointWriter::CheckpointWriter(const std::string& plumedVersion) {
  buffer.append(magic,8);
  put32(buffer,formatVersionCurrent);
  put32(buffer,byteOrderMark);
  putString(buffer,plumedVersion);
  entriesInSection.push_back(0);
}

std::string CheckpointWriter::currentPath() const {
  std::string p;
  for(const auto & s : sections) {
    p+=s+"/";
  }
  return p;
}

void CheckpointWriter::putName(const std::string& name) {
  plumed_massert(!name.empty() && name.find('/')==std::string::npos,"checkpoint: invalid entry name '"+name+"'");
  putString(buffer,name);
}

void CheckpointWriter::putShape(const std::vector<unsigned long long>& dims) {
  plumed_massert(dims.size()<256,"checkpoint: too many dimensions");
  put8(buffer,static_cast<uint8_t>(dims.size()));
  for(auto d : dims) {
    put64(buffer,d);
  }
}

void CheckpointWriter::beginSection(const std::string& name) {
  put8(buffer,tBegin);
  putName(name);
  entriesInSection.back()++;
  sections.push_back(name);
  entriesInSection.push_back(0);
}

void CheckpointWriter::endSection() {
  plumed_massert(!sections.empty(),"checkpoint: endSection without beginSection");
  put8(buffer,tEndSection);
  sections.pop_back();
  entriesInSection.pop_back();
}

unsigned CheckpointWriter::getNumberOfEntriesInSection() const {
  return entriesInSection.back();
}

void CheckpointWriter::write(const std::string& name,double value) {
  writeDoubleArray(name,&value, {});
}

void CheckpointWriter::write(const std::string& name,long long value) {
  writeIntArray(name,&value, {});
}

void CheckpointWriter::write(const std::string& name,const std::string& value) {
  put8(buffer,tString);
  putName(name);
  put64(buffer,value.size());
  buffer.append(value);
  entriesInSection.back()++;
}

void CheckpointWriter::writeIntArray(const std::string& name,const long long* data,const std::vector<unsigned long long>& dims) {
  put8(buffer,tInt);
  putName(name);
  putShape(dims);
  unsigned long long n=1;
  for(auto d : dims) {
    n*=d;
  }
  for(unsigned long long i=0; i<n; i++) {
    put64(buffer,static_cast<uint64_t>(data[i]));
  }
  entriesInSection.back()++;
}

void CheckpointWriter::writeDoubleArray(const std::string& name,const double* data,const std::vector<unsigned long long>& dims) {
  put8(buffer,tDouble);
  putName(name);
  putShape(dims);
  unsigned long long n=1;
  for(auto d : dims) {
    n*=d;
  }
  for(unsigned long long i=0; i<n; i++) {
    putDouble(buffer,data[i]);
  }
  entriesInSection.back()++;
}

std::string CheckpointWriter::finalize() const {
  plumed_massert(sections.empty(),"checkpoint: section '"+sections.back()+"' was not closed");
  std::string content=buffer;
  put8(content,tEnd);
  put32(content,crc32(content));
  return content;
}

void CheckpointWriter::writeFile(const std::string& path,const std::string& prevPath) const {
  const std::string content=finalize();
  const std::string tmp=path+".tmp";
  {
#ifndef _WIN32
    // explicit permissions (rw-r--r--, further restricted by the umask)
    std::FILE* fp=nullptr;
    const int fd=open(tmp.c_str(),O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if(fd>=0) {
      fp=fdopen(fd,"wb");
      if(!fp) {
        close(fd);
      }
    }
#else
    std::FILE* fp=fopenUTF8(tmp,"wb");
#endif
    if(!fp) {
      plumed_merror("checkpoint: cannot open "+tmp+" for writing: "+std::strerror(errno));
    }
    const size_t written=std::fwrite(content.data(),1,content.size(),fp);
    bool ok=(written==content.size()) && std::fflush(fp)==0;
#ifndef _WIN32
    if(ok) {
      ok=(fsync(fileno(fp))==0);
    }
#endif
    ok=(std::fclose(fp)==0) && ok;
    if(!ok) {
      plumed_merror("checkpoint: error writing "+tmp+": "+std::strerror(errno));
    }
  }
  std::error_code ec;
  if(!prevPath.empty() && std::filesystem::exists(pathFromUTF8(path),ec)) {
    std::filesystem::rename(pathFromUTF8(path),pathFromUTF8(prevPath),ec);
    if(ec) {
      plumed_merror("checkpoint: cannot rename "+path+" to "+prevPath+": "+ec.message());
    }
  }
  std::filesystem::rename(pathFromUTF8(tmp),pathFromUTF8(path),ec);
  if(ec) {
    plumed_merror("checkpoint: cannot rename "+tmp+" to "+path+": "+ec.message());
  }
}

// CheckpointReader

unsigned CheckpointReader::checksum(const std::string& data) {
  return crc32(data);
}

std::string CheckpointReader::readFile(const std::string& path) {
  std::ifstream in(pathFromUTF8(path),std::ios::binary);
  if(!in) {
    plumed_merror("checkpoint: cannot open "+path+" for reading");
  }
  std::ostringstream ss;
  ss<<in.rdbuf();
  return ss.str();
}

namespace {

// Check the magic, the format version (before the checksum, so that a newer
// file gets a clear message) and the CRC32 trailer of a checkpoint file.
void checkHeaderAndChecksum(const std::string& content) {
  if(content.size()<16 || std::memcmp(content.data(),magic,8)!=0) {
    plumed_merror("checkpoint: not a PLUMED checkpoint file");
  }
  Parser ph(content);
  for(unsigned i=0; i<8; i++) {
    ph.get8();
  }
  const uint32_t v=ph.get32();
  if(v>formatVersionCurrent) {
    plumed_merror("checkpoint: file format version "+std::to_string(v)+" is newer than the one supported by this PLUMED version ("+std::to_string(formatVersionCurrent)+")");
  }
  // the checksum covers everything but the last 4 bytes
  const std::string body=content.substr(0,content.size()-4);
  const std::string stored=content.substr(content.size()-4);
  Parser pt(stored);
  if(pt.get32()!=crc32(body)) {
    plumed_merror("checkpoint: checksum mismatch (corrupted or truncated file)");
  }
}

// Parse the payload of a data entry (type byte and name already consumed):
// a string, or a shaped array of int64 or float64 (ndims=0 for a scalar).
CheckpointReader::Entry parseEntry(Parser& p,uint8_t type) {
  CheckpointReader::Entry e;
  if(type==tString) {
    e.type=CheckpointReader::Type::String;
    e.str=p.getString64();
    return e;
  }
  e.type=(type==tInt) ? CheckpointReader::Type::Int : CheckpointReader::Type::Double;
  const unsigned ndims=p.get8();
  unsigned long long n=1;
  for(unsigned i=0; i<ndims; i++) {
    e.dims.push_back(p.get64());
    n*=e.dims.back();
  }
  if(type==tInt) {
    e.ints.reserve(n);
    for(unsigned long long i=0; i<n; i++) {
      e.ints.push_back(static_cast<long long>(p.get64()));
    }
  } else {
    e.doubles.reserve(n);
    for(unsigned long long i=0; i<n; i++) {
      e.doubles.push_back(p.getDouble());
    }
  }
  return e;
}

}

CheckpointReader::CheckpointReader(const std::string& content) {
  checkHeaderAndChecksum(content);
  // header: magic, format version, byte order mark, PLUMED version
  Parser p(content);
  for(unsigned i=0; i<8; i++) {
    p.get8();
  }
  formatVersion=p.get32();
  if(p.get32()!=byteOrderMark) {
    plumed_merror("checkpoint: unexpected byte order mark");
  }
  plumedVersion=p.getString32();
  // entries, keyed by their full path "section/subsection/name"
  std::vector<std::string> stack;
  auto path=[&stack]() {
    std::string s;
    for(const auto & x : stack) {
      s+=x+"/";
    }
    return s;
  };
  bool finished=false;
  while(!finished) {
    const uint8_t type=p.get8();
    switch(type) {
    case tEnd:
      finished=true;
      break;
    case tBegin:
      stack.push_back(p.getString32());
      sectionList.push_back(path());
      break;
    case tEndSection:
      if(stack.empty()) {
        plumed_merror("checkpoint: unbalanced sections");
      }
      stack.pop_back();
      break;
    case tInt:
    case tDouble:
    case tString: {
      const std::string full=path()+p.getString32();
      if(entries.count(full)) {
        plumed_merror("checkpoint: duplicate entry '"+full+"'");
      }
      entries[full]=parseEntry(p,type);
      order.push_back(full);
    }
    break;
    default:
      plumed_merror("checkpoint: unknown entry type "+std::to_string(type));
    }
  }
  if(!stack.empty()) {
    plumed_merror("checkpoint: section '"+stack.back()+"' is not closed");
  }
  // nothing after the end marker but the checksum
  if(p.position()+4!=content.size()) {
    plumed_merror("checkpoint: trailing data after the end marker");
  }
}

std::string CheckpointReader::currentPath() const {
  std::string p;
  for(const auto & s : sections) {
    p+=s+"/";
  }
  return p;
}

void CheckpointReader::beginSection(const std::string& name) {
  plumed_massert(hasSection(name),"checkpoint: section '"+currentPath()+name+"' not found");
  sections.push_back(name);
}

void CheckpointReader::endSection() {
  plumed_massert(!sections.empty(),"checkpoint: endSection without beginSection");
  sections.pop_back();
}

bool CheckpointReader::hasSection(const std::string& name) const {
  const std::string full=currentPath()+name+"/";
  for(const auto & s : sectionList) {
    if(s==full) {
      return true;
    }
  }
  return false;
}

bool CheckpointReader::has(const std::string& name) const {
  return entries.count(currentPath()+name)>0;
}

const CheckpointReader::Entry& CheckpointReader::find(const std::string& name,Type type,bool scalar) const {
  const std::string full=currentPath()+name;
  const auto it=entries.find(full);
  if(it==entries.end()) {
    plumed_merror("checkpoint: entry '"+full+"' not found");
  }
  const Entry& e=it->second;
  if(e.type!=type) {
    plumed_merror("checkpoint: entry '"+full+"' has an unexpected type");
  }
  if(scalar && !e.scalar()) {
    plumed_merror("checkpoint: entry '"+full+"' is an array, not a scalar");
  }
  return e;
}

const CheckpointReader::Entry& CheckpointReader::getEntry(const std::string& fullName) const {
  const auto it=entries.find(fullName);
  if(it==entries.end()) {
    plumed_merror("checkpoint: entry '"+fullName+"' not found");
  }
  return it->second;
}

void CheckpointReader::read(const std::string& name,double& value) const {
  value=find(name,Type::Double,true).doubles[0];
}

void CheckpointReader::read(const std::string& name,long long& value) const {
  value=find(name,Type::Int,true).ints[0];
}

void CheckpointReader::read(const std::string& name,std::string& value) const {
  value=find(name,Type::String,true).str;
}

std::vector<unsigned long long> CheckpointReader::getShape(const std::string& name) const {
  const std::string full=currentPath()+name;
  return getEntry(full).dims;
}

}

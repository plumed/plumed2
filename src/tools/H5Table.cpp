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
#include "H5Table.h"
#include "Exception.h"
#include "FileTools.h"
#include "Tools.h"
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <limits>
#include <mutex>
#include <utility>
#ifdef _WIN32
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef __PLUMED_HAS_LIBHDF5
#include <hdf5.h>
#endif

namespace PLMD {
namespace h5table {

namespace {

constexpr unsigned long long tarRecord=512;

unsigned tarChecksum(const char* h) {
  unsigned sum=0;
  for(unsigned k=0; k<tarRecord; k++) {
    sum+= (k>=148 && k<156) ? static_cast<unsigned>(' ') : static_cast<unsigned char>(h[k]);
  }
  return sum;
}

bool getOctal(const char* p,unsigned n,unsigned long long& v) {
  v=0;
  unsigned k=0;
  while(k<n && p[k]==' ') {
    k++;
  }
  const unsigned first=k;
  for(; k<n && p[k]>='0' && p[k]<='7'; k++) {
    v=v*8+(p[k]-'0');
  }
  if(k==first) {
    return false;
  }
  for(; k<n; k++) {
    if(p[k]!=' ' && p[k]!='\0') {
      return false;
    }
  }
  return true;
}

/// a valid ustar header, and the size of its entry
/// (POSIX ustar, pax and GNU formats; the size in octal, or in base 256 as GNU tar writes large sizes)
bool validTarHeader(const char* h,unsigned long long& size) {
  unsigned long long sum=0;
  if(std::memcmp(h+257,"ustar",5)!=0 || !getOctal(h+148,8,sum) || sum!=tarChecksum(h)) {
    return false;
  }
  if(static_cast<unsigned char>(h[124])&0x80) {
    size=0;
    for(unsigned k=125; k<136; k++) {
      if(size>>56) {
        return false;
      }
      size=(size<<8)|static_cast<unsigned char>(h[k]);
    }
    return (static_cast<unsigned char>(h[124])&0x7f)==0;
  }
  return getOctal(h+124,12,size);
}

bool hasArchiveExtension(const std::string& path) {
  return Tools::extension(path)=="hpa" ||
         (path.size()>7 && path.compare(path.size()-7,7,".h5.tar")==0);
}

enum class Content {Unknown,HDF5,Archive};

/// what an existing file contains, from its first bytes
Content sniff(const std::string& path) {
  // (opened with fopenUTF8, as all the files that PLUMED opens)
  unique_FILE in(fopenUTF8(path,"rb"));
  if(!in) {
    return Content::Unknown;
  }
  char h[tarRecord]= {};
  const auto n=std::fread(h,1,tarRecord,in.get());
  static const char signature[8]= {'\x89','H','D','F','\r','\n','\x1a','\n'};
  if(n>=8 && std::memcmp(h,signature,8)==0) {
    return Content::HDF5;
  }
  unsigned long long size=0;
  if(n==tarRecord && validTarHeader(h,size)) {
    return Content::Archive;
  }
  return Content::Unknown;
}

}

bool isArchive(const std::string& path) {
  // by content: HDF5 would read an archive as its first member (it looks for its
  // signature also at byte 512), so an archive is recognized before an HDF5 file
  switch(sniff(path)) {
  case Content::Archive:
    return true;
  case Content::HDF5:
    return false;
  default:
    return hasArchiveExtension(path);
  }
}

bool isArchiveName(const std::string& path) {
  return hasArchiveExtension(path);
}

bool isH5(const std::string& path) {
  return Tools::extension(path)=="h5" || hasArchiveExtension(path);
}

bool isTable(const std::string& path) {
  return isH5(path) || sniff(path)!=Content::Unknown;
}

std::string Value::text() const {
  if(type==String) {
    return s;
  }
  if(type==Int) {
    return std::to_string(i);
  }
  if(std::isnan(d)) {
    return "nan";
  }
  if(std::isinf(d)) {
    return d>0 ? "inf" : "-inf";
  }
  char buf[64];
  std::snprintf(buf,sizeof(buf),"%.17g",d);
  return buf;
}

bool Value::operator==(const Value& v) const {
  if(type!=v.type) {
    return false;
  }
  if(type==String) {
    return s==v.s;
  }
  if(type==Int) {
    return i==v.i;
  }
  return d==v.d || (std::isnan(d) && std::isnan(v.d));
}

#ifndef __PLUMED_HAS_LIBHDF5

bool available() {
  return false;
}

namespace {
[[noreturn]] void noHDF5(const std::string& path) {
  plumed_merror("file "+path+": PLUMED was compiled without libhdf5, files with extension .h5, .hpa or .h5.tar cannot be used (configure PLUMED with --enable-libhdf5)");
}
}

struct Writer::Impl {};
Writer::Writer(const std::string& path,bool,std::size_t) {
  noHDF5(path);
}
Writer::~Writer()=default;
void Writer::addRow(const std::vector<Field>&,const std::vector<Field>&) {}
void Writer::flush() {}

struct Reader::Impl {};
Reader::Reader(const std::string& path) {
  noHDF5(path);
}
Reader::~Reader()=default;
bool Reader::getline(std::string&) {
  return false;
}

void copy(const std::string& input,const std::string&,std::size_t) {
  noHDF5(input);
}

unsigned long long trimArchive(const std::string& path,double,const std::string&,bool&,bool&) {
  noHDF5(path);
}

#else

bool available() {
  return true;
}

namespace {

/// The HDF5 library is thread safe only if compiled to be so (which is not
/// possible for static libraries on Windows, as in the Python wheels), and the
/// ranks of the embedded thread-MPI are threads of the same process: all the
/// calls to HDF5 made here are serialized with this mutex.
std::recursive_mutex& hdf5Mutex() {
  static std::recursive_mutex m;
  return m;
}

using Lock=std::lock_guard<std::recursive_mutex>;

/// owns an HDF5 identifier
class Id {
  hid_t id_=-1;
  herr_t (*close_)(hid_t)=nullptr;
public:
  Id()=default;
  Id(hid_t id,herr_t (*close)(hid_t)):
    id_(id),close_(close) {}
  ~Id() {
    reset();
  }
  Id(const Id&)=delete;
  Id& operator=(const Id&)=delete;
  Id(Id&& o) noexcept:
    id_(o.id_),close_(o.close_) {
    o.id_=-1;
  }
  Id& operator=(Id&& o) noexcept {
    if(this!=&o) {
      reset();
      id_=o.id_;
      close_=o.close_;
      o.id_=-1;
    }
    return *this;
  }
  void reset() {
    if(id_>=0 && close_) {
      Lock lock(hdf5Mutex());
      close_(id_);
    }
    id_=-1;
  }
  operator hid_t() const {
    return id_;
  }
};

void check(bool ok,const std::string& path,const std::string& what) {
  if(!ok) {
    plumed_merror("HDF5 error on file "+path+": "+what);
  }
}

std::string blockName(unsigned b) {
  return "block_"+std::to_string(b);
}

Id stringType() {
  Id t(H5Tcopy(H5T_C_S1),H5Tclose);
  H5Tset_size(t,H5T_VARIABLE);
  return t;
}

void writeAttribute(hid_t obj,const std::string& name,const Value& v,const std::string& path) {
  if(H5Aexists(obj,name.c_str())>0) {
    H5Adelete(obj,name.c_str());
  }
  Id space(H5Screate(H5S_SCALAR),H5Sclose);
  herr_t err=0;
  if(v.type==Value::Double) {
    Id a(H5Acreate2(obj,name.c_str(),H5T_NATIVE_DOUBLE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(a>=0,path,"cannot create attribute "+name);
    err=H5Awrite(a,H5T_NATIVE_DOUBLE,&v.d);
  } else if(v.type==Value::Int) {
    Id a(H5Acreate2(obj,name.c_str(),H5T_NATIVE_LLONG,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(a>=0,path,"cannot create attribute "+name);
    err=H5Awrite(a,H5T_NATIVE_LLONG,&v.i);
  } else {
    Id t=stringType();
    Id a(H5Acreate2(obj,name.c_str(),t,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(a>=0,path,"cannot create attribute "+name);
    const char* c=v.s.c_str();
    err=H5Awrite(a,t,&c);
  }
  check(err>=0,path,"cannot write attribute "+name);
}

Value readAttribute(hid_t obj,const std::string& name,const std::string& path) {
  Id a(H5Aopen(obj,name.c_str(),H5P_DEFAULT),H5Aclose);
  check(a>=0,path,"cannot open attribute "+name);
  Id t(H5Aget_type(a),H5Tclose);
  Value v;
  const H5T_class_t cl=H5Tget_class(t);
  herr_t err=0;
  if(cl==H5T_FLOAT) {
    v.type=Value::Double;
    err=H5Aread(a,H5T_NATIVE_DOUBLE,&v.d);
  } else if(cl==H5T_INTEGER) {
    v.type=Value::Int;
    err=H5Aread(a,H5T_NATIVE_LLONG,&v.i);
  } else if(cl==H5T_STRING) {
    v.type=Value::String;
    if(H5Tis_variable_str(t)>0) {
      Id mt=stringType();
      char* c=nullptr;
      err=H5Aread(a,mt,&c);
      if(c) {
        v.s=c;
        H5free_memory(c);
      }
    } else {
      const std::size_t n=H5Tget_size(t);
      std::string buf(n,'\0');
      err=H5Aread(a,t,buf.data());
      // a fixed-length string is padded with null characters
      const std::size_t end=buf.find('\0');
      if(end!=std::string::npos) {
        buf.resize(end);
      }
      v.s=std::move(buf);
    }
  } else {
    check(false,path,"attribute "+name+" has an unsupported type");
  }
  check(err>=0,path,"cannot read attribute "+name);
  return v;
}

/// a column of a block, with the rows that are buffered or read
struct Column {
  std::string name;
  Value::Type type=Value::Double;
  Id dataset;
  std::vector<double> d;
  std::vector<long long> i;
  std::vector<std::string> s;
  std::size_t size() const {
    return type==Value::Double ? d.size() : type==Value::Int ? i.size() : s.size();
  }
  void clear() {
    d.clear();
    i.clear();
    s.clear();
  }
  void push(const Value& v) {
    if(type==Value::Double) {
      d.push_back(v.d);
    } else if(type==Value::Int) {
      i.push_back(v.i);
    } else {
      s.push_back(v.s);
    }
  }
  Value get(std::size_t k) const {
    Value v;
    v.type=type;
    if(type==Value::Double) {
      v.d=d[k];
    } else if(type==Value::Int) {
      v.i=i[k];
    } else {
      v.s=s[k];
    }
    return v;
  }
};

hid_t memType(Value::Type t) {
  return t==Value::Double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_LLONG;
}

unsigned long long datasetRows(hid_t ds) {
  Id space(H5Dget_space(ds),H5Sclose);
  hsize_t n=0;
  H5Sget_simple_extent_dims(space,&n,nullptr);
  return n;
}

/// create the (empty) dataset of a column in a block; with rows>0, a dataset that will have that
/// number of rows (a member of an archive): if they are few, contiguous (smaller than chunked)
void createDataset(hid_t group,Column& c,const std::string& path,hsize_t rows=0) {
  const hsize_t zero=0,unlimited=H5S_UNLIMITED,chunk=1024;
  Id space(rows>0 ? H5Screate_simple(1,&rows,&rows) : H5Screate_simple(1,&zero,&unlimited),H5Sclose);
  Id plist(H5Pcreate(H5P_DATASET_CREATE),H5Pclose);
  if(rows==0 || rows>=chunk) {
    H5Pset_chunk(plist,1,&chunk);
    if(c.type!=Value::String && H5Zfilter_avail(H5Z_FILTER_DEFLATE)>0) {
      H5Pset_shuffle(plist);
      H5Pset_deflate(plist,1);
    }
  }
  Id st=stringType();
  const hid_t ft = c.type==Value::Double ? H5T_IEEE_F64LE : c.type==Value::Int ? H5T_STD_I64LE : hid_t(st);
  c.dataset=Id(H5Dcreate2(group,c.name.c_str(),ft,space,H5P_DEFAULT,plist,H5P_DEFAULT),H5Dclose);
  check(c.dataset>=0,path,"cannot create the dataset of field "+c.name);
}

/// append the buffered rows of a column to its dataset (in a dataset of fixed size, write them at its start)
void appendRows(Column& c,const std::string& path,bool fixed=false) {
  const hsize_t m=c.size();
  if(m==0) {
    return;
  }
  const hsize_t n=fixed ? 0 : datasetRows(c.dataset);
  const hsize_t total=n+m;
  check(fixed || H5Dset_extent(c.dataset,&total)>=0,path,"cannot extend the dataset of field "+c.name);
  Id fspace(H5Dget_space(c.dataset),H5Sclose);
  H5Sselect_hyperslab(fspace,H5S_SELECT_SET,&n,nullptr,&m,nullptr);
  Id mspace(H5Screate_simple(1,&m,nullptr),H5Sclose);
  herr_t err=0;
  if(c.type==Value::String) {
    std::vector<const char*> p(m);
    for(hsize_t k=0; k<m; k++) {
      p[k]=c.s[k].c_str();
    }
    Id st=stringType();
    err=H5Dwrite(c.dataset,st,mspace,fspace,H5P_DEFAULT,p.data());
  } else {
    err=H5Dwrite(c.dataset,memType(c.type),mspace,fspace,H5P_DEFAULT,c.type==Value::Double ? static_cast<const void*>(c.d.data()) : static_cast<const void*>(c.i.data()));
  }
  check(err>=0,path,"cannot write the dataset of field "+c.name);
  c.clear();
}

/// read rows [from,from+count) of a column (all of them if count is 0)
void readRows(hid_t ds,Column& c,unsigned long long from,unsigned long long count,const std::string& path) {
  const unsigned long long n=datasetRows(ds);
  c.clear();
  if(count==0) {
    count=n-from;
  }
  if(count==0) {
    return;
  }
  const hsize_t start=from,m=count;
  Id fspace(H5Dget_space(ds),H5Sclose);
  H5Sselect_hyperslab(fspace,H5S_SELECT_SET,&start,nullptr,&m,nullptr);
  Id mspace(H5Screate_simple(1,&m,nullptr),H5Sclose);
  herr_t err=0;
  if(c.type==Value::Double) {
    c.d.resize(m);
    err=H5Dread(ds,H5T_NATIVE_DOUBLE,mspace,fspace,H5P_DEFAULT,c.d.data());
  } else if(c.type==Value::Int) {
    c.i.resize(m);
    err=H5Dread(ds,H5T_NATIVE_LLONG,mspace,fspace,H5P_DEFAULT,c.i.data());
  } else {
    std::vector<char*> p(m,nullptr);
    Id st=stringType();
    err=H5Dread(ds,st,mspace,fspace,H5P_DEFAULT,p.data());
    c.s.resize(m);
    for(hsize_t k=0; k<m; k++) {
      if(p[k]) {
        c.s[k]=p[k];
        H5free_memory(p[k]);
      }
    }
  }
  check(err>=0,path,"cannot read the dataset of field "+c.name);
}

/// the description of a block: its columns (with their types) and its constants
struct Signature {
  std::vector<std::string> names;
  std::vector<Value::Type> types;
  std::vector<Field> constants;
  bool operator==(const Signature& s) const {
    if(names!=s.names || types!=s.types || constants.size()!=s.constants.size()) {
      return false;
    }
    for(unsigned k=0; k<constants.size(); k++) {
      if(constants[k].name!=s.constants[k].name || !(constants[k].value==s.constants[k].value)) {
        return false;
      }
    }
    return true;
  }
};

std::string join(const std::vector<std::string>& v) {
  std::string r;
  for(const auto & s : v) {
    r+=(r.empty()?"":" ")+s;
  }
  return r;
}

/// read the signature of a block, and open its datasets
Signature readSignature(hid_t group,std::vector<Column>& columns,const std::string& path) {
  Signature sig;
  sig.names=Tools::getWords(readAttribute(group,"plumed_columns",path).s);
  columns.clear();
  for(const auto & n : sig.names) {
    Column c;
    c.name=n;
    c.dataset=Id(H5Dopen2(group,n.c_str(),H5P_DEFAULT),H5Dclose);
    check(c.dataset>=0,path,"cannot open the dataset of field "+n);
    Id t(H5Dget_type(c.dataset),H5Tclose);
    const H5T_class_t cl=H5Tget_class(t);
    c.type= cl==H5T_FLOAT ? Value::Double : cl==H5T_INTEGER ? Value::Int : Value::String;
    sig.types.push_back(c.type);
    columns.push_back(std::move(c));
  }
  for(const auto & n : Tools::getWords(readAttribute(group,"plumed_constants",path).s)) {
    sig.constants.push_back(Field{n,readAttribute(group,n,path)});
  }
  return sig;
}

unsigned countBlocks(hid_t file) {
  unsigned b=0;
  while(H5Lexists(file,blockName(b).c_str(),H5P_DEFAULT)>0) {
    b++;
  }
  return b;
}

void checkFormat(hid_t f,const std::string& path) {
  check(H5Aexists(f,"plumed_format")>0 && readAttribute(f,"plumed_format",path).s=="plumed-table",path,"not a table written by PLUMED");
}

Id openFile(const std::string& path,unsigned flags) {
  Id f(H5Fopen(path.c_str(),flags,H5P_DEFAULT),H5Fclose);
  check(f>=0,path,"cannot open the file");
  checkFormat(f,path);
  return f;
}

/// the attributes of the root group of a table
void writeRoot(hid_t file,const std::string& path) {
  Value v;
  v.s="plumed-table";
  writeAttribute(file,"plumed_format",v,path);
  v.type=Value::Int;
  v.i=1;
  writeAttribute(file,"plumed_format_version",v,path);
}

/// create block b with signature s, and the datasets of its columns (whose names and types are those of s)
Id createBlock(hid_t file,unsigned b,const Signature& s,std::vector<Column>& columns,const std::string& path,hsize_t rows=0) {
  Id group(H5Gcreate2(file,blockName(b).c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
  check(group>=0,path,"cannot create group "+blockName(b));
  Value v;
  v.s=join(s.names);
  writeAttribute(group,"plumed_columns",v,path);
  std::vector<std::string> cn;
  for(const auto & c : s.constants) {
    cn.push_back(c.name);
    plumed_massert(c.name.compare(0,7,"plumed_")!=0,"file "+path+": field names starting with plumed_ are reserved");
    writeAttribute(group,c.name,c.value,path);
  }
  v.s=join(cn);
  writeAttribute(group,"plumed_constants",v,path);
  for(auto & c : columns) {
    createDataset(group,c,path,rows);
  }
  return group;
}

/// empty columns with the names and types of a signature
std::vector<Column> emptyColumns(const Signature& s) {
  std::vector<Column> columns(s.names.size());
  for(unsigned k=0; k<s.names.size(); k++) {
    columns[k].name=s.names[k];
    columns[k].type=s.types[k];
  }
  return columns;
}

/// rows of a block (whose columns are open)
unsigned long long blockRows(const std::vector<Column>& columns) {
  return columns.empty() ? 0 : datasetRows(columns[0].dataset);
}

// Tables in memory (the members of an archive)

/// a new HDF5 file in memory
Id memoryFile(const std::string& path,const std::vector<char>* image=nullptr) {
  // HDF5 identifies the files in memory by name: every one gets its own
  static unsigned long long counter=0;
  const std::string name="plumed-memory-table-"+std::to_string(counter++);
  Id fapl(H5Pcreate(H5P_FILE_ACCESS),H5Pclose);
  check(fapl>=0 && H5Pset_fapl_core(fapl,1<<20,0)>=0,path,"cannot set up a file in memory");
  if(!image) {
    Id f(H5Fcreate(name.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,fapl),H5Fclose);
    check(f>=0,path,"cannot create a file in memory");
    return f;
  }
  check(H5Pset_file_image(fapl,const_cast<char*>(image->data()),image->size())>=0,path,"cannot set up a file in memory");
  Id f(H5Fopen(name.c_str(),H5F_ACC_RDONLY,fapl),H5Fclose);
  check(f>=0,path,"a member of the archive is not an HDF5 file (the archive is damaged)");
  checkFormat(f,path);
  return f;
}

/// the bytes of a file in memory
std::vector<char> fileImage(hid_t f,const std::string& path) {
  check(H5Fflush(f,H5F_SCOPE_LOCAL)>=0,path,"cannot flush a file in memory");
  const auto n=H5Fget_file_image(f,nullptr,0);
  check(n>0,path,"cannot get the image of a file in memory");
  std::vector<char> image(n);
  check(H5Fget_file_image(f,image.data(),image.size())==n,path,"cannot get the image of a file in memory");
  return image;
}

/// a complete table with a single block made of these rows (they are removed from the columns)
std::vector<char> tableImage(const Signature& s,std::vector<Column>& columns,const std::string& path) {
  Id f=memoryFile(path);
  writeRoot(f,path);
  const hsize_t rows=columns.empty() ? 0 : columns[0].size();
  plumed_assert(rows>0);
  {
    Id g=createBlock(f,0,s,columns,path,rows);
    for(auto & c : columns) {
      appendRows(c,path,true);
      c.dataset.reset();
    }
  }
  return fileImage(f,path);
}

// Archives (.hpa, also .h5.tar): uncompressed POSIX (ustar) TAR files whose members are complete tables.
// They are only appended to: a member, once written, is never modified, so that an interrupted
// run cannot damage what was written before. A member is published in two steps: its data
// behind a header of zeros (which ends the archive for the readers), then its header, each
// followed by a synchronization to the disk. An interruption can only leave a last member
// that is not published, or whose header was written in part: the readers ignore it, and the
// writer removes it when it appends to the archive.

bool getDecimal(const std::string& t,unsigned long long& v) {
  v=0;
  for(const char c : t) {
    if(c<'0' || c>'9' || v>(~0ULL)/10-9) {
      return false;
    }
    v=v*10+(c-'0');
  }
  return !t.empty();
}

unsigned long long tarPadded(unsigned long long n) {
  return (n+tarRecord-1)/tarRecord*tarRecord;
}

/// write v in octal with exactly digits digits, followed by a null character (false if it does not fit)
bool putOctal(char* p,unsigned digits,unsigned long long v) {
  for(unsigned k=digits; k>0; k--) {
    p[k-1]='0'+(v&7);
    v>>=3;
  }
  p[digits]='\0';
  return v==0;
}

/// the header of a member
std::vector<char> tarHeader(const std::string& name,unsigned long long size,const std::string& path) {
  std::vector<char> h(tarRecord,'\0');
  plumed_assert(name.size()<100);
  std::memcpy(h.data(),name.c_str(),name.size());
  const bool ok=putOctal(h.data()+100,7,0644) && putOctal(h.data()+108,7,0) && putOctal(h.data()+116,7,0)
                && putOctal(h.data()+124,11,size) && putOctal(h.data()+136,11,static_cast<unsigned long long>(std::time(nullptr)));
  check(ok,path,"a member of the archive is too large (more than 8 GiB)");
  h[156]='0';
  std::memcpy(h.data()+257,"ustar",6);
  h[263]='0';
  h[264]='0';
  putOctal(h.data()+148,6,tarChecksum(h.data()));
  h[155]=' ';
  return h;
}

/// a member: header, data and padding
std::vector<char> tarMember(const std::string& name,const std::vector<char>& data,const std::string& path) {
  std::vector<char> m=tarHeader(name,data.size(),path);
  m.insert(m.end(),data.begin(),data.end());
  m.resize(tarRecord+tarPadded(data.size()),'\0');
  return m;
}

std::string memberName(std::size_t k) {
  char buf[32];
  std::snprintf(buf,sizeof(buf),"batch_%08zu.h5",k);
  return buf;
}

/// the entries of an archive
struct Scan {
  /// header and data offsets, and size, of the members that are tables (the regular files)
  struct Member {
    unsigned long long header=0,data=0,size=0;
  };
  std::vector<Member> members;
  /// end of the last complete entry: what follows is the end-of-archive marker, or a member
  /// not published (an interrupted run)
  unsigned long long end=0;
  unsigned long long size=0;
  /// not empty if the archive is damaged (valid entries follow an invalid header)
  std::string damage;
};

void seekFile(std::FILE* f,unsigned long long offset,const std::string& path) {
#ifdef _WIN32
  check(_fseeki64(f,offset,SEEK_SET)==0,path,"cannot seek in the file");
#else
  check(fseeko(f,offset,SEEK_SET)==0,path,"cannot seek in the file");
#endif
}

/// bytes read at given offsets from a file (opened with fopenUTF8, as all the files that
/// PLUMED opens) or from memory
class Bytes {
  std::string path;
  unique_FILE f;
  const std::vector<char>* memory=nullptr;
public:
  explicit Bytes(const std::string& p):
    path(p),
    f(fopenUTF8(p,"rb")) {
    check(f!=nullptr,path,"cannot open the file");
  }
  Bytes(const std::vector<char>& m,const std::string& p):
    path(p),
    memory(&m) {}
  void read(unsigned long long offset,char* data,std::size_t n) {
    if(memory) {
      check(offset+n<=memory->size(),path,"cannot read the file");
      std::memcpy(data,memory->data()+offset,n);
      return;
    }
    seekFile(f.get(),offset,path);
    check(std::fread(data,1,n,f.get())==n,path,"cannot read the file");
  }
};

/// scan the headers of the size bytes of an archive (only the headers are read)
Scan scanEntries(Bytes& in,unsigned long long size) {
  Scan s;
  s.size=size;
  std::vector<char> h(tarRecord);
  auto readHeader=[&](unsigned long long offset) {
    in.read(offset,h.data(),tarRecord);
  };
  while(size-s.end>=tarRecord) {
    readHeader(s.end);
    if(std::all_of(h.begin(),h.end(),[](char c) {
    return c=='\0';
  })) {
      // end-of-archive marker, or a member whose header is not written yet
      break;
    }
    unsigned long long n=0;
    if(!validTarHeader(h.data(),n)) {
      // the header of the last member, written in part by an interrupted run;
      // if valid entries follow, the archive is damaged
      for(unsigned long long o=s.end+tarRecord; size-o>=tarRecord; o+=tarRecord) {
        readHeader(o);
        unsigned long long m=0;
        if(validTarHeader(h.data(),m)) {
          s.damage="the archive is damaged: invalid TAR header at byte "+std::to_string(s.end);
          break;
        }
      }
      if(s.damage.empty() && s.end==0) {
        s.damage="not a TAR file in a format that PLUMED reads (POSIX ustar or pax, GNU; not compressed)";
      }
      break;
    }
    const unsigned long long next=s.end+tarRecord+tarPadded(n);
    if(next>size) {
      // written in part
      break;
    }
    if(h[156]=='x') {
      // a pax extended header: its size record replaces the size of the next entry
      std::string records(n,'\0');
      in.read(s.end+tarRecord,records.data(),n);
      const unsigned long long header=next;
      s.end=next;
      for(std::size_t p=0; p<records.size();) {
        // "length key=value\n"
        unsigned long long length=0;
        const std::size_t sp=records.find(' ',p);
        if(sp==std::string::npos || !getDecimal(records.substr(p,sp-p),length) || length==0 || p+length>records.size()) {
          break;
        }
        const std::string record=records.substr(sp+1,p+length-sp-2);
        if(record.compare(0,5,"size=")==0 && size-header>=tarRecord) {
          unsigned long long m=0;
          if(getDecimal(record.substr(5),m)) {
            readHeader(header);
            unsigned long long ignored=0;
            if(validTarHeader(h.data(),ignored) && header+tarRecord+tarPadded(m)<=size && (h[156]=='0' || h[156]=='\0')) {
              s.members.push_back({header,header+tarRecord,m});
              s.end=header+tarRecord+tarPadded(m);
            }
          }
        }
        p+=length;
      }
      continue;
    }
    // the other entries (directories, other extended headers, ...) are skipped
    if(h[156]=='0' || h[156]=='\0') {
      s.members.push_back({s.end,s.end+tarRecord,n});
    }
    s.end=next;
  }
  return s;
}

unsigned long long fileSize(const std::string& path) {
  std::error_code ec;
  const unsigned long long n=std::filesystem::file_size(pathFromUTF8(path),ec);
  check(!ec,path,"cannot get the size of the file");
  return n;
}

Scan scanArchive(const std::string& path) {
  Bytes in(path);
  return scanEntries(in,fileSize(path));
}

void checkNotDamaged(const Scan& s,const std::string& path) {
  check(s.damage.empty(),path,s.damage);
}

/// make what has been written to the file durable (it survives a crash of the machine)
void syncFile(std::FILE* f,const std::string& path) {
  check(std::fflush(f)==0,path,"cannot write the file");
#ifdef _WIN32
  check(_commit(_fileno(f))==0,path,"cannot write the file to disk");
#else
  check(fsync(fileno(f))==0,path,"cannot write the file to disk");
#endif
}

/// File synchronization does not publish a new name or a rename on POSIX.
void syncDirectory(const std::string& path) {
#ifndef _WIN32
  auto parent=pathFromUTF8(path).parent_path();
  if(parent.empty()) {
    parent=".";
  }
  const int fd=::open(parent.c_str(),O_RDONLY);
  check(fd>=0,path,"cannot open the parent directory");
  const int status=fsync(fd);
  const int error=errno;
  ::close(fd);
  check(status==0,path,"cannot synchronize the parent directory: "+std::string(std::strerror(error)));
#else
  (void)path; // Directory changes use write-through renames below on Windows.
#endif
}

void renameSynced(const std::string& from,const std::string& to) {
#ifdef _WIN32
  check(MoveFileExW(pathFromUTF8(from).c_str(),pathFromUTF8(to).c_str(),
                    MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0,to,"cannot durably rename the file");
#else
  std::error_code ec;
  std::filesystem::rename(pathFromUTF8(from),pathFromUTF8(to),ec);
  check(!ec,to,"cannot rename the file: "+ec.message());
  syncDirectory(to);
#endif
}

void writeAt(std::FILE* f,unsigned long long offset,const char* data,std::size_t n,const std::string& path) {
  seekFile(f,offset,path);
  check(std::fwrite(data,1,n,f)==n,path,"cannot write the file");
}

// Rows removed from inside a member (trimming at a restart) are removed by replacing the
// member with a new one with the rows that are kept. So that an interruption cannot lose
// these rows, the members that replace the ones after offset are first saved in a trim
// record (path.plumed-trim, renamed when complete): then the archive is truncated at
// offset and the members appended. A trim record left by an interrupted run is applied by
// the readers (without modifying the files) and by the writer (which completes the trim).

std::string trimRecord(const std::string& path) {
  return path+".plumed-trim";
}

/// the trim record of an archive, if any: offset, and the entries that replace what follows it
bool readTrimRecord(const std::string& path,unsigned long long& offset,std::vector<char>& entries) {
  const std::string record=trimRecord(path);
  std::error_code ec;
  if(!std::filesystem::exists(pathFromUTF8(record),ec)) {
    return false;
  }
  Bytes in(record);
  std::vector<char> h(tarRecord+1,'\0');
  in.read(0,h.data(),tarRecord);
  unsigned long long n=0;
  check(std::sscanf(h.data(),"plumed-trim-record 1 %llu %llu",&offset,&n)==2 && fileSize(record)==tarRecord+n,record,"invalid trim record");
  entries.resize(n);
  in.read(tarRecord,entries.data(),n);
  return true;
}

/// complete an interrupted trim (only the writer)
void completeTrim(const std::string& path) {
  const std::string record=trimRecord(path);
  std::error_code ec;
  // a record not complete: the archive was not modified yet
  std::filesystem::remove(pathFromUTF8(record+".tmp"),ec);
  unsigned long long offset=0;
  std::vector<char> entries;
  if(!readTrimRecord(path,offset,entries)) {
    return;
  }
  check(fileSize(path)>=offset,path,"the archive is shorter than its trim record says");
  std::filesystem::resize_file(pathFromUTF8(path),offset,ec);
  check(!ec,path,"cannot truncate the file: "+ec.message());
  {
    unique_FILE f(fopenUTF8(path,"r+b"));
    check(f!=nullptr,path,"cannot open the file");
    writeAt(f.get(),offset,entries.data(),entries.size(),path);
    syncFile(f.get(),path);
  }
#ifdef _WIN32
  // A write-through rename retires the record durably. A leftover .tmp is
  // ignored by readers and removed by the next writer, never replayed.
  renameSynced(record,record+".tmp");
  std::filesystem::remove(pathFromUTF8(record+".tmp"),ec);
#else
  std::filesystem::remove(pathFromUTF8(record),ec);
#endif
  check(!ec,record,"cannot remove the file: "+ec.message());
  syncDirectory(record);
}

/// replace what follows offset with these entries
void replaceTail(const std::string& path,unsigned long long offset,const std::vector<char>& entries) {
  std::error_code ec;
  if(entries.empty()) {
    // a single operation
    std::filesystem::resize_file(pathFromUTF8(path),offset,ec);
    check(!ec,path,"cannot truncate the file: "+ec.message());
    return;
  }
  const std::string record=trimRecord(path);
  {
    unique_FILE f(fopenUTF8(record+".tmp","wb"));
    check(f!=nullptr,record,"cannot open the file");
    std::vector<char> h(tarRecord,'\0');
    std::snprintf(h.data(),h.size(),"plumed-trim-record 1 %llu %llu\n",offset,static_cast<unsigned long long>(entries.size()));
    check(std::fwrite(h.data(),1,h.size(),f.get())==h.size() && std::fwrite(entries.data(),1,entries.size(),f.get())==entries.size(),record,"cannot write the file");
    syncFile(f.get(),record);
  }
  // The recovery name must be durable before truncation can discard its data.
  renameSynced(record+".tmp",record);
  completeTrim(path);
}

/// the HDF5 files of a table: the file itself or, for an archive, its members (opened in memory)
class Files {
  std::string path;
  bool archive;
  Scan scan;
  std::unique_ptr<Bytes> in;
  /// the members of a trim record not completed, that replace those after recordOffset
  std::vector<char> recordEntries;
  struct Source {
    bool inRecord;
    unsigned long long data,size;
  };
  std::vector<Source> sources;
public:
  explicit Files(const std::string& p):
    path(p),
    archive(h5table::isArchive(p)) {
    if(!archive) {
      return;
    }
    scan=scanArchive(path);
    checkNotDamaged(scan,path);
    in=std::make_unique<Bytes>(path);
    unsigned long long offset=scan.end;
    const bool record=readTrimRecord(path,offset,recordEntries);
    for(const auto & m : scan.members) {
      if(m.header<offset) {
        sources.push_back({false,m.data,m.size});
      }
    }
    if(record) {
      Bytes rin(recordEntries,trimRecord(path));
      const Scan rs=scanEntries(rin,recordEntries.size());
      for(const auto & m : rs.members) {
        sources.push_back({true,m.data,m.size});
      }
    }
  }
  std::size_t size() const {
    return archive ? sources.size() : 1;
  }
  bool isArchive() const {
    return archive;
  }
  /// the entries of the archive (without the trim record, see trimArchive)
  const Scan& getScan() const {
    return scan;
  }
  std::vector<char> bytes(std::size_t k) {
    const auto & s=sources[k];
    if(s.inRecord) {
      return std::vector<char>(recordEntries.begin()+s.data,recordEntries.begin()+s.data+s.size);
    }
    std::vector<char> image(s.size);
    in->read(s.data,image.data(),image.size());
    return image;
  }
  Id open(std::size_t k) {
    if(!archive) {
      return openFile(path,H5F_ACC_RDONLY);
    }
    const std::vector<char> image=bytes(k);
    return memoryFile(path,&image);
  }
};

/// rows read at once, so that the memory used does not grow with the size of the table
constexpr unsigned long long readWindow=65536;

}

struct Writer::Impl {
  std::string path;
  Id file;
  Id group;
  int block=-1;
  Signature sig;
  std::vector<Column> columns;
  std::size_t buffered=0;
  std::size_t batch=4096;
  /// an archive: the buffered rows are written as new members
  bool archive=false;
  unique_FILE tar;
  /// end of the last member, and number of members
  unsigned long long end=0;
  std::size_t members=0;
  bool failed=false;
  void checkWritable() const {
    check(!failed,path,"a previous write failed; close and reopen the writer before continuing");
  }
  void writeBuffer() {
    checkWritable();
    try {
      if(archive) {
        publish();
        return;
      }
      for(auto & c : columns) {
        appendRows(c,path);
      }
      buffered=0;
    } catch(...) {
      failed=true;
      throw;
    }
  }
  /// write the buffered rows of an archive as a new member
  void publish() {
    if(buffered==0) {
      return;
    }
    const std::vector<char> image=tableImage(sig,columns,path);
    buffered=0;
    const std::vector<char> header=tarHeader(memberName(members+1),image.size(),path);
    // the data, behind a header of zeros: for the readers, the archive still ends at end
    std::vector<char> data(tarRecord+tarPadded(image.size()),'\0');
    std::memcpy(data.data()+tarRecord,image.data(),image.size());
    writeAt(tar.get(),end,data.data(),data.size(),path);
    syncFile(tar.get(),path);
    // then the header, which publishes the member
    writeAt(tar.get(),end,header.data(),header.size(),path);
    syncFile(tar.get(),path);
    if(members==0) {
      syncDirectory(path);
    }
    end+=data.size();
    members++;
  }
  void newBlock(const Signature& s) {
    writeBuffer();
    columns=emptyColumns(s);
    sig=s;
    block++;
    if(!archive) {
      group.reset();
      group=createBlock(file,block,s,columns,path);
    }
  }
};

Writer::Writer(const std::string& path,bool append,std::size_t batch):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  impl->batch=std::max<std::size_t>(batch,1);
  std::error_code ec;
  const bool exists=append && std::filesystem::exists(pathFromUTF8(path),ec);
  // an existing file is appended to as what it is, a new one is what its name says
  impl->archive= exists ? isArchive(path) : hasArchiveExtension(path);
  if(impl->archive) {
    if(exists) {
      completeTrim(path);
      // only the headers are read; what follows the last member (the end-of-archive marker,
      // or a member not published by an interrupted run) is removed
      const Scan s=scanArchive(path);
      checkNotDamaged(s,path);
      if(s.end<s.size) {
        std::filesystem::resize_file(pathFromUTF8(path),s.end,ec);
        check(!ec,path,"cannot truncate the file: "+ec.message());
      }
      impl->members=s.members.size();
      impl->end=s.end;
    }
    impl->tar.reset(fopenUTF8(path,exists ? "r+b" : "w+b"));
    check(impl->tar!=nullptr,path,"cannot open the file");
  } else {
    // a single HDF5 file is only written from scratch (h5convert): it would be modified in place
    check(!exists,path,"a single HDF5 file cannot be appended to: convert it into an archive (.hpa) with plumed h5convert");
    impl->file=Id(H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
    check(impl->file>=0,path,"cannot create the file");
    writeRoot(impl->file,path);
  }
}

Writer::~Writer() {
  Lock lock(hdf5Mutex());
  if(impl->failed) {
    return;
  }
  try {
    impl->writeBuffer();
    if(impl->archive) {
      // end-of-archive marker (removed when the archive is appended to)
      const std::vector<char> zeros(2*tarRecord,'\0');
      writeAt(impl->tar.get(),impl->end,zeros.data(),zeros.size(),impl->path);
    }
  } catch(...) {
  }
}

void Writer::addRow(const std::vector<Field>& fields,const std::vector<Field>& constants) {
  Lock lock(hdf5Mutex());
  impl->checkWritable();
  Signature s;
  for(const auto & f : fields) {
    s.names.push_back(f.name);
    s.types.push_back(f.value.type);
  }
  s.constants=constants;
  if(impl->block<0 || !(s==impl->sig)) {
    impl->newBlock(s);
  }
  for(unsigned k=0; k<fields.size(); k++) {
    impl->columns[k].push(fields[k].value);
  }
  if(++impl->buffered>=impl->batch) {
    impl->writeBuffer();
  }
}

void Writer::flush() {
  Lock lock(hdf5Mutex());
  try {
    impl->writeBuffer();
    if(!impl->archive) {
      check(H5Fflush(impl->file,H5F_SCOPE_GLOBAL)>=0,impl->path,"cannot flush the file");
    }
  } catch(...) {
    impl->failed=true;
    throw;
  }
}

struct Reader::Impl {
  std::string path;
  std::unique_ptr<Files> files;
  /// the next file (member of an archive) to read
  std::size_t next=0;
  Id file;
  unsigned nblocks=0;
  unsigned block=0;
  bool loaded=false;
  std::vector<Column> columns;
  /// rows of the current block, rows read so far, and the rows in memory (from first)
  unsigned long long nrows=0;
  unsigned long long first=0;
  unsigned long long row=0;
  std::deque<std::string> pending;
  /// the signature of the last block read: in an archive, consecutive blocks with the same
  /// signature are the same block of the table (as when the rows are written to a single file)
  Signature last;
  /// the rows of the current block from row on (at most readWindow)
  void window() {
    first=row;
    const unsigned long long count=std::min(readWindow,nrows-row);
    for(auto & c : columns) {
      readRows(c.dataset,c,first,count,path);
    }
  }
  bool load() {
    while(block>=nblocks) {
      if(next>=files->size()) {
        return false;
      }
      columns.clear();
      file=files->open(next++);
      nblocks=countBlocks(file);
      block=0;
    }
    Id g(H5Gopen2(file,blockName(block).c_str(),H5P_DEFAULT),H5Gclose);
    check(g>=0,path,"cannot open group "+blockName(block));
    const Signature sig=readSignature(g,columns,path);
    nrows=blockRows(columns);
    row=0;
    if(nrows>0) {
      window();
    }
    if(!(loaded && files->isArchive() && sig==last)) {
      pending.push_back("#! FIELDS "+join(sig.names));
      for(const auto & c : sig.constants) {
        pending.push_back("#! SET "+c.name+" "+c.value.text());
      }
    }
    last=sig;
    block++;
    loaded=true;
    return true;
  }
};

Reader::Reader(const std::string& path):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  impl->files=std::make_unique<Files>(path);
  if(!impl->files->isArchive()) {
    // a missing or invalid file is an error here
    impl->file=impl->files->open(0);
    impl->nblocks=countBlocks(impl->file);
    impl->next=1;
  }
}

Reader::~Reader() {
  Lock lock(hdf5Mutex());
  impl.reset();
}

bool Reader::getline(std::string& line) {
  Lock lock(hdf5Mutex());
  while(true) {
    if(!impl->pending.empty()) {
      line=impl->pending.front();
      impl->pending.pop_front();
      return true;
    }
    if(impl->loaded && impl->row<impl->nrows) {
      if(impl->row-impl->first>=impl->columns[0].size()) {
        impl->window();
      }
      line.clear();
      for(const auto & c : impl->columns) {
        line+=(line.empty()?"":" ")+c.get(impl->row-impl->first).text();
      }
      impl->row++;
      return true;
    }
    if(!impl->load()) {
      return false;
    }
  }
}

namespace {
/// call f(file,signature,columns,rows) for every block of every file of a table, with the columns open
template<class F>
void forEachBlock(Files& files,const std::string& path,const F& f,std::size_t from=0) {
  for(std::size_t k=from; k<files.size(); k++) {
    Id file=files.open(k);
    const unsigned nb=countBlocks(file);
    for(unsigned b=0; b<nb; b++) {
      Id g(H5Gopen2(file,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
      check(g>=0,path,"cannot open group "+blockName(b));
      std::vector<Column> cols;
      const Signature sig=readSignature(g,cols,path);
      f(k,sig,cols,blockRows(cols));
    }
  }
}

/// add rows [from,to) of a block to a writer (readWindow rows at a time)
void copyRows(Writer& out,const Signature& sig,std::vector<Column>& cols,unsigned long long from,unsigned long long to,const std::string& path) {
  std::vector<Field> row(cols.size());
  for(unsigned long long w=from; w<to; w+=readWindow) {
    const unsigned long long count=std::min(readWindow,to-w);
    for(auto & c : cols) {
      readRows(c.dataset,c,w,count,path);
    }
    for(unsigned long long r=0; r<count; r++) {
      for(unsigned j=0; j<cols.size(); j++) {
        row[j].name=cols[j].name;
        row[j].value=cols[j].get(r);
      }
      out.addRow(row,sig.constants);
    }
  }
}

/// Remove the rows of an archive from row r of member k on (rows counted over the blocks of the
/// member, 0<r<rows of member k or r=0). Only the members from k on are read: the rows removed are
/// written in the backup (if not empty), the rows of member k that are kept make new members.
unsigned long long cutArchive(const std::string& path,Files& files,std::size_t k,unsigned long long r,const std::string& backup) {
  const Scan& scan=files.getScan();
  std::unique_ptr<Writer> out;
  std::vector<char> entries;
  unsigned long long removed=0,seen=0;
  std::size_t current=files.size(),produced=0;
  forEachBlock(files,path,[&](std::size_t m,const Signature& sig,std::vector<Column>& cols,unsigned long long n) {
    if(m!=current) {
      current=m;
      seen=0;
    }
    const unsigned long long kept= m==k ? std::min(n,r>seen ? r-seen : 0) : 0;
    seen+=n;
    if(n>kept) {
      if(backup.length()>0) {
        if(!out) {
          out=std::make_unique<Writer>(backup,false);
        }
        copyRows(*out,sig,cols,kept,n,path);
      }
      removed+=n-kept;
    }
    if(kept>0) {
      for(auto & c : cols) {
        readRows(c.dataset,c,0,kept,path);
      }
      std::vector<Column> mine=emptyColumns(sig);
      for(unsigned j=0; j<cols.size(); j++) {
        std::swap(mine[j].d,cols[j].d);
        std::swap(mine[j].i,cols[j].i);
        std::swap(mine[j].s,cols[j].s);
      }
      const std::vector<char> m2=tarMember(memberName(k+1+produced++),tableImage(sig,mine,path),path);
      entries.insert(entries.end(),m2.begin(),m2.end());
    }
  },k);
  out.reset();
  replaceTail(path,scan.members[k].header,entries);
  return removed;
}

/// the time column of a block as doubles (false if the block has none)
bool timeColumn(const Signature& sig,std::vector<Column>& cols,unsigned long long n,std::vector<double>& t,const std::string& path) {
  for(unsigned j=0; j<sig.names.size(); j++) {
    if(sig.names[j]=="time" && sig.types[j]!=Value::String) {
      readRows(cols[j].dataset,cols[j],0,n,path);
      t.resize(n);
      for(unsigned long long i=0; i<n; i++) {
        t[i]= cols[j].type==Value::Double ? cols[j].d[i] : static_cast<double>(cols[j].i[i]);
      }
      cols[j].clear();
      return true;
    }
  }
  return false;
}
}

unsigned long long trimArchive(const std::string& path,double tcut,const std::string& backup,bool& hasTime,bool& increasing) {
  Lock lock(hdf5Mutex());
  completeTrim(path);
  Files files(path);
  hasTime=false;
  increasing=true;
  // the members are read from the last one, until a row with time before tcut is found:
  // the cost does not depend on the rows that are kept
  const std::size_t n=files.size();
  std::size_t k=0;
  unsigned long long r=0;
  bool found=false;
  double following=std::numeric_limits<double>::infinity();
  for(std::size_t m=n; m>0 && !found; m--) {
    Id file=files.open(m-1);
    const unsigned nb=countBlocks(file);
    std::vector<unsigned long long> start(nb+1,0);
    std::vector<std::vector<double>> times(nb);
    std::vector<bool> withTime(nb,false);
    for(unsigned b=0; b<nb; b++) {
      Id g(H5Gopen2(file,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
      check(g>=0,path,"cannot open group "+blockName(b));
      std::vector<Column> cols;
      const Signature sig=readSignature(g,cols,path);
      const unsigned long long rows=blockRows(cols);
      start[b+1]=start[b]+rows;
      withTime[b]=timeColumn(sig,cols,rows,times[b],path);
    }
    for(unsigned b=nb; b>0 && !found; b--) {
      if(!withTime[b-1]) {
        continue;
      }
      hasTime=true;
      const auto & t=times[b-1];
      for(std::size_t i=t.size(); i>0; i--) {
        if(t[i-1]<tcut) {
          found=true;
          k=m-1;
          r=start[b-1]+i;
          break;
        }
        if(t[i-1]>following) {
          increasing=false;
        }
        following=t[i-1];
      }
    }
    if(found && r==start[nb]) {
      // the cut is at the end of the member
      k++;
      r=0;
    }
  }
  if(!hasTime || !increasing) {
    return 0;
  }
  if(!found) {
    // all the rows are after tcut
    k=0;
    r=0;
  }
  if(k>=n) {
    // nothing after tcut
    return 0;
  }
  return cutArchive(path,files,k,r,backup);
}

void copy(const std::string& input,const std::string& output,std::size_t batch) {
  Lock lock(hdf5Mutex());
  std::error_code ec;
  check(!std::filesystem::equivalent(pathFromUTF8(input),pathFromUTF8(output),ec),input,
        "input and output refer to the same file");
  Files files(input);
  Writer out(output,false,batch);
  forEachBlock(files,input,[&](std::size_t,const Signature& sig,std::vector<Column>& cols,unsigned long long n) {
    copyRows(out,sig,cols,0,n,input);
  });
  out.flush();
}

#endif

}
}

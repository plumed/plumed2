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
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <utility>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef __PLUMED_HAS_LIBHDF5
#include <hdf5.h>
#endif

namespace PLMD {
namespace h5table {

bool isArchive(const std::string& path) {
  return path.size()>7 && path.compare(path.size()-7,7,".h5.tar")==0;
}

bool isH5(const std::string& path) {
  return Tools::extension(path)=="h5" || isArchive(path);
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
  plumed_merror("file "+path+": PLUMED was compiled without libhdf5, files with extension .h5 or .h5.tar cannot be used (configure PLUMED with --enable-libhdf5)");
}
}

struct Writer::Impl {};
Writer::Writer(const std::string& path,bool) {
  noHDF5(path);
}
Writer::~Writer()=default;
void Writer::addRow(const std::vector<Field>&,const std::vector<Field>&) {}
void Writer::flush() {}
std::vector<unsigned long long> Writer::rows() {
  return {};
}

struct Reader::Impl {};
Reader::Reader(const std::string& path) {
  noHDF5(path);
}
Reader::~Reader()=default;
bool Reader::getline(std::string&) {
  return false;
}

std::vector<unsigned long long> rows(const std::string& path) {
  noHDF5(path);
}

void truncate(const std::string& path,const std::vector<unsigned long long>&,const std::string&) {
  noHDF5(path);
}

void copy(const std::string& input,const std::string&) {
  noHDF5(input);
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

// Archives (.h5.tar): uncompressed POSIX (ustar) TAR files whose members are complete tables.
// They are only appended to: a member, once written, is never modified, so that an interrupted
// run cannot damage what was written before. A member written in part (the run was interrupted
// while writing it) is ignored by the readers, and removed when the archive is appended to.

constexpr unsigned long long tarRecord=512;

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

unsigned tarChecksum(const char* h) {
  unsigned sum=0;
  for(unsigned k=0; k<tarRecord; k++) {
    sum+= (k>=148 && k<156) ? static_cast<unsigned>(' ') : static_cast<unsigned char>(h[k]);
  }
  return sum;
}

/// a member of an archive: header, data and padding
std::vector<char> tarMember(const std::string& name,const std::vector<char>& data,const std::string& path) {
  std::vector<char> m(tarRecord+tarPadded(data.size()),'\0');
  char* h=m.data();
  plumed_assert(name.size()<100);
  std::memcpy(h,name.c_str(),name.size());
  bool ok=putOctal(h+100,7,0644) && putOctal(h+108,7,0) && putOctal(h+116,7,0)
          && putOctal(h+124,11,data.size()) && putOctal(h+136,11,static_cast<unsigned long long>(std::time(nullptr)));
  check(ok,path,"a member of the archive is too large (more than 8 GiB)");
  h[156]='0';
  std::memcpy(h+257,"ustar",6);
  h[263]='0';
  h[264]='0';
  putOctal(h+148,6,tarChecksum(h));
  h[155]=' ';
  std::memcpy(h+tarRecord,data.data(),data.size());
  return m;
}

std::string memberName(std::size_t k) {
  char buf[32];
  std::snprintf(buf,sizeof(buf),"batch_%08zu.h5",k);
  return buf;
}

/// what an archive contains
struct Scan {
  /// header and data offsets, and size, of the members that are tables (the regular files)
  struct Member {
    unsigned long long header=0,data=0,size=0;
  };
  std::vector<Member> members;
  /// end of the last complete member (what follows is the end-of-archive marker, or a member written in part)
  unsigned long long end=0;
  unsigned long long fileSize=0;
  /// not empty if the archive is damaged (not just interrupted while writing its last member)
  std::string damage;
};

Scan scanArchive(const std::string& path) {
  Scan s;
  std::error_code ec;
  s.fileSize=std::filesystem::file_size(pathFromUTF8(path),ec);
  check(!ec,path,"cannot get the size of the file");
  std::ifstream in(pathFromUTF8(path),std::ios::binary);
  check(static_cast<bool>(in),path,"cannot open the file");
  std::vector<char> h(tarRecord);
  // the last record can be incomplete: the run was interrupted while writing it
  while(s.fileSize-s.end>=tarRecord) {
    in.seekg(s.end);
    in.read(h.data(),tarRecord);
    check(static_cast<bool>(in),path,"cannot read the file");
    if(std::all_of(h.begin(),h.end(),[](char c) {
    return c=='\0';
  })) {
      // end-of-archive marker: only zeros can follow
      std::vector<char> rest(1<<16);
      while(in) {
        in.read(rest.data(),rest.size());
        if(!std::all_of(rest.begin(),rest.begin()+in.gcount(),[](char c) {
        return c=='\0';
      })) {
          s.damage="data after the end-of-archive marker at byte "+std::to_string(s.end);
          break;
        }
      }
      break;
    }
    unsigned long long sum=0,size=0;
    if(std::memcmp(h.data()+257,"ustar",5)!=0 || !getOctal(h.data()+148,8,sum) || sum!=tarChecksum(h.data()) || !getOctal(h.data()+124,12,size)) {
      s.damage="invalid TAR header at byte "+std::to_string(s.end);
      break;
    }
    const unsigned long long next=s.end+tarRecord+tarPadded(size);
    if(next>s.fileSize) {
      // written in part
      break;
    }
    // the other entries (directories, extended headers, ...) are skipped
    if(h[156]=='0' || h[156]=='\0') {
      s.members.push_back({s.end,s.end+tarRecord,size});
    }
    s.end=next;
  }
  return s;
}

void checkNotDamaged(const Scan& s,const std::string& path) {
  check(s.damage.empty(),path,"the archive is damaged: "+s.damage);
}

/// the HDF5 files of a table: the file itself or, for an archive, its members (opened in memory)
class Files {
  std::string path;
  bool archive;
  Scan scan;
  std::ifstream in;
public:
  explicit Files(const std::string& p):
    path(p),
    archive(h5table::isArchive(p)) {
    if(archive) {
      scan=scanArchive(path);
      checkNotDamaged(scan,path);
      in.open(pathFromUTF8(path),std::ios::binary);
      check(static_cast<bool>(in),path,"cannot open the file");
    }
  }
  std::size_t size() const {
    return archive ? scan.members.size() : 1;
  }
  bool isArchive() const {
    return archive;
  }
  const Scan& getScan() const {
    return scan;
  }
  std::vector<char> bytes(std::size_t k) {
    const auto & m=scan.members[k];
    std::vector<char> image(m.size);
    in.seekg(m.data);
    in.read(image.data(),image.size());
    check(static_cast<bool>(in),path,"cannot read a member of the archive");
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

/// make what has been written to the file durable (it survives a crash of the machine)
void syncFile(std::FILE* f,const std::string& path) {
  check(std::fflush(f)==0,path,"cannot write the file");
#ifdef _WIN32
  check(_commit(_fileno(f))==0,path,"cannot write the file to disk");
#else
  check(fsync(fileno(f))==0,path,"cannot write the file to disk");
#endif
}

}

struct Writer::Impl {
  std::string path;
  Id file;
  Id group;
  int block=-1;
  Signature sig;
  std::vector<Column> columns;
  std::size_t buffered=0;
  /// an archive: the buffered rows are written as new members
  bool archive=false;
  unique_FILE tar;
  std::size_t members=0;
  void writeBuffer() {
    if(archive) {
      publish();
      return;
    }
    for(auto & c : columns) {
      appendRows(c,path);
    }
    buffered=0;
  }
  /// write the buffered rows of an archive as a new member (with a single write, at the end of the file)
  void publish() {
    if(buffered==0) {
      return;
    }
    const std::vector<char> m=tarMember(memberName(members+1),tableImage(sig,columns,path),path);
    check(std::fwrite(m.data(),1,m.size(),tar.get())==m.size() && std::fflush(tar.get())==0,path,"cannot write the file");
    members++;
    buffered=0;
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

Writer::Writer(const std::string& path,bool append):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  impl->archive=isArchive(path);
  std::error_code ec;
  const bool exists=std::filesystem::exists(pathFromUTF8(path),ec);
  if(impl->archive) {
    if(append && exists) {
      // what follows the last complete member (the end-of-archive marker, or a member written
      // in part by an interrupted run) is removed
      const Scan s=scanArchive(path);
      checkNotDamaged(s,path);
      if(s.end<s.fileSize) {
        std::filesystem::resize_file(pathFromUTF8(path),s.end,ec);
        check(!ec,path,"cannot truncate the file: "+ec.message());
      }
      impl->members=s.members.size();
    }
    impl->tar.reset(fopenUTF8(path,append ? "ab" : "wb"));
    check(impl->tar!=nullptr,path,"cannot open the file");
  } else if(append && exists) {
    impl->file=openFile(path,H5F_ACC_RDWR);
    const unsigned nb=countBlocks(impl->file);
    if(nb>0) {
      impl->block=nb-1;
      impl->group=Id(H5Gopen2(impl->file,blockName(nb-1).c_str(),H5P_DEFAULT),H5Gclose);
      check(impl->group>=0,path,"cannot open group "+blockName(nb-1));
      impl->sig=readSignature(impl->group,impl->columns,path);
    }
  } else {
    impl->file=Id(H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
    check(impl->file>=0,path,"cannot create the file");
    writeRoot(impl->file,path);
  }
}

Writer::~Writer() {
  Lock lock(hdf5Mutex());
  try {
    impl->writeBuffer();
    if(impl->archive) {
      // end-of-archive marker (removed when the archive is appended to)
      const std::vector<char> zeros(2*tarRecord,'\0');
      std::fwrite(zeros.data(),1,zeros.size(),impl->tar.get());
    }
  } catch(...) {
  }
}

void Writer::addRow(const std::vector<Field>& fields,const std::vector<Field>& constants) {
  Lock lock(hdf5Mutex());
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
  if(++impl->buffered>=4096) {
    impl->writeBuffer();
  }
}

void Writer::flush() {
  Lock lock(hdf5Mutex());
  impl->writeBuffer();
  if(!impl->archive) {
    H5Fflush(impl->file,H5F_SCOPE_GLOBAL);
  }
}

std::vector<unsigned long long> Writer::rows() {
  Lock lock(hdf5Mutex());
  flush();
  if(impl->archive) {
    return h5table::rows(impl->path);
  }
  std::vector<unsigned long long> r;
  for(int b=0; b<=impl->block; b++) {
    Id g(H5Gopen2(impl->file,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
    std::vector<Column> cols;
    readSignature(g,cols,impl->path);
    r.push_back(blockRows(cols));
  }
  return r;
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
  unsigned long long nrows=0;
  unsigned long long row=0;
  std::deque<std::string> pending;
  /// the signature of the last block read: in an archive, consecutive blocks with the same
  /// signature are the same block of the table (as when the rows are written to a single file)
  Signature last;
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
    nrows=0;
    for(auto & c : columns) {
      readRows(c.dataset,c,0,0,path);
      nrows=c.size();
    }
    row=0;
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
      line.clear();
      for(const auto & c : impl->columns) {
        line+=(line.empty()?"":" ")+c.get(impl->row).text();
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
/// call f(member,signature,columns,rows) for every block of every file of a table, with the columns open
template<class F>
void forEachBlock(Files& files,const std::string& path,F f) {
  for(std::size_t k=0; k<files.size(); k++) {
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

/// add rows [from,to) of a block to a writer
void copyRows(Writer& out,const Signature& sig,std::vector<Column>& cols,unsigned long long from,unsigned long long to,const std::string& path) {
  if(to<=from) {
    return;
  }
  for(auto & c : cols) {
    readRows(c.dataset,c,from,to-from,path);
  }
  std::vector<Field> row(cols.size());
  for(unsigned long long r=0; r<to-from; r++) {
    for(unsigned j=0; j<cols.size(); j++) {
      row[j].name=cols[j].name;
      row[j].value=cols[j].get(r);
    }
    out.addRow(row,sig.constants);
  }
}

/// truncate an archive: only the members after the cut are rewritten, in a copy that replaces the archive
void truncateArchive(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
  // the rows to keep are the first ones of the table
  const auto current=rows(path);
  check(keep.size()<=current.size(),path,"the table has fewer blocks than expected");
  unsigned long long total=0;
  for(std::size_t b=0; b<keep.size(); b++) {
    check(keep[b]<=current[b],path,"block "+std::to_string(b)+" has fewer rows than expected");
    check(b+1==keep.size() || keep[b]==current[b],path,"only the last rows of an archive can be removed");
    total+=keep[b];
  }
  Files files(path);
  const Scan& scan=files.getScan();
  std::unique_ptr<Writer> out;
  // the first member with rows to remove, and the new members with the rows of it that are kept
  std::size_t cut=files.size();
  std::vector<char> tail;
  unsigned long long seen=0;
  forEachBlock(files,path,[&](std::size_t k,const Signature& sig,std::vector<Column>& cols,unsigned long long n) {
    const unsigned long long kept=std::min(n,total>seen ? total-seen : 0);
    seen+=n;
    if(kept==n && cut==files.size()) {
      return;
    }
    if(cut==files.size()) {
      cut=k;
    }
    // the rows removed are kept in the backup
    if(backup.length()>0 && n>kept) {
      if(!out) {
        out=std::make_unique<Writer>(backup,false);
      }
      copyRows(*out,sig,cols,kept,n,path);
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
      const auto m=tarMember(memberName(cut+1+tail.size()),tableImage(sig,mine,path),path);
      tail.insert(tail.end(),m.begin(),m.end());
    }
  });
  out.reset();
  if(cut==files.size()) {
    // nothing to remove (but what follows the last complete member)
    return;
  }
  const unsigned long long offset=scan.members[cut].header;
  std::error_code ec;
  if(tail.empty()) {
    // a single operation: the archive is either as before or truncated
    std::filesystem::resize_file(pathFromUTF8(path),offset,ec);
    check(!ec,path,"cannot truncate the file: "+ec.message());
    return;
  }
  // the members before the cut are copied (the copy keeps the permissions of the file),
  // then the copy replaces the archive: an interruption leaves the archive as it was
  const std::string tmp=path+".plumed-tmp";
  std::filesystem::copy_file(pathFromUTF8(path),pathFromUTF8(tmp),std::filesystem::copy_options::overwrite_existing,ec);
  check(!ec,path,"cannot copy the file: "+ec.message());
  std::filesystem::resize_file(pathFromUTF8(tmp),offset,ec);
  check(!ec,tmp,"cannot truncate the file: "+ec.message());
  {
    unique_FILE f(fopenUTF8(tmp,"ab"));
    check(f!=nullptr,tmp,"cannot open the file");
    check(std::fwrite(tail.data(),1,tail.size(),f.get())==tail.size(),tmp,"cannot write the file");
    syncFile(f.get(),tmp);
  }
  std::filesystem::rename(pathFromUTF8(tmp),pathFromUTF8(path),ec);
  check(!ec,path,"cannot replace the file: "+ec.message());
}
}

std::vector<unsigned long long> rows(const std::string& path) {
  Lock lock(hdf5Mutex());
  Files files(path);
  std::vector<unsigned long long> r;
  Signature last;
  forEachBlock(files,path,[&](std::size_t,const Signature& sig,std::vector<Column>&,unsigned long long n) {
    // in an archive, consecutive blocks with the same signature are the same block of the table
    if(files.isArchive() && !r.empty() && sig==last) {
      r.back()+=n;
    } else {
      r.push_back(n);
    }
    last=sig;
  });
  return r;
}

void truncate(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
  Lock lock(hdf5Mutex());
  if(isArchive(path)) {
    truncateArchive(path,keep,backup);
    return;
  }
  std::unique_ptr<Writer> out;
  {
    Id f=openFile(path,H5F_ACC_RDWR);
    const unsigned nb=countBlocks(f);
    check(keep.size()<=nb,path,"the table has fewer blocks than expected");
    for(unsigned b=0; b<nb; b++) {
      Id g(H5Gopen2(f,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
      std::vector<Column> cols;
      const Signature sig=readSignature(g,cols,path);
      const unsigned long long n=blockRows(cols);
      const unsigned long long k= b<keep.size() ? keep[b] : 0;
      check(k<=n,path,"block "+std::to_string(b)+" has fewer rows than expected");
      if(k==n && b<keep.size()) {
        continue;
      }
      // the rows removed are kept in the backup
      if(backup.length()>0 && n>k) {
        if(!out) {
          out=std::make_unique<Writer>(backup,false);
        }
        copyRows(*out,sig,cols,k,n,path);
      }
      if(b<keep.size()) {
        const hsize_t kk=k;
        for(auto & c : cols) {
          check(H5Dset_extent(c.dataset,&kk)>=0,path,"cannot shrink the dataset of field "+c.name);
        }
      } else {
        cols.clear();
        g.reset();
        check(H5Ldelete(f,blockName(b).c_str(),H5P_DEFAULT)>=0,path,"cannot remove "+blockName(b));
      }
    }
    H5Fflush(f,H5F_SCOPE_GLOBAL);
  }
  if(out) {
    out->flush();
  }
}

void copy(const std::string& input,const std::string& output) {
  Lock lock(hdf5Mutex());
  Files files(input);
  Writer out(output,false);
  forEachBlock(files,input,[&](std::size_t,const Signature& sig,std::vector<Column>& cols,unsigned long long n) {
    copyRows(out,sig,cols,0,n,input);
  });
  out.flush();
}

#endif

}
}

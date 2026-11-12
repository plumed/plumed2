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
#include "Tools.h"
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <mutex>
#include <utility>

#ifdef __PLUMED_HAS_LIBHDF5
#include <hdf5.h>
#endif

namespace PLMD {
namespace h5table {

bool isH5(const std::string& path) {
  return Tools::extension(path)=="h5";
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
  plumed_merror("file "+path+": PLUMED was compiled without libhdf5, files with extension .h5 cannot be used (configure PLUMED with --enable-libhdf5)");
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

/// create the dataset of a column in a block
void createDataset(hid_t group,Column& c,const std::string& path) {
  const hsize_t zero=0,unlimited=H5S_UNLIMITED,chunk=1024;
  Id space(H5Screate_simple(1,&zero,&unlimited),H5Sclose);
  Id plist(H5Pcreate(H5P_DATASET_CREATE),H5Pclose);
  H5Pset_chunk(plist,1,&chunk);
  if(c.type!=Value::String && H5Zfilter_avail(H5Z_FILTER_DEFLATE)>0) {
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist,1);
  }
  Id st=stringType();
  const hid_t ft = c.type==Value::Double ? H5T_IEEE_F64LE : c.type==Value::Int ? H5T_STD_I64LE : hid_t(st);
  c.dataset=Id(H5Dcreate2(group,c.name.c_str(),ft,space,H5P_DEFAULT,plist,H5P_DEFAULT),H5Dclose);
  check(c.dataset>=0,path,"cannot create the dataset of field "+c.name);
}

/// append the buffered rows of a column to its dataset
void appendRows(Column& c,const std::string& path) {
  const hsize_t m=c.size();
  if(m==0) {
    return;
  }
  const hsize_t n=datasetRows(c.dataset);
  const hsize_t total=n+m;
  check(H5Dset_extent(c.dataset,&total)>=0,path,"cannot extend the dataset of field "+c.name);
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

Id openFile(const std::string& path,unsigned flags) {
  Id f(H5Fopen(path.c_str(),flags,H5P_DEFAULT),H5Fclose);
  check(f>=0,path,"cannot open the file");
  check(H5Aexists(f,"plumed_format")>0 && readAttribute(f,"plumed_format",path).s=="plumed-table",path,"not a table written by PLUMED");
  return f;
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
  void writeBuffer() {
    for(auto & c : columns) {
      appendRows(c,path);
    }
    buffered=0;
  }
  void newBlock(const Signature& s) {
    writeBuffer();
    columns.clear();
    group.reset();
    block++;
    group=Id(H5Gcreate2(file,blockName(block).c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
    check(group>=0,path,"cannot create group "+blockName(block));
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
    for(unsigned k=0; k<s.names.size(); k++) {
      Column c;
      c.name=s.names[k];
      c.type=s.types[k];
      createDataset(group,c,path);
      columns.push_back(std::move(c));
    }
    sig=s;
  }
};

Writer::Writer(const std::string& path,bool append):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  std::error_code ec;
  if(append && std::filesystem::exists(path,ec)) {
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
    Value v;
    v.s="plumed-table";
    writeAttribute(impl->file,"plumed_format",v,path);
    v.type=Value::Int;
    v.i=1;
    writeAttribute(impl->file,"plumed_format_version",v,path);
  }
}

Writer::~Writer() {
  Lock lock(hdf5Mutex());
  try {
    impl->writeBuffer();
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
  H5Fflush(impl->file,H5F_SCOPE_GLOBAL);
}

std::vector<unsigned long long> Writer::rows() {
  Lock lock(hdf5Mutex());
  flush();
  std::vector<unsigned long long> r;
  for(int b=0; b<=impl->block; b++) {
    Id g(H5Gopen2(impl->file,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
    std::vector<Column> cols;
    readSignature(g,cols,impl->path);
    r.push_back(cols.empty() ? 0 : datasetRows(cols[0].dataset));
  }
  return r;
}

struct Reader::Impl {
  std::string path;
  Id file;
  unsigned nblocks=0;
  unsigned block=0;
  bool loaded=false;
  std::vector<Column> columns;
  unsigned long long nrows=0;
  unsigned long long row=0;
  std::deque<std::string> pending;
  bool load() {
    if(block>=nblocks) {
      return false;
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
    pending.push_back("#! FIELDS "+join(sig.names));
    for(const auto & c : sig.constants) {
      pending.push_back("#! SET "+c.name+" "+c.value.text());
    }
    block++;
    loaded=true;
    return true;
  }
};

Reader::Reader(const std::string& path):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  impl->file=openFile(path,H5F_ACC_RDONLY);
  impl->nblocks=countBlocks(impl->file);
}

Reader::~Reader()=default;

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

std::vector<unsigned long long> rows(const std::string& path) {
  Lock lock(hdf5Mutex());
  Id f=openFile(path,H5F_ACC_RDONLY);
  std::vector<unsigned long long> r;
  const unsigned nb=countBlocks(f);
  for(unsigned b=0; b<nb; b++) {
    Id g(H5Gopen2(f,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
    std::vector<Column> cols;
    readSignature(g,cols,path);
    r.push_back(cols.empty() ? 0 : datasetRows(cols[0].dataset));
  }
  return r;
}

void truncate(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
  Lock lock(hdf5Mutex());
  std::unique_ptr<Writer> out;
  {
    Id f=openFile(path,H5F_ACC_RDWR);
    const unsigned nb=countBlocks(f);
    check(keep.size()<=nb,path,"the table has fewer blocks than expected");
    for(unsigned b=0; b<nb; b++) {
      Id g(H5Gopen2(f,blockName(b).c_str(),H5P_DEFAULT),H5Gclose);
      std::vector<Column> cols;
      const Signature sig=readSignature(g,cols,path);
      const unsigned long long n=cols.empty() ? 0 : datasetRows(cols[0].dataset);
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
        for(auto & c : cols) {
          readRows(c.dataset,c,k,n-k,path);
        }
        std::vector<Field> row(cols.size());
        for(unsigned long long r=0; r<n-k; r++) {
          for(unsigned j=0; j<cols.size(); j++) {
            row[j].name=cols[j].name;
            row[j].value=cols[j].get(r);
          }
          out->addRow(row,sig.constants);
        }
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

#endif

}
}

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
#include "H5Archive.h"
#include "FileTools.h"
#include "Exception.h"
#include "Tools.h"
#include <cmath>
#include <array>
#include <cstring>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <mutex>
#include <limits>
#include <utility>

#ifdef __PLUMED_HAS_LIBHDF5
#include <hdf5.h>
#endif

namespace PLMD {
namespace h5table {

bool isH5(const std::string& path) {
  return Tools::extension(path)=="h5" || Tools::extension(path)=="hdf5" || isArchive(path);
}

bool isArchive(const std::string& path) {
  return path.size()>=7 && path.compare(path.size()-7,7,".h5.tar")==0;
}

bool isTable(const std::string& path) {
  if(Archive::recognizes(path)) {
    return true;
  }
  unique_FILE file(fopenUTF8(path,"rb"));
  if(file) {
    for(long offset=0; offset<=65536; offset=offset?offset*2:512) {
      std::array<char,8> signature= {};
      if(std::fseek(file.get(),offset,SEEK_SET)!=0 || std::fread(signature.data(),1,8,file.get())!=8) {
        break;
      }
      if(std::memcmp(signature.data(),"\211HDF\r\n\032\n",8)==0) {
        return true;
      }
    }
  }
  return false;
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

Writer::Writer(const std::string& path,bool append):Writer(path,append,4096) {}

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
Writer::Writer(const std::string& path,bool,std::size_t) {
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
Reader::Reader(const std::vector<char>&) {
  noHDF5("HDF5 image");
}
Reader::~Reader()=default;
bool Reader::getline(std::string&) {
  return false;
}

bool Reader::readRow(std::vector<Field>&,std::vector<Field>&,bool*) {
  return false;
}

std::vector<unsigned long long> rows(const std::string& path) {
  noHDF5(path);
}

std::vector<unsigned long long> keepBefore(const std::string& path,double,bool&,bool&) {
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
  check(H5Aexists(f,"plumed_format_version")>0 && readAttribute(f,"plumed_format_version",path).i==1,path,"unsupported PLUMED HDF5 table version");
  return f;
}

Id openImage(const std::vector<char>& image) {
  const std::string path="PLUMED HDF5 batch image";
  Id properties(H5Pcreate(H5P_FILE_ACCESS),H5Pclose);
  check(H5Pset_fapl_core(properties,65536,false)>=0,path,"cannot select memory driver");
  check(!image.empty() && H5Pset_file_image(properties,const_cast<char*>(image.data()),image.size())>=0,path,"cannot load memory image");
  static unsigned long long sequence=0;
  const auto name="plumed-memory-batch-"+std::to_string(++sequence);
  Id file(H5Fopen(name.c_str(),H5F_ACC_RDONLY,properties),H5Fclose);
  check(file>=0,path,"cannot open memory image");
  check(H5Aexists(file,"plumed_format")>0 && readAttribute(file,"plumed_format",path).s=="plumed-table",path,"not a PLUMED table");
  check(H5Aexists(file,"plumed_format_version")>0 && readAttribute(file,"plumed_format_version",path).i==1,path,"unsupported PLUMED HDF5 version");
  return file;
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
  std::unique_ptr<Archive> archive;
  std::unique_ptr<BatchFile> staging;
  std::unique_ptr<Writer> batch;
  std::size_t batchSize=4096,batchRecords=0,batchBytes=0;
  std::vector<unsigned long long> counts;
  bool failed=false;
  void publishBatch() {
    check(!failed,path,"a previous batch publication failed; close and reopen the archive");
    if(!batch) {
      return;
    }
    try {
      batch->flush();
      batch.reset(); // Close HDF5 handles before publishing the bytes.
      archive->add(staging->image(),counts.size()-1,batchRecords);
      staging.reset();
      batchRecords=0;
      batchBytes=0;
    } catch(...) {
      failed=true;
      throw;
    }
  }

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

Writer::Writer(const std::string& path,bool append,std::size_t batchSize):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  check(batchSize>0,path,"batch size must be positive");
  impl->batchSize=batchSize;
  if(isArchive(path)) {
    impl->archive=std::make_unique<Archive>(path,true,append);
    for(const auto& member:impl->archive->members()) {
      if(member.block>=impl->counts.size()) {
        impl->counts.push_back(0);
      }
      impl->counts.back()+=member.rows;
    }
    if(!impl->counts.empty()) {
      Id file=openImage(impl->archive->image(impl->archive->members().size()-1));
      check(countBlocks(file)==1,path,"each archive batch must contain exactly one block");
      Id group(H5Gopen2(file,"block_0",H5P_DEFAULT),H5Gclose);
      std::vector<Column> columns;
      impl->sig=readSignature(group,columns,path);
    }
    return;
  }
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
    if(impl->archive) {
      impl->publishBatch();
      impl->archive->finish();
    } else {
      impl->writeBuffer();
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
  if(impl->archive) {
    check(!impl->failed,impl->path,"a previous batch publication failed; close and reopen the archive");
    if(!impl->counts.empty() && !(impl->sig==s)) {
      impl->publishBatch();
    }
    if(!impl->batch) {
      impl->staging=std::make_unique<BatchFile>();
      impl->batch=std::make_unique<Writer>(impl->staging->path,false);
    }
    impl->batch->addRow(fields,constants);
    if(impl->counts.empty() || !(impl->sig==s)) {
      impl->counts.push_back(0);
    }
    ++impl->counts.back();
    impl->sig=s;
    ++impl->batchRecords;
    for(const auto& f:fields) {
      impl->batchBytes+=f.value.type==Value::String?f.value.s.size():8;
    }
    if(impl->batchRecords>=impl->batchSize || impl->batchBytes>=64*1024*1024) {
      impl->publishBatch();
    }
    return;
  }
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
  if(impl->archive) {
    impl->publishBatch();
    return;
  }
  impl->writeBuffer();
  check(H5Fflush(impl->file,H5F_SCOPE_GLOBAL)>=0,impl->path,"cannot flush HDF5 file");
}

std::vector<unsigned long long> Writer::rows() {
  Lock lock(hdf5Mutex());
  flush();
  if(impl->archive) {
    return impl->counts;
  }
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
  std::unique_ptr<Archive> archive;
  std::size_t member=0;
  unsigned nblocks=0,block=0;
  bool loaded=false,haveSignature=false,newBlock=false;
  std::uint64_t archiveBlock=0;
  Signature sig;
  std::vector<Column> columns;
  unsigned long long nrows=0,row=0,chunkStart=0;
  std::deque<std::string> pending;
  bool load() {
    while(block>=nblocks) {
      if(!archive || member>=archive->members().size()) {
        return false;
      }
      columns.clear();
      file=openImage(archive->image(member++));
      nblocks=countBlocks(file);
      check(nblocks==1,path,"each TAR batch must contain exactly one HDF5 block");
      block=0;
    }
    Id group(H5Gopen2(file,blockName(block++).c_str(),H5P_DEFAULT),H5Gclose);
    check(group>=0,path,"cannot open table block");
    const auto next=readSignature(group,columns,path);
    if(archive) {
      const auto logical=archive->members()[member-1].block;
      if(haveSignature && logical==archiveBlock) {
        check(sig==next,path,"schema/constants change within a logical archive block");
      }
      newBlock=!haveSignature || logical!=archiveBlock;
      archiveBlock=logical;
    } else {
      newBlock=true;
    }
    sig=next;
    haveSignature=true;
    nrows=columns.empty()?0:datasetRows(columns[0].dataset);
    for(const auto& column:columns) {
      check(datasetRows(column.dataset)==nrows,path,"columns have inconsistent lengths");
    }
    if(archive) {
      check(nrows==archive->members()[member-1].rows,path,"batch row count disagrees with TAR header");
    }
    row=0;
    chunkStart=0;
    if(newBlock) {
      pending.push_back("#! FIELDS "+join(sig.names));
      for(const auto& c:sig.constants) {
        pending.push_back("#! SET "+c.name+" "+c.value.text());
      }
    }
    loaded=true;
    return true;
  }
  void data(std::vector<Field>& fields,std::vector<Field>& constants) {
    if(columns[0].size()==0 || row-chunkStart>=columns[0].size()) {
      chunkStart=row;
      for(auto& c:columns) {
        readRows(c.dataset,c,row,std::min<unsigned long long>(4096,nrows-row),path);
      }
    }
    fields.resize(columns.size());
    for(std::size_t k=0; k<columns.size(); ++k) {
      fields[k]=Field{columns[k].name,columns[k].get(row-chunkStart)};
    }
    constants=sig.constants;
    ++row;
  }
};

Reader::Reader(const std::string& path):impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  if(Archive::recognizes(path)) {
    impl->archive=std::make_unique<Archive>(path);
  } else {
    impl->file=openFile(path,H5F_ACC_RDONLY);
    impl->nblocks=countBlocks(impl->file);
  }
}
Reader::Reader(const std::vector<char>& image):impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path="HDF5 batch image";
  impl->file=openImage(image);
  impl->nblocks=countBlocks(impl->file);
}
Reader::~Reader()=default;

bool Reader::readRow(std::vector<Field>& fields,std::vector<Field>& constants,bool* newBlock) {
  Lock lock(hdf5Mutex());
  while(!impl->loaded || impl->row>=impl->nrows) {
    if(!impl->load()) {
      return false;
    }
  }
  impl->pending.clear();
  if(newBlock) {
    *newBlock=impl->newBlock;
  }
  impl->newBlock=false;
  impl->data(fields,constants);
  return true;
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
      std::vector<Field> fields,constants;
      impl->data(fields,constants);
      line.clear();
      for(const auto& f:fields) {
        line+=(line.empty()?"":" ")+f.value.text();
      }
      return true;
    }
    if(!impl->load()) {
      return false;
    }
  }
}

std::vector<unsigned long long> rows(const std::string& path) {
  Lock lock(hdf5Mutex());
  if(Archive::recognizes(path)) {
    Archive archive(path);
    std::vector<unsigned long long> result;
    for(const auto& member:archive.members()) {
      if(member.block>=result.size()) {
        result.push_back(0);
      }
      result.back()+=member.rows;
    }
    return result;
  }
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

std::vector<unsigned long long> keepBefore(const std::string& path,double cut,bool& hasTime,bool& increasing) {
  Lock lock(hdf5Mutex());
  Archive archive(path);
  std::vector<unsigned long long> keep;
  for(const auto& member:archive.members()) {
    if(member.block>=keep.size()) {
      keep.push_back(0);
    }
    keep.back()+=member.rows;
  }
  hasTime=false;
  increasing=true;
  double next=std::numeric_limits<double>::infinity();
  for(std::size_t k=archive.members().size(); k>0; --k) {
    const auto& member=archive.members()[k-1];
    Reader reader(archive.image(k-1));
    std::vector<Field> fields,constants;
    std::vector<double> times;
    while(reader.readRow(fields,constants)) {
      bool found=false;
      for(const auto& field:fields) {
        if(field.name=="time") {
          double time=0;
          if(field.value.type==Value::Double) {
            time=field.value.d;
          } else if(field.value.type==Value::Int) {
            time=static_cast<double>(field.value.i);
          } else {
            check(Tools::convertNoexcept(field.value.s,time),path,"invalid time field in batch");
          }
          times.push_back(time);
          found=true;
          hasTime=true;
          break;
        }
      }
      if(!found) {
        hasTime=false;
        return {};
      }
    }
    check(times.size()==member.rows,path,"batch row count disagrees with TAR header");
    for(std::size_t j=times.size(); j>0; --j) {
      const double time=times[j-1];
      if(time<cut) {
        return keep;
      }
      if(time>next) {
        increasing=false;
      }
      next=time;
      --keep[member.block];
    }
    while(!keep.empty() && keep.back()==0) {
      keep.pop_back();
    }
  }
  return keep;
}

void truncate(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
  Lock lock(hdf5Mutex());
  if(Archive::recognizes(path)) {
    Archive archive(path,true,true);
    std::vector<unsigned long long> current;
    for(const auto& member:archive.members()) {
      if(member.block>=current.size()) {
        current.push_back(0);
      }
      current.back()+=member.rows;
    }
    check(keep.size()<=current.size(),path,"fewer blocks than expected");
    for(std::size_t b=0; b<keep.size(); ++b) {
      check(keep[b]<=current[b],path,"fewer rows than expected");
      if(b+1<keep.size()) {
        check(keep[b]==current[b],path,"archive trimming requires a retained prefix");
      }
    }
    if(keep==current) {
      return;
    }
    std::vector<unsigned long long> before(current.size(),0);
    std::size_t first=0;
    unsigned long long retained=0;
    for(; first<archive.members().size(); ++first) {
      const auto& member=archive.members()[first];
      const auto wanted=member.block<keep.size()?keep[member.block]:0;
      if(before[member.block]+member.rows>wanted) {
        retained=wanted-before[member.block];
        break;
      }
      before[member.block]+=member.rows;
    }
    check(first<archive.members().size(),path,"cannot locate trim boundary");
    std::vector<char> replacement;
    std::unique_ptr<Archive> saved;
    if(!backup.empty()) {
      saved=std::make_unique<Archive>(backup,true,false);
    }
    const auto cutBlock=archive.members()[first].block;
    if(retained) {
      BatchFile keptFile,removedFile;
      {
        Writer kept(keptFile.path,false);
        std::unique_ptr<Writer> removed;
        if(saved) {
          removed=std::make_unique<Writer>(removedFile.path,false);
        }
        Reader reader(archive.image(first));
        std::vector<Field> fields,constants;
        unsigned long long row=0;
        while(reader.readRow(fields,constants)) {
          if(row++<retained) {
            kept.addRow(fields,constants);
          } else if(removed) {
            removed->addRow(fields,constants);
          }
        }
        check(row==archive.members()[first].rows,path,"batch row count disagrees with TAR header");
        kept.flush();
        if(removed) {
          removed->flush();
        }
      }
      replacement=keptFile.image();
      if(saved) {
        saved->add(removedFile.image(),0,archive.members()[first].rows-retained);
      }
    }
    if(saved) {
      for(std::size_t k=first+(retained?1:0); k<archive.members().size(); ++k) {
        const auto& member=archive.members()[k];
        saved->add(archive.image(k),member.block-cutBlock,member.rows);
      }
      saved->finish(); // Preserve the discarded tail before changing the source.
    }
    archive.replaceTail(first,replacement,cutBlock,retained);
    archive.finish();
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

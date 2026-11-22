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
#include "FileTools.h"
#include "Exception.h"
#include "Tools.h"
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <mutex>
#include <utility>
#include <cerrno>
#include <cstring>
#include <random>

#ifdef __PLUMED_HAS_LIBHDF5
#include <hdf5.h>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/stat.h>
#ifdef __APPLE__
#include <sys/clonefile.h>
#elif defined(__linux__)
#include <sys/ioctl.h>
#include <linux/fs.h>
#endif
#endif
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

bool atomicOutput() {
  std::string value;
  if(!getenvUTF8("PLUMED_H5_ATOMIC",value) || value=="0") {
    return false;
  }
  plumed_massert(value=="1","PLUMED_H5_ATOMIC must be 0 or 1");
  return true;
}

/// OS locks are independent of HDF5 handles: published snapshots are protected
/// against in-place writers while remaining readable by ordinary HDF5 clients.
class FileLock {
#ifdef _WIN32
  HANDLE handle=INVALID_HANDLE_VALUE;
#else
  int descriptor=-1;
#endif
public:
  FileLock(const std::string& path,bool exclusive,bool create=false,bool writable=false) {
#ifdef _WIN32
    handle=CreateFileW(pathFromUTF8(path).c_str(),
                       GENERIC_READ|((create || writable)?GENERIC_WRITE:0),
                       FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,nullptr,
                       create?OPEN_ALWAYS:OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    check(handle!=INVALID_HANDLE_VALUE,path,"cannot open publication lock");
    OVERLAPPED position= {};
    if(!LockFileEx(handle,LOCKFILE_FAIL_IMMEDIATELY|(exclusive?LOCKFILE_EXCLUSIVE_LOCK:0),
                   0,1,0,&position)) {
      CloseHandle(handle);
      handle=INVALID_HANDLE_VALUE;
      check(false,path,"file is already open for writing or cannot be locked");
    }
#else
    int flags=(create?O_CREAT:0)|((create || writable)?O_RDWR:O_RDONLY);
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
    descriptor=::open(path.c_str(),flags,0666);
    check(descriptor>=0,path,"cannot open publication lock: "+std::string(std::strerror(errno)));
    if(flock(descriptor,(exclusive?LOCK_EX:LOCK_SH)|LOCK_NB)!=0) {
      ::close(descriptor);
      descriptor=-1;
      check(false,path,"file is already open for writing or cannot be locked");
    }
#endif
  }
  FileLock(const FileLock&)=delete;
  FileLock& operator=(const FileLock&)=delete;
  ~FileLock() {
#ifdef _WIN32
    if(handle!=INVALID_HANDLE_VALUE) {
      CloseHandle(handle);
    }
#else
    if(descriptor>=0) {
      ::close(descriptor);
    }
#endif
  }
  void sync(const std::string& path) {
#ifdef _WIN32
    check(FlushFileBuffers(handle)!=0,path,"cannot synchronize snapshot data");
#else
    check(fsync(descriptor)==0,path,"cannot synchronize snapshot data");
#ifdef F_FULLFSYNC
    if(fcntl(descriptor,F_FULLFSYNC)!=0) {
      check(errno==EINVAL || errno==ENOTSUP,path,"cannot synchronize device buffers");
    }
#endif
#endif
  }
};

/// A private work file and an immutable, atomically published v1 table.
class Snapshot {
  std::string path;
  std::string pending;
  std::unique_ptr<FileLock> owner;
  std::unique_ptr<FileLock> published;
public:
  explicit Snapshot(const std::string& original) {
    std::error_code ec;
    const auto resolved=std::filesystem::weakly_canonical(pathFromUTF8(original),ec);
    check(!ec,original,"cannot resolve snapshot path");
    path=pathToUTF8(resolved);
    owner=std::make_unique<FileLock>(path+".plumed-lock",true,true);
    const bool exists=std::filesystem::exists(resolved,ec);
    check(!ec,original,"cannot inspect snapshot path: "+ec.message());
    if(exists) {
      check(std::filesystem::is_regular_file(resolved,ec) && !ec,original,
            "atomic output requires a regular file");
    }
    if(exists) {
      published=std::make_unique<FileLock>(path,false);
    }
  }
  ~Snapshot() {
    if(!pending.empty()) {
      std::error_code ec;
      std::filesystem::remove(pathFromUTF8(pending),ec);
    }
  }
  std::pair<std::string,bool> prepare(bool append) {
    plumed_assert(pending.empty());
    const auto target=pathFromUTF8(path);
    const auto directory=target.parent_path();
    // Reserve a unique name with exclusive creation; never reuse a work file
    // left by a crashed process.
    std::random_device random;
    for(unsigned attempt=0; attempt<128; ++attempt) {
      const auto temporary=directory/pathFromUTF8("."+pathToUTF8(target.filename())+".pending."+
                           std::to_string(random())+"."+std::to_string(random()));
      const std::string candidate=pathToUTF8(temporary);
#ifdef _WIN32
      HANDLE file=CreateFileW(temporary.c_str(),GENERIC_READ|GENERIC_WRITE,
                              FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,
                              nullptr,CREATE_NEW,FILE_ATTRIBUTE_NORMAL,nullptr);
      if(file!=INVALID_HANDLE_VALUE) {
        CloseHandle(file);
        pending=candidate;
        break;
      }
      check(GetLastError()==ERROR_FILE_EXISTS,path,"cannot create private snapshot");
#else
      const int file=::open(candidate.c_str(),O_CREAT|O_EXCL|O_RDWR,0666);
      if(file>=0) {
        ::close(file);
        pending=candidate;
        break;
      }
      check(errno==EEXIST,path,"cannot create private snapshot");
#endif
    }
    check(!pending.empty(),path,"cannot reserve a private snapshot name");
    std::error_code ec;
    const bool copy=append && std::filesystem::exists(target,ec);
    check(!ec,path,"cannot inspect published snapshot: "+ec.message());
    if(copy) {
      const auto permissions=std::filesystem::status(target,ec).permissions();
      check(!ec,path,"cannot inspect published snapshot permissions");
      std::filesystem::permissions(pathFromUTF8(pending),permissions,ec);
      check(!ec,path,"cannot preserve snapshot permissions");
      bool cloned=false;
#ifdef __APPLE__
      std::filesystem::remove(pathFromUTF8(pending));
      cloned=clonefile(path.c_str(),pending.c_str(),0)==0;
#elif defined(__linux__) && defined(FICLONE)
      const int source=::open(path.c_str(),O_RDONLY);
      const int destination=::open(pending.c_str(),O_RDWR);
      if(source>=0 && destination>=0) {
        cloned=ioctl(destination,FICLONE,source)==0;
      }
      if(source>=0) {
        ::close(source);
      }
      if(destination>=0) {
        ::close(destination);
      }
#endif
      if(!cloned) {
        std::filesystem::copy_file(target,pathFromUTF8(pending),
                                   std::filesystem::copy_options::overwrite_existing,ec);
        check(!ec,path,"cannot copy the published snapshot: "+ec.message());
      }
    }
    return {pending,copy};
  }
  void publish() {
    plumed_assert(!pending.empty());
    auto next=std::make_unique<FileLock>(pending,false,false,true);
    next->sync(path);
#ifdef _WIN32
    const auto source=pathFromUTF8(pending);
    const auto destination=pathFromUTF8(path).native();
    HANDLE file=CreateFileW(source.c_str(),DELETE|GENERIC_READ|GENERIC_WRITE,
                            FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,nullptr,
                            OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    check(file!=INVALID_HANDLE_VALUE,path,"cannot open snapshot for publication");
    const DWORD bytes=static_cast<DWORD>(destination.size()*sizeof(wchar_t));
    std::vector<char> storage(sizeof(FILE_RENAME_INFO)+bytes);
    auto* rename=reinterpret_cast<FILE_RENAME_INFO*>(storage.data());
    // FileRenameInfoEx + POSIX semantics preserve existing readers of the old
    // inode, including HDF5 clients which do not request delete sharing.
    const DWORD flags=0x00000001|0x00000002;
    std::memcpy(storage.data(),&flags,sizeof(flags));
    rename->RootDirectory=nullptr;
    rename->FileNameLength=bytes;
    std::memcpy(rename->FileName,destination.data(),bytes);
    const bool renamed=SetFileInformationByHandle(file,static_cast<FILE_INFO_BY_HANDLE_CLASS>(22),
                       rename,static_cast<DWORD>(storage.size()))!=0;
    const bool synced=renamed && FlushFileBuffers(file)!=0;
    CloseHandle(file);
    check(renamed,path,"cannot atomically publish snapshot (POSIX rename semantics required)");
    pending.clear();
    published=std::move(next);
    check(synced,path,"cannot synchronize published snapshot");
#else
    check(renameUTF8(pending,path)==0,path,"cannot atomically publish snapshot");
    pending.clear();
    published=std::move(next);
    const std::string directory=pathToUTF8(pathFromUTF8(path).parent_path());
    const int descriptor=::open(directory.c_str(),O_RDONLY);
    check(descriptor>=0,path,"cannot open snapshot directory for synchronization");
    const bool synced=fsync(descriptor)==0;
    ::close(descriptor);
    check(synced,path,"cannot synchronize snapshot directory");
#endif
  }
};

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
  bool atomic=false;
  bool append=false;
  bool committed=false;
  std::unique_ptr<Snapshot> snapshot;
  std::unique_ptr<Writer> staging;
  void prepareSnapshot() {
    if(!staging) {
      const auto work=snapshot->prepare(append || committed);
      staging.reset(new Writer(work.first,work.second,false));
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

Writer::Writer(const std::string& path,bool append):
  Writer(path,append,atomicOutput()) {
}

Writer::Writer(const std::string& path,bool append,bool atomic):
  impl(std::make_unique<Impl>()) {
  Lock lock(hdf5Mutex());
  impl->path=path;
  if(atomic) {
    impl->atomic=true;
    impl->append=append;
    impl->snapshot=std::make_unique<Snapshot>(path);
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
    if(impl->atomic) {
      flush();
    } else {
      impl->writeBuffer();
    }
  } catch(...) {
  }
}

void Writer::addRow(const std::vector<Field>& fields,const std::vector<Field>& constants) {
  Lock lock(hdf5Mutex());
  if(impl->atomic) {
    impl->prepareSnapshot();
    impl->staging->addRow(fields,constants);
    if(++impl->buffered>=4096) {
      flush();
    }
    return;
  }
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
  if(impl->atomic) {
    if(!impl->staging && impl->committed) {
      return;
    }
    impl->prepareSnapshot();
    impl->staging->flush();
    impl->staging.reset();
    impl->snapshot->publish();
    impl->committed=true;
    impl->buffered=0;
    return;
  }
  impl->writeBuffer();
  check(H5Fflush(impl->file,H5F_SCOPE_GLOBAL)>=0,impl->path,"cannot flush the file");
}

std::vector<unsigned long long> Writer::rows() {
  Lock lock(hdf5Mutex());
  flush();
  if(impl->atomic) {
    return h5table::rows(impl->path);
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

namespace {
void truncateInPlace(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
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
    check(H5Fflush(f,H5F_SCOPE_GLOBAL)>=0,path,"cannot flush the trimmed file");
  }
  if(out) {
    out->flush();
  }
}
}

void truncate(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup) {
  Lock lock(hdf5Mutex());
  if(atomicOutput()) {
    Snapshot snapshot(path);
    const auto work=snapshot.prepare(true);
    check(work.second,path,"cannot trim a missing table");
    truncateInPlace(work.first,keep,backup);
    snapshot.publish();
  } else {
    truncateInPlace(path,keep,backup);
  }
}

#endif

}
}

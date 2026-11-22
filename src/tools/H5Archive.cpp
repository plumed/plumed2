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
#include "H5Archive.h"
#include "Exception.h"
#include "FileTools.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <random>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace PLMD {
namespace h5table {
namespace {
void require(bool ok,const std::string& path,const std::string& message) {
  plumed_massert(ok,"HDF5 archive "+path+": "+message);
}

// Unbuffered OS I/O: a completed member is synchronized before another starts.
// Readers share writes/deletion on Windows and never acquire the writer lock.
class File {
  std::string path;
#ifdef _WIN32
  HANDLE handle=INVALID_HANDLE_VALUE;
#else
  int descriptor=-1;
#endif
public:
  File(const std::string& name,bool writable,bool create=false,bool exclusive=false):path(name) {
#ifdef _WIN32
    handle=CreateFileW(pathFromUTF8(path).c_str(),GENERIC_READ|(writable?GENERIC_WRITE:0),
                       FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,nullptr,
                       create?OPEN_ALWAYS:OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    require(handle!=INVALID_HANDLE_VALUE,path,"cannot open file");
    if(exclusive) {
      OVERLAPPED position= {};
      if(!LockFileEx(handle,LOCKFILE_EXCLUSIVE_LOCK|LOCKFILE_FAIL_IMMEDIATELY,0,1,0,&position)) {
        CloseHandle(handle);
        handle=INVALID_HANDLE_VALUE;
        require(false,path,"another writer owns this archive");
      }
    }
#else
    descriptor=::open(path.c_str(),(writable?O_RDWR:O_RDONLY)|(create?O_CREAT:0),0666);
    require(descriptor>=0,path,"cannot open file: "+std::string(std::strerror(errno)));
    if(exclusive && flock(descriptor,LOCK_EX|LOCK_NB)!=0) {
      ::close(descriptor);
      descriptor=-1;
      require(false,path,"another writer owns this archive");
    }
#endif
  }
  ~File() {
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
  void seek(std::uint64_t offset) const {
#ifdef _WIN32
    LARGE_INTEGER position;
    position.QuadPart=static_cast<LONGLONG>(offset);
    require(offset<=static_cast<std::uint64_t>(std::numeric_limits<LONGLONG>::max()) &&
            SetFilePointerEx(handle,position,nullptr,FILE_BEGIN),path,"cannot seek");
#else
    require(offset<=static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()) &&
            lseek(descriptor,static_cast<off_t>(offset),SEEK_SET)>=0,path,"cannot seek");
#endif
  }
  std::uint64_t size() const {
#ifdef _WIN32
    LARGE_INTEGER n;
    require(GetFileSizeEx(handle,&n),path,"cannot obtain size");
    return static_cast<std::uint64_t>(n.QuadPart);
#else
    struct stat st;
    require(fstat(descriptor,&st)==0,path,"cannot obtain size");
    return static_cast<std::uint64_t>(st.st_size);
#endif
  }
  void read(char* data,std::size_t count) const {
    while(count) {
#ifdef _WIN32
      DWORD n=0;
      require(ReadFile(handle,data,static_cast<DWORD>(std::min<std::size_t>(count,1U<<30)),&n,nullptr) && n,path,"cannot read member");
#else
      const ssize_t n=::read(descriptor,data,std::min<std::size_t>(count,1U<<30));
      if(n<0 && errno==EINTR) {
        continue;
      }
      require(n>0,path,"cannot read member");
#endif
      data+=n;
      count-=n;
    }
  }
  void write(const char* data,std::size_t count) {
    while(count) {
#ifdef _WIN32
      DWORD n=0;
      require(WriteFile(handle,data,static_cast<DWORD>(std::min<std::size_t>(count,1U<<30)),&n,nullptr) && n,path,"cannot write member");
#else
      const ssize_t n=::write(descriptor,data,std::min<std::size_t>(count,1U<<30));
      if(n<0 && errno==EINTR) {
        continue;
      }
      require(n>0,path,"cannot write member: "+std::string(std::strerror(errno)));
#endif
      data+=n;
      count-=n;
    }
  }
  void truncate(std::uint64_t offset) {
    seek(offset);
#ifdef _WIN32
    require(SetEndOfFile(handle),path,"cannot remove incomplete tail");
#else
    require(ftruncate(descriptor,static_cast<off_t>(offset))==0,path,"cannot remove incomplete tail");
#endif
  }
  void sync() {
#ifdef _WIN32
    require(FlushFileBuffers(handle),path,"cannot synchronize file");
#else
    require(fsync(descriptor)==0,path,"cannot synchronize file");
#ifdef F_FULLFSYNC
    if(fcntl(descriptor,F_FULLFSYNC)!=0) {
      require(errno==EINVAL || errno==ENOTSUP,path,"cannot synchronize device buffers");
    }
#endif
#endif
  }
};
void syncDirectory(const std::string& path) {
#ifndef _WIN32
  auto parent=pathFromUTF8(path).parent_path();
  if(parent.empty()) {
    parent=".";
  }
  File directory(pathToUTF8(parent),false);
  directory.sync();
#else
  (void)path; // Windows publication uses a synchronized handle / WRITE_THROUGH.
#endif
}
std::uint64_t octal(const char* data,std::size_t size,const std::string& path) {
  std::uint64_t value=0;
  bool ended=false;
  for(std::size_t i=0; i<size; ++i) {
    const char c=data[i];
    if(c==' ' || c=='\0') {
      if(value) {
        ended=true;
      }
    } else {
      require(!ended && c>='0' && c<='7' && value<=(std::numeric_limits<std::uint64_t>::max()-7)/8,path,"invalid TAR number");
      value=value*8+static_cast<unsigned>(c-'0');
    }
  }
  return value;
}
std::string nameOf(const char* data,std::size_t length) {
  return std::string(data,std::find(data,data+length,'\0'));
}
unsigned checksum(const std::array<char,512>& header) {
  unsigned sum=0;
  for(unsigned i=0; i<header.size(); ++i) {
    sum+=i>=148 && i<156 ? 32 : static_cast<unsigned char>(header[i]);
  }
  return sum;
}
bool zero(const std::array<char,512>& header) {
  return std::all_of(header.begin(),header.end(),[](char c) {
    return c==0;
  });
}
}

struct Archive::Impl {
  std::string path;
  std::unique_ptr<File> owner;
  std::unique_ptr<File> file;
  std::unique_ptr<File> recovery;
  std::vector<Member> members;
  std::uint64_t end=0;
  bool writing=false,finished=false,failed=false;
};

bool Archive::recognizes(const std::string& requested) {
  std::error_code error;
  const auto canonical=std::filesystem::weakly_canonical(pathFromUTF8(requested),error);
  if(error) {
    return false;
  }
  const auto path=pathToUTF8(canonical);
  unique_FILE file(fopenUTF8(path,"rb"));
  if(!file) {
    return false;
  }
  std::array<char,512> header= {};
  const auto n=std::fread(header.data(),1,header.size(),file.get());
  if(n==512 && std::memcmp(header.data()+257,"ustar",5)==0) {
    return true;
  }
  if(n==512 && zero(header)) {
    std::array<char,512> footer= {};
    if(std::fread(footer.data(),1,footer.size(),file.get())==512 && zero(footer)) {
      return true; // Standard empty TAR, including one trimmed to zero rows.
    }
  }
  unique_FILE recovery(fopenUTF8(path+".plumed-trim","rb"));
  std::array<char,18> magic= {};
  return recovery && std::fread(magic.data(),1,17,recovery.get())==17 &&
         std::memcmp(magic.data(),"plumed-h5-trim-v1 ",17)==0;
}

Archive::Archive(const std::string& requested,bool writing,bool append):impl(std::make_unique<Impl>()) {
  const auto path=pathToUTF8(std::filesystem::weakly_canonical(pathFromUTF8(requested)));
  impl->path=path;
  impl->writing=writing;
  if(writing) {
    impl->owner=std::make_unique<File>(path+".plumed-lock",true,true,true);
  }
  impl->file=std::make_unique<File>(path,writing,writing);
  if(writing && !append) {
    impl->file->truncate(0);
  }
  auto length=impl->file->size(); // A live reader sees a finite prefix.
  std::uint64_t recoveryCount=0,recoverySize=0,recoveryBlock=0,recoveryRows=0;
  const auto recoveryPath=path+".plumed-trim";
  std::error_code error;
  if(std::filesystem::exists(pathFromUTF8(recoveryPath),error) && (!writing || append)) {
    impl->recovery=std::make_unique<File>(recoveryPath,false);
    require(impl->recovery->size()>=512,path,"incomplete trim recovery record");
    std::array<char,513> header= {};
    impl->recovery->read(header.data(),512);
    unsigned long long prefix=0,count=0,size=0,block=0,rows=0;
    require(std::sscanf(header.data(),"plumed-h5-trim-v1 %llu %llu %llu %llu %llu",&prefix,&count,&size,&block,&rows)==5,path,"invalid trim recovery record");
    require(prefix<=length && prefix%4096==0 && size==impl->recovery->size()-512,path,"inconsistent trim recovery record");
    length=prefix;
    recoveryCount=count;
    recoverySize=size;
    recoveryBlock=block;
    recoveryRows=rows;
  } else if(writing && !append) {
    std::filesystem::remove(pathFromUTF8(recoveryPath),error);
  }
  bool incomplete=false;
  std::uint64_t cursor=0;
  while(length-cursor>=512) {
    impl->file->seek(cursor);
    std::array<char,512> header;
    impl->file->read(header.data(),header.size());
    if(zero(header)) {
      break;
    }
    bool valid=false;
    try {
      valid=std::memcmp(header.data()+257,"ustar",5)==0 &&
            octal(header.data()+148,8,path)==checksum(header);
    } catch(const std::exception&) {
      // A torn publication header must not hide earlier committed members.
    }
    if(!valid) {
      require(!impl->members.empty(),path,"not a supported USTAR HDF5 archive: invalid first header/checksum");
      std::fprintf(stderr,"PLUMED: ignoring damaged/unpublished HDF5 archive tail in %s\n",path.c_str());
      break;
    }
    require(header[156]=='0' || header[156]=='\0',path,"only regular numbered HDF5 members are supported");
    require(nameOf(header.data()+265,32)=="plumed-h5-v1",path,"unsupported PLUMED archive version/profile");
    require(nameOf(header.data(),100)=="batch_"+std::to_string(impl->members.size()+1)+".h5",path,"members must be consecutively numbered batch_N.h5");
    require(header[345]=='\0',path,"TAR prefixes are not supported");
    const auto size=octal(header.data()+124,12,path);
    const auto rounded=(size+511)/512*512;
    require(size && rounded>=size,path,"invalid TAR member size");
    if(rounded>length-cursor-512) {
      incomplete=true;
      break; // Last header/data/padding has not been completely published.
    }
    unsigned long long block=0,rows=0;
    const auto metadata=nameOf(header.data()+297,32);
    char trailing=0;
    require(std::sscanf(metadata.c_str(),"%llu:%llu%c",&block,&rows,&trailing)==2 && rows>0,path,"invalid batch block/row metadata");
    require((impl->members.empty() && block==0) || (!impl->members.empty() &&
            (block==impl->members.back().block || block==impl->members.back().block+1)),path,"invalid logical block sequence");
    impl->members.push_back(Member{cursor+512,size,false,block,rows});
    cursor+=512+rounded;
    impl->end=cursor;
  }
  if(incomplete || (cursor<length && length-cursor<512)) {
    std::fprintf(stderr,"PLUMED: ignoring incomplete HDF5 archive tail in %s\n",path.c_str());
  }
  if(impl->recovery) {
    require(impl->end==length && impl->members.size()==recoveryCount,path,"trim prefix does not match recovery record");
    if(writing) {
      auto image=std::vector<char>(static_cast<std::size_t>(recoverySize));
      impl->recovery->seek(512);
      impl->recovery->read(image.data(),image.size());
      impl->file->truncate(impl->end);
      if(!image.empty()) {
        add(image,recoveryBlock,recoveryRows);
      }
      impl->file->sync();
      std::filesystem::remove(pathFromUTF8(recoveryPath),error);
      require(!error,path,"cannot remove completed trim recovery record");
      syncDirectory(path);
      impl->recovery.reset();
    } else if(recoverySize) {
      impl->members.push_back(Member{512,recoverySize,true,recoveryBlock,recoveryRows});
    }
  }
  if(writing) {
    require(impl->end%4096==0,path,"cannot append an archive without aligned batch boundaries");
    impl->file->truncate(impl->end);
  }
}
Archive::~Archive() {
  if(impl->writing && !impl->failed) {
    try {
      finish();
    } catch(...) {
    }
  }
}
const std::vector<Archive::Member>& Archive::members() const {
  return impl->members;
}
std::vector<char> Archive::image(std::size_t index) const {
  const auto member=impl->members.at(index);
  require(member.size<=std::numeric_limits<std::size_t>::max(),impl->path,"member too large for memory");
  std::vector<char> result(static_cast<std::size_t>(member.size));
  const auto& file=member.recovery?impl->recovery:impl->file;
  file->seek(member.offset);
  file->read(result.data(),result.size());
  return result;
}
void Archive::add(const std::vector<char>& image,std::uint64_t block,std::uint64_t rows) {
  require(impl->writing && !impl->finished && !impl->failed,impl->path,"archive is not writable");
  try {
    // Align the end to 4096 bytes. The extra zeros belong to the HDF5 member:
    // they are legal trailing bytes, so ordinary TAR extraction remains valid.
    const std::uint64_t size=(static_cast<std::uint64_t>(image.size())+512+4095)/4096*4096-512;
    require(size>=image.size() && size<=077777777777ULL,impl->path,"batch exceeds USTAR size limit; use smaller batches");
    require(impl->members.size()<std::numeric_limits<std::uint64_t>::max(),impl->path,"batch sequence overflow");
    require(rows>0 && (impl->members.empty()?block==0:(block==impl->members.back().block || block==impl->members.back().block+1)),impl->path,"invalid batch block/row metadata");
    const auto metadata=std::to_string(block)+":"+std::to_string(rows);
    require(metadata.size()<32,impl->path,"batch metadata exceeds USTAR group-name field");
    std::array<char,512> header= {};
    const auto name="batch_"+std::to_string(impl->members.size()+1)+".h5";
    std::memcpy(header.data(),name.data(),name.size());
    std::snprintf(header.data()+100,8,"%07o",0644);
    std::snprintf(header.data()+108,8,"%07o",0);
    std::snprintf(header.data()+116,8,"%07o",0);
    std::snprintf(header.data()+124,12,"%011llo",static_cast<unsigned long long>(size));
    std::snprintf(header.data()+136,12,"%011o",0);
    header[156]='0';
    std::memcpy(header.data()+257,"ustar",5);
    std::memcpy(header.data()+263,"00",2);
    std::memcpy(header.data()+265,"plumed-h5-v1",12);
    std::memcpy(header.data()+297,metadata.data(),metadata.size());
    std::snprintf(header.data()+148,7,"%06o",checksum(header));
    header[155]=' ';
    impl->file->seek(impl->end);
    const std::array<char,512> unpublished= {};
    impl->file->write(unpublished.data(),unpublished.size());
    impl->file->write(image.data(),image.size());
    std::array<char,4096> padding= {};
    impl->file->write(padding.data(),static_cast<std::size_t>(size-image.size()));
    impl->file->sync(); // The payload is durable before its header publishes it.
    impl->file->seek(impl->end);
    impl->file->write(header.data(),header.size());
    impl->file->sync();
    if(impl->members.empty()) {
      syncDirectory(impl->path);
    }
    impl->members.push_back(Member{impl->end+512,size,false,block,rows});
    impl->end+=512+size;
  } catch(...) {
    impl->failed=true; // Never finish or append over a failed publication.
    throw;
  }
}
void Archive::replaceTail(std::size_t first,const std::vector<char>& image,std::uint64_t block,std::uint64_t rows) {
  require(impl->writing && !impl->failed && !impl->finished && first<=impl->members.size(),impl->path,"invalid tail replacement");
  const std::uint64_t offset=first?impl->members[first-1].offset+(impl->members[first-1].size+511)/512*512:0;
  require(offset%4096==0,impl->path,"trim requires aligned batches");
  auto parent=pathFromUTF8(impl->path).parent_path();
  BatchFile pending(parent.empty()?".":pathToUTF8(parent));
  {
    File record(pending.path,true,true);
    std::array<char,512> header= {};
    std::snprintf(header.data(),header.size(),"plumed-h5-trim-v1 %llu %llu %llu %llu %llu\n",
                  static_cast<unsigned long long>(offset),static_cast<unsigned long long>(first),
                  static_cast<unsigned long long>(image.size()),static_cast<unsigned long long>(block),static_cast<unsigned long long>(rows));
    record.write(header.data(),header.size());
    record.write(image.data(),image.size());
    record.sync();
  }
  const auto recoveryPath=impl->path+".plumed-trim";
  std::error_code permissionsError;
  const auto permissions=std::filesystem::status(pathFromUTF8(impl->path),permissionsError).permissions();
  require(!permissionsError,impl->path,"cannot inspect archive permissions");
  std::filesystem::permissions(pathFromUTF8(pending.path),permissions,permissionsError);
  require(!permissionsError,impl->path,"cannot protect trim recovery record");
  pending.publish(recoveryPath); // Durable intent before touching the old tail.
  try {
    impl->members.resize(first);
    impl->end=offset;
    impl->file->truncate(offset);
    if(!image.empty()) {
      add(image,block,rows);
    }
    impl->file->sync();
    std::error_code error;
    std::filesystem::remove(pathFromUTF8(recoveryPath),error);
    require(!error,impl->path,"cannot remove completed trim recovery record");
    syncDirectory(impl->path);
  } catch(...) {
    impl->failed=true; // Recovery record must remain for the next reader/writer.
    throw;
  }
}

void Archive::finish() {
  if(impl->writing && !impl->finished) {
    require(!impl->failed,impl->path,"previous publication failed");
    const std::array<char,1024> footer= {};
    impl->file->seek(impl->end);
    impl->file->write(footer.data(),footer.size());
    impl->file->sync();
    syncDirectory(impl->path);
    impl->finished=true;
  }
}

BatchFile::BatchFile(const std::string& parent) {
  const auto base=parent.empty()?std::filesystem::temp_directory_path():pathFromUTF8(parent);
  std::random_device random;
  for(unsigned attempt=0; attempt<100; ++attempt) {
    const auto candidate=base/("plumed-h5-"+std::to_string(random())+"-"+std::to_string(random()));
    std::error_code error;
    if(std::filesystem::create_directory(candidate,error)) {
      directory=pathToUTF8(candidate);
#ifndef _WIN32
      if(chmod(candidate.c_str(),0700)!=0) {
        std::filesystem::remove(candidate,error);
        require(false,directory,"cannot protect staging directory");
      }
#endif
      path=pathToUTF8(candidate/"batch.h5");
      return;
    }
    require(!error || error==std::errc::file_exists,pathToUTF8(base),"cannot create batch staging directory: "+error.message());
  }
  require(false,pathToUTF8(base),"cannot allocate private batch staging directory");
}
BatchFile::~BatchFile() {
  std::error_code error;
  std::filesystem::remove_all(pathFromUTF8(directory),error);
}
std::vector<char> BatchFile::image() const {
  File file(path,false);
  const auto size=file.size();
  require(size<=std::numeric_limits<std::size_t>::max(),path,"batch too large");
  std::vector<char> result(static_cast<std::size_t>(size));
  file.read(result.data(),result.size());
  return result;
}
void BatchFile::publish(const std::string& requested) {
  const auto destination=pathToUTF8(std::filesystem::weakly_canonical(pathFromUTF8(requested)));
  std::unique_ptr<File> writerLock;
  if(destination.size()>=7 && destination.compare(destination.size()-7,7,".h5.tar")==0) {
    writerLock=std::make_unique<File>(destination+".plumed-lock",true,true,true);
    require(!std::filesystem::exists(pathFromUTF8(destination+".plumed-trim")),destination,
            "recover the pending trim before replacing this archive");
  }
  std::error_code error;
  const auto status=std::filesystem::status(pathFromUTF8(destination),error);
  if(!error && std::filesystem::exists(status)) {
    require(std::filesystem::is_regular_file(status),destination,"destination must be a regular file");
    std::filesystem::permissions(pathFromUTF8(path),status.permissions(),error);
    require(!error,path,"cannot preserve destination permissions");
  }
  {
    File file(path,true);
    file.sync();
  }
#ifdef _WIN32
  require(MoveFileExW(pathFromUTF8(path).c_str(),pathFromUTF8(destination).c_str(),
                      MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH),destination,"cannot publish converted/trimmed archive");
#else
  require(renameUTF8(path,destination)==0,destination,"cannot publish converted/trimmed archive");
  syncDirectory(destination);
#endif
}
}
}

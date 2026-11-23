#include "plumed/tools/H5Table.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#ifndef _WIN32
#include <csignal>
#include <sys/resource.h>
#include <unistd.h>
#endif

using namespace PLMD::h5table;

void require(bool ok,const std::string& message) {
  if(!ok) {
    throw std::runtime_error(message);
  }
}
template<class F> void rejects(F action) {
  bool failed=false;
  try {
    action();
  } catch(const std::exception&) {
    failed=true;
  }
  require(failed,"operation unexpectedly succeeded");
}
void row(Writer& writer,unsigned i) {
  Value value;
  value.type=Value::Double;
  value.d=i;
  writer.addRow({{"time",value}}, {});
}
std::string bytes(const std::string& path) {
  std::ifstream file(path,std::ios::binary);
  return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}
unsigned count(const std::string& path) {
  Reader reader(path);
  std::string line;
  unsigned n=0;
  while(reader.getline(line)) {
    if(!line.empty() && line[0]!='#') {
      n++;
    }
  }
  return n;
}
void create(const std::string& path,unsigned n) {
  Writer writer(path,false,10);
  for(unsigned i=0; i<n; i++) {
    row(writer,i);
  }
  writer.flush();
}
void ioFailure(const std::string& path) {
#ifndef _WIN32
  {
    Writer writer(path,false);
    for(unsigned i=0; i<10; i++) {
      row(writer,i);
    }
    writer.flush();
    // A real short write, not a mocked return value. Restore the resource limit
    // before retrying: the writer must remember the failure after space returns.
    struct Limit {
      struct rlimit previous;
      using Handler=void (*)(int);
      Handler handler;
      Limit(const std::string& path) {
        require(getrlimit(RLIMIT_FSIZE,&previous)==0,"getrlimit failed");
        handler=std::signal(SIGXFSZ,SIG_IGN);
        struct rlimit limit=previous;
        limit.rlim_cur=std::filesystem::file_size(path)+1024;
        require(setrlimit(RLIMIT_FSIZE,&limit)==0,"setrlimit failed");
      }
      void restore() {
        require(setrlimit(RLIMIT_FSIZE,&previous)==0,"restoring rlimit failed");
      }
      ~Limit() {
        setrlimit(RLIMIT_FSIZE,&previous);
        std::signal(SIGXFSZ,handler);
      }
    } limit(path);
    for(unsigned i=10; i<20; i++) {
      row(writer,i);
    }
    rejects([&] { writer.flush(); });
    limit.restore();
    rejects([&] { writer.flush(); });
    rejects([&] { row(writer,20); });
  }
  require(count(path)==10,"short write lost the published prefix");
  {
    Writer writer(path,true);
    row(writer,20);
    writer.flush();
  }
  require(count(path)==11,"reopening did not recover the published prefix");
#else
  (void)path;
  std::cerr<<"short-write resource-limit probe is POSIX only\n";
#endif
}
void directoryFailure(const std::string& path) {
  const auto parent=std::filesystem::path(path).parent_path();
  std::filesystem::create_directories(parent);
  create(path,30);
  const auto original=bytes(path);
  bool hasTime=false,increasing=false;
#ifndef _WIN32
  if(geteuid()!=0) {
    const auto permissions=std::filesystem::status(parent).permissions();
    // Files can still be opened by name, created, and renamed in this directory,
    // but opening the directory itself for synchronization must fail.
    std::filesystem::permissions(parent,std::filesystem::perms::owner_write|std::filesystem::perms::owner_exec);
    bool failed=false;
    try {
      trimArchive(path,25,"",hasTime,increasing);
    } catch(const std::exception&) {
      failed=true;
    }
    std::filesystem::permissions(parent,permissions);
    require(failed,"trim did not report directory synchronization failure");
    require(bytes(path)==original,"archive changed before its recovery name was durable");
    require(std::filesystem::exists(path+".plumed-trim"),"missing published recovery record");
    require(count(path)==25,"reader did not apply the pending trim record");
    {
      Writer writer(path,true);
      writer.flush();
    }
    require(!std::filesystem::exists(path+".plumed-trim"),"completed trim record was not retired");
    require(count(path)==25,"writer did not replay the pending trim");
    return;
  }
  std::cerr<<"directory permission fault requires an unprivileged POSIX user\n";
#endif
  require(trimArchive(path,25,"",hasTime,increasing)==5,"boundary trim failed");
  require(count(path)==25,"wrong number of retained rows");
}

int main(int argc,char** argv) {
  if(argc==1) {
    return 0;
  }
  try {
    require(argc==3,"expected mode and path");
    const std::string mode=argv[1],path=argv[2];
    if(mode=="create") {
      create(path,30);
    } else if(mode=="create25") {
      create(path,25);
    } else if(mode=="copy") {
      const auto original=bytes(path);
      rejects([&] { copy(path,path); });
      require(bytes(path)==original,"typed copy modified its input");
    } else if(mode=="iofailure") {
      ioFailure(path);
    } else if(mode=="directory") {
      directoryFailure(path);
    } else if(mode=="resume") {
      Writer writer(path,true);
      row(writer,30);
      writer.flush();
    } else if(mode=="read") {
      std::cout<<count(path)<<"\n";
    } else {
      throw std::runtime_error("unknown mode");
    }
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<"\n";
    return 1;
  }
}

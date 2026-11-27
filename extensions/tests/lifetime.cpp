#define __PLUMED_WRAPPER_LINK_RUNTIME 1
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_HAS_DLOPEN 1
#include "plumed/wrapper/Plumed.h"
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef LIFETIME_DIRECT_CLI
#include "plumed/core/CLToolMain.h"
#endif

namespace {
using PLMD::Plumed;
bool kernelCanUnload=false;

void require(bool condition,const char* message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

bool resident(const char* path) {
  // Balance the probe's reference; this must not itself retain the extension.
  void* handle=dlopen(path,RTLD_NOW|RTLD_LOCAL|RTLD_NOLOAD);
  if(handle) {
    dlclose(handle);
  }
  return handle!=nullptr;
}

void* kernelSymbol(const char* path) {
  void* handle=dlopen(path,RTLD_NOW|RTLD_LOCAL);
  require(handle,"Cannot open the kernel for the residency probe");
  void* symbol=dlsym(handle,"plumed_plumedmain_create");
  dlclose(handle);
  require(symbol,"Cannot find the kernel entry point");
  return symbol;
}

void checkKernel(void* symbol) {
  Dl_info info{};
  if(kernelCanUnload) {
    require(!dladdr(symbol,&info),"Extension lifetime prevented kernel unloading");
  }
}

std::unique_ptr<Plumed> object(const char* kernel) {
  auto p=std::make_unique<Plumed>(Plumed::dlopen(kernel));
  require(p->valid(),"Cannot create PLUMED");
  int atoms=0;
  p->cmd("setNatoms",&atoms);
  p->cmd("setLogFile","/dev/null");
  p->cmd("init");
  p->cmd("readInputLine","probe: LIFETIME_PROBE");
  return p;
}

void overlap(const char* kernel,const char* library) {
  auto a=object(kernel), b=object(kernel);
  void* symbol=kernelSymbol(kernel);
  require(resident(library),"Autoload did not open the extension");
  a.reset();
  require(resident(library),"Extension unloaded while another client was alive");
  int present=0;
  b->cmd("checkAction LIFETIME_PROBE",&present);
  require(present,"Surviving client lost the extension registry");
  b.reset();
  require(!resident(library),"Extension retained after the last client");
  checkKernel(symbol);
}
}

int main(int argc,char** argv) try {
  require(argc==4 || argc==6,"Expected mode, kernel, library, [hidden, restored]");
  const std::string mode=argv[1];
  const char* kernel=argv[2];
  const char* library=argv[3];
  // Distinguish the autoload lifetime from native restrictions (macOS TLS,
  // GNU-unique symbols, etc.). Extension unloading is always checked strictly.
  setenv("PLUMED_AUTOLOAD","0",1);
  plumed baseline=plumed_create_dlopen(kernel);
  require(plumed_valid(baseline),"Cannot create the baseline client");
  void* symbol=kernelSymbol(kernel);
  plumed_finalize(baseline);
  Dl_info info{};
  kernelCanUnload=!dladdr(symbol,&info);
  std::printf("Kernel unload baseline: %s\n",kernelCanUnload?"yes":"native retention");
  setenv("PLUMED_AUTOLOAD","1",1);
  if(mode=="failure") {
    require(argc==6,"Failure mode needs the missing library paths");
    // The catalog opens valid modules before encountering the missing final one.
    // Failed PlumedMain construction must close those modules, then permit retry.
    plumed failed=plumed_create_dlopen(kernel);
    const bool valid=plumed_valid(failed);
    plumed_finalize(failed);
    require(!valid,"The missing extension did not fail construction");
    require(!resident(library),"Failed construction retained an earlier extension");
    require(std::rename(argv[4],argv[5])==0,"Cannot restore the missing extension");
    overlap(kernel,library);
  } else if(mode=="overlap") {
    for(int cycle=0; cycle<3; ++cycle) {
      overlap(kernel,library);
    }
  } else if(mode=="threads") {
    // Keep one client alive so every worker shares one autoload lifetime.
    auto anchor=object(kernel);
    symbol=kernelSymbol(kernel);
    std::vector<std::thread> threads;
    for(int i=0; i<4; ++i) {
      threads.emplace_back([&] {
        for(int cycle=0; cycle<5; ++cycle) {
          auto p=object(kernel);
        }
      });
    }
    for(auto& thread:threads) {
      thread.join();
    }
    require(resident(library),"Workers released the live anchor's extension");
    anchor.reset();
    require(!resident(library),"Concurrent clients retained the extension");
    checkKernel(symbol);
  } else if(mode=="disabled") {
    setenv("PLUMED_AUTOLOAD","0",1);
    auto disabled=std::make_unique<Plumed>(Plumed::dlopen(kernel));
    require(disabled->valid(),"Cannot create an autoload-disabled client");
    require(!resident(library),"Autoload opt-out opened the extension");
    setenv("PLUMED_AUTOLOAD","1",1);
    auto enabled=object(kernel);
    symbol=kernelSymbol(kernel);
    enabled.reset();
    require(resident(library),"Shared registry outlived its owning handles");
    disabled.reset();
    require(!resident(library),"Autoload-disabled client retained shared handles");
    checkKernel(symbol);
#ifdef LIFETIME_DIRECT_CLI
  } else if(mode=="cli") {
    PLMD::CLToolMain cltool;
    cltool.cmd("setArgvLine","plumed modules --has lifetime_root -q");
    int result=-1;
    cltool.cmd("run",&result);
    require(result==0,"Direct CLToolMain invocation failed");
    require(!resident(library),"Direct CLI retained the extension");
#endif
  } else {
    throw std::runtime_error("Unknown lifetime test mode");
  }
  std::puts("Module lifetime passed");
} catch(const std::exception& error) {
  std::fprintf(stderr,"%s\n",error.what());
  return 1;
}

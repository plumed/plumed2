#include "plumed/core/ActionRegister.h"
#include "plumed/core/ModuleLoader.h"
#include "plumed/tools/DLLoader.h"
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace PLMD;

struct ContextUnwind:public std::exception {};

void check(bool pamm,bool adjmat) {
  if(actionRegister().check("PAMM")!=pamm || actionRegister().check("CONTACT_MATRIX")!=adjmat) {
    throw std::runtime_error("Incorrect CLI context action visibility");
  }
  DLLoader inherited;
  ModuleLoader::inheritContext(inherited);
  if(ModuleLoader::available("pamm",inherited)!=pamm || ModuleLoader::available("adjmat",inherited)!=adjmat) {
    throw std::runtime_error("Incorrect inherited module context");
  }
}

void nested(const DLLoader& a,const DLLoader& b) {
  check(false,false);
  {
    ModuleLoader::Context outer(a);
    check(true,true);
    try {
      ModuleLoader::Context inner(b);
      check(false,true);
      throw ContextUnwind();
    } catch(const ContextUnwind&) {
      // The outer context must be restored during exception unwinding.
    }
    check(true,true);
    {
      DLLoader empty;
      ModuleLoader::Context inner(empty);
      check(false,false);
    }
    check(true,true);
  }
  check(false,false);
}

int main() {
  DLLoader a,b;
  ModuleLoader::load("pamm",a);
  ModuleLoader::load("adjmat",b);
  // The original loader can expire while another owner retains its libraries.
  // Context metadata must stay valid independently of that stack object.
  std::unique_ptr<ModuleLoader::Context> snapshot;
  {
    DLLoader temporary;
    ModuleLoader::load("pamm",temporary);
    snapshot=std::make_unique<ModuleLoader::Context>(temporary);
  }
  check(true,true);
  snapshot.reset();
  check(false,false);
  nested(a,b);
  std::vector<std::thread> threads;
  for(unsigned i=0; i<4; ++i) {
    threads.emplace_back([&a,&b] {
      for(unsigned repeat=0; repeat<3; ++repeat) {
        nested(a,b);
      }
    });
  }
  for(auto& thread:threads) {
    thread.join();
  }
}

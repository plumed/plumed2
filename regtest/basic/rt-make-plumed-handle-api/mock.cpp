#include <cstring>

namespace {
struct State {
  int api=1;
  int counts[4] {};
};
void* create() {
  return new State;
}
void finalize(void* ptr) {
  delete static_cast<State*>(ptr);
}
void command(void* ptr,const char* key,const void* value) {
  auto& state=*static_cast<State*>(ptr);
  if(std::strcmp(key,"getApiVersion")==0) {
    *static_cast<int*>(const_cast<void*>(value))=state.api;
  } else if(std::strcmp(key,"setMockApi")==0) {
    state.api=*static_cast<const int*>(value);
  } else if(std::strcmp(key,"getMockCounts")==0) {
    std::memcpy(const_cast<void*>(value),state.counts,sizeof(state.counts));
  } else if(std::strcmp(key,"performCalcNoUpdate")==0 && state.api>=4) {
    ++state.counts[0];
  } else if(std::strcmp(key,"performCalcNoForces")==0 && state.api>=9) {
    ++state.counts[1];
  } else if(std::strcmp(key,"setExchangeStep")==0 && state.api>=13) {
    ++state.counts[2];
  } else {
    ++state.counts[3]; // Any unsupported dispatch is a regression.
  }
}
struct Functions {
  void* (*create)();
  void (*cmd)(void*,const char*,const void*);
  void (*finalize)(void*);
};
struct Table {
  int version;
  Functions functions;
};
}
extern "C" {
// A version-one table exercises the loaded-kernel compatibility path without
// allowing the real kernel's dependency symbols to replace the mock functions.
  __attribute__((visibility("default"))) Table plumed_symbol_table= {1,{create,command,finalize}};
}

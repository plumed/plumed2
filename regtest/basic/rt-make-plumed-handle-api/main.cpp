#include "plumed/tools/PlumedHandle.h"
#include <array>
#include <fstream>
#include <stdexcept>

int main() {
  for(const int api : {
        1,3,4,8,9,10,13
      }) {
    auto p=PLMD::PlumedHandle::dlopen("./mock." __PLUMED_SOEXT);
    p.cmd("setMockApi",&api);
    auto expect=[&](auto call,bool allowed) {
      bool success=false;
      try {
        call();
        success=true;
      } catch(const std::exception&) {
      }
      if(success!=allowed) {
        throw std::runtime_error("wrong API availability guard");
      }
    };
    expect([&]() {
      p.setExchangeStep(false);
    },true);
    expect([&]() {
      p.setExchangeStep(true);
    },api>=13);
    expect([&]() {
      p.performCalcNoUpdate(false);
    },api>=4);
    expect([&]() {
      p.performCalcNoUpdate(true);
    },api>=9);
    std::array<int,4> counts{};
    p.cmd("getMockCounts",counts.data());
    if(counts!=std::array<int,4> {{api>=4,api>=9,api>=13 ? 2 : 0,0}}) {
      throw std::runtime_error("unsupported or missing guest command");
    }
  }
  std::ofstream("output")<<"Loaded guests: API 1/3/4/8/9/10/13 guards and dispatch OK\n";
}

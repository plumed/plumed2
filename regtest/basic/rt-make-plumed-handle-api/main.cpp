#include "plumed/tools/PlumedHandle.h"
#include <array>
#include <fstream>
#include <stdexcept>
#include <utility>

int main() {
  for(const int api : {
        1,3,4,8,9,10,13
      }) {
    auto p=PLMD::PlumedHandle::dlopen("./mock." __PLUMED_SOEXT);
    p.cmd("setMockApi",&api);
    int initialQueries=0;
    p.cmd("getMockApiQueries",&initialQueries);
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
    int queries=0;
    p.cmd("getMockApiQueries",&queries);
    if(queries!=initialQueries+1) {
      throw std::runtime_error("loaded guest API version was not cached");
    }
    auto moved=std::move(p);
    moved.setExchangeStep(false);
    auto assigned=PLMD::PlumedHandle::dlopen("./mock." __PLUMED_SOEXT);
    const int otherApi=api>=13 ? 4 : 13;
    assigned.cmd("setMockApi",&otherApi);
    assigned.setExchangeStep(false); // Fill a different guest's cache first.
    assigned=std::move(moved);
    expect([&]() {
      assigned.setExchangeStep(true);
    },api>=13);
    expect([&]() {
      assigned.performCalcNoUpdate(true);
    },api>=9);
    assigned.cmd("getMockApiQueries",&queries);
    assigned.cmd("getMockCounts",counts.data());
    if(queries!=initialQueries+1 || counts[3]!=0) {
      throw std::runtime_error("moving the handle lost or mixed guest API capabilities");
    }
  }
  std::ofstream("output")<<"Loaded guests: API 1/3/4/8/9/10/13 guards and dispatch OK\n";
}

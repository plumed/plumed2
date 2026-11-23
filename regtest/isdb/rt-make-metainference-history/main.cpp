#include "plumed/isdb/MetainferenceCore.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

int main() {
  using PLMD::isdb::SigmaMeanHistory;
  for(unsigned capacity : {
        1U,2U,3U,9U
      }) {
    for(unsigned size=1; size<=capacity; ++size) {
      for(double candidate : {
            0.,4.,16.,100.
          }) {
        SigmaMeanHistory history(capacity);
        for(unsigned i=0; i<size; ++i) {
          history.push(100.-i*i);
        }
        const double expected=history.maxAfterPush(candidate);
        history.commitPush(candidate);
        if(std::abs(history.max()-expected)>1.e-14 ||
            std::abs(history.maxAfterPush(candidate)-expected)>1.e-14) {
          throw std::runtime_error("committed sample changed its own replay window");
        }
        const double second=history.maxAfterPush(2.);
        history.commitPush(2.);
        if(std::abs(history.max()-second)>1.e-14) {
          throw std::runtime_error("same-step commit did not replace its predecessor");
        }
        history.beginStep();
        const double next=history.maxAfterPush(3.);
        history.commitPush(3.);
        if(std::abs(history.max()-next)>1.e-14) {
          throw std::runtime_error("new step did not use committed history");
        }
      }
    }
    SigmaMeanHistory raised(capacity);
    raised.push(4.);
    raised.commitRaiseFirst(16.);
    if(std::abs(raised.maxAfterRaiseFirst(9.)-9.)>1.e-14) {
      throw std::runtime_error("raise-first replay lost its original sample");
    }
  }
  std::ofstream("output")<<"History replay, partial windows, repeated commits and step boundaries: OK\n";
}

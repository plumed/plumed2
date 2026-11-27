#include "plumed/isdb/MetainferenceCore.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

int main() {
  using PLMD::isdb::SigmaMeanHistory;
  // Compare arbitrary replacements against a materialized history, including
  // decreases of a unique/duplicate maximum and rotated, partially filled rings.
  for(unsigned capacity : {
        1U,2U,3U,9U
      }) {
    for(unsigned size=1; size<=capacity; ++size) {
      SigmaMeanHistory replaced(capacity);
      std::vector<double> samples;
      for(unsigned i=0; i<size; ++i) {
        const double sample=i%2==0 ? 9. : 4.;
        replaced.push(sample);
        samples.push_back(sample);
      }
      for(unsigned step=0; step<40; ++step) {
        for(double candidate : {
              0.,4.,9.,16.,2.
            }) {
          replaced.setFirst(candidate);
          samples.front()=candidate;
          if(replaced.values()!=samples ||
              replaced.max()!=*std::max_element(samples.begin(),samples.end())) {
            throw std::runtime_error("replacement maximum disagrees with stored history");
          }
        }
        const double sample=step%3==0 ? 9. : 4.;
        replaced.push(sample);
        if(samples.size()==capacity) {
          samples.erase(samples.begin());
        }
        samples.push_back(sample);
        if(replaced.values()!=samples ||
            replaced.max()!=*std::max_element(samples.begin(),samples.end())) {
          throw std::runtime_error("rotation changed replacement history or its maximum");
        }
      }
    }
  }
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
  // Cumulative (AVERAGING=0) estimates must retain the maximum on every
  // ordinary step, with a mixture of increases and decreases across data.
  std::vector<SigmaMeanHistory> cumulative(15, SigmaMeanHistory(1));
  std::vector<double> expected(15,1.e-8);
  for(auto& history : cumulative) {
    history.push(1.e-8);
  }
  for(unsigned step=0; step<40; ++step) {
    for(unsigned datum=0; datum<cumulative.size(); ++datum) {
      auto& history=cumulative[datum];
      history.beginStep();
      const double candidate=std::pow(0.001*(1+(step*(datum+1))%17),2);
      expected[datum]=std::max(expected[datum],candidate);
      if(history.maxAfterRaiseFirst(candidate)!=expected[datum]) {
        throw std::runtime_error("prospective cumulative maximum lost a past sample");
      }
      history.commitRaiseFirst(candidate);
      if(history.get(0)!=expected[datum] || history.max()!=expected[datum]) {
        throw std::runtime_error("committed cumulative maximum disagrees with its stored entry");
      }
    }
  }
  std::ofstream("output")<<"History replay, partial windows, repeated commits and step boundaries: OK\n";
}

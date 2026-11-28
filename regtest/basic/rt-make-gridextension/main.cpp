#include "plumed/tools/GridExtension.h"
#include "plumed/tools/GridCoordinatesObject.h"
#include "testUtils.h"
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace PLMD;
using Mode=GridExtension::Mode;

template <typename T>
std::string str(const std::vector<T>& v) {
  std::string s;
  for(unsigned i=0; i<v.size(); ++i) {
    s+=(i>0?" ":"")+std::to_string(v[i]);
  }
  return s;
}

void print(tee& out, const std::string& name, const GridExtension::Plan& p) {
  out << name << ": changed " << p.changed << " lower " << str(p.lower) << " upper " << str(p.upper)
      << " npoints " << p.npoints << " toolarge " << p.toolarge
      << " nonfinite " << (p.nonfinite==std::string::npos ? std::string("none") : std::to_string(p.nonfinite)) << "\n";
}

void print(tee& out, const std::string& name, const GridCoordinatesObject& g) {
  out << name << ": min";
  for(const auto& s : g.getMin()) {
    out << " " << s;
  }
  out << " max";
  for(const auto& s : g.getMax()) {
    out << " " << s;
  }
  out << " points " << str(g.getNbin(true)) << " total " << g.getNumberOfPoints() << "\n";
}

int main() {
  tee out("output");
  out << "checkSizeLimit(0.5): " << GridExtension::checkSizeLimit(0.5) << "\n";
  out << "checkSizeLimit(1e10): " << GridExtension::checkSizeLimit(1e10) << "\n";
  out << "checkSizeLimit(1e8): '" << GridExtension::checkSizeLimit(1e8) << "'\n";

  GridExtension e;
  out << "raised " << e.setSizeLimit(1000,{11,11}) << " limit " << e.getSizeLimit() << "\n";
  // nothing needed beyond the grid (101 points: indices 0..100)
  print(out,"inside",e.plan({101},{Mode::extend},{0},{100}));
  // two points missing below: growth to a quarter of the points (26)
  print(out,"below",e.plan({101},{Mode::extend},{-2},{50}));
  // 40 points missing above: more than the growth
  print(out,"above",e.plan({101},{Mode::extend},{10},{140}));
  // periodic (fixed) direction and an extended one
  print(out,"fixed",e.plan({20,11},{Mode::fixed,Mode::extend},{-5,-1},{30,10}));
  // the growth would exceed the limit: only the missing points
  print(out,"fallback",e.plan({30,30},{Mode::extend,Mode::extend},{-1,0},{29,29}));
  // too large even without growth
  print(out,"toolarge",e.plan({30,30},{Mode::extend,Mode::extend},{-10,0},{29,29}));
  // placement: exactly the needed range, also fewer points than now
  print(out,"place",e.plan({2},{Mode::place},{3},{9}));
  print(out,"place-shift",e.plan({2},{Mode::place},{-7},{-3}));
  print(out,"nonfinite",e.plan({11,11},{Mode::extend,Mode::extend},{0,std::numeric_limits<double>::quiet_NaN()},{10,10}));
  print(out,"far",e.plan({11},{Mode::extend},{-1e16},{10}));

  // a grid on a lattice, extended as GridBase::resize
  GridCoordinatesObject g;
  g.setup("flat",{false,true},0,0.0);
  g.setLattice({-1.0,-3.14159265358979},{0.1,3.14159265358979/10},{20,20},{"-1","-pi"},{"1","pi"});
  print(out,"lattice",g);
  out << "lattice coordinate of 0.25: " << g.getLatticeCoordinate(0,0.25) << "\n";
  GridCoordinatesObject h(g);
  h.extend({5,0},{-2,0});
  print(out,"extended",h);
  out << "same lattice " << h.onSameLattice(g) << " shift " << str(h.getShiftFrom(g)) << " back " << str(g.getShiftFrom(h)) << "\n";
  GridCoordinatesObject k(h);
  k.extend({-3,0},{7,0});
  out << "same lattice " << k.onSameLattice(g) << " shift " << str(k.getShiftFrom(g)) << "\n";
  GridCoordinatesObject other;
  other.setup("flat",{false,true},0,0.0);
  other.setLattice({-1.0,-3.14159265358979},{0.1,3.14159265358979/10},{20,20},{"-1","-pi"},{"1","pi"});
  out << "a new lattice with the same numbers: same lattice " << other.onSameLattice(g) << "\n";
  other.setLattice({-1.05,-3.14159265358979},{0.1,3.14159265358979/10},{20,20});
  out << "a shifted lattice: same lattice " << other.onSameLattice(g) << "\n";
  return 0;
}

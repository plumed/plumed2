#include "plumed/core/Value.h"
#include <fstream>
int main() {
  PLMD::Value scalar(nullptr,"scalar",true);
  scalar.resizeDerivatives(2);
  scalar.setDerivative(0,2);
  scalar.addDerivative(0,3);
  scalar.setDerivative(1,7);
  bool ok=scalar.getDerivative(0)==5 && scalar.getDerivative(1)==7;
  scalar.set(3);
  scalar.setConstant();
  scalar.addDerivative(0,4);
  scalar.setDerivative(1,8);
  ok=ok && scalar.get()==3 && scalar.getNumberOfDerivatives()==0;
  // Constant grids retain spatial derivatives for interpolation.
  PLMD::Value grid(nullptr,"grid",true,{2});
  grid.setConstant();
  grid.setDerivative(0,2);
  grid.addDerivative(0,3);
  ok=ok && grid.getGridDerivative(0,0)==5;
  std::ofstream("output") << (ok ? "PASS" : "FAIL") << "\n";
  return ok ? 0 : 1;
}

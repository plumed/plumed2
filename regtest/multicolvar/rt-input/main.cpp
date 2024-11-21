#include "plumed/colvar/MultiColvarTemplate.h"
#include "plumed/tools/Vector.h"

#include <fstream>
#include <iostream>
#include <variant>
#include <vector>

using namespace PLMD;
using colvar::multiColvars::Input;

template <typename T> void pv(std::ofstream &stream, const std::vector<T> &v) {
  // todo: put this is a common header
  std::string div = "";
  stream << "(" << v.size() << ") ";
  for (const auto &t : v) {
    stream << div << t;
    div = ", ";
  }
}

void print(std::ofstream &stream, const std::vector<double> &masses,
           const std::vector<double> &charges,
           const std::vector<Vector> positions) {
  stream << "Masses: ";
  pv(stream, masses);
  stream << "\nCharges: ";
  pv(stream, charges);
  stream << "\nPositions: ";
  pv(stream, positions);
  stream << "\n";
}

void inputTest(std::ofstream &stream, Input input) {
  try {
    const auto &masses = input.masses();
    stream << " Masses: ";
    pv(stream, masses);
    stream << "\n";
  } catch (std::bad_variant_access &) {
    stream << " There are no masses\n";
  }

  try {
    const auto &charges = input.charges();
    stream << " Charges: ";
    pv(stream, charges);
    stream << "\n";
  } catch (std::bad_variant_access &) {
    stream << " There are no charges\n";
  }
  try {
    const auto &positions = input.positions();
    stream << " Positions: ";
    pv(stream, positions);
    stream << "\n";
  } catch (std::bad_variant_access &) {
    stream << " There are no positions\n";
  }
  stream << std::endl;
}

void workOnMass(Input input) {
  auto &masses = input.var_masses();
  for (auto &m : masses) {
    m *= 2;
  }
}

void workOnCharges(Input input) {
  auto &charges = input.var_charges();
  for (auto &c : charges) {
    c *= 2;
  }
}

void workOnPositions(Input input) {
  auto &pos = input.var_positions();
  for (auto &p : pos) {
    p = p * 2.0;
  }
}

class testData {
  using vd = std::vector<double>;
  using vv = std::vector<Vector>;
  vv p = {{0, 0, 0}, {1, 1, 1}};
  vd c = {-1, 1};
  vd m = {3, 4};

public:
  testData() = default;
  vd &masses() { return m; }
  vd &charges() { return c; }
  vv &positions() { return p; }
  const vd &c_masses() const { return m; }
  const vd &c_charges() const { return c; }
  const vv &c_positions() const { return p; }
};

class test {
  using vd = std::vector<double>;
  using vv = std::vector<Vector>;

public:
  test() = default;
  void dotestWithReference(std::ofstream &output) {
    testData td;
    output << "*dotestWithReference:" << std::endl;
    output << "*Passing full input:\n";
    inputTest(output, Input(td.masses(), td.charges(), td.positions()));
    output << "*Passing only positions:\n";
    inputTest(output, Input().positions(td.positions()));
    output << "*Passing only masses:\n";
    inputTest(output, Input().masses(td.masses()));
    output << "*Passing only charges:\n";
    inputTest(output, Input().charges(td.charges()));
  }
  void dotestWithConstReference(std::ofstream &output) {
    testData td;
    output << "*dotestWithConstReference:" << std::endl;
    output << "*Passing full input:\n";
    inputTest(output, Input(td.c_masses(), td.c_charges(), td.c_positions()));
    output << "*Passing only positions:\n";
    inputTest(output, Input().positions(td.c_positions()));
    output << "*Passing only masses:\n";
    inputTest(output, Input().masses(td.c_masses()));
    output << "*Passing only charges:\n";
    inputTest(output, Input().charges(td.c_charges()));
  }

  void testUsingInputAsBuffer_masses(std::ofstream &output) {
    testData td;
    output << "*testUsingInputAsBuffer_masses:" << std::endl;
    output << "*Passing constant_reference:\n";
    try {
      workOnMass(Input(td.c_masses(), td.c_charges(), td.c_positions()));
    } catch (std::bad_variant_access &) {
      output << " should throw\n";
    }
    output << "*Passing reference:\n";
    output << " Initial values\n Masses:";
    pv(output, td.masses());
    workOnMass(Input(td.masses(), td.charges(), td.positions()));
    output << "\n";
    // this extra line is to remember that something has changed!!!
    output << " Values after call\n Masses:";
    pv(output, td.masses());
    output << std::endl;
  }
  void testUsingInputAsBuffer_charges(std::ofstream &output) {
    testData td;
    output << "*testUsingInputAsBuffer_charges:" << std::endl;
    output << "*Passing constant_reference:\n";
    try {
      workOnCharges(Input(td.c_masses(), td.c_charges(), td.c_positions()));
    } catch (std::bad_variant_access &) {
      output << " should throw\n";
    }
    output << "*Passing reference:\n";
    output << " Initial values\n Charges:";
    pv(output, td.c_charges());
    workOnCharges(Input(td.masses(), td.charges(), td.positions()));
    output << "\n";
    // this extra line is to remember that something has changed!!!
    output << " Values after call\n Charges:";
    pv(output, td.c_charges());
    output << std::endl;
  }
  void testUsingInputAsBuffer_positions(std::ofstream &output) {
    testData td;
    output << "*testUsingInputAsBuffer_positions:" << std::endl;
    output << "*Passing constant_reference:\n";
    try {
      workOnPositions(Input(td.c_masses(), td.c_charges(), td.c_positions()));
    } catch (std::bad_variant_access &) {
      output << " should throw\n";
    }
    output << "*Passing reference:\n";
    output << " Initial values\n Positions:";
    pv(output, td.positions());
    workOnPositions(Input(td.masses(), td.charges(), td.positions()));
    output << "\n";
    // this extra line is to remember that something has changed!!!
    output << " Values after call\n Positions:";
    pv(output, td.positions());
    output << std::endl;
  }
};

int main() {
  std::ofstream output("output");

  test t;
  t.dotestWithReference(output);
  t.dotestWithConstReference(output);
  t.testUsingInputAsBuffer_masses(output);
  t.testUsingInputAsBuffer_charges(output);
  t.testUsingInputAsBuffer_positions(output);
  return 0;
}

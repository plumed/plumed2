#include "plumed/lepton/Lepton.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PLMD::lepton;

static bool same(double actual, double expected) {
  if (std::isnan(expected)) return std::isnan(actual);
  if (std::isinf(expected)) return actual == expected;
  return std::abs(actual - expected) <= 1e-12 * (1 + std::abs(expected));
}

int main() {
  if (!useAsmJit()) throw std::runtime_error("JIT was not enabled");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<std::string> expressions = {
    "step(x)", "delta(x)", "nandelta(x)", "select(x,2,3)",
    "1/(1+x^6)", "x^2+x^4+x^6", "x^(-2)+x^(-4)",
    "atan2(x,2)", "acot(x)", "sech(x)", "erf(x)", "sin(x)+exp(x)"
  };
  for (const auto& text : expressions) {
    auto parsed = Parser::parse(text);
    auto compiled = parsed.createCompiledExpression();
    auto derivative = parsed.differentiate("x").optimize();
    auto compiledDerivative = derivative.createCompiledExpression();
    for (double x : {-2.0, -0.5, -0.0, 0.0, 0.5, 2.0, nan, inf, -inf}) {
      std::map<std::string, double> variables{{"x", x}};
      if (compiled.getVariables().count("x")) compiled.getVariableReference("x") = x;
      if (compiledDerivative.getVariables().count("x")) compiledDerivative.getVariableReference("x") = x;
      if (!same(compiled.evaluate(), parsed.evaluate(variables)))
        throw std::runtime_error("JIT value mismatch: " + text + " x=" + std::to_string(x));
      if (!same(compiledDerivative.evaluate(), derivative.evaluate(variables)))
        throw std::runtime_error("JIT derivative mismatch: " + text + " x=" + std::to_string(x) + " actual=" + std::to_string(compiledDerivative.evaluate()) + " expected=" + std::to_string(derivative.evaluate(variables)));
    }
  }
  // Repeated rebinding regenerates code, including in copied expressions.
  auto parsed = Parser::parse("1/(1+x^6)");
  auto compiled = parsed.createCompiledExpression();
  for (int i = 0; i < 50; ++i) {
    double x = 0.1 * i;
    std::map<std::string, double*> locations{{"x", &x}};
    compiled.setVariableLocations(locations);
    auto copy = compiled;
    copy.setVariableLocations(locations);
    if (!same(compiled.evaluate(), 1/(1+std::pow(x, 6))) || !same(copy.evaluate(), compiled.evaluate()))
      throw std::runtime_error("Rebinding or copying failed");
  }
  std::ofstream("output") << "Lepton JIT values, derivatives, signed zero, NaN, and rebinding passed\n";
}

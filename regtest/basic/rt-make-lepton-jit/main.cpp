#include "plumed/lepton/Lepton.h"
#include <cmath>
#include <cstdlib>
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
  const char* jitSetting = std::getenv("PLUMED_USE_ASMJIT");
  const bool expectJit = jitSetting && std::string(jitSetting) == "yes";
  if (useAsmJit() != expectJit) throw std::runtime_error("Unexpected JIT setting");
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
  // Check the derivative convention against independent expectations, including
  // signed zero. Comparing two evaluators alone cannot detect a symbolic error.
  struct DerivativeCase {
    std::string expression;
    std::vector<double> first;
    std::vector<double> second;
  };
  const std::vector<DerivativeCase> derivativeCases = {
    {"step(x)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"step(-x)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"step(y)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"delta(x)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"delta(y)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"abs(x)", {-1, 1, 1, 1}, {0, 0, 0, 0}},
    {"x^2*step(x)", {0, 0, 0, 2}, {0, 2, 2, 2}},
    {"x^2*step(-x)", {-2, 0, 0, 0}, {2, 2, 2, 0}},
    {"min(x,0)", {1, 0, 0, 0}, {0, 0, 0, 0}},
    {"max(x,0)", {0, 1, 1, 1}, {0, 0, 0, 0}},
    {"floor(x)", {0, 0, 0, 0}, {0, 0, 0, 0}},
    {"ceil(x)", {0, 0, 0, 0}, {0, 0, 0, 0}}
  };
  const std::vector<double> points = {-1.0, -0.0, 0.0, 1.0};
  for (const auto& test : derivativeCases) {
    auto derivative = Parser::parse(test.expression);
    for (unsigned order = 1; order <= 2; ++order) {
      derivative = derivative.differentiate("x").optimize();
      auto compiled = derivative.createCompiledExpression();
      const auto& expected = order == 1 ? test.first : test.second;
      for (unsigned i = 0; i < points.size(); ++i) {
        const std::map<std::string, double> variables{{"x", points[i]}, {"y", 0.0}};
        for (const auto& name : compiled.getVariables())
          compiled.getVariableReference(name) = variables.at(name);
        if (!same(derivative.evaluate(variables), expected[i]) || !same(compiled.evaluate(), expected[i]))
          throw std::runtime_error("Wrong derivative: " + test.expression + " order=" + std::to_string(order) + " x=" + std::to_string(points[i]));
      }
    }
  }
  // Derivative changes must preserve explicit PLUMED step/delta/nandelta values.
  for (const auto& test : std::map<std::string, double>{{"step(x)", 1.0}, {"delta(x)", inf}, {"nandelta(x)", nan}}) {
    auto parsed = Parser::parse(test.first);
    auto compiled = parsed.createCompiledExpression();
    for (double x : {-0.0, 0.0}) {
      compiled.getVariableReference("x") = x;
      if (!same(parsed.evaluate({{"x", x}}), test.second) || !same(compiled.evaluate(), test.second))
        throw std::runtime_error("Changed value at zero: " + test.first);
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
  std::ofstream("output") << "Lepton values, derivatives, signed zero, NaN, and rebinding passed\n";
}

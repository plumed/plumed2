#include "plumed/tools/Tools.h"
#include <fstream>
#include <string>

using namespace PLMD;

int main() {
  std::ofstream ofs("output");
  for(const auto& in: {
        R"(file\.ext)",
        R"(file.ext)",
        R"([a-zA-Z0-9]\.ext)",
        R"(\[squared\])",
        R"(\{graph\})",
        R"(\(round\))",
        R"(\\slash)",
        R"(\\slashes\\)",
        R"(anyspace\s)",
        R"(\\\^various\$\.)",
        R"(\|\?)",
        R"(\*\+)",
        R"( \[with extra space\] )",
        R"<<(lab: PRODUCT ARG=(lab_diag\.vals-[0-9]))<<"
      }) {
    std::string out = Tools::convertRegexForJson(in);
    ofs <<"\""<< in<< "\" -> \"" << out <<"\"\n";
  }
  return 0;
}

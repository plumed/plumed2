#include "core/CLTool.h"
#include "core/CLToolRegister.h"

namespace PLMD {
class AliasTool : public CLTool {
public:
  static void registerKeywords(Keywords& keys) {
    CLTool::registerKeywords(keys);
    keys.add("compulsory","--input/-i","input string");
    keys.add("compulsory","--number/-n","7","numeric value");
    keys.add("optional","--values/-v","vector values");
    keys.addFlag("--flag/-f",false,"boolean flag");
  }
  explicit AliasTool(const CLToolOptions& options): CLTool(options) {
    inputdata=inputType::commandline;
  }
  int main(FILE*,FILE* out,Communicator&) override {
    std::string input;
    int number;
    bool flag;
    std::vector<int> values;
    parse("--input",input);
    parse("--number",number);
    parseFlag("--flag",flag);
    parseVector("--values",values);
    std::fprintf(out,"input=%s number=%d flag=%d",input.c_str(),number,int(flag));
    for(int value : values) {
      std::fprintf(out," value=%d",value);
    }
    std::fprintf(out,"\n");
    // Forwarded driver input must use canonical names, without alias records.
    for(const auto& entry : getInputData()) {
      plumed_massert(entry.first.find('/')==std::string::npos,"noncanonical input key");
      plumed_massert(entry.first.size()>2,"short input key");
    }
    return 0;
  }
};
PLUMED_REGISTER_CLTOOL(AliasTool,"alias_tool")
}

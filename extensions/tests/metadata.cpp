#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/CLTool.h"
#include "core/CLToolRegister.h"

namespace PLMD {
class MetadataAction:public Action {
public:
  static void registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
    keys.add("optional","VALUE","metadata fixture value");
  }
  explicit MetadataAction(const ActionOptions& options):Action(options) {
    checkRead();
  }
  void calculate() override {}
  void apply() override {}
};
PLUMED_REGISTER_ACTION(MetadataAction,"METADATA_ACTION")

class MetadataTool:public CLTool {
public:
  static void registerKeywords(Keywords& keys) {
    CLTool::registerKeywords(keys);
    keys.addFlag("--probe",false,"metadata fixture flag");
  }
  explicit MetadataTool(const CLToolOptions& options):CLTool(options) {
    inputdata=inputType::commandline;
  }
  int main(FILE*,FILE*,Communicator&) override {
    return 0;
  }
};
PLUMED_REGISTER_CLTOOL(MetadataTool,"metadata_tool")
}

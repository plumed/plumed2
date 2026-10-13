#include "core/Action.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/OFile.h"
#include <exception>
#include <string>

namespace PLMD {

// Checks that multi_sim_comm is marked as a stand-in exactly on the ranks
// that were not given the inter-replica communicator, and that using it there
// raises an error in debug builds only.
class TestStandIn : public Action {
public:
  static void registerKeywords( Keywords& keys );
  explicit TestStandIn(const ActionOptions&);
  void calculate() override {}
  void apply() override {}
};

PLUMED_REGISTER_ACTION(TestStandIn,"TEST_STAND_IN")

void TestStandIn::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
}

TestStandIn::TestStandIn(const ActionOptions&ao):
  Action(ao) {
// driver gives the inter-replica communicator to the first rank of each replica only
  const bool expectStandIn = getNumberOfReplicas()>1 && comm.Get_rank()!=0;
#ifdef NDEBUG
  const bool expectThrow = false;
#else
  const bool expectThrow = expectStandIn;
#endif
  bool threw=false;
  try {
    multi_sim_comm.Get_size();
  } catch(std::exception&e) {
    // catch std::exception rather than ExceptionDebug: in a library loaded
    // with LOAD, a class defined only in a header can have a type identity
    // distinct from the one used where it is thrown
    threw=std::string(e.what()).find("stand-in")!=std::string::npos;
  }
  int failures=0;
  if(multi_sim_comm.isStandIn()!=expectStandIn) {
    failures++;
  }
  if(threw!=expectThrow) {
    failures++;
  }
  comm.Sum(failures);
  OFile of;
  of.link(*this);
  of.open("stand-in.log");
  of.printf("replicas %d ranks per replica %d failures %d\n", getNumberOfReplicas(), comm.Get_size(), failures);
  of.close();
}

}

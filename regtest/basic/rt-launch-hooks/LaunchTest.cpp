#include "core/Colvar.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include <cstdio>
#include <vector>

// Test of Action::launch() and cancelLaunch(), see plumed.dat.
// The value is the sum of the x coordinates of the atoms (LAUNCH_TEST) or of the arguments (LAUNCH_TEST_ARG).

namespace PLMD {

// every call is written in the file "order"
static void writeOrder(const Action& a,const char* what) {
  static std::FILE* fp=std::fopen("order","w");
  std::fprintf(fp,"%ld %s %s\n",static_cast<long>(a.getStep()),what,a.getLabel().c_str());
  std::fflush(fp);
}

class LaunchTest : public Colvar {
  double launchedSum=0.0;
  double sum() const {
    double s=0.0;
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      s+=getPosition(i)[0];
    }
    return s;
  }
  void write(const char* what) {
    writeOrder(*this,what);
  }
public:
  static void registerKeywords( Keywords& keys ) {
    Colvar::registerKeywords( keys );
    keys.add("atoms","ATOMS","the atoms");
    keys.setValueDescription("scalar","the sum of the x coordinates of the atoms");
  }
  explicit LaunchTest(const ActionOptions&ao):
    PLUMED_COLVAR_INIT(ao) {
    std::vector<AtomNumber> atoms;
    parseAtomList("ATOMS",atoms);
    checkRead();
    addValueWithDerivatives();
    setNotPeriodic();
    requestAtoms(atoms);
    setLaunchable();
  }
  void launch() override {
    write("launch");
    launchedSum=sum();
  }
  void cancelLaunch() override {
    write("cancelLaunch");
  }
  void calculate() override {
    const bool launched=takeLaunch();
    write(launched ? "calculate (launched)" : "calculate");
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      setAtomsDerivatives(i,Vector(1,0,0));
    }
    setBoxDerivativesNoPbc();
    setValue(launched ? launchedSum : sum());
  }
};

PLUMED_REGISTER_ACTION(LaunchTest,"LAUNCH_TEST")

class LaunchTestArg :
  public ActionWithValue,
  public ActionWithArguments {
  double launchedSum=0.0;
  double sum() const {
    double s=0.0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      s+=getArgument(i);
    }
    return s;
  }
public:
  static void registerKeywords( Keywords& keys ) {
    Action::registerKeywords( keys );
    ActionWithValue::registerKeywords( keys );
    ActionWithArguments::registerKeywords( keys );
    keys.addInputKeyword("compulsory","ARG","scalar","the arguments");
    keys.setValueDescription("scalar","the sum of the arguments");
  }
  explicit LaunchTestArg(const ActionOptions&ao):
    Action(ao),
    ActionWithValue(ao),
    ActionWithArguments(ao) {
    addValueWithDerivatives();
    setNotPeriodic();
    setLaunchable();
  }
  unsigned getNumberOfDerivatives() override {
    return getNumberOfArguments();
  }
  void launch() override {
    writeOrder(*this,"launch");
    launchedSum=sum();
  }
  void calculate() override {
    const bool launched=takeLaunch();
    writeOrder(*this,launched ? "calculate (launched)" : "calculate");
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      getPntrToValue()->addDerivative(i,1.0);
    }
    setValue(launched ? launchedSum : sum());
  }
  void apply() override {
    if(checkForForces()) {
      unsigned start=0;
      addForcesOnArguments(0,getForcesToApply(),start);
    }
  }
};

PLUMED_REGISTER_ACTION(LaunchTestArg,"LAUNCH_TEST_ARG")

}

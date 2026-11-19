#include "core/Colvar.h"
#include "core/ActionRegister.h"
#include <cstdio>
#include <vector>

// Test of Action::launch(), cancelLaunch() and launchApply(), see plumed.dat.
// The value is the sum of the x coordinates of the atoms.

namespace PLMD {

class LaunchTest : public Colvar {
  double launchedSum=0.0;
  double launchedForce=0.0;
  double sum() const {
    double s=0.0;
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      s+=getPosition(i)[0];
    }
    return s;
  }
  void write(const char* what) {
    static std::FILE* fp=std::fopen("order","w");
    std::fprintf(fp,"%ld %s %s\n",static_cast<long>(getStep()),what,getLabel().c_str());
    std::fflush(fp);
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
    setLaunchable( true, true );
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
  void launchApply() override {
    write("launchApply");
    launchedForce=getPntrToValue()->getForce();
  }
  void apply() override {
    const bool launched=takeApplyLaunch();
    write(launched ? "apply (launched)" : "apply");
    if(launched && launchedForce!=getPntrToValue()->getForce()) {
      error("the force changed after launchApply()");
    }
    Colvar::apply();
  }
};

PLUMED_REGISTER_ACTION(LaunchTest,"LAUNCH_TEST")

}

// A bias written for PLUMED 2.10, where Bias::registerKeywords reserved ARG and actions called keys.use("ARG").
// The base class now adds ARG itself: keys.use("ARG") must still work (it stopped with "the ARG keyword is not reserved").
#include "bias/Bias.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace rtloadkeysuse {

class HarmonicOld : public bias::Bias {
  double kappa;
public:
  explicit HarmonicOld(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(HarmonicOld,"HARMONIC_OLD")

void HarmonicOld::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","KAPPA","the force constant");
}

HarmonicOld::HarmonicOld(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  kappa(0.0) {
  parse("KAPPA",kappa);
  checkRead();
}

void HarmonicOld::calculate() {
  const double x=getArgument(0);
  setBias(0.5*kappa*x*x);
  setOutputForce(0,-kappa*x);
}

}
}

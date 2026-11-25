#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/OFile.h"

namespace PLMD {

class TestDepreciatedAtoms : public Action {
public:
  static void registerKeywords( Keywords& keys );
  explicit TestDepreciatedAtoms(const ActionOptions&);
  void calculate() {}
  void apply() {}
};

PLUMED_REGISTER_ACTION(TestDepreciatedAtoms,"TEST_DEPRECIATED_ATOMS")

void TestDepreciatedAtoms::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
}

TestDepreciatedAtoms::TestDepreciatedAtoms(const ActionOptions&ao):
Action(ao)
{
  OFile of; of.link(*this);
  of.open("depreciated.log");
  of.printf("Boltzmann's constant is %f \n", plumed.getAtoms().getKBoltzmann());
  of.printf("Number of atoms is %d \n", plumed.getAtoms().getNatoms());
  of.printf("Are we using natural units %d \n", plumed.getAtoms().usingNaturalUnits());
  of.printf("So kBT is %f \n", plumed.getAtoms().getKbT());
  of.close();
}

}

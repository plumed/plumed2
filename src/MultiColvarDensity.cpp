#include "MultiColvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC MCOLVAR DENSITY
/*
Calculate functions of the density of atoms as a function of the box.  This allows one to calculate
density gradients, number of atoms in half the box and so on.

\par Examples 

The following example calculates the number of atoms in one half of the simulation box. 

\verbatim
DENSITY SPECIES=1-100 SUBCELL=(XLOWER=0.0 XUPPER=0.5) LABEL=d1
PRINT ARG=d1.* FILE=colvar1 FMT=%8.4f

\endverbatim

*/
//+ENDPLUMEDOC


class MultiColvarDensity : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarDensity(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
  void getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv );
  /// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives(){ plumed_assert(0); };
  bool isPeriodic(const unsigned nn){ return false; }
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"DENSITY")

void MultiColvarDensity::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  // Note we don't parallelize this as it would be stupid
  keys.use("SPECIES"); keys.remove("AVERAGE"); keys.remove("LESS_THAN"); 
  keys.remove("MIN"); keys.remove("MORE_THAN"); keys.remove("HISTOGRAM");
  keys.remove("WITHIN");
  // Use density keywords
  keys.use("SUBCELL"); keys.use("GRADIENT");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  int nat; readAtoms( nat ); 
  requestDistribution();
  // And check everything has been read in correctly
  checkRead(); 
}

double MultiColvarDensity::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
  return 1.0;
}

void MultiColvarDensity::getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv ){
   cpos=pos[0]; deriv[0]=Tensor::identity();
}

}


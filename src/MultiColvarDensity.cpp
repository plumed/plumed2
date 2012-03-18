#include "MultiColvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC MCOLVAR DENSITY
/**
Calculate functions of the density of atoms as a function of the box.  This allows one to calculate
density gradients, number of atoms in half the box and so on.

\par Examples 

*/
//+ENDPLUMEDOC


class MultiColvarDensity : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarDensity(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
  void getCentralAtom( const std::vector<Vector>& pos, std::vector<Value>& pos);
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"DENSITY")

void MultiColvarDensity::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  // Note we don't parallelize this as it would be stupid
  keys.use("SPECIES"); keys.remove("AVERAGE"); keys.remove("LESS_THAN"); 
  keys.remove("MIN"); keys.remove("MORE_THAN"); keys.remove("HISTOGRAM");
  keys.remove("WITHIN");
  // Use density keywords
  keys.use("SUBCELL"); 
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  int nat; readAtoms( nat ); 
  // And check everything has been read in correctly
  checkRead(); 
}

double MultiColvarDensity::compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
  return 1.0;
}

void MultiColvarDensity::getCentralAtom( const std::vector<Vector>& pos, std::vector<Value>& cpos){
  plumed_assert( cpos.size()==3 );
  Vector fracp; fracp=getPbc().realToScaled(pos[0]);
  Vector ff,cc;
  cpos[0].set(fracp[0]); 
  ff.clear(); ff[0]=1.0; cc=getPbc().realToScaled(ff);
  for(unsigned i=0;i<3;++i) cpos[0].addDerivative( i, cc[i] );  
  cpos[1].set(fracp[1]); 
  ff.clear(); ff[1]=1.0; cc=getPbc().realToScaled(ff);
  for(unsigned i=0;i<3;++i) cpos[1].addDerivative( i, cc[i] );
  cpos[2].set(fracp[2]);
  ff.clear(); ff[2]=1.0; cc=getPbc().realToScaled(ff);
  for(unsigned i=0;i<3;++i) cpos[2].addDerivative( i, cc[i] );
}

}


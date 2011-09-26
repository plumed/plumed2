#include "Colvar.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/**
Calculate the RMSD with respect to a reference structure

\par Syntax
\verbatim
RMSD REFERENCE=file.pdb
\endverbatim

...

*/
//+ENDPLUMEDOC
   
class ColvarRMSD : public Colvar {
private:
  RMSD rmsd;
  vector<Vector> pos;
public:
  ColvarRMSD(const ActionOptions&);
  double calcFunction( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
Colvar(ao)
{
  allowKeyword("ATOMS");
  registerKeyword( 1, "REFERENCE", "the reference structure we are calculating the rmsd from");
  std::vector<double> domain( 2, 0.0 );
  int natoms=-1; readActionColvar( natoms, domain );
  pos.resize(natoms); 

  string reference;
  parse("REFERENCE",reference);
  log.printf("  reference from file %s\n",reference.c_str());
  checkRead();

  PDB pdb;
  pdb.read(reference,0.1/plumed.getAtoms().getUnits().length);
  if( pdb.size()!=natoms ) error("number of atoms in reference pdb file " + reference + " does not match number of atoms in input");
  rmsd.setFromPDB(pdb);
}

double ColvarRMSD::calcFunction( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( pos.size()==indexes.size() && derivatives.size()==indexes.size() );
  for(unsigned i=0;i<pos.size();++i) pos[i]=getPositions( indexes[i] );
  double r=rmsd.calculate( pos, derivatives );
  virial.clear(); for(unsigned i=0;i<pos.size();i++) virial=virial+( -1.0*Tensor(pos[i],derivatives[i]) );
  return r;
}

}




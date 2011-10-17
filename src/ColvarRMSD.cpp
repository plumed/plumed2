#include "ColvarDistinguishable.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/**
Calculate RMSD distances between user defined sets of atoms and a reference structure.

\par Example

The following calculates the RMSD distance between the atoms listed and the positions of the 10 atoms in the reference file file.pdb
\verbatim
RMSD ATOMS=1-10 REFERENCE=file.pdb
\endverbatim

The following claculates the RMSD distace between the two sets of listed atoms and the 10 atoms in the reference file file.pdb
\verbatim
RMSD ATOMS=1-10 ATOMS=11-20 REFERENCE=file.pdb
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarRMSD : public ColvarDistinguishable {
private:
  RMSD rmsd;
  vector<Vector> pos;
public:
  ColvarRMSD(const ActionOptions&);
  double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
ColvarDistinguishable(ao)
{
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

double ColvarRMSD::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( pos.size()==indexes.size() && derivatives.size()==indexes.size() );
  for(unsigned i=0;i<pos.size();++i) pos[i]=getPositions( indexes[i] );
  double r=rmsd.calculate( pos, derivatives, virial );
  return r;
}

}




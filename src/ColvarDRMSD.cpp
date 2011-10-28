#include "Matrix.h"
#include "ColvarDistinguishable.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "DRMSD.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DRMSD
/**
Calculate DRMSD distances between user defined sets of atoms and a reference structure.

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
   
class ColvarDRMSD : public ColvarDistinguishable {
private:
  std::vector<Vector> pos;
  DRMSD drmsd;
public:
  ColvarDRMSD(const ActionOptions&);
  double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarDRMSD,"DRMSD")

ColvarDRMSD::ColvarDRMSD(const ActionOptions&ao):
ColvarDistinguishable(ao)
{
  registerKeyword( 1, "REFERENCE", "a file containing the reference structure we are calculating the drmsd from");
  registerKeyword( 1, "BOND_CUTOFF", "specifies that all bonds less than this value are chemical bonds and so their separation should be constant.  The DRMSD will thus not include any distances from the reference structure that are less than this value.");
  std::vector<double> domain( 2, 0.0 );
  int natoms=-1; readActionColvar( natoms, domain );
  pos.resize(natoms);

  string reference; parse("REFERENCE",reference);
  log.printf("  reference structure from file %s\n",reference.c_str());
  double bondc; parse("BOND_CUTOFF",bondc);
  log.printf("  all distances in refernce less than %f are bonds and will be ignored\n",bondc);
  checkRead();

  PDB pdb; pdb.read(reference,0.1/plumed.getAtoms().getUnits().length);
  if( pdb.size()!=natoms ) error("number of atoms in reference pdb file " + reference + " does not match number of atoms in input");
  drmsd.setFromPDB(bondc, pdb);
}

double ColvarDRMSD::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( pos.size()==indexes.size() && derivatives.size()==indexes.size() );

  for(unsigned i=0;i<pos.size();++i) pos[i]=getPositions( indexes[i] );
  double r=drmsd.calculate( pos, derivatives, virial );
  return r;
}

}




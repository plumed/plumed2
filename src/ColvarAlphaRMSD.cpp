#include "ColvarDistinguishable.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"
#include "DRMSD.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR ALPHARMSD
/**
Calculate the rmsd/drmsd distance of portions of the protein backbone from an alpha helix.  By default the distances from
the alpha helix of all the portions of the backbone defined in the input are calculated.  In other words this keyword is only
 equivalent to the ALPHARMSD keyword in plumed 1.0 if you use it in conjunction with the 
LESS_THAN colvar modifier described in \ref Colvar.

\par Example

The following calculates the number of alpha helix portions in the backbone between residues 1 and 10 that are within 2.0 of
an alpha helix

\verbatim
PROTEIN_TOPOLOGY FILE=topol.pdb
ALPHARMSD BACKBONE=1-10 LESS_THAN=2.0
\endverbatim

(See also \ref PROTEIN_TOPOLOGY )

*/
//+ENDPLUMEDOC

class ColvarAlphaRMSD : public ColvarDistinguishable {
private:
  bool used;
  RMSD rmsd;
  DRMSD drmsd;
  vector<Vector> pos;
public:
  ColvarAlphaRMSD(const ActionOptions&);
  double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarAlphaRMSD,"ALPHARMSD")

ColvarAlphaRMSD::ColvarAlphaRMSD(const ActionOptions&ao):
ColvarDistinguishable(ao),
used(false)
{
  registerKeyword( 2, "BACKBONE", "(strongly recommended) build the secondary structure from the backbone atoms in the specified residues");
  allowKeyword("BACKBONE");
  registerKeyword( 0, "DRMSD", "calcuate the distance using a DRMSD metric rather than the RMSD metric");


  std::vector< std::vector<unsigned> > backbone; 
  if( readBackboneAtoms( "protein", backbone ) ){
      unsigned n, nres, nprevious=0; std::vector<unsigned> atom_list(30);
      for(unsigned i=0;i<backbone.size();++i){
         if( backbone[i].size()<30 ) error("segment of backbone defined is not long enough to form an alpha helix. Each backbone fragment must contain a minimum of 6 residues");
         nres=backbone[i].size()/5; assert( backbone[i].size()%5==0 ); n=0;
         for(unsigned ires=0;ires<nres-5;ires++){
           n=0; for(unsigned iat=ires*5;iat<(ires+6)*5;iat++){ atom_list[n] = backbone[i][ nprevious + iat ]; n++; }
           assert(n==30); addIndexes( nprevious + 5*ires, atom_list );
         }
         nprevious+=backbone[i].size();
      }
  } 
  std::vector<double> domain( 2, 0.0 );
  readActionColvar( 30, domain ); pos.resize(30);
  parseFlag("DRMSD",used);
  checkRead();

  // Built the reference structure
  std::vector<Vector> reference(30);
  reference[0] = Vector( 0.733,  0.519,  5.298 ); // N    i
  reference[1] = Vector( 1.763,  0.810,  4.301 ); // CA
  reference[2] = Vector( 3.166,  0.543,  4.881 ); // CB
  reference[3] = Vector( 1.527, -0.045,  3.053 ); // C
  reference[4] = Vector( 1.646,  0.436,  1.928 ); // O
  reference[5] = Vector( 1.180, -1.312,  3.254 ); // N    i+1
  reference[6] = Vector( 0.924, -2.203,  2.126 ); // CA
  reference[7] = Vector( 0.650, -3.626,  2.626 ); // CB
  reference[8] = Vector(-0.239, -1.711,  1.261 ); // C
  reference[9] = Vector(-0.190, -1.815,  0.032 ); // O
  reference[10] = Vector(-1.280, -1.172,  1.891 ); // N    i+2
  reference[11] = Vector(-2.416, -0.661,  1.127 ); // CA
  reference[12] = Vector(-3.548, -0.217,  2.056 ); // CB
  reference[13] = Vector(-1.964,  0.529,  0.276 ); // C
  reference[14] = Vector(-2.364,  0.659, -0.880 ); // O
  reference[15] = Vector(-1.130,  1.391,  0.856 ); // N    i+3
  reference[16] = Vector(-0.620,  2.565,  0.148 ); // CA
  reference[17] = Vector( 0.228,  3.439,  1.077 ); // CB
  reference[18] = Vector( 0.231,  2.129, -1.032 ); // C
  reference[19] = Vector( 0.179,  2.733, -2.099 ); // O
  reference[20] = Vector( 1.028,  1.084, -0.833 ); // N    i+4
  reference[21] = Vector( 1.872,  0.593, -1.919 ); // CA
  reference[22] = Vector( 2.850, -0.462, -1.397 ); // CB
  reference[23] = Vector( 1.020,  0.020, -3.049 ); // C
  reference[24] = Vector( 1.317,  0.227, -4.224 ); // O
  reference[25] = Vector(-0.051, -0.684, -2.696 ); // N    i+5
  reference[26] = Vector(-0.927, -1.261, -3.713 ); // CA
  reference[27] = Vector(-1.933, -2.219, -3.074 ); // CB
  reference[28] = Vector(-1.663, -0.171, -4.475 ); // C
  reference[29] = Vector(-1.916, -0.296, -5.673 ); // O
  // The above coordinates are in angstroms - we convert them to plumed units (nm)
  // and into whetver unit the user is using
  Tools::scale( reference, 0.1*plumed.getAtoms().getUnits().length );    

  if(used){
     log.printf("  using DRMSD metric to measure distances from structure\n");
     drmsd.setReference( 0.17, reference );
  } else {
     log.printf("  using RMSD metric to measure distances from structure\n"); 
     rmsd.setReference( reference );
  }
}

double ColvarAlphaRMSD::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( pos.size()==indexes.size() && derivatives.size()==indexes.size() );
  for(unsigned i=0;i<pos.size();++i) pos[i]=getPositions( indexes[i] );

  double r;
  if (used) r=drmsd.calculate( pos, derivatives, virial );
  else r=rmsd.calculate( pos, derivatives, virial );
  return r;
}

}


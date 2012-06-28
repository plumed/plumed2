#include "MultiColvarSecondaryStructureRMSD.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

//+PLUMEDOC MCOLVAR ALPHARMSD
/**

\par Examples

*/
//+ENDPLUMEDOC

class ColvarAlphaRMSD : public MultiColvarSecondaryStructureRMSD {
public:
  static void registerKeywords( Keywords& keys );
  ColvarAlphaRMSD(const ActionOptions&);
}; 

PLUMED_REGISTER_ACTION(ColvarAlphaRMSD,"ALPHARMSD")

void ColvarAlphaRMSD::registerKeywords( Keywords& keys ){
  MultiColvarSecondaryStructureRMSD::registerKeywords( keys );
}

ColvarAlphaRMSD::ColvarAlphaRMSD(const ActionOptions&ao):
Action(ao),
MultiColvarSecondaryStructureRMSD(ao)
{
  // read in the backbone atoms
  std::vector<std::string> backnames(5); std::vector<unsigned> chains;
  backnames[0]="N"; backnames[1]="CA"; backnames[2]="CB"; backnames[3]="C"; backnames[4]="O";
  readBackboneAtoms( backnames, chains);

  // This constructs all conceivable sections of alpha helix in the backbone of the chains
  unsigned nres, n, nprevious=0; std::vector<unsigned> nlist(30);
  for(unsigned i=0;i<chains.size();++i){
     if( chains[i]<30 ) error("segment of backbone defined is not long enough to form an alpha helix. Each backbone fragment must contain a minimum of 6 residues");
     nres=chains[i]/5; plumed_assert( chains[i]%5==0 ); n=0;
     for(unsigned ires=0;ires<nres-5;ires++){
       unsigned accum=nprevious + 5*ires; 
       for(unsigned k=0;k<30;++k) nlist[k] = accum+k;
       addColvar( nlist );
     }
     nprevious+=chains[i];
  }

  // Build the reference structure ( in angstroms )
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
  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( reference, 0.17/atoms.getUnits().length, 0.1/atoms.getUnits().length ); 
}

}

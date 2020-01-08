/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "SecondaryStructureRMSD.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace secondarystructure {

//+PLUMEDOC COLVAR ALPHARMSD
/*
Probe the alpha helical content of a protein structure.

Any chain of six contiguous residues in a protein chain can form an alpha helix. This
colvar thus generates the set of all possible six residue sections and calculates
the RMSD distance between the configuration in which the residues find themselves
and an idealized alpha helical structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures.

This colvar is based on the following reference \cite pietrucci09jctc.  The authors of
this paper use the set of distances from the alpha helix configurations to measure
the number of segments that have an alpha helical configuration. This is done by calculating
the following sum of functions of the rmsd distances:

\f[
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
\f]

where the sum runs over all possible segments of alpha helix.  By default the
NN, MM and D_0 parameters are set equal to those used in \cite pietrucci09jctc.  The R_0
parameter must be set by the user - the value used in \cite pietrucci09jctc was 0.08 nm.

If you change the function in the above sum you can calculate quantities such as the average
distance from a purely the alpha helical configuration or the distance between the set of
residues that is closest to an alpha helix and the reference configuration. To do these sorts of
calculations you can use the AVERAGE and MIN keywords. In addition you can use the LESS_THAN
keyword if you would like to change the form of the switching function. If you use any of these
options you no longer need to specify NN, R_0, MM and D_0.

Please be aware that for codes like gromacs you must ensure that plumed
reconstructs the chains involved in your CV when you calculate this CV using
anything other than TYPE=DRMSD.  For more details as to how to do this see \ref WHOLEMOLECULES.

\par Examples

The following input calculates the number of six residue segments of
protein that are in an alpha helical configuration.

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
alpha: ALPHARMSD RESIDUES=all
\endplumedfile

Here the same is done use RMSD instead of DRMSD

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
alpha: ALPHARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1
\endplumedfile

*/
//+ENDPLUMEDOC

class AlphaRMSD : public SecondaryStructureRMSD {
public:
  static void registerKeywords( Keywords& keys );
  explicit AlphaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AlphaRMSD,"ALPHARMSD")

void AlphaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureRMSD::registerKeywords( keys );
}

AlphaRMSD::AlphaRMSD(const ActionOptions&ao):
  Action(ao),
  SecondaryStructureRMSD(ao)
{
  // read in the backbone atoms
  std::vector<unsigned> chains; readBackboneAtoms( "protein", chains);

  // This constructs all conceivable sections of alpha helix in the backbone of the chains
  unsigned nprevious=0; std::vector<unsigned> nlist(30);
  for(unsigned i=0; i<chains.size(); ++i) {
    if( chains[i]<30 ) error("segment of backbone defined is not long enough to form an alpha helix. Each backbone fragment must contain a minimum of 6 residues");
    unsigned nres=chains[i]/5;
    if( chains[i]%5!=0 ) error("backbone segment received does not contain a multiple of five residues");
    for(unsigned ires=0; ires<nres-5; ires++) {
      unsigned accum=nprevious + 5*ires;
      for(unsigned k=0; k<30; ++k) nlist[k] = accum+k;
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
  setSecondaryStructure( reference, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );
}

}
}

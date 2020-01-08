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

//+PLUMEDOC COLVAR PARABETARMSD
/*
Probe the parallel beta sheet content of your protein structure.

Two protein segments containing three contiguous residues can form a parallel beta sheet.
Although if the two segments are part of the same protein chain they must be separated by
a minimum of 3 residues to make room for the turn. This colvar thus generates the set of
all possible six residue sections that could conceivably form a parallel beta sheet
and calculates the RMSD distance between the configuration in which the residues find themselves
and an idealized parallel beta sheet structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures.

This colvar is based on the following reference \cite pietrucci09jctc.  The authors of
this paper use the set of distances from the parallel beta sheet configurations to measure
the number of segments whose configuration resembles a parallel beta sheet. This is done by calculating
the following sum of functions of the rmsd distances:

\f[
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
\f]

where the sum runs over all possible segments of parallel beta sheet.  By default the
NN, MM and D_0 parameters are set equal to those used in \cite pietrucci09jctc.  The R_0
parameter must be set by the user - the value used in \cite pietrucci09jctc was 0.08 nm.

If you change the function in the above sum you can calculate quantities such as the average
distance from a structure composed of only parallel beta sheets or the distance between the set of
residues that is closest to a parallel beta sheet and the reference configuration. To do these sorts of
calculations you can use the AVERAGE and MIN keywords. In addition you can use the LESS_THAN
keyword if you would like to change the form of the switching function. If you use any of these
options you no longer need to specify NN, R_0, MM and D_0.

Please be aware that for codes like gromacs you must ensure that plumed
reconstructs the chains involved in your CV when you calculate this CV using
anything other than TYPE=DRMSD.  For more details as to how to do this see \ref WHOLEMOLECULES.

\par Examples

The following input calculates the number of six residue segments of
protein that are in an parallel beta sheet configuration.

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=beta.pdb
pb: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1
\endplumedfile

Here the same is done use RMSD instead of DRMSD

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
hh: PARABETARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1  STRANDS_CUTOFF=1
\endplumedfile

*/
//+ENDPLUMEDOC

class ParabetaRMSD : public SecondaryStructureRMSD {
public:
  static void registerKeywords( Keywords& keys );
  explicit ParabetaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ParabetaRMSD,"PARABETARMSD")

void ParabetaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureRMSD::registerKeywords( keys );
  keys.add("compulsory","STYLE","all","Parallel beta sheets can either form in a single chain or from a pair of chains. If STYLE=all all "
           "chain configuration with the appropriate geometry are counted.  If STYLE=inter "
           "only sheet-like configurations involving two chains are counted, while if STYLE=intra "
           "only sheet-like configurations involving a single chain are counted");
  keys.use("STRANDS_CUTOFF");
}

ParabetaRMSD::ParabetaRMSD(const ActionOptions&ao):
  Action(ao),
  SecondaryStructureRMSD(ao)
{
  // read in the backbone atoms
  std::vector<unsigned> chains; readBackboneAtoms( "protein", chains );

  bool intra_chain(false), inter_chain(false);
  std::string style; parse("STYLE",style);
  if( style=="all" ) {
    intra_chain=true; inter_chain=true;
  } else if( style=="inter") {
    intra_chain=false; inter_chain=true;
  } else if( style=="intra") {
    intra_chain=true; inter_chain=false;
  } else {
    error( style + " is not a valid directive for the STYLE keyword");
  }

  // Align the atoms based on the positions of these two atoms
  setAtomsFromStrands( 6, 21 );

  // This constructs all conceivable sections of antibeta sheet in the backbone of the chains
  if( intra_chain ) {
    unsigned nprevious=0; std::vector<unsigned> nlist(30);
    for(unsigned i=0; i<chains.size(); ++i) {
      if( chains[i]<40 ) error("segment of backbone is not long enough to form an antiparallel beta hairpin. Each backbone fragment must contain a minimum of 8 residues");
      // Loop over all possible triples in each 8 residue segment of protein
      unsigned nres=chains[i]/5;
      if( chains[i]%5!=0 ) error("backbone segment received does not contain a multiple of five residues");
      for(unsigned ires=0; ires<nres-8; ires++) {
        for(unsigned jres=ires+6; jres<nres-2; jres++) {
          for(unsigned k=0; k<15; ++k) {
            nlist[k]=nprevious + ires*5+k;
            nlist[k+15]=nprevious + jres*5+k;
          }
          addColvar( nlist );
        }
      }
      nprevious+=chains[i];
    }
  }
  // This constructs all conceivable sections of antibeta sheet that form between chains
  if( inter_chain ) {
    if( chains.size()==1 && style!="all" ) error("there is only one chain defined so cannot use inter_chain option");
    std::vector<unsigned> nlist(30);
    for(unsigned ichain=1; ichain<chains.size(); ++ichain) {
      unsigned iprev=0; for(unsigned i=0; i<ichain; ++i) iprev+=chains[i];
      unsigned inres=chains[ichain]/5;
      if( chains[ichain]%5!=0 ) error("backbone segment received does not contain a multiple of five residues");
      for(unsigned ires=0; ires<inres-2; ++ires) {
        for(unsigned jchain=0; jchain<ichain; ++jchain) {
          unsigned jprev=0; for(unsigned i=0; i<jchain; ++i) jprev+=chains[i];
          unsigned jnres=chains[jchain]/5;
          if( chains[jchain]%5!=0 ) error("backbone segment received does not contain a multiple of five residues");
          for(unsigned jres=0; jres<jnres-2; ++jres) {
            for(unsigned k=0; k<15; ++k) {
              nlist[k]=iprev + ires*5+k;
              nlist[k+15]=jprev + jres*5+k;
            }
            addColvar( nlist );
          }
        }
      }
    }
  }

  // Build the reference structure ( in angstroms )
  std::vector<Vector> reference(30);
  reference[0]=Vector( 1.244, -4.620, -2.127); // N    i
  reference[1]=Vector(-0.016, -4.500, -1.395); // CA
  reference[2]=Vector( 0.105, -5.089,  0.024); // CB
  reference[3]=Vector(-0.287, -3.000, -1.301); // C
  reference[4]=Vector( 0.550, -2.245, -0.822); // O
  reference[5]=Vector(-1.445, -2.551, -1.779); // N    i+1
  reference[6]=Vector(-1.752, -1.130, -1.677); // CA
  reference[7]=Vector(-2.113, -0.550, -3.059); // CB
  reference[8]=Vector(-2.906, -0.961, -0.689); // C
  reference[9]=Vector(-3.867, -1.738, -0.695); // O
  reference[10]=Vector(-2.774,  0.034,  0.190); // N    i+2
  reference[11]=Vector(-3.788,  0.331,  1.201); // CA
  reference[12]=Vector(-3.188,  0.300,  2.624); // CB
  reference[13]=Vector(-4.294,  1.743,  0.937); // C
  reference[14]=Vector(-3.503,  2.671,  0.821); // O
  reference[15]=Vector( 4.746, -2.363,  0.188); // N    j
  reference[16]=Vector( 3.427, -1.839,  0.545); // CA
  reference[17]=Vector( 3.135, -1.958,  2.074); // CB
  reference[18]=Vector( 3.346, -0.365,  0.181); // C
  reference[19]=Vector( 4.237,  0.412,  0.521); // O
  reference[20]=Vector( 2.261,  0.013, -0.487); // N    j+1
  reference[21]=Vector( 2.024,  1.401, -0.875); // CA
  reference[22]=Vector( 1.489,  1.514, -2.313); // CB
  reference[23]=Vector( 0.914,  1.902,  0.044); // C
  reference[24]=Vector(-0.173,  1.330,  0.052); // O
  reference[25]=Vector( 1.202,  2.940,  0.828); // N    j+2
  reference[26]=Vector( 0.190,  3.507,  1.718); // CA
  reference[27]=Vector( 0.772,  3.801,  3.104); // CB
  reference[28]=Vector(-0.229,  4.791,  1.038); // C
  reference[29]=Vector( 0.523,  5.771,  0.996); // O
  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( reference, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );

  reference[0]=Vector(-1.439, -5.122, -1.144); // N    i
  reference[1]=Vector(-0.816, -3.803, -1.013); // CA
  reference[2]=Vector( 0.099, -3.509, -2.206); // CB
  reference[3]=Vector(-1.928, -2.770, -0.952); // C
  reference[4]=Vector(-2.991, -2.970, -1.551); // O
  reference[5]=Vector(-1.698, -1.687, -0.215); // N    i+1
  reference[6]=Vector(-2.681, -0.613, -0.143); // CA
  reference[7]=Vector(-3.323, -0.477,  1.267); // CB
  reference[8]=Vector(-1.984,  0.681, -0.574); // C
  reference[9]=Vector(-0.807,  0.921, -0.273); // O
  reference[10]=Vector(-2.716,  1.492, -1.329); // N    i+2
  reference[11]=Vector(-2.196,  2.731, -1.883); // CA
  reference[12]=Vector(-2.263,  2.692, -3.418); // CB
  reference[13]=Vector(-2.989,  3.949, -1.433); // C
  reference[14]=Vector(-4.214,  3.989, -1.583); // O
  reference[15]=Vector( 2.464, -4.352,  2.149); // N    j
  reference[16]=Vector( 3.078, -3.170,  1.541); // CA
  reference[17]=Vector( 3.398, -3.415,  0.060); // CB
  reference[18]=Vector( 2.080, -2.021,  1.639); // C
  reference[19]=Vector( 0.938, -2.178,  1.225); // O
  reference[20]=Vector( 2.525, -0.886,  2.183); // N    j+1
  reference[21]=Vector( 1.692,  0.303,  2.346); // CA
  reference[22]=Vector( 1.541,  0.665,  3.842); // CB
  reference[23]=Vector( 2.420,  1.410,  1.608); // C
  reference[24]=Vector( 3.567,  1.733,  1.937); // O
  reference[25]=Vector( 1.758,  1.976,  0.600); // N    j+2
  reference[26]=Vector( 2.373,  2.987, -0.238); // CA
  reference[27]=Vector( 2.367,  2.527, -1.720); // CB
  reference[28]=Vector( 1.684,  4.331, -0.148); // C
  reference[29]=Vector( 0.486,  4.430, -0.415); // O
  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( reference, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );
}

}
}

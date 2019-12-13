/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

//+PLUMEDOC COLVAR ANTIBETARMSD
/*
Probe the antiparallel beta sheet content of your protein structure.

Two protein segments containing three contiguous residues can form an antiparallel beta sheet.
Although if the two segments are part of the same protein chain they must be separated by
a minimum of 2 residues to make room for the turn. This colvar thus generates the set of
all possible six residue sections that could conceivably form an antiparallel beta sheet
and calculates the RMSD distance between the configuration in which the residues find themselves
and an idealized antiparallel beta sheet structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures.

This colvar is based on the following reference \cite pietrucci09jctc.  The authors of
this paper use the set of distances from the anti parallel beta sheet configurations to measure
the number of segments that have an configuration that resembles an anti parallel beta sheet. This is done by calculating
the following sum of functions of the rmsd distances:

\f[
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
\f]

where the sum runs over all possible segments of antiparallel beta sheet.  By default the
NN, MM and D_0 parameters are set equal to those used in \cite pietrucci09jctc.  The R_0
parameter must be set by the user - the value used in \cite pietrucci09jctc was 0.08 nm.

If you change the function in the above sum you can calculate quantities such as the average
distance from a purely configuration composed of pure anti-parallel beta sheets or the distance between the set of
residues that is closest to an anti-parallel beta sheet and the reference configuration. To do these sorts of
calculations you can use the AVERAGE and MIN keywords. In addition you can use the LESS_THAN
keyword if you would like to change the form of the switching function. If you use any of these
options you no longer need to specify NN, R_0, MM and D_0.

Please be aware that for codes like gromacs you must ensure that plumed
reconstructs the chains involved in your CV when you calculate this CV using
anything other than TYPE=DRMSD.  For more details as to how to do this see \ref WHOLEMOLECULES.

\par Examples

The following input calculates the number of six residue segments of
protein that are in an antiparallel beta sheet configuration.

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=beta.pdb
ab: ANTIBETARMSD RESIDUES=all STRANDS_CUTOFF=1
\endplumedfile

Here the same is done use RMSD instead of DRMSD

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
hh: ANTIBETARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1  STRANDS_CUTOFF=1
\endplumedfile
*/
//+ENDPLUMEDOC

class AntibetaRMSD : public SecondaryStructureRMSD {
public:
  static void registerKeywords( Keywords& keys );
  explicit AntibetaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AntibetaRMSD,"ANTIBETARMSD")

void AntibetaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureRMSD::registerKeywords( keys );
  keys.add("compulsory","STYLE","all","Antiparallel beta sheets can either form in a single chain or from a pair of chains. If STYLE=all all "
           "chain configuration with the appropriate geometry are counted.  If STYLE=inter "
           "only sheet-like configurations involving two chains are counted, while if STYLE=intra "
           "only sheet-like configurations involving a single chain are counted");
  keys.use("STRANDS_CUTOFF");
}

AntibetaRMSD::AntibetaRMSD(const ActionOptions&ao):
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
      for(unsigned ires=0; ires<nres-7; ires++) {
        for(unsigned jres=ires+7; jres<nres; jres++) {
          for(unsigned k=0; k<15; ++k) {
            nlist[k]=nprevious + ires*5+k;
            nlist[k+15]=nprevious + (jres-2)*5+k;
          }
          addColvar( nlist );
        }
      }
      nprevious+=chains[i];
    }
  }
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
              nlist[k]=iprev+ ires*5+k;
              nlist[k+15]=jprev+ jres*5+k;
            }
            addColvar( nlist );
          }
        }
      }
    }
  }

  // Build the reference structure ( in angstroms )
  std::vector<Vector> reference(30);
  reference[0]=Vector( 2.263, -3.795,  1.722); // N    i
  reference[1]=Vector( 2.493, -2.426,  2.263); // CA
  reference[2]=Vector( 3.847, -1.838,  1.761); // CB
  reference[3]=Vector( 1.301, -1.517,  1.921); // C
  reference[4]=Vector( 0.852, -1.504,  0.739); // O
  reference[5]=Vector( 0.818, -0.738,  2.917); // N    i+1
  reference[6]=Vector(-0.299,  0.243,  2.748); // CA
  reference[7]=Vector(-1.421, -0.076,  3.757); // CB
  reference[8]=Vector( 0.273,  1.680,  2.854); // C
  reference[9]=Vector( 0.902,  1.993,  3.888); // O
  reference[10]=Vector( 0.119,  2.532,  1.813); // N    i+2
  reference[11]=Vector( 0.683,  3.916,  1.680); // CA
  reference[12]=Vector( 1.580,  3.940,  0.395); // CB
  reference[13]=Vector(-0.394,  5.011,  1.630); // C
  reference[14]=Vector(-1.459,  4.814,  0.982); // O
  reference[15]=Vector(-2.962,  3.559, -1.359); // N    j-2
  reference[16]=Vector(-2.439,  2.526, -2.287); // CA
  reference[17]=Vector(-1.189,  3.006, -3.087); // CB
  reference[18]=Vector(-2.081,  1.231, -1.520); // C
  reference[19]=Vector(-1.524,  1.324, -0.409); // O
  reference[20]=Vector(-2.326,  0.037, -2.095); // N    j-1
  reference[21]=Vector(-1.858, -1.269, -1.554); // CA
  reference[22]=Vector(-3.053, -2.199, -1.291); // CB
  reference[23]=Vector(-0.869, -1.949, -2.512); // C
  reference[24]=Vector(-1.255, -2.070, -3.710); // O
  reference[25]=Vector( 0.326, -2.363, -2.072); // N    j
  reference[26]=Vector( 1.405, -2.992, -2.872); // CA
  reference[27]=Vector( 2.699, -2.129, -2.917); // CB
  reference[28]=Vector( 1.745, -4.399, -2.330); // C
  reference[29]=Vector( 1.899, -4.545, -1.102); // O

  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( reference, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );
}

}
}

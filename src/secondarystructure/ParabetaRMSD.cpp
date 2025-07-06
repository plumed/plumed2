/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "SecondaryStructureBase.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace secondarystructure {

//+PLUMEDOC COLVAR PARABETARMSD
/*
Probe the parallel beta sheet content of your protein structure.

Two protein segments containing three contiguous residues can form a parallel beta sheet.
Although if the two segments are part of the same protein chain they must be separated by
a minimum of 3 residues to make room for the turn. This colvar thus generates the set of
all possible six residue sections that could conceivably form a parallel beta sheet
and calculates the [DRMSD](DRMSD.md) or [RMSD](RMSD.md) distance between the configuration in which the residues find themselves
and an idealized parallel beta sheet structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures.

This colvar is based on the reference in the bibliography below.  The authors of
this paper use the set of distances from the parallel beta sheet configurations to measure
the number of segments whose configuration resembles a parallel beta sheet. This is done by calculating
the following sum of functions of the rmsd distances:

$$
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
$$

where the sum runs over all possible segments of parallel beta sheet.  By default the
NN, MM and D_0 parameters are set equal to those used in the paper cited below.  The R_0
parameter must be set by the user - the value used in the paper cited below was 0.08 nm.

If you change the function in the above sum you can calculate quantities such as the average
distance from a structure composed of only parallel beta sheets or the distance between the set of
residues that is closest to a parallel beta sheet and the reference configuration. To do these sorts of
calculations you can use the AVERAGE and MIN keywords. In addition you can use the LESS_THAN
keyword if you would like to change the form of the switching function. If you use any of these
options you no longer need to specify NN, R_0, MM and D_0.

The following input calculates the number of six residue segments of
protein that are in an parallel beta sheet configuration.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
pb: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1 R_0=0.1
PRINT ARG=pb FILE=colvar
```

Here the same is done use [RMSD](RMSD.md) instead of [DRMSD](DRMSD.md), which is normally faster.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
WHOLEMOLECULES ENTITY0=1-100
hh: PARABETARMSD RESIDUES=all TYPE=OPTIMAL LESS_THAN={RATIONAL R_0=0.1 NN=8 MM=12}  STRANDS_CUTOFF=1
PRINT ARG=hh.lessthan FILE=colvar
```

__YOUR CALCULATION WILL BE MUCH FASTER IF YOU USE THE `STRANDS_CUTOFF` KEYWORD.__  As you can see from the
expanded version of the inputs above this keyword reduces the computational cost of the calculation by
avoiding calculations of the RMSD values for segments that have the two strands of the beta sheet further apart
than a cutoff.

## Periodic boundary conditions

You can turn off periodic boundary conditions by using the NOPBC flag as shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
pb: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1 R_0=0.1 NOPBC
PRINT ARG=pb FILE=colvar
```

If you are using [DRMSD](DRMSD.md) to measure distances and you __don't__ use the NOPBC flag then
all distances in the instaneous structure are evaluated in a way that takes the periodic boundary conditions
into account. Using the NOPBC flag turns off this treatment.

If you are using [RMSD](RMSD.md) to measure distances and you __don't__ use the NOPBC flag the the instaneous positions of
each segment for which the RMSD is computed is reconstructed using the procedure outlined in the documentation for [WHOLEMOLECULES](WHOLEMOLECULES.md)
before the RMSD is computed. Using the NOPBC flag turns off this treatment.

*/
//+ENDPLUMEDOC

class ParabetaRMSD : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit ParabetaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ParabetaRMSD,"PARABETARMSD")

void ParabetaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureBase<Vector>::registerKeywords( keys );
  keys.remove("ATOMS");
  keys.remove("SEGMENT");
  keys.remove("STRUCTURE");
  keys.remove("MASK");
  keys.remove("ALIGN_STRANDS");
  keys.setValueDescription("scalar/vector","if LESS_THAN is present the RMSD distance between each residue and the ideal parallel beta sheet.  If LESS_THAN is not present the number of residue segments where the structure is similar to a parallel beta sheet");
  keys.add("compulsory","STYLE","all","Parallel beta sheets can either form in a single chain or from a pair of chains. If STYLE=all all "
           "chain configuration with the appropriate geometry are counted.  If STYLE=inter "
           "only sheet-like configurations involving two chains are counted, while if STYLE=intra "
           "only sheet-like configurations involving a single chain are counted");
  keys.add("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart then the calculation "
           "of the actual RMSD is skipped as the structure is very far from being beta-sheet like. "
           "This keyword speeds up the calculation enormously when you are using the LESS_THAN option. "
           "However, if you are using some other option, then this cannot be used");
  keys.needsAction("LOWEST");
  keys.needsAction("DISTANCE");
  keys.needsAction("CUSTOM");
}

ParabetaRMSD::ParabetaRMSD(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the input and create a string that describes how to compute the less than
  std::string ltmap;
  bool uselessthan=SecondaryStructureBase<Vector>::readShortcutWords( ltmap, this );
  // read in the backbone atoms
  std::vector<unsigned> chains;
  std::vector<std::string> all_atoms;
  SecondaryStructureBase<Vector>::readBackboneAtoms( this, plumed, "protein", chains, all_atoms );

  bool intra_chain(false), inter_chain(false);
  std::string style;
  parse("STYLE",style);
  if( Tools::caseInSensStringCompare(style, "all") ) {
    intra_chain=true;
    inter_chain=true;
  } else if( Tools::caseInSensStringCompare(style, "inter") ) {
    intra_chain=false;
    inter_chain=true;
  } else if( Tools::caseInSensStringCompare(style, "intra") ) {
    intra_chain=true;
    inter_chain=false;
  } else {
    error( style + " is not a valid directive for the STYLE keyword");
  }

  double strands_cutoff=-1.0;
  parse("STRANDS_CUTOFF",strands_cutoff);
  std::string scutoff_action;
  if( strands_cutoff>0 ) {
    scutoff_action=getShortcutLabel() + "_cut_dists: DISTANCE ";
  }

  // This constructs all conceivable sections of antibeta sheet in the backbone of the chains
  std::string seglist;
  unsigned jjkk=1;
  if( intra_chain ) {
    unsigned nprevious=0;
    std::vector<unsigned> nlist(30);
    for(unsigned i=0; i<chains.size(); ++i) {
      if( chains[i]<40 ) {
        error("segment of backbone is not long enough to form an antiparallel beta hairpin. Each backbone fragment must contain a minimum of 8 residues");
      }
      // Loop over all possible triples in each 8 residue segment of protein
      unsigned nres=chains[i]/5;
      if( chains[i]%5!=0 ) {
        error("backbone segment received does not contain a multiple of five residues");
      }
      for(unsigned ires=0; ires<nres-8; ires++) {
        for(unsigned jres=ires+6; jres<nres-2; jres++) {
          for(unsigned k=0; k<15; ++k) {
            nlist[k]=nprevious + ires*5+k;
            nlist[k+15]=nprevious + jres*5+k;
          }
          std::string nlstr, num;
          Tools::convert( nlist[0], nlstr );
          Tools::convert(jjkk, num);
          jjkk++;
          seglist += " SEGMENT" + num + "=" + nlstr;
          for(unsigned kk=1; kk<nlist.size(); ++kk ) {
            Tools::convert( nlist[kk], nlstr );
            seglist += "," + nlstr;
          }
          if( strands_cutoff>0 ) {
            scutoff_action += " ATOMS" + num + "=" + all_atoms[nlist[6]] + "," + all_atoms[nlist[21]];
          }
        }
      }
      nprevious+=chains[i];
    }
  }
  // This constructs all conceivable sections of antibeta sheet that form between chains
  if( inter_chain ) {
    if( chains.size()==1 && !Tools::caseInSensStringCompare(style, "all") ) {
      error("there is only one chain defined so cannot use inter_chain option");
    }
    std::vector<unsigned> nlist(30);
    for(unsigned ichain=1; ichain<chains.size(); ++ichain) {
      unsigned iprev=0;
      for(unsigned i=0; i<ichain; ++i) {
        iprev+=chains[i];
      }
      unsigned inres=chains[ichain]/5;
      if( chains[ichain]%5!=0 ) {
        error("backbone segment received does not contain a multiple of five residues");
      }
      for(unsigned ires=0; ires<inres-2; ++ires) {
        for(unsigned jchain=0; jchain<ichain; ++jchain) {
          unsigned jprev=0;
          for(unsigned i=0; i<jchain; ++i) {
            jprev+=chains[i];
          }
          unsigned jnres=chains[jchain]/5;
          if( chains[jchain]%5!=0 ) {
            error("backbone segment received does not contain a multiple of five residues");
          }
          for(unsigned jres=0; jres<jnres-2; ++jres) {
            for(unsigned k=0; k<15; ++k) {
              nlist[k]=iprev + ires*5+k;
              nlist[k+15]=jprev + jres*5+k;
            }
            std::string nlstr, num;
            Tools::convert( nlist[0], nlstr );
            Tools::convert(jjkk, num);
            jjkk++;
            seglist += " SEGMENT" + num + "=" + nlstr;
            for(unsigned kk=1; kk<nlist.size(); ++kk ) {
              Tools::convert( nlist[kk], nlstr );
              seglist += "," + nlstr;
            }
            if( strands_cutoff>0 ) {
              scutoff_action += " ATOMS" + num + "=" + all_atoms[nlist[6]] + "," + all_atoms[nlist[21]];
            }
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
  std::string ref0, ref1, ref2;
  Tools::convert(  reference[0][0], ref0 );
  Tools::convert(  reference[0][1], ref1 );
  Tools::convert(  reference[0][2], ref2 );
  std::string structure=" STRUCTURE1=" + ref0 + "," + ref1 + "," + ref2;
  for(unsigned i=1; i<30; ++i) {
    for(unsigned k=0; k<3; ++k) {
      Tools::convert( reference[i][k], ref0 );
      structure += "," + ref0;
    }
  }

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
  Tools::convert(  reference[0][0], ref0 );
  Tools::convert(  reference[0][1], ref1 );
  Tools::convert(  reference[0][2], ref2 );
  structure +=" STRUCTURE2=" + ref0 + "," + ref1 + "," + ref2;
  for(unsigned i=1; i<30; ++i) {
    for(unsigned k=0; k<3; ++k) {
      Tools::convert( reference[i][k], ref0 );
      structure += "," + ref0;
    }
  }

  std::string nopbcstr="";
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  if( nopbc ) {
    nopbcstr = " NOPBC";
  }
  std::string usegpustr="";
  {
    bool usegpu;
    parseFlag("USEGPU",usegpu);
    if( usegpu ) {
      usegpustr = " USEGPU";
    }
  }
  std::string type;
  parse("TYPE",type);
  std::string lab = getShortcutLabel() + "_low";
  if( uselessthan ) {
    lab = getShortcutLabel();
  }
  std::string atoms="ATOMS=" + all_atoms[0];
  for(unsigned i=1; i<all_atoms.size(); ++i) {
    atoms += "," + all_atoms[i];
  }
  std::string inputLine =  getShortcutLabel() + "_both:";
  if( type=="DRMSD" ) {
    inputLine+=" SECONDARY_STRUCTURE_DRMSD BONDLENGTH=0.17";
  } else {
    inputLine+=" SECONDARY_STRUCTURE_RMSD TYPE=" +type;
  }
  inputLine+=" ALIGN_STRANDS " + seglist + structure
             + " " + atoms + nopbcstr + usegpustr;

  if( strands_cutoff>0 ) {
    readInputLine( scutoff_action );
    std::string str_cut;
    Tools::convert( strands_cutoff, str_cut );
    readInputLine( getShortcutLabel() + "_cut: CUSTOM ARG=" + getShortcutLabel()
                   + "_cut_dists FUNC=step(" + str_cut + "-x) PERIODIC=NO");
    inputLine+=" MASK=" + getShortcutLabel() + "_cut";
    readInputLine(inputLine);
    if( ltmap.length()>0 ) {
      // Create the lowest line
      readInputLine( lab + ": LOWEST ARG=" + getShortcutLabel() + "_both.struct-1,"
                     + getShortcutLabel() + "_both.struct-2 MASK=" + getShortcutLabel() + "_cut");
      // Create the less than object
      readInputLine( getShortcutLabel() + "_ltu: LESS_THAN ARG=" + lab
                     + " SWITCH={" + ltmap  +"} MASK=" + getShortcutLabel() + "_cut");
      // Multiply by the strands cutoff
      readInputLine( getShortcutLabel()
                     + "_lt: CUSTOM ARG=" + getShortcutLabel() + "_ltu,"
                     + getShortcutLabel() + "_cut"
                     " FUNC=x*y PERIODIC=NO");
    }
  } else {
    readInputLine(inputLine);
    if( ltmap.length()>0 ) {
      // Create the lowest line
      readInputLine( lab + ": LOWEST ARG=" + getShortcutLabel() + "_both.struct-1,"
                     + getShortcutLabel() + "_both.struct-2" );
      // Create the less than object
      readInputLine( getShortcutLabel() + "_lt: LESS_THAN ARG=" + lab + " SWITCH={" + ltmap  +"}");
    }
  }
  // Create the less than object
  if( ltmap.length()>0 ) {
    if( uselessthan ) {
      readInputLine( getShortcutLabel() + "_lessthan: SUM ARG=" + getShortcutLabel() + "_lt PERIODIC=NO");
    } else {
      readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_lt PERIODIC=NO");
    }
  }
}

}
}

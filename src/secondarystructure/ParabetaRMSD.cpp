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
a minimum of three residues to make room for the turn. This colvar thus generates the set of
all possible six residue sections (residues $i$, $i+1$, $i+2$ and $j$, $j+1$, $j+2$) that could conceivably form a parallel beta sheet
and calculates the [DRMSD](DRMSD.md) or [RMSD](RMSD.md) distance between the configuration in which the residues find themselves
and an idealized parallel beta sheet structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures. The backbone atoms N, CA, CB, C and O
of each residue are used (for glycine the HA2 atom is used in place of CB), so each section contains 30 atoms.

Two reference structures are used for the parallel beta sheet, and for each segment the distance is the lowest of the
two distances (the vector `label_low` in the expanded input). The vectors of the distances from the two references are the
components `struct-1` and `struct-2` of the action labelled `label_both` in the expanded input.

This colvar is based on the paper cited below.  The authors of
this paper use the set of distances from the parallel beta sheet configurations to measure
the number of segments whose configuration resembles a parallel beta sheet. This is done by calculating
the following sum of functions of the rmsd distances:

$$
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
$$

where the sum runs over all possible segments of parallel beta sheet.  By default the
NN, MM, D_0 and R_0 parameters are set equal to those used in the paper cited below, which are
$n=8$, $m=12$, $d_0=0$ and $r_0=0.08$ nm. The value of R_0 is given in the length units of [UNITS](UNITS.md);
the default is converted, so it is always 0.08 nm. As in the other actions that take R_0, NN, MM and D_0, the function above, $f(r)$,
is stretched so that it goes to zero at $d_{max}=d_0+r_0\,10^{5/(m-n)}$, where $f(d_{max})=10^{-5}$ (1.42 nm with the default
parameters): each segment contributes $\left(f(r_i)-f(d_{max})\right)/\left(1-f(d_{max})\right)$ if $r_i<d_{max}$ and zero otherwise.
The value of $d_{max}$ is the D_MAX of the switching function in the expanded input. The segments with a distance from the
reference larger than $d_{max}$, whose strands are far apart, are then not computed (see below).
To use the function without stretching, use LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12}.

If you would like to use a different switching function you can use the LESS_THAN keyword, as in the second example below,
where D_MAX is given as for the default function, so that the pairs of strands that are too far apart can be skipped also with this function.
In this case the label of the action is the vector of the distances between each segment and the ideal parallel beta sheet, and
the number of segments that are similar to a parallel beta sheet is the `lessthan` component.

The following input calculates the number of six residue segments of
protein that are in a parallel beta sheet configuration.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
pb: PARABETARMSD RESIDUES=all R_0=0.1
PRINT ARG=pb FILE=colvar
```

Here the same is done using [RMSD](RMSD.md) instead of [DRMSD](DRMSD.md), and with the LESS_THAN keyword. The cost of the two methods is similar.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
hh: PARABETARMSD RESIDUES=all TYPE=OPTIMAL LESS_THAN={RATIONAL R_0=0.1 NN=8 MM=12 D_MAX=1.78}
PRINT ARG=hh.lessthan FILE=colvar
```

## Strands cutoff

The segments that have the two strands of the beta sheet further apart than `STRANDS_CUTOFF` are not computed, and they contribute
zero to the number of segments that are similar to a beta sheet. The distance between the two strands is the distance between the
CA atoms of their central residues. The default cutoff is 1 nm (the cutoff is in the length units of [UNITS](UNITS.md), the default
is converted). This makes the calculation much faster, as the number of pairs of strands grows with the square of the number of residues.
Use `STRANDS_CUTOFF=0` to compute all the pairs: the number of segments then also includes the small contributions of the many segments
whose strands are far apart (with 18 peptides of 8 residues about 1 more segment). If you use LESS_THAN together with STRANDS_CUTOFF,
notice that the distances of the skipped segments are zero in the vector that is the value of the action, so this vector
should not be used to find the segment that is closest to a beta sheet.

## How the number of segments is computed

Without `USEGPU`, the number of segments that are similar to a sheet (the value of the action, or its `lessthan` component when
LESS_THAN is used) is computed by a single action, PARABETARMSD_SCALAR, which loops over the pairs of strands, skips those that are further
apart than `STRANDS_CUTOFF` and adds the derivatives of each segment directly to the atoms. When the switching function goes to zero at
a finite distance (as the default one, or with D_MAX in LESS_THAN), it also skips the pairs of strands that are so far apart that the
distance of the segment from the reference is certainly larger than that distance (this does not change the result, and makes
STRANDS_CUTOFF=0 cheaper). The vectors that appear in the expanded input
(the distances of the strands, the distances from the reference structure and their switching functions) are only computed when they are
used elsewhere in the input, e.g. when the vector that is the value of the action with LESS_THAN is printed. The result is the same.
With `USEGPU` and TYPE=DRMSD the same calculation is done on the GPU by PARABETARMSD_SCALAR_ACC, from the openaccPTM plugin: one GPU thread per pair of strands within `STRANDS_CUTOFF`. With `USEGPU` and the other types the number of segments is the sum of the vector of the switching functions.

## Chains and residues

The STYLE keyword selects whether the two strands of a segment are taken from the same chain (STYLE=intra), from different
chains (STYLE=inter) or both (STYLE=all, the default). With STYLE=intra and STYLE=all each chain must contain at least 8 residues,
otherwise PLUMED stops with an error; for systems made of short peptides (for example in studies of aggregation) use STYLE=inter.

## Periodic boundary conditions

Unless the NOPBC flag is used, each six residue section is made whole before computing its distance from the reference.
The two strands of three residues are made whole separately: the atoms N, CA, CB, C and O of each strand are taken in this order,
and each atom is moved to the periodic image that is closest to the previous one, as is done by [WHOLEMOLECULES](WHOLEMOLECULES.md).
The second strand is then moved to the periodic image where the CA atom of its central residue is closest to the CA atom of the
central residue of the first strand (this is the ALIGN_STRANDS flag of [SECONDARY_STRUCTURE_RMSD](SECONDARY_STRUCTURE_RMSD.md)).
This is done in the same way for [DRMSD](DRMSD.md) and [RMSD](RMSD.md), so there is no need to use [WHOLEMOLECULES](WHOLEMOLECULES.md)
before this action. The distance used with STRANDS_CUTOFF is also computed with the periodic boundary conditions.
You can turn off this treatment by using the NOPBC flag as shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
pb: PARABETARMSD RESIDUES=all R_0=0.1 NOPBC
PRINT ARG=pb FILE=colvar
```

*/
//+ENDPLUMEDOC

class ParabetaRMSD : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit ParabetaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ParabetaRMSD,"PARABETARMSD")

void ParabetaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureInput::registerKeywords( keys );
  keys.remove("ATOMS");
  keys.remove("SEGMENT");
  keys.remove("STRUCTURE");
  keys.remove("MASK");
  keys.remove("ALIGN_STRANDS");
  keys.setValueDescription("scalar/vector","if LESS_THAN is present the vector of the distances between each segment of six residues and the ideal parallel beta sheet. If LESS_THAN is not present the number of segments where the structure is similar to a parallel beta sheet");
  keys.add("compulsory","STYLE","all","Parallel beta sheets can either form in a single chain or from a pair of chains. If STYLE=all all "
           "chain configuration with the appropriate geometry are counted.  If STYLE=inter "
           "only sheet-like configurations involving two chains are counted, while if STYLE=intra "
           "only sheet-like configurations involving a single chain are counted");
  keys.add("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart than this distance (between the CA atoms of their central residues) "
           "the calculation of the distance from the reference is skipped, as the structure is very far from being beta-sheet like, and the segment does not "
           "contribute to the number of segments that are similar to a beta sheet. The default is 1 nm (the value is in the length units of UNITS, "
           "the default is converted). Use STRANDS_CUTOFF=0 to compute all the pairs of strands");
  keys.needsAction("LOWEST");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_SCALAR_ACC");
  keys.needsAction("DISTANCE");
  keys.needsAction("CUSTOM");
}

ParabetaRMSD::ParabetaRMSD(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the input and create a string that describes how to compute the less than
  std::string ltmap;
  bool uselessthan=SecondaryStructureInput::readShortcutWords( ltmap, this );
  // read in the backbone atoms
  std::vector<unsigned> chains;
  std::vector<std::string> all_atoms;
  SecondaryStructureInput::readBackboneAtoms( this, plumed, "protein", chains, all_atoms );

  std::string style;
  parse("STYLE",style);
  // All the conceivable pairs of strands in the backbone of the chains
  std::vector<std::pair<unsigned,unsigned> > pairs;
  SecondaryStructureInput::getSheetPairs( this, chains, true, style, pairs );

  // 1 nm by default, in the length units of UNITS; zero or negative to use all the pairs of strands
  double strands_cutoff=1.0/plumed.getUnits().getLength();
  parse("STRANDS_CUTOFF",strands_cutoff);
  std::string scutoff_action;
  if( strands_cutoff>0 ) {
    scutoff_action=getShortcutLabel() + "_cut_dists: DISTANCE ";
  }
  std::string seglist;
  std::vector<unsigned> nlist(30);
  for(unsigned k=0; k<pairs.size(); ++k) {
    for(unsigned kk=0; kk<15; ++kk) {
      nlist[kk]=pairs[k].first+kk;
      nlist[kk+15]=pairs[k].second+kk;
    }
    SecondaryStructureInput::addSegment( seglist, k+1, nlist );
    if( strands_cutoff>0 ) {
      // the central CA atoms of the two strands
      scutoff_action += " ATOMS" + std::to_string(k+1) + "=" + all_atoms[nlist[6]] + "," + all_atoms[nlist[21]];
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
  std::string structure;
  SecondaryStructureInput::addStructure( structure, 1, reference );

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
  SecondaryStructureInput::addStructure( structure, 2, reference );

  std::string type;
  bool nopbc=false, usegpu=false;
  std::string distanceAction = SecondaryStructureInput::readDistanceAction( this, type, nopbc, usegpu );
  std::string lab = getShortcutLabel() + "_low";
  if( uselessthan ) {
    lab = getShortcutLabel();
  }
  std::string atoms="ATOMS=" + all_atoms[0];
  for(unsigned i=1; i<all_atoms.size(); ++i) {
    atoms += "," + all_atoms[i];
  }
  std::string inputLine = getShortcutLabel() + "_both: " + distanceAction + " ALIGN_STRANDS" + seglist + structure + " " + atoms;

  if( strands_cutoff>0 ) {
    // The distance used for the cutoff has to treat the periodic boundary conditions as the rest of the calculation
    readInputLine( scutoff_action + (nopbc ? " NOPBC" : "") );
    std::string str_cut;
    Tools::convert( strands_cutoff, str_cut );
    readInputLine( getShortcutLabel() + "_cut: CUSTOM ARG=" + getShortcutLabel()
                   + "_cut_dists FUNC=step(" + str_cut + "-x) PERIODIC=NO");
    inputLine+=" MASK=" + getShortcutLabel() + "_cut";
    readInputLine(inputLine);
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
  } else {
    readInputLine(inputLine);
    // Create the lowest line
    readInputLine( lab + ": LOWEST ARG=" + getShortcutLabel() + "_both.struct-1,"
                   + getShortcutLabel() + "_both.struct-2" );
    // Create the less than object
    readInputLine( getShortcutLabel() + "_lt: LESS_THAN ARG=" + lab + " SWITCH={" + ltmap  +"}");
  }
  if( usegpu && type!="DRMSD" ) {
    // there is no GPU version of the single action with the RMSD
    SecondaryStructureInput::addSum( this, uselessthan );
  } else {
    // The number of segments that are similar to a sheet is computed by a single action: without the vectors above,
    // which are only computed if they are used elsewhere in the input. With USEGPU this is its GPU version
    std::string chainsizes=std::to_string(chains[0]);
    for(unsigned i=1; i<chains.size(); ++i) {
      chainsizes += "," + std::to_string(chains[i]);
    }
    std::string sumline = getShortcutLabel() + (uselessthan ? "_lessthan" : "") + ": PARABETARMSD_SCALAR" + (usegpu ? "_ACC " : " ") + atoms
                          + " CHAIN_SIZES=" + chainsizes + " STYLE=" + style + " TYPE=" + type + " SWITCH={" + ltmap + "}" + structure;
    if( strands_cutoff>0 ) {
      std::string str_cut;
      Tools::convert( strands_cutoff, str_cut );
      sumline += " STRANDS_CUTOFF=" + str_cut;
    }
    if( nopbc ) {
      sumline += " NOPBC";
    }
    readInputLine( sumline );
  }
}

}
}

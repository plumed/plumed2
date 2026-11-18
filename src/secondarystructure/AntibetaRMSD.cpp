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

//+PLUMEDOC COLVAR ANTIBETARMSD
/*
Probe the antiparallel beta sheet content of your protein structure.

Two protein segments containing three contiguous residues can form an antiparallel beta sheet.
Although if the two segments are part of the same protein chain they must be separated by
a minimum of two residues to make room for the turn. This colvar thus generates the set of
all possible six residue sections (residues $i$, $i+1$, $i+2$ and $j-2$, $j-1$, $j$) that could conceivably form an antiparallel beta sheet
and calculates the [DRMSD](DRMSD.md) or [RMSD](RMSD.md) distance between the configuration in which the residues find themselves
and an idealized antiparallel beta sheet structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures. The backbone atoms N, CA, CB, C and O
of each residue are used (for glycine the HA2 atom is used in place of CB), so each section contains 30 atoms.

This colvar is based on the paper cited below.  The authors of
this paper use the set of distances from the antiparallel beta sheet configurations to measure
the number of segments whose configuration resembles an antiparallel beta sheet. This is done by calculating
the following sum of functions of the rmsd distances:

$$
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
$$

where the sum runs over all possible segments of antiparallel beta sheet.  By default the
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
In this case the label of the action is the vector of the distances between each segment and the ideal antiparallel beta sheet, and
the number of segments that are similar to an antiparallel beta sheet is the `lessthan` component.

The following input calculates the number of six residue segments of
protein that are in an antiparallel beta sheet configuration.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
ab: ANTIBETARMSD RESIDUES=all R_0=0.1
PRINT ARG=ab FILE=colvar
```

Here the same is done using [RMSD](RMSD.md) instead of [DRMSD](DRMSD.md), and with the LESS_THAN keyword. The cost of the two methods is similar.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
hh: ANTIBETARMSD RESIDUES=all TYPE=OPTIMAL LESS_THAN={RATIONAL R_0=0.1 NN=8 MM=12 D_MAX=1.78}
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
LESS_THAN is used) is computed by a single action, ANTIBETARMSD_SCALAR, which loops over the pairs of strands, skips those that are further
apart than `STRANDS_CUTOFF` and adds the derivatives of each segment directly to the atoms. When the switching function goes to zero at
a finite distance (as the default one, or with D_MAX in LESS_THAN), it also skips the pairs of strands that are so far apart that the
distance of the segment from the reference is certainly larger than that distance (this does not change the result, and makes
STRANDS_CUTOFF=0 cheaper). The vectors that appear in the expanded input
(the distances of the strands, the distances from the reference structure and their switching functions) are only computed when they are
used elsewhere in the input, e.g. when the vector that is the value of the action with LESS_THAN is printed. The result is the same.
With `USEGPU` and TYPE=DRMSD the same calculation is done on the GPU by ANTIBETARMSD_SCALAR_ACC, from the openaccPTM plugin: one GPU thread per pair of strands within `STRANDS_CUTOFF`. With `USEGPU` and the other types the number of segments is the sum of the vector of the switching functions.

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
ab: ANTIBETARMSD RESIDUES=all R_0=0.1 NOPBC
PRINT ARG=ab FILE=colvar
```

*/
//+ENDPLUMEDOC

class AntibetaRMSD : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit AntibetaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AntibetaRMSD,"ANTIBETARMSD")

void AntibetaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureInput::registerKeywords( keys );
  keys.remove("ATOMS");
  keys.remove("SEGMENT");
  keys.remove("STRUCTURE");
  keys.remove("MASK");
  keys.remove("ALIGN_STRANDS");
  keys.setValueDescription("scalar/vector","if LESS_THAN is present the vector of the distances between each segment of six residues and the ideal antiparallel beta sheet. If LESS_THAN is not present the number of segments where the structure is similar to an antiparallel beta sheet");
  keys.add("compulsory","STYLE","all","Antiparallel beta sheets can either form in a single chain or from a pair of chains. If STYLE=all all "
           "chain configuration with the appropriate geometry are counted.  If STYLE=inter "
           "only sheet-like configurations involving two chains are counted, while if STYLE=intra "
           "only sheet-like configurations involving a single chain are counted");
  keys.add("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart than this distance (between the CA atoms of their central residues) "
           "the calculation of the distance from the reference is skipped, as the structure is very far from being beta-sheet like, and the segment does not "
           "contribute to the number of segments that are similar to a beta sheet. The default is 1 nm (the value is in the length units of UNITS, "
           "the default is converted). Use STRANDS_CUTOFF=0 to compute all the pairs of strands");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_SCALAR_ACC");
  keys.needsAction("DISTANCE");
  keys.needsAction("CUSTOM");
}

AntibetaRMSD::AntibetaRMSD(const ActionOptions&ao):
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
  SecondaryStructureInput::getSheetPairs( this, chains, false, style, pairs );

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
  std::string structure;
  SecondaryStructureInput::addStructure( structure, 1, reference );

  std::string type;
  bool nopbc=false, usegpu=false;
  std::string distanceAction = SecondaryStructureInput::readDistanceAction( this, type, nopbc, usegpu );
  std::string lab = getShortcutLabel() + "_rmsd";
  if( uselessthan ) {
    lab = getShortcutLabel();
  }
  std::string atoms="ATOMS=" + all_atoms[0];
  for(unsigned i=1; i<all_atoms.size(); ++i) {
    atoms += "," + all_atoms[i];
  }
  std::string inputLine = lab + ": " + distanceAction + " ALIGN_STRANDS" + seglist + structure + " " + atoms;

  if( strands_cutoff>0 ) {
    // The distance used for the cutoff has to treat the periodic boundary conditions as the rest of the calculation
    readInputLine( scutoff_action + (nopbc ? " NOPBC" : "") );
    std::string str_cut;
    Tools::convert( strands_cutoff, str_cut );
    readInputLine( getShortcutLabel() + "_cut: CUSTOM ARG=" + getShortcutLabel()
                   + "_cut_dists FUNC=step(" + str_cut + "-x) PERIODIC=NO");
    inputLine+=" MASK=" + getShortcutLabel() + "_cut";
    readInputLine(inputLine);
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
    std::string sumline = getShortcutLabel() + (uselessthan ? "_lessthan" : "") + ": ANTIBETARMSD_SCALAR" + (usegpu ? "_ACC " : " ") + atoms
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

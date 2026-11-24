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

//+PLUMEDOC COLVAR ALPHARMSD
/*
Probe the alpha helical content of a protein structure.

Any chain of six contiguous residues in a protein chain can form an alpha helix. This
colvar thus generates the set of all possible six residue sections and calculates
the [DRMSD](DRMSD.md) or [RMSD](RMSD.md) distance between the configuration in which the residues find themselves
and an idealized alpha helical structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures. The backbone atoms N, CA, CB, C and O
of each residue are used (for glycine the HA2 atom is used in place of CB), so each section contains 30 atoms.

In the original paper on this method (see below) the authors used the set of distances from the alpha helix configurations to measure
the number of segments that have an alpha helical configuration. This is done by calculating
the following sum of functions of the rmsd distances:

$$
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
$$

where the sum runs over all possible segments of alpha helix.  By default the
NN, MM, D_0 and R_0 parameters are set equal to those used in the paper cited below, which are
$n=8$, $m=12$, $d_0=0$ and $r_0=0.08$ nm. The value of R_0 is given in the length units of [UNITS](UNITS.md);
the default is converted, so it is always 0.08 nm. As in the other actions that take R_0, NN, MM and D_0, the function above, $f(r)$,
is stretched so that it goes to zero at $d_{max}=d_0+r_0\,10^{5/(m-n)}$, where $f(d_{max})=10^{-5}$ (1.42 nm with the default
parameters): each segment contributes $\left(f(r_i)-f(d_{max})\right)/\left(1-f(d_{max})\right)$ if $r_i<d_{max}$ and zero otherwise.
The value of $d_{max}$ is the D_MAX of the switching function in the expanded input. To use the function without stretching,
use LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12}.

If you would like to use a different switching function you can use the LESS_THAN keyword, as in the second example below,
where D_MAX is given as for the default function.
In this case the label of the action is the vector of the distances between each segment and the ideal alpha helix, and
the number of segments that are similar to an alpha helix is the `lessthan` component. The vector can be given in input to other
actions: for example [LOWEST](LOWEST.md) gives the distance of the segment that is closest to an alpha helix.

The following input calculates the number of six residue segments of
protein that are in an alpha helical configuration.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
alpha: ALPHARMSD RESIDUES=all R_0=0.1
PRINT ARG=alpha FILE=colvar
```

Here the same is done using [RMSD](RMSD.md) instead of [DRMSD](DRMSD.md), and with the LESS_THAN keyword. The cost of the two methods is similar.
The distance of the segment that is closest to an alpha helix is also computed.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
alpha: ALPHARMSD RESIDUES=all TYPE=OPTIMAL LESS_THAN={RATIONAL R_0=0.1 NN=8 MM=12 D_MAX=1.78}
closest: LOWEST ARG=alpha
PRINT ARG=alpha.lessthan,closest FILE=colvar
```

## Periodic boundary conditions

Unless the NOPBC flag is used, each six residue section is made whole before computing its distance from the reference: the
atoms N, CA, CB, C and O of the six residues are taken in this order, and each atom is moved to the periodic image that is closest
to the previous one, as is done by [WHOLEMOLECULES](WHOLEMOLECULES.md). This is done in the same way for [DRMSD](DRMSD.md) and
[RMSD](RMSD.md), so there is no need to use [WHOLEMOLECULES](WHOLEMOLECULES.md) before this action.
You can turn off this treatment by using the NOPBC flag as shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
alpha: ALPHARMSD RESIDUES=all R_0=0.1 NOPBC
PRINT ARG=alpha FILE=colvar
```

## How the number of segments is computed

The number of segments that are similar to an alpha helix (the value of the action, or its `lessthan` component when
LESS_THAN is used) is computed by a single action, ALPHARMSD_SCALAR, which loops over the segments and adds the derivatives
of each segment directly to the atoms, so that the segments are computed once per step. The vectors that appear in the expanded
input (the distances from the reference structure and their switching functions) are only computed when they are used
elsewhere in the input, e.g. when the vector that is the value of the action with LESS_THAN is printed. The result is the same.
With `USEGPU` and TYPE=DRMSD the same calculation is done on the GPU by ALPHARMSD_SCALAR_ACC, from the openaccPTM plugin.
With `USEGPU` and the other types the number of segments is the sum of the vector of the switching functions.

*/
//+ENDPLUMEDOC

class AlphaRMSD : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit AlphaRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AlphaRMSD,"ALPHARMSD")

void AlphaRMSD::registerKeywords( Keywords& keys ) {
  SecondaryStructureInput::registerKeywords( keys );
  keys.remove("ATOMS");
  keys.remove("SEGMENT");
  keys.remove("STRUCTURE");
  keys.remove("MASK");
  keys.remove("ALIGN_STRANDS");
  keys.setValueDescription("scalar/vector","if LESS_THAN is present the vector of the distances between each segment of six residues and the ideal alpha helix. If LESS_THAN is not present the number of segments where the structure is similar to an alpha helix");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_SCALAR_ACC");
}

AlphaRMSD::AlphaRMSD(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the input and create a string that describes how to compute the less than
  std::string ltmap;
  bool uselessthan=SecondaryStructureInput::readShortcutWords( ltmap, this );
  // read in the backbone atoms
  std::vector<unsigned> chains;
  std::vector<std::string> all_atoms;
  SecondaryStructureInput::readBackboneAtoms( this, plumed, "protein", chains, all_atoms );

  // This constructs all conceivable sections of alpha helix in the backbone of the chains
  unsigned nprevious=0, segno=1;
  std::string seglist;
  std::vector<unsigned> nlist(30);
  for(unsigned i=0; i<chains.size(); ++i) {
    if( chains[i]<30 ) {
      error("segment of backbone defined is not long enough to form an alpha helix. Each backbone fragment must contain a minimum of 6 residues");
    }
    unsigned nres=chains[i]/5;
    if( chains[i]%5!=0 ) {
      error("backbone segment received does not contain a multiple of five residues");
    }
    for(unsigned ires=0; ires<nres-5; ires++) {
      for(unsigned k=0; k<30; ++k) {
        nlist[k] = nprevious + 5*ires + k;
      }
      SecondaryStructureInput::addSegment( seglist, segno, nlist );
      segno++;
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
  readInputLine( lab + ": " + distanceAction + seglist + structure + " " + atoms );
  // Create the less than object
  readInputLine( getShortcutLabel() + "_lt: LESS_THAN ARG=" + lab + " SWITCH={" + ltmap  +"}");
  if( usegpu && type!="DRMSD" ) {
    // there is no GPU version of the single action with the RMSD
    SecondaryStructureInput::addSum( this, uselessthan );
  } else {
    // The number of segments that are similar to a helix is computed by a single action: the vectors above are
    // only computed if they are used elsewhere in the input. With USEGPU this is its GPU version
    std::string chainsizes=std::to_string(chains[0]);
    for(unsigned i=1; i<chains.size(); ++i) {
      chainsizes += "," + std::to_string(chains[i]);
    }
    std::string sumline = getShortcutLabel() + (uselessthan ? "_lessthan" : "") + ": ALPHARMSD_SCALAR" + (usegpu ? "_ACC " : " ") + atoms
                          + " CHAIN_SIZES=" + chainsizes + " TYPE=" + type + " SWITCH={" + ltmap + "}" + structure;
    if( nopbc ) {
      sumline += " NOPBC";
    }
    readInputLine( sumline );
  }
}

}
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Rangsiman Ketkaew

   The xsprint module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The xsprint module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MATRIXF XSPRINT
/*
Calculate eXtended SPRINT (xSPRINT) topological variables from adjacency matrices over multiple shells.

The xSPRINT coordinates, introduced in \cite ketkaew2022, extend the standard [SPRINT](SPRINT.md) coordinates
by incorporating structural information from multiple concentric spherical shells around each atom.
For a system divided into $N$ shells, the xSPRINT coordinate of atom $i$ is defined as:

$$
S_i^x = \begin{cases}
  S_i, & \text{if } r_{ij} < r_0 \\
  \frac{1}{\omega} \sum_{n=1}^{N} r_n\, S_i^{(n)}, & \text{otherwise}
\end{cases}
$$

where $S_i^{(n)}$ is the SPRINT coordinate computed from the adjacency matrix of the $n$-th shell,
$r_n = n\,r_0$ is the shell radius, $r_0$ is the user-defined cutoff for the first shell, and
$\omega = \left(\sum_{n=1}^{N} r_n\right)^2$ is a normalization weight.  Adjacency matrix elements
from inner shells are zeroed out before computing each shell's SPRINT to avoid redundancy.

When `NSHELLS=1` the xSPRINT coordinate reduces to the standard SPRINT coordinate.

The following example computes the xSPRINT coordinates for a 7-atom cluster using 3 shells
with a shell radius of 0.1:

```plumed
xs: XSPRINT GROUP=1-7 SWITCH={RATIONAL R_0=0.1} R0=0.1 NSHELLS=3
PRINT ARG=xs.* FILE=colvar
```

Note that `R0` defines the xSPRINT shell radius ($r_0$ in the equation) and is independent
of `R_0` inside `SWITCH`, which controls the shape of the switching function.  For shell $n$,
the switching function `R_0` is automatically set to $n \times$ `R0`.

This example computes xSPRINT with 2 shells for two groups of atoms:

```plumed
xs: XSPRINT ...
  GROUP1=1-7 GROUP2=8-14
  SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}
  SWITCH12={RATIONAL R_0=2.2 NN=6 MM=12}
  SWITCH22={RATIONAL R_0=2.2 NN=6 MM=12}
  R0=2.0
  NSHELLS=2
...

PRINT ARG=xs.* FILE=colvar
```

You can also provide a pre-built matrix (for a single-shell case this is equivalent to SPRINT):

```plumed
b: BRIDGE_MATRIX GROUP=1-7 BRIDGING_ATOMS=8-100 SWITCH={RATIONAL R_0=0.2}
xs: XSPRINT MATRIX=b NSHELLS=1
PRINT ARG=xs.* FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace xsprint {

class Xsprint : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Xsprint(const ActionOptions&);
};

// Function to replace all R_0=<value> entries in a string with R_0=<newval>
static std::string replaceR0(const std::string& input, double newval) {
  std::string result = input;
  std::string needle = "R_0=";
  std::string newvalstr;
  Tools::convert(newval, newvalstr);
  size_t pos = 0;
  while((pos = result.find(needle, pos)) != std::string::npos) {
    size_t valstart = pos + needle.length();
    size_t valend = result.find_first_of(" \t},", valstart);
    if(valend == std::string::npos) valend = result.length();
    result = result.substr(0, valstart) + newvalstr + result.substr(valend);
    pos = valstart + newvalstr.length();
  }
  return result;
}

PLUMED_REGISTER_ACTION(Xsprint,"XSPRINT")

void Xsprint::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","MATRIX","the matrix that you would like to perform xSPRINT on (used only when NSHELLS=1)");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms");
  keys.add("compulsory","NSHELLS","1","the number of concentric shells to use for the xSPRINT calculation");
  keys.add("optional","R0","the radius of the first shell for xSPRINT (shell n has radius n*R0). Required when NSHELLS>1. This is independent of the R_0 in SWITCH which controls the switching function shape.");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("SPRINT");
  keys.needsAction("CUSTOM");
  keys.needsAction("SELECT_COMPONENTS");
  keys.needsAction("SORT");
  keys.needsAction("COMBINE");
  keys.addOutputComponent("coord","default","scalar","the xsprint coordinates");
  keys.addDOI("10.1021/acs.jpclett.1c04004");
}

Xsprint::Xsprint(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  const std::string label = getShortcutLabel();

  unsigned nshells = 1;
  parse("NSHELLS", nshells);
  std::string matinp;
  parse("MATRIX", matinp);
  if( !matinp.empty() && nshells > 1 ) {
    plumed_merror("MATRIX keyword can only be used with NSHELLS=1");
  }
  double r0 = 0.0;
  parse("R0", r0);
  if( nshells > 1 && r0 <= 0.0 ) {
    plumed_merror("R0 keyword is required and must be positive when NSHELLS > 1");
  }

  // NSHELLS == 1 : use existing SPRINT
  if( nshells == 1 ) {
    if( !matinp.empty() ) {
      readInputLine( label + "_sprint: SPRINT MATRIX=" + matinp );
    } else {
      readInputLine( label + "_sprint: SPRINT " + convertInputLineToString() );
    }
    // Forward SPRINT's coord-N outputs as our own coord-N outputs
    for(unsigned i = 0;; ++i) {
      std::string inum;
      Tools::convert( i, inum );
      ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>(
        label + "_sprint_coord-" + inum );
      if( !av ) break;
      readInputLine( label + "_coord-" + inum + ": COMBINE ARG=" +
        label + "_sprint_coord-" + inum + " PERIODIC=NO" );
    }
    return;
  }

  // NSHELLS > 1 : compute xSPRINT using SPRINT per shell
  std::string baseargs = convertInputLineToString();

  // Shell 1: create CONTACT_MATRIX with R_0=r0 and compute SPRINT
  std::string shell1args = replaceR0(baseargs, r0);
  readInputLine( label + "_shell1_mat: CONTACT_MATRIX " + shell1args );
  readInputLine( label + "_shell1_sprint: SPRINT MATRIX=" + label + "_shell1_mat" );

  // Determine group sizes from shell-1 matrix (this is needed for final sorting)
  std::string mat1 = label + "_shell1_mat";
  std::vector<unsigned> nin_group;
  for(unsigned i = 1;; ++i) {
    std::string inum;
    Tools::convert( i, inum );
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( mat1 + inum + inum );
    if( !av ) break;
    nin_group.push_back( (av->copyOutput(0))->getShape()[0] );
  }
  if( nin_group.empty() ) {
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( mat1 );
    nin_group.push_back( (av->copyOutput(0))->getShape()[0] );
  }

  // For shells 2..N: create CONTACT_MATRIX with R_0=s*r0, 
  // subtract inner shell, compute SPRINT, respectively
  for(unsigned s = 2; s <= nshells; ++s) {
    std::string snum, sprev;
    Tools::convert( s, snum );
    Tools::convert( s - 1, sprev );

    double shell_r0 = static_cast<double>(s) * r0;
    std::string shellargs = replaceR0(baseargs, shell_r0);
    readInputLine( label + "_shell" + snum + "_mat: CONTACT_MATRIX " + shellargs );

    // Subtract inner-shell matrix to isolate contacts in this shell only
    readInputLine( label + "_shell" + snum + "_diff: CUSTOM ARG=" +
                   label + "_shell" + snum + "_mat," +
                   label + "_shell" + sprev + "_mat FUNC=x-y PERIODIC=NO" );

    readInputLine( label + "_shell" + snum + "_sprint: SPRINT MATRIX=" +
                   label + "_shell" + snum + "_diff" );
  }

  // Combine per-shell SPRINT vectors with weights r_n / omega
  // where r_n = n*r0, omega = (sum_{n=1}^{N} r_n)^2
  double sum_rn = r0 * static_cast<double>(nshells) * static_cast<double>(nshells + 1) / 2.0;
  double omega = sum_rn * sum_rn;

  std::string combine_args, coefficients;
  for(unsigned s = 1; s <= nshells; ++s) {
    std::string snum;
    Tools::convert( s, snum );
    if( s > 1 ) {
      combine_args += ",";
      coefficients += ",";
    }
    combine_args += label + "_shell" + snum + "_sprint_sp";
    std::string wstr;
    Tools::convert( static_cast<double>(s) * r0 / omega, wstr );
    coefficients += wstr;
  }

  readInputLine( label + "_xsp: COMBINE ARG=" + combine_args +
                 " COEFFICIENTS=" + coefficients + " PERIODIC=NO" );

  // Sort xSPRINT coordinates per group
  unsigned k = 0, kk = 0;
  for(unsigned j = 0; j < nin_group.size(); ++j) {
    std::string jnum, knum;
    Tools::convert( j+1, jnum );
    Tools::convert( k+1, knum );
    k++;
    std::string sort_act = label + "_selection" + jnum + ": SELECT_COMPONENTS ARG=" + label + "_xsp COMPONENTS=" + knum;
    for(unsigned n = 1; n < nin_group[j]; ++n) {
      Tools::convert( k+1, knum );
      sort_act += "," + knum;
      k++;
    }
    readInputLine( sort_act );
    readInputLine( label + jnum + ": SORT ARG=" + label + "_selection" + jnum );
    for(unsigned n = 0; n < nin_group[j]; ++n) {
      std::string nnum;
      Tools::convert( kk, knum );
      Tools::convert( n+1, nnum );
      kk++;
      readInputLine( label + "_coord-" + knum + ": COMBINE ARG=" + label + jnum + "." + nnum + " PERIODIC=NO" );
    }
  }
}

}
}

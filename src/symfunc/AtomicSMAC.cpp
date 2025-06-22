/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR ATOMIC_SMAC
/*
Calculate the atomic smac CV

This shortcut offers another example of a [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html).
This action was inspired by [SMAC](SMAC.md) and the distribution of angles in the first coordination sphere around an atom to determine if the
environment around the atom is ordered. For atom $i$ the symmetry function is calculated using:

$$
s_i = [ 1 - \gamma(c_i) ] \frac{ \sum_j \sum{k \ne j} \sigma(r_{ij})\sigma(r_{ik}) \sum_n G\left( \frac{ \theta_{jik} - \phi_n}{b_n}\right) }{\sum{k \ne j} \sigma(r_{ij})\sigma(r_{ik})} \qquad \textrm{where} \qquad c_i = \sum_j \sigma(r_{ij})
$$

In this expression $r_{ij}$ is the distance between atom $i$ and atom $j$ and $\sigma$ is a switching function that acts upon this distance.   $c_i$ is thus the number of atoms that are within
a certain cutoff of atom $i$ and $\gamma$ is another switching function that acts upon this quantity.  This switching function ensures that the symmetry function is zero for atoms that are
regions where the density is low.  $\theta_{jik}$ is the angle between the vector connecting atoms $i$ and $j$ and the vector connecting atoms $i$ and $k$.  This angle is the argument for the
set of Gaussian kernel functions, $G$, that are centered on $\phi_n$ and that have bandwidths of $b_n$.  The function above is thus determining if the angles between the bonds in the first coordination
sphere around atom $i$ are similar to the $\phi_n$ values that have been specified by the user or not.

The following example demonstrates how this symmetry function can be used in practise.

```plumed
smac: ATOMIC_SMAC SPECIES=1-64 KERNEL1={GAUSSIAN CENTER=pi/2 SIGMA=1.0} KERNEL2={GAUSSIAN CENTER=pi/4 SIGMA=1.0} SWITCH_COORD={EXP R_0=4.0} SWITCH={RATIONAL R_0=2.0 D_0=2.0} SUM
PRINT ARG=smac.* FILE=colvar
```

The input above would calculate 64 instances of $s_i$ using the formula above.  In each of these two Gaussian Kernels are used in the sum over $n$.  The parameters for the switching function
$\sigma$ are specified using the SWITCH keyword, while the parameters for $\gamma$ are specified using SWITCH_COORD.  As you can see if you expand the shortcut in the input above, the 64 values
for $s_i$ are stored in a vector.  All the elements of this vector are then added together to produce the single quantity that is output in the colvar file.

Notice that you can also use different groups of atoms by using the SPECIESA and SPECIESB keywords as shown below:

```plumed
smac: ATOMIC_SMAC SPECIESA=1-64 SPECIESB=65-128 KERNEL1={GAUSSIAN CENTER=pi/2 SIGMA=1.0} KERNEL2={GAUSSIAN CENTER=pi/4 SIGMA=1.0} SWITCH_COORD={EXP R_0=4.0} SWITCH={RATIONAL R_0=2.0 D_0=2.0} SUM
PRINT ARG=smac.* FILE=colvar
```

In this input one $s_i$ value is evaluated for each of the 64 atoms specified in the SPECIESA keyword.  In evaluating each of these $s_i$ values using the sum in the expression above runs over the 64 atoms
that were specified using SPECIESB.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class AtomicSMAC : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit AtomicSMAC(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AtomicSMAC,"ATOMIC_SMAC")

void AtomicSMAC::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms-3","SPECIES","the list of atoms for which the symmetry function is being calculated and the atoms that can be in the environments");
  keys.add("atoms-4","SPECIESA","the list of atoms for which the symmetry function is being calculated.  This keyword must be used in conjunction with SPECIESB, which specifies the atoms that are in the environment.");
  keys.add("atoms-4","SPECIESB","the list of atoms that can be in the environments of each of the atoms for which the symmetry function is being calculated.  This keyword must be used in conjunction with SPECIESA, which specifies the atoms for which the symmetry function is being calculated.");
  keys.add("optional","SWITCH","the switching function that it used in the construction of the contact matrix");
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.add("optional","SWITCH_COORD","This keyword is used to define the coordination switching function.");
  keys.reset_style("KERNEL","optional");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("GSYMFUNC_THREEBODY");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

AtomicSMAC::AtomicSMAC(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create the matrices
  std::string sw_input;
  parse("SWITCH",sw_input);
  std::string sp_lab, sp_laba;
  parse("SPECIES",sp_lab);
  parse("SPECIESA",sp_laba);
  std::string cmap_input = getShortcutLabel() + "_cmap: CONTACT_MATRIX";
  if( sp_lab.length()>0 ) {
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_lab + " COMPONENTS SWITCH={" + sw_input + "}");
  } else if( sp_laba.length()>0 ) {
    std::string sp_labb;
    parse("SPECIESB",sp_labb);
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + sp_laba + " GROUPB=" + sp_labb + " COMPONENTS SWITCH={" + sw_input + "}");
  }
  // Now need the Gaussians
  std::string mykernels;
  for(unsigned i=1;; ++i) {
    std::string kstr_inpt, istr, kern_str;
    Tools::convert( i, istr );
    if( !parseNumbered("KERNEL",i,kstr_inpt ) ) {
      break;
    }
    std::vector<std::string> words = Tools::getWords(kstr_inpt);
    if( words[0]=="GAUSSIAN" ) {
      kern_str="gaussian";
    } else {
      error("unknown kernel type");
    }
    std::string center, var;
    Tools::parse(words,"CENTER",center);
    Tools::parse(words,"SIGMA",var);
    if( mykernels.length()==0 ) {
      mykernels = "exp(-(ajik-" + center + ")^2/(2*" + var + "*" + var + "))";
    } else {
      mykernels = mykernels + "+exp(-(ajik-" + center + ")^2/(2*" + var + "*" + var + "))";
    }
  }
  // Hard coded switching function on minimum distance here -- should be improved
  readInputLine( getShortcutLabel() + "_ksum: GSYMFUNC_THREEBODY WEIGHT=" + getShortcutLabel() + "_cmap.w " +
                 "ARG=" + getShortcutLabel() + "_cmap.x," + getShortcutLabel() + "_cmap.y," + getShortcutLabel() + "_cmap.z"
                 " FUNCTION1={FUNC=" + mykernels + " LABEL=n} FUNCTION2={FUNC=1 LABEL=d}" );
  // And just the sum of the coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmap");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cmap.w," + getShortcutLabel() + "_ones");
  // And the transformed switching functions
  std::string swcoord_str;
  parse("SWITCH_COORD",swcoord_str);
  readInputLine( getShortcutLabel() + "_mtdenom: MORE_THAN ARG=" + getShortcutLabel() + "_denom SWITCH={" + swcoord_str +"}");
  // And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_smac: CUSTOM ARG=" + getShortcutLabel() + "_ksum.n," + getShortcutLabel() + "_mtdenom," + getShortcutLabel() + "_ksum.d FUNC=x*y/z PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_smac", "", this );
}

}
}

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
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR SMAC
/*
Calculate the SMAC order parameter for a set of molecules

This shortcut action provides a quantity that can be thought of as a [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html) for molecules.
You can thus use it to detect whether the molecules in a sample are ordered in the way that they would expect to be ordered in the crystal. If you are using this CV you normally use
the [COM](COM.md) command to define centers of mass for your molecules and the [DISTANCE](DISTANCE.md) or [PLANE](PLANE.md) command to define each molecules orientation.  You can thus calculate
the relative orientation molecules $i$ and $j$ by calculating the torsional angle, $\theta_{ij}$, between the two molecular orientations around the vector $r_{ij}$ connecting the centers of
mass of the two molecules. The value of the molecular order parameter for molecule $i$ is thus calculated as:

$$
s_i = frac{1-\gamma(c_i)}{c_i} \sum_j \sigma(|r_{ij}|) \sum_k K\left( \frac{\theta_{ij} - \phi_k}{b} \right) \qquad \textrm{where} \qquad c_i = \sum_j \sigma(|r_{ij}|)
$$

In this expression $\sigma$ is a switching function that acts on the distance between the centers of mass of molecules $i$ and $j$. $c_i$ is thus the number of molecules that are within
a certain cutoff of molecule $i$ and $\gamma$ is another switching function that acts upon this quantity. This switching function ensures that the symmetry function is zero for atoms that are
regions where the density of molecules are low.  $K$ is then a kernel function with bandwidth $b$ that is centered at $\phi_k$.  The function above is thus only large if molecule $i$ is surrounded
by molecules whose relative orientations are as the user has requested by specifying $\phi_k$ parameters.

The following example illustrates how the SMAC order parameter in PLUMED is used:

```plumed
m3: DISTANCES ...
   ATOMS1=9,10 LOCATION1=9
   ATOMS2=89,90 LOCATION2=89
   ATOMS3=473,474 LOCATION3=473
   ATOMS4=1161,1162 LOCATION4=1161
   ATOMS5=1521,1522 LOCATION5=1521
   ATOMS6=1593,1594 LOCATION6=1593
   ATOMS7=1601,1602 LOCATION7=1601
   ATOMS8=2201,2202 LOCATION8=2201
   COMPONENTS
...

s2: SMAC SPECIES=m3 KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480} KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480} SWITCH={RATIONAL R_0=0.6} MORE_THAN={RATIONAL R_0=0.7} SWITCH_COORD={EXP R_0=4}

PRINT ARG=s2_morethan FILE=colvar
```

Here the orientations of the molecules are specified by calculating the vectors connecting pairs of atoms in the molecules.  The LOCATION keywords in the distance command are used to specify the
positions of the molecules from which the $r_{ij}$ vectors in the above expression are calculated.  The sum over $k$ in the above expression has two terms corresponding to the two Gaussian kernels
that have been specified using KERNEL keywords.  The SWITCH keyword has been used to specify the parameters of the switching function $\sigma$, while the SWITCH_COORD keyword has been used to specify
the parameters of the switching function $\gamma$.

A vector of 8 values for $s_i$ is calculated.  The elements of this vector are transformed by a switching function that is one if the $s_i$ value is larger than 0.7.  The eight elements of the resulting vector
of transformed $s_i$ are then added together to determine the final quantity that is output to the colvar file.

Incidentally, the authors who designed the SMAC symmetry function have forgotten what the letters in this acronym stand for.

## Working with two types of molecule

You can use two different sets of vectors in the input to this action as shown below:

```plumed
m3: DISTANCES ...
   ATOMS1=9,10 LOCATION1=9
   ATOMS2=89,90 LOCATION2=89
   ATOMS3=473,474 LOCATION3=473
   ATOMS4=1161,1162 LOCATION4=1161
   COMPONENTS
...

m4: DISTANCES ...
   ATOMS1=1521,1522 LOCATION1=1521
   ATOMS2=1593,1594 LOCATION2=1593
   ATOMS3=1601,1602 LOCATION3=1601
   ATOMS4=2201,2202 LOCATION4=2201
   COMPONENTS
...

s2: SMAC SPECIESA=m3 SPECIESB=m4 KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480} KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480} SWITCH={RATIONAL R_0=0.6} MORE_THAN={RATIONAL R_0=0.7} SWITCH_COORD={EXP R_0=4}

PRINT ARG=s2_morethan FILE=colvar
```

In this input the $s_i$ values defined using the equation above is evaluated for the four vectors that are specified in the DISTANCE action with label m3.  In evaluating the sum in that expression
we run over the four vectors that are specified in the action with label `m4` above.


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class SMAC : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit SMAC(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SMAC,"SMAC")

void SMAC::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","SPECIES","the label of the DISTANCES or PLANES action that computes the orientations of the molecules for which you would like to compute SMAC.  The coordination sphere contains the same list of molecules if you use this keyword.");
  keys.add("optional","SPECIESA","the label of the DISTANCES or PLANES action that computes the orientations of the molecules for which you would like to compute SMAC.  This keyword must be used with SPECIESB.");
  keys.add("optional","SPECIESB","the label of the DISTANCES or PLANES action that computes the orientations of the molecules that should be considered as part of the coordination sphere.  This keyword must be used with SPECIESA.");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.add("optional","SWITCH_COORD","This keyword is used to define the coordination switching function.");
  keys.linkActionInDocs("SWITCH_COORD","LESS_THAN");
  keys.reset_style("KERNEL","optional");
  keys.setValueDescription("vector","the value of the smac parameter for each of the input molecules");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("VSTACK");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("TORSIONS_MATRIX");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("MORE_THAN");
  keys.addDOI("10.1016/j.ces.2014.08.032");
  keys.addDOI("10.1021/acs.jctc.6b01073");
}

SMAC::SMAC(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create the matrices
  std::string sw_input;
  parse("SWITCH",sw_input);
  std::string sp_lab, sp_laba;
  parse("SPECIES",sp_lab);
  parse("SPECIESA",sp_laba);
  if( sp_lab.length()>0 ) {
    readInputLine( getShortcutLabel() + "_vecs: VSTACK ARG=" + sp_lab + ".x," + sp_lab + ".y," + sp_lab + ".z" );
    readInputLine( getShortcutLabel() + "_vecsT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecs" );
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_lab + " SWITCH={" + sw_input + "}");
    readInputLine( getShortcutLabel() + "_tpmat: TORSIONS_MATRIX ARG=" + getShortcutLabel() + "_vecs," + getShortcutLabel() + "_vecsT POSITIONS1=" + sp_lab + " POSITIONS2=" + sp_lab + " MASK=" + getShortcutLabel() + "_cmap");
  } else if( sp_laba.length()>0 ) {
    std::string sp_labb;
    parse("SPECIESB",sp_labb);
    readInputLine( getShortcutLabel() + "_vecsa: VSTACK ARG=" + sp_laba + ".x," + sp_laba + ".y," + sp_laba + ".z" );
    readInputLine( getShortcutLabel() + "_vecsb: VSTACK ARG=" + sp_labb + ".x," + sp_labb + ".y," + sp_labb + ".z" );
    readInputLine( getShortcutLabel() + "_vecsbT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecsb" );
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + sp_laba + " GROUPB=" + sp_labb + " SWITCH={" + sw_input + "}");
    readInputLine( getShortcutLabel() + "_tpmat: TORSIONS_MATRIX ARG=" + getShortcutLabel() + "_vecsa," + getShortcutLabel() + "_vecsbT POSITIONS1=" + sp_laba + " POSITIONS2=" + sp_labb + " MASK=" + getShortcutLabel() + "_cmap");
  }
  // Now need the Gaussians
  std::string kmap_input= getShortcutLabel() + "_ksum: COMBINE MASK=" + getShortcutLabel() + "_cmap PERIODIC=NO";
  for(unsigned i=1;; ++i) {
    std::string kstr_inpt, istr;
    Tools::convert( i, istr );
    if( !parseNumbered("KERNEL",i,kstr_inpt ) ) {
      break;
    }
    std::vector<std::string> words = Tools::getWords(kstr_inpt);
    std::string center, var;
    Tools::parse(words,"CENTER",center);
    Tools::parse(words,"SIGMA",var);
    double nsig;
    Tools::convert( var, nsig );
    std::string coeff;
    Tools::convert( 1/(nsig*nsig), coeff );
    readInputLine( getShortcutLabel() + "_kf" + istr + "_r2: COMBINE MASK=" + getShortcutLabel() + "_cmap PERIODIC=NO ARG=" + getShortcutLabel() + "_tpmat COEFFICIENTS=" + coeff + " PARAMETERS=" + center + " POWERS=2");
    if( words[0]=="GAUSSIAN" ) {
      readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM MASK=" + getShortcutLabel() + "_cmap PERIODIC=NO FUNC=exp(-x/2) ARG=" +  getShortcutLabel() + "_kf" + istr + "_r2" );
    } else if( words[0]=="TRIANGULAR" ) {
      readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM MASK=" + getShortcutLabel() + "_cmap PERIODIC=NO FUNC=step(1-sqrt(x))*(1-sqrt(x)) ARG=" + getShortcutLabel() + "_kf" + istr + "_r2" );
    } else {
      readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM PERIODIC=NO MASK=" + getShortcutLabel() + "_cmap FUNC=" + words[0] + " ARG=" + getShortcutLabel() + "_kf" + istr + "_r2" );
    }
    if( i==1 ) {
      kmap_input += " ARG=" + getShortcutLabel() + "_kf" + istr;
    } else {
      kmap_input += "," + getShortcutLabel() + "_kf" + istr;
    }
  }
  readInputLine( kmap_input );
  // Now create the product matrix
  readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() + "_ksum FUNC=x*y PERIODIC=NO");
  // Now the sum of coordination numbers times the switching functions
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmap");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_prod," + getShortcutLabel() + "_ones");
  // And just the sum of the coordination numbers
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() + "_ones");
  // And the transformed switching functions
  std::string swcoord_str;
  parse("SWITCH_COORD",swcoord_str);
  readInputLine( getShortcutLabel() + "_mtdenom: MORE_THAN ARG=" + getShortcutLabel() + "_denom SWITCH={" + swcoord_str +"}");
// And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_smac: CUSTOM ARG=" + getShortcutLabel() + "," + getShortcutLabel() + "_mtdenom," + getShortcutLabel() + "_denom FUNC=(x*y)/z PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_smac", "", this );
}

}
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR HEXACTIC_PARAMETER
/*
Calculate the hexatic order parameter

This [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html) can be used to understand phase transitions in two dimensional systems.
The symmetry function for atom $k$ is calculated using:

$$
s_k = \left| \frac{\sum_j \sigma(r_{kj}) e^{6i\theta_j} }{\sum_j \sigma(r_{kj})} \right|
$$

In this expression, the sum run over all the atoms of interest and $r_{kj}$ is the distance between atom $k$ and atom $j$ and $\sigma$ is
a switching function that acts upon this distance.  $\theta_j$ is the angle between either the $x$, $y$ or $z$ axis and the bond connecting
atom $k$ and atom $j$.  This angle is multiplied by the imaginary number $i$ - the square root of minus one.  In the code, we thus calculate
$e^{i\theta_j}$ as follows:

$$
e^{i\theta_j) = \frac{x_{kj}}{r_{kj}} + i \frac{y_{kj}}{r_{kj}}
$$

We then take the 6th power of this complex number directly before compupting the magnitude by multiplying the result by its complex conjugate.
Notice, furthermore, that we can replace $x_{kj}$ or $y_{kj}$ with $z_{kj}$ by using PLANE=xz or PLANE=yz in place of PLANE=xy.

An example that shows how you can use this shortcut is shown below:

```plumed
hex: HEXACTIC_PARAMETER SPECIES=1-400 PLANE=xy SWITCH={RATIONAL D_0=1.4 R_0=0.2} MEAN
PRINT ARG=hex.mean FILE=colvar
```

As you can see if you expand the shortcut above, this input calculates the quantity defined in the equation above for the 400 atoms in the simulated system and stores them in a vector.
The elements of this vector are then added together so the mean value can be computed.

In papers where symmetry functions similar to this one have been used a switching function is not employed. The sums over $j$ in the expression above are replaced by sums over the
six nearest neighbours to each atom.  If you would like to calculate this quantity using PLUMED you can use an input like this:

```plumed
dmat: DISTANCE_MATRIX GROUP=1-400 CUTOFF=3.0 COMPONENTS
neigh: NEIGHBORS ARG=dmat.w NLOWEST=6
harm: CYLINDRICAL_HARMONIC DEGREE=6 ARG=dmat.x,dmat.y
rprod: CUSTOM ARG=neigh,harm.rm FUNC=x*y PERIODIC=NO
iprod: CUSTOM ARG=neigh,harm.im FUNC=x*y PERIODIC=NO
hex2_ones: ONES SIZE=400
hex2_denom: MATRIX_VECTOR_PRODUCT ARG=neigh,hex2_ones
harm_rm: MATRIX_VECTOR_PRODUCT ARG=rprod,hex2_ones
harm_im: MATRIX_VECTOR_PRODUCT ARG=iprod,hex2_ones
hex2_rmn: CUSTOM ARG=harm_rm,hex2_denom FUNC=x/y PERIODIC=NO
hex2_imn: CUSTOM ARG=harm_im,hex2_denom FUNC=x/y PERIODIC=NO
DUMPATOMS ATOMS=1-400 ARG=hex2_rmn,hex2_imn,hex2_denom FILE=hexparam.xyz
```

This input outputs the values of the order parameters for all the atoms to an extended xyz file .

> ![CAUTION]
> Virial is not working currently


*/
//+ENDPLUMEDOC


class HexacticParameter : public ActionShortcut {
private:
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  explicit HexacticParameter(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(HexacticParameter,"HEXACTIC_PARAMETER")

void HexacticParameter::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.add("compulsory","PLANE","the plane to use when calculating the value of the order parameter should be xy, xz or yz");
  keys.setValueDescription("matrix","the value of the cylindrical harmonic for each bond vector specified");
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","scalar","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","scalar","the norm of the mean vector");
  keys.needsAction("CYLINDRICAL_HARMONIC_MATRIX");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("MEAN");
  keys.needsAction("SUM");
  keys.needsAction("COMBINE");
  keys.addDOI("https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.215701");
}

HexacticParameter::HexacticParameter( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this );
  std::string myplane;
  parse("PLANE",myplane);
  if( myplane=="xy" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="xz" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.z," + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="yz" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z," + getShortcutLabel() + "_mat.w" );
  } else {
    error("invalid input for plane -- should be xy, xz or yz");
  }
  // And coordination number
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_rm: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + ".rm," + getShortcutLabel() + "_ones");
  readInputLine( getShortcutLabel() + "_im: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + ".im," + getShortcutLabel() + "_ones");
  // Input for denominator (coord)
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat.w," + getShortcutLabel() + "_ones");
  // Divide real part by coordination numbers
  readInputLine( getShortcutLabel() + "_rmn: CUSTOM ARG=" + getShortcutLabel() + "_rm," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  // Devide imaginary part by coordination number
  readInputLine( getShortcutLabel() + "_imn: CUSTOM ARG=" + getShortcutLabel() + "_im," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");

  // If we are doing VMEAN determine sum of vector components
  bool do_vmean;
  parseFlag("VMEAN",do_vmean);
  if( do_vmean ) {
    // Real part
    readInputLine( getShortcutLabel() + "_rms: MEAN ARG=" + getShortcutLabel() + "_rmn PERIODIC=NO");
    // Imaginary part
    readInputLine( getShortcutLabel() + "_ims: MEAN ARG=" + getShortcutLabel() + "_imn PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", "ms" );
  }
  bool do_vsum;
  parseFlag("VSUM",do_vsum);
  if( do_vsum ) {
    // Real part
    readInputLine( getShortcutLabel() + "_rmz: SUM ARG=" + getShortcutLabel() + "_rmn PERIODIC=NO");
    // Imaginary part
    readInputLine( getShortcutLabel() + "_imz: SUM ARG=" + getShortcutLabel() + "_imn PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", "mz" );
  }

  // Now calculate the total length of the vector
  createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_norm", "mn" );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_norm", "", this );
}

void HexacticParameter::createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab ) {
  readInputLine( olab + "2: COMBINE PERIODIC=NO ARG=" + ilab + "_r" + vlab + "," + ilab + "_i" + vlab + " POWERS=2,2" );
  readInputLine( olab + ": CUSTOM ARG=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}


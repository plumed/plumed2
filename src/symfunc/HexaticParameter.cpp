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
e^{i\theta_j} = \frac{x_{kj}}{r_{kj}} + i \frac{y_{kj}}{r_{kj}}
$$

We then take the 6th power of this complex number directly before compupting the magnitude by multiplying the result by its complex conjugate.
Notice, furthermore, that we can replace $x_{kj}$ or $y_{kj}$ with $z_{kj}$ by using PLANE=xz or PLANE=yz in place of PLANE=xy.

An example that shows how you can use this shortcut is shown below:

```plumed
hex: HEXACTIC_PARAMETER SPECIES=1-400 PLANE=xy R_0=0.2 D_0=1.4 NN=6 MM=12
hex_mean: MEAN ARG=hex_norm PERIODIC=NO
PRINT ARG=hex_mean FILE=colvar
```

As you can see if you expand the shortcut above, this input calculates the quantity defined in the equation above for the 400 atoms in the simulated system and stores them in a vector.
The elements of this vector are then added together so the mean value can be computed.

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:


```plumed
hex: HEXACTIC_PARAMETER SPECIES=1-400 PLANE=xy SWITCH={RATIONAL D_0=1.4 R_0=0.2 D_MAX=3.0}
hex_mean: MEAN ARG=hex_norm PERIODIC=NO
PRINT ARG=hex_mean FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## The VMEAN and VSUM options

If you expand the inputs in the previous section you will see that the value `hex_norm` that we are calculating the average of is a vector that contains the magnitude of the
complex number defined in the equation above.  If you would like to add up the vector of complex numbers directly (or take an average of the vector) __before__ computing the
magnitude of the complex number you can use the VSUM and VMEAN flags as sillstrated below:

```plumed
hex: HEXACTIC_PARAMETER SPECIES=1-400 PLANE=xy SWITCH={RATIONAL D_0=1.4 R_0=0.2} VMEAN VSUM
PRINT ARG=hex_vmean,hex_vsum FILE=colvar
```

If you expand the shortcut in the input above you will see that this input adds the complex numbers together directly before calculating the modulus.  This contrasts with the approach
in the inputs above, which computes a vector that contains the square moduli for each of the individual hexatic order parameters and then computes the average from this vector.

## Using two types of atom

If you would like to calculate the hexatic parameters using the directors of the bonds connecting the atoms in GROUPA to the atoms in GROUPB you can use an input like the one
shown below:

```plumed
hex: HEXACTIC_PARAMETER SPECIESA=1-200 SPECIESB=201-400 PLANE=xy SWITCH={RATIONAL D_0=1.4 R_0=0.2 D_MAX=3.0}
hex_mean: MEAN ARG=hex_norm PERIODIC=NO
PRINT ARG=hex_mean FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the hexatic parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the tetrahedral parameter of the atoms
hex: HEXACTIC_PARAMETER ...
   SPECIES=1-400 MASK=sphere PLANE=xy
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply hexatic parameters by sphere vector
prod_rm: CUSTOM ARG=hex_rm,sphere FUNC=x*y PERIODIC=NO
prod_im: CUSTOM ARG=hex_im,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer_rm: SUM ARG=prod_rm PERIODIC=NO
numer_im: SUM ARG=prod_im PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av_rm: CUSTOM ARG=numer_rm,denom FUNC=x/y PERIODIC=NO
av_im: CUSTOM ARG=numer_im,denom FUNC=x/y PERIODIC=NO
# Take the square modulus
av: CUSTOM ARG=av_rm,av_im FUNC=x*x+y*y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the (complex) hexatic parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.  The square modulus of this average is then output.

## Using nearest neighbours

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

!!! warning "Broken virial"

    Virial is not working currently


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
  keys.addDOI("10.1103/PhysRevLett.99.215701");
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


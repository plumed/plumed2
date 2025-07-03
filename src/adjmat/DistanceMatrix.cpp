/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX DISTANCE_MATRIX
/*
Calculate a matrix of distances between atoms.

To calculate the matrix of distances between every distinct pair of atoms in a single group you can use the following command:

```plumed
d1: DISTANCE_MATRIX GROUP=1-7
```

If you would like to calculate the matrix of distances between the atoms in two different groups of atoms you can use the following command:

```plumed
d2: DISTANCE_MATRIX GROUPA=1-7 GROUPB=8-20
```

For both these inputs the distances between atoms are calculated in a way that takes the periodic boundary conditions into account. If you want to
ignore the periodic boundaries when calculating distances you use the NOPBC flag as shown below:

```plumed
d3: DISTANCE_MATRIX GROUP=1-7 NOPBC
```

Once you have calculated your distance matrix in this way you can do many of the operations that were discussed for [CONTACT_MATRIX](CONTACT_MATRIX.md) with the output.
For example, you can use the COMPONENTS flag to calcuate the $x$, $y$ and $z$ components of the vectors connecting the atoms in your two groups by using
an input like that shown below:

```plumed
d1: DISTANCE_MATRIX GROUP=1-7 COMPONENTS
```

## Optimisation details

If for some reaon, you only want to calculate the distances if they are less than a certain CUTOFF you can add the cutoff keyword as follows:

```plumed
d3: DISTANCE_MATRIX GROUP=1-7 CUTOFF=1.0
```

This command will only store those distances that are less than 1 nm.  Notice, however, that the cutoff implemented here is __not__ a continuous function.
You should thus be careful when commands such as this one above to ensure that any quantities that are forced have continuous derivatives.  If you use
the CUTOFF keyword, however, many of the features that are used to optimise [CONTACT_MATRIX](CONTACT_MATRIX.md) are used for this action.

Also notice that you can use MASK to calculate a subset of the rows in the DISTANCE matrix as is done in the following example:

```plumed
# The atoms that are of interest
ow: GROUP ATOMS=1-16500
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=ow CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# The distance matrix
dmap: DISTANCE_MATRIX COMPONENTS GROUP=ow CUTOFF=1.0 MASK=sphere
# Find the four nearest neighbors
acv_neigh: NEIGHBORS ARG=dmap.w NLOWEST=4 MASK=sphere
# Compute a function for the atoms that are in the first coordination sphere
acv_g8: GSYMFUNC_THREEBODY ...
  WEIGHT=acv_neigh ARG=dmap.x,dmap.y,dmap.z
  FUNCTION1={FUNC=(cos(ajik)+1/3)^2 LABEL=g8}
  MASK=sphere
...
# Now compute the value of the function above for those atoms that are in the
# sphere of interest
acv: CUSTOM ARG=acv_g8.g8,sphere FUNC=y*(1-(3*x/8)) PERIODIC=NO
# And now compute the final average
acv_sum: SUM ARG=acv PERIODIC=NO
acv_norm: SUM ARG=sphere PERIODIC=NO
mean: CUSTOM ARG=acv_sum,acv_norm FUNC=x/y PERIODIC=NO
PRINT ARG=mean FILE=colvar
```

This input calculates the average value for measure of tetrahedral order that is introduced in the documentation for the [TETRA_ANGULAR](TETRA_ANGULAR.md) shortcut
for those atom that are within a sphere that is centered on the point $(2.5,2.5,2.5)$.

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class DistanceMatrix {
public:
  double cutoff;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<DistanceMatrix>* action );
  static void calculateWeight( const DistanceMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
};

typedef AdjacencyMatrixBase<DistanceMatrix> dmap;
PLUMED_REGISTER_ACTION(dmap,"DISTANCE_MATRIX")

void DistanceMatrix::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
}

void DistanceMatrix::parseInput( AdjacencyMatrixBase<DistanceMatrix>* action ) {
  // And set the link cell cutoff
  action->log.printf("  weight is distance between atoms \n");
  action->parse("CUTOFF",cutoff);
  if( cutoff<0 ) {
    action->setLinkCellCutoff( true, std::numeric_limits<double>::max() );
  } else {
    action->log.printf("  ignoring distances that are larger than %f \n", cutoff);
    action->setLinkCellCutoff( true, cutoff );
  }
}

void DistanceMatrix::calculateWeight( const DistanceMatrix& data,
                                      const AdjacencyMatrixInput& input,
                                      MatrixOutput output ) {
  output.val[0] = input.pos.modulo();
  if( data.cutoff<0 || output.val[0]<data.cutoff ) {
    double invd = 1.0/output.val[0];
    Vector v = (-invd)*input.pos;
    output.deriv[0] = v[0];
    output.deriv[1] = v[1];
    output.deriv[2] = v[2];
    output.deriv[3] = -v[0];
    output.deriv[4] = -v[1];
    output.deriv[5] = -v[2];

    output.assignOuterProduct(6,v,input.pos);

  }
}

}
}


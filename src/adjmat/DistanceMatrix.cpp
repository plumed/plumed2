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

Using a CUTOFF ensures that PLUMED can use the link cell technique that is described in the documentation for the [CONTACT_MATRIX](CONTACT_MATRIX.md) action to optimise the calculation.
Using this technique ensures that many of the distance calculations are avoided. However, this __does not__ mean that PLUMED will not evaluate and store the distances between pairs of
atoms that are more than the cutoff apart. The derivatives for such pairs are not evaluated but __the distances are stored__.

You can see how to work around this strange implementation detail in the following example input.  Lets suppose that we want to calculate the average distances
between atoms 1-10 and all the atoms that are within 1 nm of them.  To do this we would use an input similar to the one shown below:

```plumed
d5: DISTANCE_MATRIX GROUPA=1-10 GROUPB=1-250 CUTOFF=1.0
# Apply a switching function to determine the elements in the matrix d5
# where the distance is less than the cutoff
cut: CUSTOM ARG=d5 FUNC=step(1-x) PERIODIC=NO
# Taking the element-wise product in the next command gives us a matrix
# where every element that is greater than the cutoff is zero.
d5cut: CUSTOM ARG=d5,cut FUNC=x*y PERIODIC=NO
# We can now calculate the average distances by multiplying these matrices by
# a vector of ones and thus summing the rows.
ones: ONES SIZE=250
totdist: MATRIX_VECTOR_PRODUCT ARG=d5cut,ones
ndist: MATRIX_VECTOR_PRODUCT ARG=cut,ones
average: CUSTOM ARG=totdist,ndist FUNC=x/y PERIODIC=NO
DUMPATOMS ATOMS=1-10 ARG=average FILE=avdist.xyz
```

__In short, if you use DISTANCE_MATRIX and CUTOFF and what to ignore distances that are larger than the CUTOFF you need to use an additional [CUSTOM](CUSTOM.md) command later in the input.__
You should thus only use this combination of action and keyword it you are certain it is necessary. Normally, you are far better using the [CONTACT_MATRIX](CONTACT_MATRIX.md) command in place
of the DISTANCE_MATRIX command. To obtain a result similar to the one above using this command you would use the following input:

```plumed
cmap: CONTACT_MATRIX GROUPA=1-10 GROUPB=1-250 SWITCH={RATIONAL R_0=0.5 D_MAX=1.0 NN=6 MM=12} COMPONENTS
# Evaluate the distances for all pairs of atoms that are within D_MAX of each other
dmat: CUSTOM ARG=cmap.x,cmap.y,cmap.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
cdist: CUSTOM ARG=cmap.w,dmat FUNC=x*y PERIODIC=NO
ones: ONES SIZE=250
totdist: MATRIX_VECTOR_PRODUCT ARG=cdist,ones
ndist: MATRIX_VECTOR_PRODUCT ARG=cmap.w,ones
average: CUSTOM ARG=totdist,ndist FUNC=x/y PERIODIC=NO
DUMPATOMS ATOMS=1-10 ARG=average FILE=avdist.xyz
``

The advantage when using this input is the the final vector `average` that is evaluated here is a continous function. You can thus evaluate derivatives for it and use it as input in a biasing
method.

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
  keys.add("compulsory","CUTOFF","-1","use a link cells algorithm with this cutoff to optimise the calculation - distances (but not derivatives) for atoms that are further apart than this cutoff and that in the same link cells will still be computed");
}

void DistanceMatrix::parseInput( AdjacencyMatrixBase<DistanceMatrix>* action ) {
  // And set the link cell cutoff
  action->log.printf("  weight is distance between atoms \n");
  action->parse("CUTOFF",cutoff);
  if( cutoff<0 ) {
    action->setLinkCellCutoff( true, std::numeric_limits<double>::max() );
  } else {
    action->log.printf("  using link cells with cutoff %f to optimise calculation \n", cutoff);
    action->warning("some distances that are greater than the cutoff will still be evaluated. If you want to ignore them you will need to use a further CUSTOM action in your PLUMED input as detailed in the manual");
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


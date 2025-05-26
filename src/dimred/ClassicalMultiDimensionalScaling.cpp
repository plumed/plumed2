/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC DIMRED CLASSICAL_MDS
/*
Create a low-dimensional projection of a trajectory using the classical multidimensional
 scaling algorithm.

Multidimensional scaling (MDS) is similar to what is done when you make a map. You start with distances
between London, Belfast, Paris and Dublin and then you try to arrange points on a piece of paper so that the (suitably scaled)
distances between the points in your map representing each of those cities are related to the true distances between the cities.
Stating this more mathematically MDS endeavors to find an [isometry](http://en.wikipedia.org/wiki/Isometry)
between points distributed in a high-dimensional space and a set of points distributed in a low-dimensional plane.
In other words, if we have $M$ $D$-dimensional points, $\mathbf{X}$,
and we can calculate dissimilarities between pairs them, $D_{ij}$, we can, with an MDS calculation, try to create $M$ projections,
$\mathbf{x}$, of the high dimensionality points in a $d$-dimensional linear space by trying to arrange the projections so that the
Euclidean distances between pairs of them, $d_{ij}$, resemble the dissimilarities between the high dimensional points.  In short we minimize:

$$
\chi^2 = \sum_{i \ne j} \left( D_{ij} - d_{ij} \right)^2
$$

where $D_{ij}$ is the distance between point $X^{i}$ and point $X^{j}$ and $d_{ij}$ is the distance between the projection
of $X^{i}$, $x^i$, and the projection of $X^{j}$, $x^j$.  A tutorial on this approach can be used to analyze simulations
can be found in [this tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html) and in the following [short video](https://www.youtube.com/watch?v=ofC2qz0_9_A&feature=youtu.be).

## Examples

The following command instructs plumed to construct a classical multidimensional scaling projection of a trajectory.
The RMSD distance between atoms 1-256 have moved is used to measure the distances in the high-dimensional space.

```plumed
ff: COLLECT_FRAMES ATOMS=1-256
mds: CLASSICAL_MDS ARG=ff NLOW_DIM=2
weights: CUSTOM ARG=ff.logweights FUNC=exp(x) PERIODIC=NO
DUMPPDB ATOM_INDICES=1-256 ATOMS=ff_data ARG=mds,weights FILE=embed.pdb
```

By contrast the following input instructs PLUMED to calculate the distances between atoms by taking the differences in five torsional angles.
The MDS algorithm is then used to arrange a set of points in a low dimensional space in a way that reproduces these differences.

```plumed
phi1: TORSION ATOMS=1,2,3,4
phi2: TORSION ATOMS=5,6,7,8
phi3: TORSION ATOMS=9,10,11,12
phi4: TORSION ATOMS=13,14,15,16
phi5: TORSION ATOMS=17,18,19,20

angles: COLLECT_FRAMES ARG=phi1,phi2,phi3,phi4,phi5
mds: CLASSICAL_MDS ARG=angles NLOW_DIM=2
weights: CUSTOM ARG=angles.logweights FUNC=exp(x) PERIODIC=NO
DUMPVECTOR ARG=mds,weights FILE=list_embed
```

The following section is for people who are interested in how this method works in detail. A solid understanding of this material is
not necessary to use MDS.

## Method of optimization

The stress function can be minimized using the conjugate gradients or steepest descent optimization algorithms that are implemented in [ARRANGE_POINTS](ARRANGE_POINTS.md).
However, it is more common to do this minimization using a technique known as classical scaling.  Classical scaling works by
recognizing that each of the distances $D_{ij}$ in the above sum can be written as:

$$
D_{ij}^2 = \sum_{\alpha} (X^i_\alpha - X^j_\alpha)^2 = \sum_\alpha (X^i_\alpha)^2 + (X^j_\alpha)^2 - 2X^i_\alpha X^j_\alpha
$$

We can use this expression and matrix algebra to calculate multiple distances at once.  For instance if we have three points,
$\mathbf{X}$, we can write distances between them as:

$$
\begin{aligned}
D^2(\mathbf{X}) &= \left[ \begin{array}{ccc}
0 & d_{12}^2 & d_{13}^2 \\
d_{12}^2 & 0 & d_{23}^2 \\
d_{13}^2 & d_{23}^2 & 0
\end{array}\right] \\
&=
\sum_\alpha \left[ \begin{array}{ccc}
(X^1_\alpha)^2 & (X^1_\alpha)^2 & (X^1_\alpha)^2 \\
(X^2_\alpha)^2 & (X^2_\alpha)^2 & (X^2_\alpha)^2 \\
(X^3_\alpha)^2 & (X^3_\alpha)^2 & (X^3_\alpha)^2 \\
\end{array}\right]
 + \sum_\alpha \left[ \begin{array}{ccc}
(X^1_\alpha)^2 & (X^2_\alpha)^2 & (X^3_\alpha)^2 \\
(X^1_\alpha)^2 & (X^2_\alpha)^2 & (X^3_\alpha)^2 \\
(X^1_\alpha)^2 & (X^2_\alpha)^2 & (X^3_\alpha)^2 \\
\end{array}\right]
- 2 \sum_\alpha \left[ \begin{array}{ccc}
X^1_\alpha X^1_\alpha & X^1_\alpha X^2_\alpha & X^1_\alpha X^3_\alpha \\
X^2_\alpha X^1_\alpha & X^2_\alpha X^2_\alpha & X^2_\alpha X^3_\alpha \\
X^1_\alpha X^3_\alpha & X^3_\alpha X^2_\alpha & X^3_\alpha X^3_\alpha
\end{array}\right] \nonumber \\
&= \mathbf{c 1^T} + \mathbf{1 c^T} - 2 \sum_\alpha \mathbf{x}_a \mathbf{x}^T_a =  \mathbf{c 1^T} + \mathbf{1 c^T} - 2\mathbf{X X^T}
\end{aligned}
$$

This last equation can be extended to situations when we have more than three points.  In it $\mathbf{X}$ is a matrix that has
one high-dimensional point on each of its rows and $\mathbf{X^T}$ is its transpose.  $\mathbf{1}$ is an $M \times 1$ vector
of ones and $\mathbf{c}$ is a vector with components given by:

$$
c_i = \sum_\alpha (x_\alpha^i)^2
$$

These quantities are the diagonal elements of $\mathbf{X X^T}$, which is a dot product or Gram Matrix that contains the
dot product of the vector $X_i$ with the vector $X_j$ in element $i,j$.

In classical scaling we introduce a centering matrix $\mathbf{J}$ that is given by:

$$
\mathbf{J} = \mathbf{I} - \frac{1}{M} \mathbf{11^T}
$$

where $\mathbf{I}$ is the identity.  Multiplying the equations above from the front and back by this matrix and a factor of a $-\frac{1}{2}$ gives:

$$
\begin{aligned}
 -\frac{1}{2} \mathbf{J} \mathbf{D}^2(\mathbf{X}) \mathbf{J} &= -\frac{1}{2}\mathbf{J}( \mathbf{c 1^T} + \mathbf{1 c^T} - 2\mathbf{X X^T})\mathbf{J} \\
 &= -\frac{1}{2}\mathbf{J c 1^T J} - \frac{1}{2} \mathbf{J 1 c^T J} + \frac{1}{2} \mathbf{J}(2\mathbf{X X^T})\mathbf{J} \\
 &= \mathbf{ J X X^T J } = \mathbf{X X^T } \label{eqn:scaling}
\end{aligned}
$$

The fist two terms in this expression disappear because $\mathbf{1^T J}=\mathbf{J 1} =\mathbf{0}$, where $\mathbf{0}$
is a matrix containing all zeros.  In the final step meanwhile we use the fact that the matrix of squared distances will not
change when we translate all the points.  We can thus assume that the mean value, $\mu$, for each of the components, $\alpha$:

$$
\mu_\alpha = \frac{1}{M} \sum_{i=1}^N \mathbf{X}^i_\alpha
$$

is equal to 0 so the columns of $\mathbf{X}$ add up to 0.  This in turn means that each of the columns of
$\mathbf{X X^T}$ adds up to zero, which is what allows us to write $\mathbf{ J X X^T J } = \mathbf{X X^T }$.

The matrix of squared distances is symmetric and positive-definite we can thus use the spectral decomposition to decompose it as:

$$
\Phi= \mathbf{V} \Lambda \mathbf{V}^T
$$

Furthermore, because the matrix we are diagonalizing, $\mathbf{X X^T}$, is the product of a matrix and its transpose
we can use this decomposition to write:

$$
\mathbf{X} =\mathbf{V} \Lambda^\frac{1}{2}
$$

Much as in PCA there are generally a small number of large eigenvalues in $\Lambda$ and many small eigenvalues.
We can safely use only the large eigenvalues and their corresponding eigenvectors to express the relationship between
the coordinates $\mathbf{X}$.  This gives us our set of low-dimensional projections.

This derivation makes a number of assumptions about the how the low dimensional points should best be arranged to minimize
the stress. If you use an interactive optimization algorithm such as SMACOF you may thus be able to find a better
(lower-stress) projection of the points.  For more details on the assumptions made
see [this website](http://quest4rigor.com/tag/multidimensional-scaling/).
*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class ClassicalMultiDimensionalScaling : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit ClassicalMultiDimensionalScaling( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(ClassicalMultiDimensionalScaling,"CLASSICAL_MDS")

void ClassicalMultiDimensionalScaling::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that you would like to make the histogram for");
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.setValueDescription("matrix","the low dimensional projections for the input data points");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("DISSIMILARITIES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("VSTACK");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("DIAGONALIZE");
}

ClassicalMultiDimensionalScaling::ClassicalMultiDimensionalScaling( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Find the argument name
  std::string argn;
  parse("ARG",argn);
  std::string dissimilarities="";
  ActionShortcut* as = plumed.getActionSet().getShortcutActionWithLabel( argn );
  if( !as ) {
    error("found no action with name " + argn );
  }
  if( as->getName()!="COLLECT_FRAMES" ) {
    if( as->getName().find("LANDMARK_SELECT")==std::string::npos ) {
      error("found no COLLECT_FRAMES or LANDMARK_SELECT action with label " + argn );
    } else {
      ActionWithValue* dissims = plumed.getActionSet().selectWithLabel<ActionWithValue*>( argn + "_sqrdissims");
      if( dissims ) {
        dissimilarities = argn + "_sqrdissims";
      }
    }
  }
  if( dissimilarities.length()==0 ) {
    dissimilarities = getShortcutLabel() + "_mat";
    // Transpose matrix of stored data values
    readInputLine( argn + "_dataT: TRANSPOSE ARG=" + argn + "_data");
    // Calculate the dissimilarity matrix
    readInputLine( getShortcutLabel() + "_mat: DISSIMILARITIES SQUARED ARG=" + argn + "_data," + argn + "_dataT");
  }
  // Center the matrix
  // Step 1: calculate the sum of the rows and duplicate them into a matrix
  readInputLine( getShortcutLabel() + "_rsums: MATRIX_VECTOR_PRODUCT ARG=" + dissimilarities + "," + argn + "_ones" );
  readInputLine( getShortcutLabel() + "_nones: SUM ARG=" + argn + "_ones PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_rsumsn: CUSTOM ARG=" + getShortcutLabel() + "_rsums," + getShortcutLabel() + "_nones FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_rsummat: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_rsumsn," + argn + "_ones");
  // Step 2: Multiply matrix by -0.5 and subtract row sums
  readInputLine( getShortcutLabel() + "_int: CUSTOM ARG=" + getShortcutLabel() + "_rsummat," + dissimilarities + " FUNC=-0.5*y+0.5*x PERIODIC=NO");
  // Step 3: Calculate column sums for new matrix and duplicate them into a matrix
  readInputLine( getShortcutLabel() + "_intT: TRANSPOSE ARG=" + getShortcutLabel() + "_int");
  readInputLine( getShortcutLabel() + "_csums: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_intT," + argn + "_ones" );
  readInputLine( getShortcutLabel() + "_csumsn: CUSTOM ARG=" + getShortcutLabel() + "_csums," + getShortcutLabel() + "_nones FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_csummat: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_csumsn," + argn + "_ones");
  // Step 4: subtract the column sums
  readInputLine( getShortcutLabel() + "_cmat: CUSTOM ARG=" + getShortcutLabel() + "_csummat," + getShortcutLabel() + "_intT FUNC=y-x PERIODIC=NO");
  // And generate the multidimensional scaling projection
  unsigned ndim;
  parse("NLOW_DIM",ndim);
  std::string vecstr="1";
  for(unsigned i=1; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    vecstr += "," + num;
  }
  readInputLine( getShortcutLabel() + "_eig: DIAGONALIZE ARG=" + getShortcutLabel() + "_cmat VECTORS=" + vecstr );
  for(unsigned i=0; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    readInputLine( getShortcutLabel() + "-" +  num + ": CUSTOM ARG=" + getShortcutLabel() + "_eig.vecs-" + num + "," + getShortcutLabel() + "_eig.vals-" + num + " FUNC=sqrt(y)*x PERIODIC=NO");
  }
  std::string eigvec_args = " ARG=" + getShortcutLabel() + "-1";
  // The final output is a stack of all the low dimensional coordinates
  for(unsigned i=1; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    eigvec_args += "," + getShortcutLabel() + "-" + num;
  }
  readInputLine( getShortcutLabel() + ": VSTACK" + eigvec_args );
}

}
}

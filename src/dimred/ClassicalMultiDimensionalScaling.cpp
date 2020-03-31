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
#include "ClassicalMultiDimensionalScaling.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/AverageBase.h"

//+PLUMEDOC DIMRED CLASSICAL_MDS
/*
Create a low-dimensional projection of a trajectory using the classical multidimensional
 scaling algorithm.

Multidimensional scaling (MDS) is similar to what is done when you make a map. You start with distances
between London, Belfast, Paris and Dublin and then you try to arrange points on a piece of paper so that the (suitably scaled)
distances between the points in your map representing each of those cities are related to the true distances between the cities.
Stating this more mathematically MDS endeavors to find an <a href="http://en.wikipedia.org/wiki/Isometry">isometry</a>
between points distributed in a high-dimensional space and a set of points distributed in a low-dimensional plane.
In other words, if we have \f$M\f$ \f$D\f$-dimensional points, \f$\mathbf{X}\f$,
and we can calculate dissimilarities between pairs them, \f$D_{ij}\f$, we can, with an MDS calculation, try to create \f$M\f$ projections,
\f$\mathbf{x}\f$, of the high dimensionality points in a \f$d\f$-dimensional linear space by trying to arrange the projections so that the
Euclidean distances between pairs of them, \f$d_{ij}\f$, resemble the dissimilarities between the high dimensional points.  In short we minimize:

\f[
\chi^2 = \sum_{i \ne j} \left( D_{ij} - d_{ij} \right)^2
\f]

where \f$D_{ij}\f$ is the distance between point \f$X^{i}\f$ and point \f$X^{j}\f$ and \f$d_{ij}\f$ is the distance between the projection
of \f$X^{i}\f$, \f$x^i\f$, and the projection of \f$X^{j}\f$, \f$x^j\f$.  A tutorial on this approach can be used to analyze simulations
can be found in the tutorial \ref lugano-5 and in the following <a href="https://www.youtube.com/watch?v=ofC2qz0_9_A&feature=youtu.be" > short video.</a>

\par Examples

The following command instructs plumed to construct a classical multidimensional scaling projection of a trajectory.
The RMSD distance between atoms 1-256 have moved is used to measure the distances in the high-dimensional space.

\plumedfile
data: COLLECT_FRAMES ATOMS=1-256
mat: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data
mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=mat NLOW_DIM=2
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=mds FILE=rmsd-embed
\endplumedfile

The following section is for people who are interested in how this method works in detail. A solid understanding of this material is
not necessary to use MDS.

\section dim-sec Method of optimization

The stress function can be minimized using a standard optimization algorithm such as conjugate gradients or steepest descent.
However, it is more common to do this minimization using a technique known as classical scaling.  Classical scaling works by
recognizing that each of the distances $D_{ij}$ in the above sum can be written as:

\f[
D_{ij}^2 = \sum_{\alpha} (X^i_\alpha - X^j_\alpha)^2 = \sum_\alpha (X^i_\alpha)^2 + (X^j_\alpha)^2 - 2X^i_\alpha X^j_\alpha
\f]

We can use this expression and matrix algebra to calculate multiple distances at once.  For instance if we have three points,
\f$\mathbf{X}\f$, we can write distances between them as:

\f{eqnarray*}{
D^2(\mathbf{X}) &=& \left[ \begin{array}{ccc}
0 & d_{12}^2 & d_{13}^2 \\
d_{12}^2 & 0 & d_{23}^2 \\
d_{13}^2 & d_{23}^2 & 0
\end{array}\right] \\
&=&
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
&=& \mathbf{c 1^T} + \mathbf{1 c^T} - 2 \sum_\alpha \mathbf{x}_a \mathbf{x}^T_a =  \mathbf{c 1^T} + \mathbf{1 c^T} - 2\mathbf{X X^T}
\f}

This last equation can be extended to situations when we have more than three points.  In it \f$\mathbf{X}\f$ is a matrix that has
one high-dimensional point on each of its rows and \f$\mathbf{X^T}\f$ is its transpose.  \f$\mathbf{1}\f$ is an \f$M \times 1\f$ vector
of ones and \f$\mathbf{c}\f$ is a vector with components given by:

\f[
c_i = \sum_\alpha (x_\alpha^i)^2
\f]

These quantities are the diagonal elements of \f$\mathbf{X X^T}\f$, which is a dot product or Gram Matrix that contains the
dot product of the vector \f$X_i\f$ with the vector \f$X_j\f$ in element \f$i,j\f$.

In classical scaling we introduce a centering matrix \f$\mathbf{J}\f$ that is given by:

\f[
\mathbf{J} = \mathbf{I} - \frac{1}{M} \mathbf{11^T}
\f]

where \f$\mathbf{I}\f$ is the identity.  Multiplying the equations above from the front and back by this matrix and a factor of a \f$-\frac{1}{2}\f$ gives:

\f{eqnarray*}{
 -\frac{1}{2} \mathbf{J} \mathbf{D}^2(\mathbf{X}) \mathbf{J} &=& -\frac{1}{2}\mathbf{J}( \mathbf{c 1^T} + \mathbf{1 c^T} - 2\mathbf{X X^T})\mathbf{J} \\
 &=& -\frac{1}{2}\mathbf{J c 1^T J} - \frac{1}{2} \mathbf{J 1 c^T J} + \frac{1}{2} \mathbf{J}(2\mathbf{X X^T})\mathbf{J} \\
 &=& \mathbf{ J X X^T J } = \mathbf{X X^T } \label{eqn:scaling}
\f}

The fist two terms in this expression disappear because \f$\mathbf{1^T J}=\mathbf{J 1} =\mathbf{0}\f$, where \f$\mathbf{0}\f$
is a matrix containing all zeros.  In the final step meanwhile we use the fact that the matrix of squared distances will not
change when we translate all the points.  We can thus assume that the mean value, \f$\mu\f$, for each of the components, \f$\alpha\f$:
\f[
\mu_\alpha = \frac{1}{M} \sum_{i=1}^N \mathbf{X}^i_\alpha
\f]
is equal to 0 so the columns of \f$\mathbf{X}\f$ add up to 0.  This in turn means that each of the columns of
\f$\mathbf{X X^T}\f$ adds up to zero, which is what allows us to write \f$\mathbf{ J X X^T J } = \mathbf{X X^T }\f$.

The matrix of squared distances is symmetric and positive-definite we can thus use the spectral decomposition to decompose it as:

\f[
\Phi= \mathbf{V} \Lambda \mathbf{V}^T
\f]

Furthermore, because the matrix we are diagonalizing, \f$\mathbf{X X^T}\f$, is the product of a matrix and its transpose
we can use this decomposition to write:

\f[
\mathbf{X} =\mathbf{V} \Lambda^\frac{1}{2}
\f]

Much as in PCA there are generally a small number of large eigenvalues in \f$\Lambda\f$ and many small eigenvalues.
We can safely use only the large eigenvalues and their corresponding eigenvectors to express the relationship between
the coordinates \f$\mathbf{X}\f$.  This gives us our set of low-dimensional projections.

This derivation makes a number of assumptions about the how the low dimensional points should best be arranged to minimize
the stress. If you use an interactive optimization algorithm such as SMACOF you may thus be able to find a better
(lower-stress) projection of the points.  For more details on the assumptions made
see <a href="http://quest4rigor.com/tag/multidimensional-scaling/"> this website.</a>
*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

PLUMED_REGISTER_ACTION(ClassicalMultiDimensionalScaling,"CLASSICAL_MDS")

void ClassicalMultiDimensionalScaling::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that you would like to make the histogram for");
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
}

ClassicalMultiDimensionalScaling::ClassicalMultiDimensionalScaling( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string arg; parse("ARG",arg); unsigned anum=1, natoms3=0; std::string argstr; 
  AverageBase* mydata = plumed.getActionSet().selectWithLabel<AverageBase*>(arg); bool isatoms=true;
  if( !mydata ) error("input to PCA should be a COLLECT_FRAMES or COLLECT_REPLICAS object");
  for(unsigned i=0;i<mydata->getNumberOfComponents();++i) {
      std::string thislab = mydata->copyOutput(i)->getName();
      if( thislab.find(".logweights")==std::string::npos ) {
          std::string num; Tools::convert( anum, num );
          argstr += " GROUP" + num + "=" + thislab; anum++;
          if( thislab.find(".pos")!=std::string::npos ) natoms3++;
          else isatoms=false;
      }
  }
  // Calculate the dissimilarity matrix
  if( isatoms ) {
      std::string nat_str;
      readInputLine( getShortcutLabel() + "_matful: DISSIMILARITIES SQUARED " + argstr ); Tools::convert( natoms3/3, nat_str );
      readInputLine( getShortcutLabel() + "_mat: MATHEVAL PERIODIC=NO ARG1=" + getShortcutLabel() + "_matful FUNC=x/" + nat_str ); 
  } else readInputLine( getShortcutLabel() + "_mat: DISSIMILARITIES SQUARED " + argstr ); 
  // And generate the multidimensional scaling projection 
  unsigned ndim; parse("NLOW_DIM",ndim); createMDSProjection( getShortcutLabel(), getShortcutLabel() + "_mat", ndim, this ); 
}

void ClassicalMultiDimensionalScaling::createMDSProjection( const std::string& lab, const std::string& matname, const unsigned& ndim, ActionShortcut* myact ) { 
  // Center the dissimilarity matrix
  myact->readInputLine( lab + "_cmat: CENTER_MATRIX ARG=" + matname );
  // Diagonalize the centered dissimilarity matrix
  std::string vecstr="1"; for(unsigned i=1;i<ndim;++i){ std::string num; Tools::convert( i+1, num ); vecstr += "," + num; }
  myact->readInputLine( lab + "_eig: DIAGONALIZE ARG=" + lab + "_cmat VECTORS=" + vecstr );
  // And calculate the mds projections
  for(unsigned i=0;i<ndim;++i) {
      std::string num; Tools::convert( i+1, num );
      myact->readInputLine( lab + "-" +  num + ": MATHEVAL ARG1=" + lab + "_eig.vals-" + num + " ARG2=" + lab + "_eig.vecs-" + num + " FUNC=sqrt(x)*y PERIODIC=NO");
  }
}

}
}

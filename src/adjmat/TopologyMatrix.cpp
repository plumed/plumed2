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
#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"
#include "tools/HistogramBead.h"
#include "tools/Matrix.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

//+PLUMEDOC MATRIX TOPOLOGY_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if they are connected topologically

The functionality in this action was developed as part of a project that attempted to study
the nucleation of bubbles.  The idea was to develop a criterion that would allow one to determine
if two gas atoms $i$ and $j$ are part of the same bubble or not.  This criterion would then be used
to construct a adjacency matrix, which could be used in the same way that [CONTACT_MATRIX](CONTACT_MATRIX.md) is used in other
methods.

The criterion that was developed to determine whether atom $i$ and $j$ are connected in this way works by
determining if the density within a cylinder that is centered on the vector connecting atoms $i$ and $j$ is
less than a certain threshold value.  To make this determination we first determine whether any given atom $k$
is within the cylinder centered on the vector connecting atoms $i$ and $j$ by using the following expression

$$
f(\mathbf{r}_{ik}, \mathbf{r}_{ij}) = s_1\left( \sqrt{ |\mathbf{r}_{ij}|^2 - \left( \frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{ik}}{|\mathbf{r}_{ij} |} \right)^2} \right)s_2\left( -\frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{ik}}{|\mathbf{r}_{ij} |} \right) s_2\left( \frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{ik}}{|\mathbf{r}_{ij} |} - |\mathbf{r}_{ij}| \right)
$$

In this expression $s_1$ and $s_2$ are switching functions, while $\mathbf{r}_{lm}$ is used to indicate the vector connecting atoms $l$ and $m$.

We then calculate the density for a grid of $M$ points along the vector connecting atom $i$ and atom $j$ using and find the maximum density on this grid using:

$$
\rho_{ij} = \max_{m \in M} \left[ \frac{M}{d_\textrm{max}} \right] \sum_k f(r_{ik}, r_{ij}) \int_{(m-1)d_{\textrm{max}}/M}^{ md_{\textrm{max}} /M } \textrm{d}x \quad K\left( \frac{x - r_{ks} \cdot r_{ij} }{ | r_{ks} | }\right)
$$

where $d_\textrm{max}$ is the `D_MAX` parameter of the switching function $s_3$ that appears in the next equation, $K$ is a kernal function and $s$ is used to represent a point in space that is $d_\textrm{max}$ from atom $j$ along the vector connecting atom $j$ to atom $i$.

The final value that is stored in the $i, j$ element of the output matrix is:

$$
T_{ij} = s_3(|\mathbf{r}_{ij}|)s_4(\rho_{ij})
$$

where $s_3$ and $s_4$ are switching functions.

We ended up abandoning this method and the project (if you want drive bubble formation you are probably better off adding a bias on the volume of the simulation cell).
However, if you would like to try this method an example input that uses this action would look like this:

```plumed
mat: TOPOLOGY_MATRIX ...
    GROUP=1-85 BACKGROUND_ATOMS=86-210
    BIN_SIZE=1.02 SIGMA=0.17 KERNEL=triangular
    CYLINDER_SWITCH={RATIONAL R_0=0.5 D_MAX=1.0}
    SWITCH={RATIONAL D_0=30 R_0=0.5 D_MAX=32}
    RADIUS={RATIONAL D_0=0.375 R_0=0.1 D_MAX=0.43}
    DENSITY_THRESHOLD={RATIONAL R_0=0.1 D_MAX=0.5}
...
```

The switching functions $s_1$, $s_2$, $s_3$ and $s_4$ are specified using the `RADIUS`, `CYLINDER_SWITCH`, `SWITCH` and `DENSITY_THRESHOLD` keywords respectively.
We loop over the atoms in the group specified using the `BACKGROUND_ATOMS` keyword when looping over $k$ in the formulas above.  An $85 \times 85$ matrix is output
from the method as we are determining the connectivity between the atoms specified via the `GROUP` keyword.

The above example assumes that you want to calculate the connectivity within a single group of atoms.  If you would calculate connectivity between two different groups
of atoms you use the GROUPA and GROUPB keywords as shown below:

```plumed
mat: TOPOLOGY_MATRIX ...
    GROUPA=1-20 GROUPB=21-85 BACKGROUND_ATOMS=86-210
    BIN_SIZE=1.02 SIGMA=0.17 KERNEL=triangular
    CYLINDER_SWITCH={RATIONAL R_0=0.5 D_MAX=1.0}
    SWITCH={RATIONAL D_0=30 R_0=0.5 D_MAX=32}
    RADIUS={RATIONAL D_0=0.375 R_0=0.1 D_MAX=0.43}
    DENSITY_THRESHOLD={RATIONAL R_0=0.1 D_MAX=0.5}
    COMPONENTS NOPBC
...
```

Notice that we have also added the NOPBC and COMPONENTS keywords in this input. The action above thus outputs four matrices with the labels
`mat.w`, `mat.x`, `mat.y` and `mat.z.`  The matrix with the label `mat.w` is the adjacency matrix
that would be output if you had not added the COMPONENTS flag. The $i,j$ component of the matrices `mat.x`, `mat.y` and `mat.z` contain the $x$, $y$ and $z$
components of the vector connecting atoms $i$ and $k$. Importantly, however, the components of these vectors are only stored in `mat.x`, `mat.y` and `mat.z`
if the elements of `mat.w` are non-zero. Using the COMPONENTS flag in this way ensures that you can use HBOND_MATRIX in tandem with many of the functionalities
that are part of the [symfunc module](module_symfunc.md).

The NOPBC flag, meanwhile, ensures that all distances are calculated in a way that __does not__ take the periodic boundary conditions into account. By default,
distances are calculated in a way that takes periodic boundary conditions into account.

## The MASK keyword

You use the MASK keyword with TOPOLOGY_MATRIX in the same way that is used in [CONTACT_MATRIX](CONTACT_MATRIX.md).  This keyword thus expects a vector in input,
which tells TOPOLOGY_MATRIX that it is safe to not calculate certain rows of the output matrix.  An example where this keyword is used is shown below:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-85 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculates cooordination numbers
cmap: TOPOLOGY_MATRIX ...
  GROUP=1-85 BACKGROUND_ATOMS=86-210
  BIN_SIZE=1.02 SIGMA=0.17 KERNEL=triangular
  CYLINDER_SWITCH={RATIONAL R_0=0.5 D_MAX=1.0}
  SWITCH={RATIONAL D_0=30 R_0=0.5 D_MAX=32}
  RADIUS={RATIONAL D_0=0.375 R_0=0.1 D_MAX=0.43}
  DENSITY_THRESHOLD={RATIONAL R_0=0.1 D_MAX=0.5}
  MASK=sphere
...
ones: ONES SIZE=85
cc: MATRIX_VECTOR_PRODUCT ARG=cmap,ones
# Multiply coordination numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=prod,sphere FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average number of topological connections there are for each of the atoms that are within a spherical region
that is centered on the point $(2.5,2.5,2.5)$.


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class TopologyMatrix {
public:
  double sigma;
  HistogramBead::KernelType kerneltype;
/// The maximum number of bins that will be used
/// This is calculated based on the dmax of the switching functions
  unsigned maxbins;
/// The volume of the cells
  double cell_volume;
  double binw_mat;
/// switching functions
  SwitchingFunction switchingFunction;
  SwitchingFunction cylinder_sw;
  SwitchingFunction low_sf;
  SwitchingFunction threshold_switch;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<TopologyMatrix>* action );
  static void calculateWeight( const TopologyMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
};

typedef AdjacencyMatrixBase<TopologyMatrix> tmap;
PLUMED_REGISTER_ACTION(tmap,"TOPOLOGY_MATRIX")

void TopologyMatrix::registerKeywords( Keywords& keys ) {
  keys.add("atoms","BACKGROUND_ATOMS","the list of atoms that should be considered as part of the background density");
  keys.add("compulsory","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("compulsory","RADIUS","swtiching function that acts upon the radial distance of the atom from the center of the cylinder");
  keys.linkActionInDocs("RADIUS","LESS_THAN");
  keys.add("compulsory","CYLINDER_SWITCH","a switching function on ( r_ij . r_ik - 1 )/r_ij");
  keys.linkActionInDocs("CYLINDER_SWITCH","LESS_THAN");
  keys.add("compulsory","BIN_SIZE","the size to use for the bins");
  keys.add("compulsory","DENSITY_THRESHOLD","a switching function that acts upon the maximum density in the cylinder");
  keys.linkActionInDocs("DENSITY_THRESHOLD","LESS_THAN");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
}

void TopologyMatrix::parseInput( AdjacencyMatrixBase<TopologyMatrix>* action ) {
  std::string errors;

  std::string sfinput;
  action->parse("SWITCH",sfinput);
  if( sfinput.length()==0 ) {
    action->error("could not find SWITCH keyword");
  }
  switchingFunction.set(sfinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading SWITCH keyword : " + errors );
  }

  std::string lowsfinput;
  action->parse("CYLINDER_SWITCH",lowsfinput);
  if( lowsfinput.length()==0 ) {
    action->error("could not find CYLINDER_SWITCH keyword");
  }
  low_sf.set(lowsfinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading CYLINDER_SWITCH keyword : " + errors );
  }

  std::string cyinput;
  action->parse("RADIUS",cyinput);
  if( cyinput.length()==0 ) {
    action->error("could not find RADIUS keyword");
  }
  cylinder_sw.set(cyinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading RADIUS keyword : " + errors );
  }

  std::string thresholdinput;
  action->parse("DENSITY_THRESHOLD",thresholdinput);
  if( thresholdinput.length()==0 ) {
    action->error("could not find DENSITY_THRESHOLD keyword");
  }
  threshold_switch.set(thresholdinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading DENSITY_THRESHOLD keyword : " + errors );
  }
  // Read in stuff for grid
  action->parse("SIGMA",sigma);
  std::string mykerneltype;
  action->parse("KERNEL",mykerneltype);
  kerneltype = HistogramBead::getKernelType(mykerneltype);
  action->parse("BIN_SIZE",binw_mat);

  // Set the link cell cutoff
  action->setLinkCellCutoff( true, switchingFunction.get_dmax(),
                             std::numeric_limits<double>::max() );
  // Set the number of bins
  maxbins = std::floor( switchingFunction.get_dmax() / binw_mat ) + 1;
  // Set the cell volume
  double r=cylinder_sw.get_d0() + cylinder_sw.get_r0();
  cell_volume=binw_mat*PLMD::pi*r*r;
}

void TopologyMatrix::calculateWeight( const TopologyMatrix& data,
                                      const AdjacencyMatrixInput& input,
                                      MatrixOutput output ) {
  // Compute switching function on distance between atoms
  Vector distance = input.pos;
  double len2 = distance.modulo2();
  if( len2>data.switchingFunction.get_dmax2() ) {
    return;
  }
  double dfuncl;
  double  sw = data.switchingFunction.calculateSqr( len2, dfuncl );

  // Now run through all sea atoms
  HistogramBead bead( data.kerneltype, 0.0, data.binw_mat, data.sigma  );
  Vector g1derivf,g2derivf,lderivf;
  Tensor vir;
  double binlength = data.maxbins * data.binw_mat;
  std::vector<double> tvals( data.maxbins, 0 );
  Matrix<double> tvals_derivs( data.maxbins, 6 + 3*input.natoms + 9 );
  tvals_derivs = 0;
  // tvals.resize( data.maxbins, 6 + 3*input.natoms + 9, 0 );
  for(unsigned i=0; i<input.natoms; ++i) {
    // Position of sea atom (this will be the origin)
    Vector d2(input.extra_positions[i][0],
              input.extra_positions[i][1],
              input.extra_positions[i][2]);
    // Vector connecting sea atom and first in bond taking pbc into account
    Vector d20 = input.pbc->distance( d2, Vector(0,0,0) );
    // Vector connecting sea atom and second in bond taking pbc into account
    Vector d21 = input.pbc->distance( d2, input.pos );
    // Now length of bond modulus and so on -- no pbc here as we want sea atom in middle
    Vector d1 = delta( d20, d21 );
    double d1_len = d1.modulo();
    d1 = d1 / d1_len;
    // Switching function on distance between nodes
    if( d1_len>data.switchingFunction.get_dmax() ) {
      continue ;
    }
    // Ensure that the center of the bins are on the center of the bond connecting the two atoms
    double start2atom = 0.5*(binlength-d1_len);
    Vector dstart = d20 - start2atom*d1;
    // Now calculate projection of axis of cylinder
    double proj=dotProduct(-dstart,d1);
    // Calculate length of vector connecting start of cylinder to first atom
    // And now work out projection on vector connecting start and end of cylinder
    double proj_between = proj - start2atom;
    // This tells us if we are outside the end of the cylinder
    double excess = proj_between - d1_len;
    // Return if we are outside of the cylinder as calculated based on excess
    if( excess>data.low_sf.get_dmax() || -proj_between>data.low_sf.get_dmax() ) {
      continue;
    }
    // Calculate the excess swiching functions
    double edf1;
    double eval1 = data.low_sf.calculate( excess, edf1 );
    double edf2;
    double eval2 = data.low_sf.calculate( -proj_between, edf2 );
    // Calculate the projection on the perpendicular distance from the center of the tube
    double cm = dstart.modulo2() - proj*proj;

    // Now calculate the density in the cylinder
    if( cm>0 && cm<data.cylinder_sw.get_dmax2() ) {
      double dfuncr;
      double val = data.cylinder_sw.calculateSqr( cm, dfuncr );
      Vector dc1, dc2, dc3, dd1, dd2, dd3, de1, de2, de3;
      if( !input.noderiv ) {
        Tensor d1_a1;
        // Derivative of director connecting atom1 - atom2 wrt the position of atom 1
        d1_a1(0,0) = ( -(d1[1]*d1[1]+d1[2]*d1[2])/d1_len );   // dx/dx
        d1_a1(0,1) = (  d1[0]*d1[1]/d1_len );                 // dx/dy
        d1_a1(0,2) = (  d1[0]*d1[2]/d1_len );                 // dx/dz
        d1_a1(1,0) = (  d1[1]*d1[0]/d1_len );                 // dy/dx
        d1_a1(1,1) = ( -(d1[0]*d1[0]+d1[2]*d1[2])/d1_len );   // dy/dy
        d1_a1(1,2) = (  d1[1]*d1[2]/d1_len );
        d1_a1(2,0) = (  d1[2]*d1[0]/d1_len );
        d1_a1(2,1) = (  d1[2]*d1[1]/d1_len );
        d1_a1(2,2) = ( -(d1[1]*d1[1]+d1[0]*d1[0])/d1_len );

        // Calculate derivatives of dot product
        dd1 = matmul(-dstart, d1_a1) - 0.5*d1;
        dd2 = matmul(-dstart, -d1_a1) - 0.5*d1;
        dd3 = d1;

        // Calculate derivatives of cross product
        Vector der( -0.5*binlength*matmul( d1_a1,dstart ) );
        dc1 = dfuncr*( 0.5*dstart + der - proj*dd1 );
        dc2 = dfuncr*( 0.5*dstart - der - proj*dd2 );
        dc3 = dfuncr*( -dstart - proj*dd3 );

        // Calculate derivatives of excess
        de1 = eval2*edf1*excess*(dd1 + 0.5*d1 )
              + eval1*edf2*proj_between*(dd1 - 0.5*d1);
        de2 = eval2*edf1*excess*(dd2 - 0.5*d1 )
              + eval1*edf2*proj_between*(dd2 + 0.5*d1);
        de3 = ( eval2*edf1*excess + eval1*edf2*proj_between )*dd3;
      }
      for(unsigned bin=0; bin<data.maxbins; ++bin) {
        bead.set( bin*data.binw_mat, (bin+1)*data.binw_mat, data.sigma );
        if( proj<(bin*data.binw_mat-bead.getCutoff())
            || proj>data.binw_mat*(bin+1)+bead.getCutoff() ) {
          continue;
        }
        double der;
        double contr=bead.calculateWithCutoff( proj, der ) / data.cell_volume;
        der /= data.cell_volume;
        tvals[bin] += contr*val*eval1*eval2;

        if( !input.noderiv ) {
          g1derivf=contr*eval1*eval2*dc1 + val*eval1*eval2*der*dd1 + contr*val*de1;
          tvals_derivs[bin][0] += g1derivf[0];
          tvals_derivs[bin][1] += g1derivf[1];
          tvals_derivs[bin][2] += g1derivf[2];
          g2derivf=contr*eval1*eval2*dc2 + val*eval1*eval2*der*dd2 + contr*val*de2;
          tvals_derivs[bin][3] += g2derivf[0];
          tvals_derivs[bin][4] += g2derivf[1];
          tvals_derivs[bin][5] += g2derivf[2];
          lderivf=contr*eval1*eval2*dc3 + val*eval1*eval2*der*dd3 + contr*val*de3;
          tvals_derivs[bin][6 + 3*i+0] += lderivf[0];
          tvals_derivs[bin][6 + 3*i+1] += lderivf[1];
          tvals_derivs[bin][6 + 3*i+2] += lderivf[2];
          // Virial
          vir = - Tensor( d20, g1derivf ) - Tensor( d21, g2derivf );
          unsigned nbase = 6 + 3*input.natoms;
          tvals_derivs[bin][nbase+0] += vir(0,0);
          tvals_derivs[bin][nbase+1] += vir(0,1);
          tvals_derivs[bin][nbase+2] += vir(0,2);
          tvals_derivs[bin][nbase+3] += vir(1,0);
          tvals_derivs[bin][nbase+4] += vir(1,1);
          tvals_derivs[bin][nbase+5] += vir(1,2);
          tvals_derivs[bin][nbase+6] += vir(2,0);
          tvals_derivs[bin][nbase+7] += vir(2,1);
          tvals_derivs[bin][nbase+8] += vir(2,2);
        }
      }
    }
  }
  // Find maximum density
  double max = tvals[0];
  unsigned vout = 0;
  for(unsigned i=1; i<data.maxbins; ++i) {
    if( tvals[i]>max ) {
      max=tvals[i];
      vout=i;
    }
  }
  // Transform the density
  double df;
  double tsw = data.threshold_switch.calculate( max, df );
  if( fabs(sw*tsw)<epsilon ) {
    return;
  }

  if( !input.noderiv ) {
    Vector ddd = tsw*dfuncl*distance;
    output.deriv[0] = sw*df*max*tvals_derivs[vout][0] - ddd[0];
    output.deriv[1] = sw*df*max*tvals_derivs[vout][1] - ddd[1];
    output.deriv[2] = sw*df*max*tvals_derivs[vout][2] - ddd[2];
    output.deriv[3] = sw*df*max*tvals_derivs[vout][3] + ddd[0];
    output.deriv[4] = sw*df*max*tvals_derivs[vout][4] + ddd[1];
    output.deriv[5] = sw*df*max*tvals_derivs[vout][5] + ddd[2];
    for(unsigned i=0; i<input.natoms; ++i) {
      output.deriv[6 + 3*i + 0] = sw*df*max*tvals_derivs[vout][6 + 3*i + 0];
      output.deriv[6 + 3*i + 1] = sw*df*max*tvals_derivs[vout][6 + 3*i + 1];
      output.deriv[6 + 3*i + 2] = sw*df*max*tvals_derivs[vout][6 + 3*i + 2];
    }
    unsigned nbase = 6 + 3*input.natoms;
    Tensor vird(ddd,distance);
    output.deriv[nbase + 0] = sw*df*max*tvals_derivs[vout][nbase+0] - vird(0,0);
    output.deriv[nbase + 1] = sw*df*max*tvals_derivs[vout][nbase+1] - vird(0,1);
    output.deriv[nbase + 2] = sw*df*max*tvals_derivs[vout][nbase+2] - vird(0,2);
    output.deriv[nbase + 3] = sw*df*max*tvals_derivs[vout][nbase+3] - vird(1,0);
    output.deriv[nbase + 4] = sw*df*max*tvals_derivs[vout][nbase+4] - vird(1,1);
    output.deriv[nbase + 5] = sw*df*max*tvals_derivs[vout][nbase+5] - vird(1,2);
    output.deriv[nbase + 6] = sw*df*max*tvals_derivs[vout][nbase+6] - vird(2,0);
    output.deriv[nbase + 7] = sw*df*max*tvals_derivs[vout][nbase+7] - vird(2,1);
    output.deriv[nbase + 8] = sw*df*max*tvals_derivs[vout][nbase+8] - vird(2,2);
  }
  output.val[0] = sw*tsw;
}

}
}

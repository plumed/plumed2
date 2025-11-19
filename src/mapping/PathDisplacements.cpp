/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"
#include "PathProjectionCalculator.h"

//+PLUMEDOC ANALYSIS AVERAGE_PATH_DISPLACEMENT
/*
Accumulate the distances between the reference frames in the paths and the configurations visited

This action is used in the shortcut for [ADAPTIVE_PATH](ADAPTIVE_PATH.md) as it is helps us to adapt and refine
the defintition of a path CV based on the sampling that is observed in the trajectory.

By way of reminder path CVs were introduced by Branduardi _et al._ in the first paper cited in the bibliography below.
In PLUMED this method is implemented in [PATH](PATH.md) and works by definition a position along the path as:

$$
s = \frac{ \sum_{i=1}^N i \exp( -\lambda R[X - X_i] ) }{ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) }
$$

while the distance from the path (z) is measured using:

$$
z = -\frac{1}{\lambda} \ln\left[ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) \right]
$$

In these expressions $N$ high-dimensional frames ($X_i$) are used to describe the path in the high-dimensional
space. The two expressions above are then functions of the distances from each of the high-dimensional frames $R[X - X_i]$.

When a simulator uses this method of path CVs the assumption is that when the system undergoes the transition from state A
to state B it passes through a narrow (non-linear) tube that is centered on the $s$ coordinate that is defiend above in terms the
various milestones that are used to define the path. In other words, when the system transitions from A to be B it does not pass along the $s$
coordinate. If, however, we average the vector connecting each position the system passes through to the nearest
point on the path that underpins the definition of the $s$ we should get zero as the $s$ coordinate is centered on the path connecting the two states.

With that theory in mind this action can be used to collect the average displacement between the points that trajectories are passing through and the path.
This action can thus be used to update the definitions of the milestones that are used in the definition of a [PATH](PATH.md) or [GPATH](GPATH.md) coordinate
so that the PATH more clearly passes through the center of the path that the system has traversed along in passing from state A to state B.  These path displacements
are accumulated as follows.  You first determine the vectors, $\mathbf{v}_1$ and $\mathbf{v}_3$ that connect the instaneous position to the closest and second closest milestones on the
current path.  You can then compute the following:

$$
\delta = \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_2|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_3|^2) } + \mathbf{v}_1\cdot\mathbf{v}_2 }{2|\mathbf{v}_2|^2} - \frac{1}{2}
$$

where $\mathbf{v}_2$ is the vector connecting the closest and second closest milestone on your path. Displacements for these two milestones are then computed as:

$$
d_1 = (1 + \delta) \left(\mathbf{v}_1 - \delta \mathbf{v}_2\right) \qquad d_2 = -\delta \left(\mathbf{v}_1 - \delta \mathbf{v}_2\right)
$$

These displacement vectors for the two closest nodes are computed on every accumulate step and are then averaged with the factors of $(1 + \delta)$ and $-\delta$ in these expressions
serving as weights that are used when normalising these averages.

## Examples

Suppose that you are interested in a transition between two states A and B and that you have representative structures for these two states in the files `start.pdb` and `end.pdb`.  You can
use [pathtools](pathtools.md) to create an initial (linear) path connecting these two states as follows:

```plumed
plumed pathtools --start start.pdb --end end.pdb --nframes 20 --metric OPTIMAL --out path.pdb
```

This path is likely not even close to the transition pathway that the system takes between these two states as the assumption that the path is linear is pretty severe.  We can nevertheless use
a path generated by such a command in a [MOVINGRESTRAINT](MOVINGRESTRAINT.md) command like the one shown below:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1: PATH REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL LAMBDA=69087
mv: MOVINGRESTRAINT ...
   ARG=p1.s KAPPA0=1000
   STEP0=0 AT0=1
   STEP1=1000000 AT1=42
   STEP2=2000000 AT2=1
   STEP3=3000000 AT3=42
   STEP4=4000000 AT4=1
   STEP5=5000000 AT5=42
   STEP6=6000000 AT6=1
   STEP7=7000000 AT7=42
   STEP8=8000000 AT8=1
...
PRINT ARG=p1.s,p1.z FILE=colvar
```

This input hopefully (you need to check by looking at your trajectory) drives our system between the two states of interest 8 times - four of these transitions are from A to B and four are from B to A.
You can then analyze the trajectories that are generated using [driver](driver.md) and reparameterize a new path that passes more closely through the sampled trajectories using:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
rmsd: RMSD REFERENCE=regtest/trajectories/path_msd/all.pdb DISPLACEMENT TYPE=OPTIMAL
# Accumulate the average displacement between the reference path and the trajectories that have sampled the transition
disp: AVERAGE_PATH_DISPLACEMENT ARG=rmsd.disp STRIDE=1 METRIC={RMSD DISPLACEMENT TYPE=OPTIMAL ALIGN=1,1,1,1,1,1,1,1,1,1,1,1,1 DISPLACE=1,1,1,1,1,1,1,1,1,1,1,1,1} METRIC_COMPONENT=disp REFERENCE=rmsd_ref
# Now displace the original path by the accumulated displacement and reparameterize so that all frames are equally spaced
REPARAMETERIZE_PATH DISPLACE_FRAMES=disp FIXED=1,42 METRIC={RMSD DISPLACEMENT TYPE=OPTIMAL ALIGN=1,1,1,1,1,1,1,1,1,1,1,1,1 DISPLACE=1,1,1,1,1,1,1,1,1,1,1,1,1} METRIC_COMPONENT=disp REFERENCE=rmsd_ref
# And output the final reparameterized path at the end of the simulation
DUMPPDB DESCRIPTION=PATH STRIDE=0 FILE=outpatb.pdb ATOMS=rmsd_ref ATOM_INDICES=1-13
```

You can then try the same calculation with the reparameterized path before eventually using a method such as [METAD](METAD.md) to obtain the free energy as a function of your $s$ path coordinate.

Notice that the METRIC keyword appears here as the calculation involves calculating the vectors connecting the milestones on your path.  This keyword operates in the same way that is described in
the documentation for [GEOMETRIC_PATH](GEOMETRIC_PATH.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class PathDisplacements : public ActionWithValue, public ActionPilot, public ActionWithArguments {
private:
  bool clearnextstep;
  unsigned clearstride;
  double fadefact;
  std::vector<double> wsum, displace_v;
  Matrix<double> displacements;
  PathProjectionCalculator path_projector;
public:
  static void registerKeywords( Keywords& keys );
  explicit PathDisplacements(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void clearDerivatives( const bool& force=false ) {}
  void calculate() {}
  void apply() {}
  void update();
};

PLUMED_REGISTER_ACTION(PathDisplacements,"AVERAGE_PATH_DISPLACEMENT")

void PathDisplacements::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  PathProjectionCalculator::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("compulsory","STRIDE","1","the frequency with which the average displacements should be collected and added to the average displacements");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50 percent in the average. This option may increase convergence by allowing to forget the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.setValueDescription("matrix","matrix containing the average displacement between the trajectory and each of the landmarks that makes up the path");
  keys.addDOI("10.1063/1.2432340");
  keys.addDOI("10.1063/1.3660208");
}

PathDisplacements::PathDisplacements(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  clearnextstep(false),
  path_projector(this) {
  // Read in clear instructions
  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) {
      error("CLEAR parameter must be a multiple of STRIDE");
    }
    log.printf("  clearing average every %u steps \n",clearstride);
  }
  double halflife;
  parse("HALFLIFE",halflife);
  log.printf("  weight of contribution to frame halves every %f steps \n",halflife);
  if( halflife<0 ) {
    fadefact=1.0;
  } else {
    fadefact = exp( -0.693147180559945 / static_cast<double>(halflife) );
  }
  // Now create the weights vector and displacements matrix
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];
  wsum.resize( nrows );
  displacements.resize( nrows, ncols );
  for(unsigned i=0; i<nrows; ++i) {
    wsum[i]=0;
    for(unsigned j=0; j<ncols; ++j) {
      displacements(i,j)=0;
    }
  }
  // Add bibliography
  log<<"  Bibliography "<<plumed.cite("Diaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
  // And create a value to hold the displacements
  std::vector<std::size_t> shape(2);
  shape[0]=nrows;
  shape[1]=ncols;
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
}

unsigned PathDisplacements::getNumberOfDerivatives() {
  return 0;
}

void PathDisplacements::update() {
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];

  if( clearnextstep ) {
    unsigned k=0;
    for(unsigned i=0; i<nrows; ++i) {
      for(unsigned j=0; j<ncols; ++j) {
        displacements(i,j)=0;
        getPntrToComponent(0)->set(k,0);
        k++;
      }
    }
    clearnextstep=false;
  }

  unsigned k=0, iclose1=0, iclose2=0;
  double v1v1=0, v3v3=0;
  for(unsigned i=0; i<nrows; ++i) {
    double dist = 0;
    for(unsigned j=0; j<ncols; ++j) {
      double tmp = getPntrToArgument(0)->get(k);
      dist += tmp*tmp;
      k++;
    }
    if( i==0 ) {
      v1v1 = dist;
      iclose1 = 0;
    } else if( dist<v1v1 ) {
      v3v3=v1v1;
      v1v1=dist;
      iclose2=iclose1;
      iclose1=i;
    } else if( i==1 ) {
      v3v3=dist;
      iclose2=1;
    } else if( dist<v3v3 ) {
      v3v3=dist;
      iclose2=i;
    }
  }
  // And find third closest point
  int isign = iclose1 - iclose2;
  if( isign>1 ) {
    isign=1;
  } else if( isign<-1 ) {
    isign=-1;
  }
  int iclose3 = iclose1 + isign;
  unsigned ifrom=iclose1, ito=iclose3;
  if( iclose3<0 || static_cast<unsigned>(iclose3)>=nrows ) {
    ifrom=iclose2;
    ito=iclose1;
  }

  // Calculate the dot product of v1 with v2
  path_projector.getDisplaceVector( ifrom, ito, displace_v );
  double v2v2=0, v1v2=0;
  unsigned kclose1 = iclose1*ncols;
  for(unsigned i=0; i<displace_v.size(); ++i) {
    v2v2 += displace_v[i]*displace_v[i];
    v1v2 += displace_v[i]*getPntrToArgument(0)->get(kclose1+i);
  }

  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double weight2 = -1.* dx;
  double weight1 = 1.0 + dx;
  if( weight1>1.0 ) {
    weight1=1.0;
    weight2=0.0;
  } else if( weight2>1.0 ) {
    weight1=0.0;
    weight2=1.0;
  }

  // Accumulate displacements for path
  for(unsigned i=0; i<ncols; ++i) {
    double displace = getPntrToArgument(0)->get(kclose1+i) - dx*displace_v[i];
    displacements(iclose1,i) += weight1 * displace;
    displacements(iclose2,i) += weight2 * displace;
  }

  // Update weight accumulators
  wsum[iclose1] *= fadefact;
  wsum[iclose2] *= fadefact;
  wsum[iclose1] += weight1;
  wsum[iclose2] += weight2;

  // Update numbers in values
  if( wsum[iclose1] > epsilon ) {
    for(unsigned i=0; i<ncols; ++i) {
      getPntrToComponent(0)->set( kclose1+i, displacements(iclose1,i) / wsum[iclose1] );
    }
  }
  if( wsum[iclose2] > epsilon ) {
    unsigned kclose2 = iclose2*ncols;
    for(unsigned i=0; i<ncols; ++i) {
      getPntrToComponent(0)->set( kclose2+i, displacements(iclose2,i) / wsum[iclose2] );
    }
  }

  // Clear if required
  if( (getStep()>0 && clearstride>0 && getStep()%clearstride==0) ) {
    clearnextstep=true;
  }
}

}
}

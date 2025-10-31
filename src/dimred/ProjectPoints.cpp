/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "core/ActionRegister.h"
#include "tools/ConjugateGradient.h"
#include "tools/SwitchingFunction.h"
#include "tools/OpenMP.h"
#include "tools/Random.h"

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED PROJECT_POINTS
/*
Find the projection of a point in a low dimensional space by matching the (transformed) distance between it and a series of reference configurations that were input

This action and [ARRANGE_POINTS](ARRANGE_POINTS.md) are the workhorses for the implementation of [SKETCHMAP](SKETCHMAP.md) that is provided within PLUMED.
PROJECT_POINTS allows you to provide the low dimensional coordinate $y_\textrm{min}$ at which the following stress function is minimised:

$$
\chi(y) = \sum_{i=1}^N w_i [ D(X_i, Y) - d(x_i,y) ]^2
$$

where $Y$ is a set of coordinates in some high dimensional space, $X_i$ is a set of coordinates for one of $N$ landmark points in this high-dimensional space,
$x_i$ is a projection for $X_i$ that we have found in some lower dimensional space and $w_i$ is a weight.  The $D$ indicates that we are calculating the dissimilarity between the point $X_i$ and
$Y$, while $d$ represents the distance between the point $x_i$ and $y$.  In minimising the expression above we are thus finding the point $y$ at which the distances between $y$ and
each of the projections, $x_i$, of the $N$ landmark points most closely resembles the dissimiarities between $Y$ and the $N$ landmark points in the high-dimensional points.

The example input below illustrates how you can use PROJECT_POINTS to find the projection of a high-dimensional point in practice.

```plumed
# The coordinates of the landmarks in the high dimensional space
d1_ref: CONSTANT VALUES=1.0,2.0,1.5,2.1
d2_ref: CONSTANT VALUES=0.5,0.7,0.2,1.3
d3_ref: CONSTANT VALUES=3.1,2.0,1.5,0.5

# The weights of the landmark
weights: CONSTANT VALUES=1,1,1,1

# The projections of the landmarks in the low dimensional space
proj1_ref: CONSTANT VALUES=0.5,0.8,0.2,0.4
proj2_ref: CONSTANT VALUES=0.2,0.9,0.3,0.7

# Calcuate the instantaneous values of the three distances
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6

# Calculate the distances between the instananeous points and the current positions
ed: EUCLIDEAN_DISTANCE SQUARED ARG2=d1,d2,d3 ARG1=d1_ref,d2_ref,d3_ref

# And generate the projection
proj: PROJECT_POINTS ARG=proj1_ref,proj2_ref TARGET1=ed WEIGHTS1=weights

# And output the projection to a file
PRINT ARG=proj.* FILE=colvar
```

In this example, we use three distances to define the high dimensional coordinates and have four landmarks points.

## Projecting multiple coordinates at once

The input to the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md) shortcut in the example in the previous section consisted of three
scalar-valued quantities.  The two components output by project points are thus also scalars. By contrast, in the input below four
three-dimenional vectors are input to the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md) shortcut. The components output by proj and thus
four-dimensional vectors.

```plumed
# The coordinates of the landmarks in the high dimensional space
d1_ref: CONSTANT VALUES=1.0,2.0,1.5,2.1
d2_ref: CONSTANT VALUES=0.5,0.7,0.2,1.3
d3_ref: CONSTANT VALUES=3.1,2.0,1.5,0.5

# The weights of the landmark
weights: CONSTANT VALUES=1,1,1,1

# The projections of the landmarks in the low dimensional space
proj1_ref: CONSTANT VALUES=0.5,0.8,0.2,0.4
proj2_ref: CONSTANT VALUES=0.2,0.9,0.3,0.7

# Calcuate the instantaneous values of the distances
d1: DISTANCE ATOMS1=1,2 ATOMS2=7,8   ATOMS3=13,14 ATOMS4=19,20
d2: DISTANCE ATOMS1=3,4 ATOMS2=9,10  ATOMS3=15,16 ATOMS4=21,22
d3: DISTANCE ATOMS1=5,6 ATOMS2=11,12 ATOMS3=17,18 ATOMS4=23,24

# Calculate the distances between the instananeous points and the current positions
ed: EUCLIDEAN_DISTANCE SQUARED ARG1=d1,d2,d3 ARG2=d1_ref,d2_ref,d3_ref

# And generate the projection
proj: PROJECT_POINTS ARG=proj1_ref,proj2_ref TARGET1=ed WEIGHTS1=weights

# And output the projection to a file
PRINT ARG=proj.* FILE=colvar
```

## Using RMSD distances

One can use [RMSD](RMSD.md) distances as the dissimilarities rather than distances in some space of arguments as is illustrated below:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/allv.pdb

# This action reads in the landmarks in the high dimensional space and calculates the
# distances from the instantaneous configuration
rmsd: RMSD SQUARED TYPE=OPTIMAL REFERENCE=regtest/trajectories/path_msd/allv.pdb

# The weights of the landmarks
weights: ONES SIZE=42

# The projections of the landmarks in the low dimensional space
X: PDB2CONSTANT ARG=X NOARGS REFERENCE=regtest/trajectories/path_msd/allv.pdb
Y: PDB2CONSTANT ARG=Y NOARGS REFERENCE=regtest/trajectories/path_msd/allv.pdb

# Generate the projection of the instantaneous coordinates
proj: PROJECT_POINTS ARG=X,Y TARGET1=rmsd WEIGHTS1=weights

# And output the projection to a file
PRINT ARG=proj.* FILE=colvar
```

For ths input there are 42 landmark points and dissimilarities are computed by computing the RMSD distance between the 13 atoms in
each of landmark coordinates and the instaneous positions of those 13 atoms.

## Using transformed distances

In [SKETCHMAP](SKETCHMAP.md) the stress function that is minimised is not the one given above.  Instead of seeking to generate a projection,
$y_\textrm{min}$, which is at a point where the distances between it and each projection the landmarks is the same as the dissimilarities between
the high-dimensional coordinate of the point and the high-dimensional landmarks, the dissimilarities and distances are transformed by functions as illustrated below:

$$
\chi(y) = \sum_{i=1}^N w_i [ F[D(X_i, Y)] - f[d(x_i,y)] ]^2
$$

The two functions $F$ and $f$ in this expression are usually different as you can see in the input below:

```plumed
# The coordinates of the landmarks in the high dimensional space
d1_ref: CONSTANT VALUES=1.0,2.0,1.5,2.1
d2_ref: CONSTANT VALUES=0.5,0.7,0.2,1.3
d3_ref: CONSTANT VALUES=3.1,2.0,1.5,0.5

# The weights of the landmark
weights: CONSTANT VALUES=1,1,1,1

# The projections of the landmarks in the low dimensional space
proj1_ref: CONSTANT VALUES=0.5,0.8,0.2,0.4
proj2_ref: CONSTANT VALUES=0.2,0.9,0.3,0.7

# Calcuate the instantaneous values of the distances
d1: DISTANCE ATOMS1=1,2 ATOMS2=7,8   ATOMS3=13,14 ATOMS4=19,20
d2: DISTANCE ATOMS1=3,4 ATOMS2=9,10  ATOMS3=15,16 ATOMS4=21,22
d3: DISTANCE ATOMS1=5,6 ATOMS2=11,12 ATOMS3=17,18 ATOMS4=23,24

# Calculate the distances between the instananeous points and the current positions
ed: EUCLIDEAN_DISTANCE SQUARED ARG1=d1,d2,d3 ARG2=d1_ref,d2_ref,d3_ref

# Transform the dissimilarities by applying the funciton F
fed: MORE_THAN ARG=ed SQUARED SWITCH={SMAP R_0=4 A=3 B=2}

# And generate the projection
proj: PROJECT_POINTS ARG=proj1_ref,proj2_ref TARGET1=fed FUNC1={SMAP R_0=4 A=1 B=2} WEIGHTS1=weights

# And output the projection to a file
PRINT ARG=proj.* FILE=colvar
```

In the input above the function, $F$, that is applied on the dissimilarities is implemented using a [MORE_THAN](MORE_THAN.md) action. The function, $f$,
that is applied on the distances in the low-dimensional space is specified using the `FUNC` keyword that is input to PROJECT_POINTS.

## Using multiple targets

At its most complex this action allows you to minimise a stress function such as the one below:

$$
\chi(y) = \sum_{i=1}^N \sum_{j=1}^M w_{ij} [ F_j[D(X_i, Y)] - f_j[d(x_i,y)] ]^2
$$

The input below shows how this can be implemted within PLUMED:

```plumed
# The coordinates of the landmarks in the high dimensional space
d1_ref: CONSTANT VALUES=1.0,2.0,1.5,2.1
d2_ref: CONSTANT VALUES=0.5,0.7,0.2,1.3
d3_ref: CONSTANT VALUES=3.1,2.0,1.5,0.5

# The weights of the landmark
weights: CONSTANT VALUES=1,1,1,1
w1: CUSTOM ARG=weights FUNC=0.3*x PERIODIC=NO
w2: CUSTOM ARG=weights FUNC=(1-0.3)*x PERIODIC=NO

# The projections of the landmarks in the low dimensional space
proj1_ref: CONSTANT VALUES=0.5,0.8,0.2,0.4
proj2_ref: CONSTANT VALUES=0.2,0.9,0.3,0.7

# Calcuate the instantaneous values of the distances
d1: DISTANCE ATOMS1=1,2 ATOMS2=7,8   ATOMS3=13,14 ATOMS4=19,20
d2: DISTANCE ATOMS1=3,4 ATOMS2=9,10  ATOMS3=15,16 ATOMS4=21,22
d3: DISTANCE ATOMS1=5,6 ATOMS2=11,12 ATOMS3=17,18 ATOMS4=23,24

# Calculate the distances between the instananeous points and the current positions
ed: EUCLIDEAN_DISTANCE SQUARED ARG1=d1,d2,d3 ARG2=d1_ref,d2_ref,d3_ref
# Transform the dissimilarities by applying the funciton F
fed: MORE_THAN ARG=ed SQUARED SWITCH={SMAP R_0=4 A=3 B=2}

# And generate the projection
proj: PROJECT_POINTS ...
  ARG=proj1_ref,proj2_ref
  TARGET1=ed WEIGHTS1=w1 FUNC1={CUSTOM FUNC=1-sqrt(x2) R_0=1.0}
  TARGET2=fed WEIGHTS2=w2 FUNC2={SMAP R_0=4 A=1 B=2}
...

# And output the projection to a file
PRINT ARG=proj.* FILE=colvar
```

Here the sum over $M$ in the expression above has two terms. In the first of these terms $F_1$ is the identity so the
input for `TARGET1` is the output from [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md). $f_1$ is similarly the identity.  To
implement the identity here we use the input to `FUNC1` shown above.  The input to this function is the input for one of
the switching functions described in the documentation for [LESS_THAN](LESS_THAN.md). What we compute for the transformed
distance is $1-s(d)$ where $s(d)$ is the switching function that is specified in input.  Consequently, applying the
function `1-sqrt(x2)` returns the distance.

The second term in our sum over $M$ in the input above has the dissimilarities and distances transformed by the functions that
we introduced in the previous section.

*/
//+ENDPLUMEDOC

class ProjectPoints;

class ProjectPointsInput {
public:
  double cgtol;
  ProjectPoints* action;
};

class ProjectPoints : public ActionWithVector {
public:
  using input_type = ProjectPointsInput;
  using PTM = ParallelTaskManager<ProjectPoints>;
private:
  unsigned dimout;
  mutable std::vector<unsigned> rowstart;
  std::vector<SwitchingFunction> switchingFunction;
  ConjugateGradient<ProjectPoints> myminimiser;
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  ProjectPoints( const ActionOptions& );
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void prepare() override ;
  static void performTask( std::size_t task_index, const ProjectPointsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(ProjectPoints,"PROJECT_POINTS")

void ProjectPoints::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","vector","the projections of the landmark points");
  keys.addInputKeyword("numbered","TARGET","vector/matrix","the matrix of target quantities that you would like to match");
  keys.add("numbered","FUNC","a function that is applied on the distances between the points in the low dimensional space");
  keys.addInputKeyword("numbered","WEIGHTS","vector","the matrix with the weights of the target quantities");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
  keys.addOutputComponent("coord","default","scalar/vector","the coordinates of the points in the low dimensional space");
  PTM::registerKeywords( keys );
}


ProjectPoints::ProjectPoints( const ActionOptions& ao ) :
  Action(ao),
  ActionWithVector(ao),
  rowstart(OpenMP::getNumThreads()),
  myminimiser(this),
  taskmanager(this) {
  dimout = getNumberOfArguments();
  unsigned nvals=getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( nvals!=getPntrToArgument(i)->getNumberOfValues() ) {
      error("mismatch between numbers of projections");
    }
  }
  std::vector<Value*> args( getArguments() ), target, weights;
  std::string sfd, errors;
  unsigned ntoproj=0;
  // Read in target "distances" and target weights
  for(unsigned i=1;; ++i) {
    target.resize(0);
    if( !parseArgumentList("TARGET",i,target) ) {
      break;
    }
    std::string inum;
    Tools::convert( i, inum );
    if( target.size()!=1 ) {
      error("should only be one value in input to TARGET" + inum );
    }
    if( (target[0]->getRank()!=1 && target[0]->getRank()!=2) || target[0]->hasDerivatives() ) {
      error("input to TARGET" + inum + " keyword should be a vector/matrix");
    }
    if( target[0]->getShape()[0]!=nvals ) {
      error("number of rows in target matrix should match number of input coordinates");
    }
    if( i==1 && target[0]->getRank()==1 ) {
      ntoproj = 1;
    } else if( ntoproj==1 && target[0]->getRank()!=1 ) {
      error("mismatch between numbers of target distances");
    } else if( i==1 ) {
      ntoproj = target[0]->getShape()[1];
    } else if( target[0]->getRank()>1 && ntoproj!=target[0]->getShape()[1] ) {
      error("mismatch between numbers of target distances");
    }
    if( !parseArgumentList("WEIGHTS",i,weights) ) {
      error("missing WEIGHTS" + inum + " keyword in input");
    }
    if( weights.size()!=1 ) {
      error("should only be one value in input to WEIGHTS" + inum );
    }
    if( weights[0]->getRank()!=1 || weights[0]->hasDerivatives() ) {
      error("input to WEIGHTS" + inum + " keyword should be a vector");
    }
    if( weights[0]->getShape()[0]!=nvals ) {
      error("number of weights should match number of input coordinates");
    }
    args.push_back( target[0] );
    args.push_back( weights[0] );
    bool has_sf = parseNumbered("FUNC",i,sfd);
    switchingFunction.push_back( SwitchingFunction() );
    if( !has_sf ) {
      switchingFunction[i-1].set( "CUSTOM FUNC=1-sqrt(x2) R_0=1.0", errors );
    } else {
      switchingFunction[i-1].set( sfd, errors );
      if( errors.length()!=0 ) {
        error("problem reading switching function description " + errors);
      }
    }
    log.printf("  %sth term seeks to match tranformed distances with those in matrix %s \n", inum.c_str(), target[0]->getName().c_str() );
    log.printf("  in %sth term distances are transformed by 1-switching function with r_0=%s \n", inum.c_str(), switchingFunction[i-1].description().c_str() );
    log.printf("  in %sth term weights of matrix elements in stress function are given by %s \n", inum.c_str(), weights[0]->getName().c_str() );
  }
  std::vector<std::size_t> shape(1);
  shape[0]=ntoproj;
  if( ntoproj==1 ) {
    shape.resize(0);
  }
  for(unsigned i=0; i<dimout; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    addComponent( "coord-" + num, shape );
    componentIsNotPeriodic( "coord-" + num );
  }
  // Create a list of tasks to perform
  double cgtol;
  parse("CGTOL",cgtol);
  log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  requestArguments( args );
  checkRead();

  // Setup parallel task manager
  ProjectPointsInput input;
  input.cgtol=cgtol;
  input.action=this;
  if( ntoproj!=1 ) {
    taskmanager.setupParallelTaskManager( 0, 0 );
  }
  taskmanager.setActionInput( input );
}

void ProjectPoints::prepare() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    return;
  }

  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(dimout)->getShape()[0];
  for(unsigned i=0; i<dimout; ++i) {
    if( getPntrToComponent(i)->getShape()[0]!=shape[0] ) {
      getPntrToComponent(i)->setShape(shape);
    }
  }
}

double ProjectPoints::calculateStress( const std::vector<double>& pp, std::vector<double>& der ) {
  unsigned nmatrices = ( getNumberOfArguments() - dimout ) / 2;
  double stress=0;

  unsigned t=OpenMP::getThreadNum();
  std::vector<double> dtmp( pp.size() );
  unsigned nland = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<nland; ++i) {
    // Calculate distance in low dimensional space
    double dd2 = 0;
    for(unsigned k=0; k<pp.size(); ++k) {
      dtmp[k] = pp[k] - getPntrToArgument(k)->get(i);
      dd2 += dtmp[k]*dtmp[k];
    }

    for(unsigned k=0; k<nmatrices; ++k ) {
      // Now do transformations and calculate differences
      double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
      // Get the weight for this connection
      double weight = getPntrToArgument( dimout + 2*k + 1 )->get( i );
      // Get the difference for the connection
      double fdiff = fd - getPntrToArgument( dimout + 2*k )->get( rowstart[t]+i );
      // Calculate derivatives
      double pref = -2.*weight*fdiff*df;
      for(unsigned n=0; n<pp.size(); ++n) {
        der[n]+=pref*dtmp[n];
      }
      // Accumulate the total stress
      stress += weight*fdiff*fdiff;
    }
  }
  return stress;
}

void ProjectPoints::performTask( std::size_t task_index, const ProjectPointsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  // I doubt we are ever going to implement this on the GPU so I think we can leave this declaration here
  std::vector<double> point( input.ncomponents );
  std::size_t nland = input.shapedata[0];
  std::size_t base = task_index;
  if( input.ranks[input.ncomponents]==2 ) {
    auto myargh = ArgumentBookeepingHolder::create( input.ncomponents, input );
    base = task_index*myargh.shape[1];
  }
  unsigned closest=0;
  double mindist = input.inputdata[input.argstarts[input.ncomponents] + base];
  for(unsigned i=1; i<nland; ++i) {
    double dist = input.inputdata[input.argstarts[input.ncomponents] + base+i];
    if( dist<mindist ) {
      mindist=dist;
      closest=i;
    }
  }
  // Put the initial guess near to the closest landmark  -- may wish to use grid here again Sandip??
  Random random;
  random.setSeed(-1234);
  for(unsigned j=0; j<input.ncomponents; ++j) {
    point[j] = input.inputdata[input.argstarts[j] + closest] + (random.RandU01() - 0.5)*0.01;
  }
  // And do the optimisation
  actiondata.action->rowstart[OpenMP::getThreadNum()]=task_index;
  if( input.ranks[input.ncomponents]==2 ) {
    auto myargh=ArgumentBookeepingHolder::create( input.ncomponents, input );
    actiondata.action->rowstart[OpenMP::getThreadNum()] = task_index*myargh.shape[1];
  }
  actiondata.action->myminimiser.minimise( actiondata.cgtol, point, &ProjectPoints::calculateStress );
  for(unsigned i=0; i<input.ncomponents; ++i) {
    output.values[i] = point[i];
  }
}

void ProjectPoints::calculate() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    auto myinput =ParallelActionsInput::create( getPbc() );
    myinput.noderiv = true;
    myinput.ncomponents = getNumberOfComponents();
    std::vector<double> input_buffer;
    getInputData( input_buffer );
    myinput.dataSize = input_buffer.size();
    myinput.inputdata = input_buffer.data();
    ArgumentsBookkeeping abk;
    abk.setupArguments( this );
    myinput.setupArguments( abk );
    std::vector<double> buffer;
    std::vector<double> derivatives, point( getNumberOfComponents() );
    auto output = ParallelActionsOutput::create( myinput.ncomponents, point.data(), 0, derivatives.data(), 0, buffer.data() );
    performTask( 0, taskmanager.getActionInput(), myinput, output );
    for(unsigned i=0; i<point.size(); ++i) {
      getPntrToComponent(i)->set(point[i]);
    }
  } else {
    taskmanager.runAllTasks();
  }
}

}
}

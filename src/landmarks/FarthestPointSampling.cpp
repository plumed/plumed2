/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "matrixtools/MatrixOperationBase.h"
#include "core/ActionRegister.h"
#include "tools/Random.h"

//+PLUMEDOC LANDMARKS FARTHEST_POINT_SAMPLING
/*
Select a set of landmarks using farthest point sampling.

This action is used within the [LANDMARK_SELECT_FPS](LANDMARK_SELECT_FPS.md) shortcut, which performs farthest point sampling.
Farthest point sampling is a method of selecting a subset of input coordinates that works by selecting a first point at random.
The remaining points are then selected by taking the unselected point in the input data set that is the furthest
from all the points that have been selected thus far.

This particular action is designed to be used in conjunction with [SELECT_WITH_MASK](SELECT_WITH_MASK.md) in the same way that
[CREATE_MASK](CREATE_MASK.md) is designed to be be used with that action.  This action takes an $N\times N$ matrix of dissimilarities
in input and outputs and $N$ dimensional vector whose elements are ones and zeros. As shown in the example input below, this output
vector is passed to a [SELECT_WITH_MASK](SELECT_WITH_MASK.md) using the MASK keyword.

The example below thus shows how this action is used in the [LANDMARK_SELECT_FPS](LANDMARK_SELECT_FPS.md) shortcut to select 100 landmark
frames from the trajectory using farthest point sampling.

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL

# We now compute the dissimilarities between these frames
cc_dataT: TRANSPOSE ARG=cc_data
dd: DISSIMILARITIES ARG=cc_data,cc_dataT

# Create a mask that will be used to create our landmark data
mask: FARTHEST_POINT_SAMPLING ARG=dd NZEROS=100

# These are the landmark points
landmarks: SELECT_WITH_MASK ARG=cc_data ROW_MASK=mask

# Output the landmarks to a file
DUMPPDB ATOMS=landmarks ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb
```

This only saves the coordinates of the landmark points.  If you look at the shortcuts in the documentation for [LANDMARK_SELECT_FPS](LANDMARK_SELECT_FPS.md)
you can see how the shortcut also saves information on the dissimilarities between the landmarks, the dissimilarities between the landmarks and all the other
points and the weights of the landmarks, which are determined using a [VORONOI](VORONOI.md) analysis.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace landmarks {

class FarthestPointSampling : public matrixtools::MatrixOperationBase {
private:
  unsigned seed;
  unsigned nlandmarks;
public:
  static void registerKeywords( Keywords& keys );
  explicit FarthestPointSampling( const ActionOptions& ao );
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void prepare() override ;
  void calculate() override ;
  void apply() override {}
  double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const override {
    plumed_merror("this should not be called");
  }
};

PLUMED_REGISTER_ACTION(FarthestPointSampling,"FARTHEST_POINT_SAMPLING")

void FarthestPointSampling::registerKeywords( Keywords& keys ) {
  matrixtools::MatrixOperationBase::registerKeywords( keys );
  keys.add("compulsory","NZEROS","the number of landmark points that you want to select");
  keys.add("compulsory","SEED","1234","a random number seed");
  keys.setValueDescription("vector","a vector which has as many elements as there are rows in the input matrix of dissimilarities. NZEROS of the elements in this vector are equal to one, the rest of the elements are equal to zero.  The nodes that have elements equal to one are the NZEROS points that are farthest appart according to the input dissimilarities");
}

FarthestPointSampling::FarthestPointSampling( const ActionOptions& ao ):
  Action(ao),
  MatrixOperationBase(ao) {
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) {
    error("input to this argument should be a square matrix of dissimilarities");
  }
  parse("NZEROS",nlandmarks);
  parse("SEED",seed);
  log.printf("  selecting %d landmark points \n", nlandmarks );

  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
}

void FarthestPointSampling::prepare() {
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]!=getPntrToArgument(0)->getShape()[0] ) {
    std::vector<std::size_t> shape(1);
    shape[0] = getPntrToArgument(0)->getShape()[0];
    myval->setShape(shape);
  }
  for(unsigned i=0; i<nlandmarks; ++i) {
    myval->set( i, 0.0 );
  }
  for(unsigned i=nlandmarks; i<myval->getShape()[0]; ++i) {
    myval->set( i, 1.0 );
  }
}

void FarthestPointSampling::calculate() {
  Value* myval=getPntrToComponent(0);
  unsigned npoints = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<npoints; ++i) {
    myval->set( i, 1.0 );
  }
  std::vector<unsigned> landmarks( nlandmarks );

  // Select first point at random
  Random random;
  random.setSeed(-seed);
  double rand=random.RandU01();
  landmarks[0] = std::floor( npoints*rand );
  myval->set( landmarks[0], 0 );

  // Now find distance to all other points (N.B. We can use squared distances here for speed)
  Matrix<double> distances( nlandmarks, npoints );
  Value* myarg = getPntrToArgument(0);
  for(unsigned i=0; i<npoints; ++i) {
    distances(0,i) = myarg->get( landmarks[0]*npoints + i );
  }

  // Now find all other landmarks
  for(unsigned i=1; i<nlandmarks; ++i) {
    // Find point that has the largest minimum distance from the landmarks selected thus far
    double maxd=0;
    for(unsigned j=0; j<npoints; ++j) {
      double mind=distances(0,j);
      for(unsigned k=1; k<i; ++k) {
        if( distances(k,j)<mind ) {
          mind=distances(k,j);
        }
      }
      if( mind>maxd ) {
        maxd=mind;
        landmarks[i]=j;
      }
    }
    myval->set( landmarks[i], 0 );
    for(unsigned k=0; k<npoints; ++k) {
      distances(i,k) = myarg->get( landmarks[i]*npoints + k );
    }
  }
}

}
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "PathProjectionCalculator.h"

//+PLUMEDOC COLVAR GEOMETRIC_PATH
/*
Distance along and from a path calculated using geometric formulas

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class GeometricPath : public ActionWithVector {
private:
  PathProjectionCalculator path_projector;
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPath(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override {
    plumed_error();
  }
};

PLUMED_REGISTER_ACTION(GeometricPath,"GEOMETRIC_PATH")

void GeometricPath::registerKeywords(Keywords& keys) {
  ActionWithVector::registerKeywords(keys);
  PathProjectionCalculator::registerKeywords(keys);
  keys.addInputKeyword("compulsory","PROPERTY","vector","the label of a value that contains the coordinates we are projecting these points onto");
  keys.addOutputComponent("s","default","scalar","the position on the path");
  keys.addOutputComponent("z","default","scalar","the distance from the path");
  keys.needsAction("GEOMETRIC_PATH");
  keys.needsAction("PDB2CONSTANT");
}

GeometricPath::GeometricPath(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  path_projector(this) {
  plumed_assert( !actionInChain() );
  // Get the coordinates in the low dimensional space
  std::vector<std::string> pcoord;
  parseVector("PROPERTY", pcoord );
  std::vector<Value*> theprop;
  ActionWithArguments::interpretArgumentList( pcoord, plumed.getActionSet(), this, theprop );
  if( theprop.size()!=1 ) {
    error("did not find property to project on");
  }
  if( theprop[0]->getNumberOfValues()!=getPntrToArgument(0)->getShape()[0] ) {
    error("mismatch between number of frames and property of interest");
  }
  log.printf("  projecting onto : %s \n", theprop[0]->getName().c_str() );
  std::vector<Value*> args( getArguments() );
  args.push_back( theprop[0] );
  requestArguments( args );
  // Create the values to store the output
  addComponentWithDerivatives("s");
  componentIsNotPeriodic("s");
  addComponentWithDerivatives("z");
  componentIsNotPeriodic("z");
}

unsigned GeometricPath::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getShape()[0]*getPntrToArgument(0)->getShape()[1] + getPntrToArgument(1)->getShape()[0];
}

void GeometricPath::calculate() {
  unsigned k=0, iclose1=0, iclose2=0;
  double v1v1=0, v3v3=0;
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];
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
  if( iclose3<0 || iclose3>=nrows ) {
    ifrom=iclose2;
    ito=iclose1;
  }

  // And calculate projection of vector connecting current point to closest frame on vector connecting nearest two frames
  std::vector<double> displace;
  path_projector.getDisplaceVector( ifrom, ito, displace );
  double v2v2=0, v1v2=0;
  k=ncols*iclose1;
  for(unsigned i=0; i<displace.size(); ++i) {
    v2v2 += displace[i]*displace[i];
    v1v2 += displace[i]*getPntrToArgument(0)->get(k+i);
  }

  // This computes s value
  double spacing = getPntrToArgument(1)->get(iclose1) - getPntrToArgument(1)->get(iclose2);
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double path_s = getPntrToArgument(1)->get(iclose1) + spacing * dx;
  Value* sp = getPntrToComponent(0);
  sp->set( path_s );
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<ncols; ++i) {
      sp->addDerivative( ncols*iclose1 + i, 0.5*spacing*(v1v2*displace[i]/v2v2 - getPntrToArgument(0)->get(ncols*iclose1 + i))/root + 0.5*spacing*displace[i]/v2v2 );
      sp->addDerivative( ncols*iclose2 + i, 0.5*spacing*getPntrToArgument(0)->get(ncols*iclose2 + i)/root );
    }
  }

  // This computes z value
  path_projector.getDisplaceVector( iclose2, iclose1, displace );
  double v4v4=0, proj=0;
  k=ncols*iclose1;
  for(unsigned i=0; i<displace.size(); ++i) {
    v4v4 += displace[i]*displace[i];
    proj += displace[i]*getPntrToArgument(0)->get(k+i);
  }
  double path_z = v1v1 + dx*dx*v4v4 - 2*dx*proj;
  path_z = sqrt(path_z);
  Value* zp = getPntrToComponent(1);
  zp->set( path_z );
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<ncols; ++i) {
      zp->addDerivative( ncols*iclose1 + i, (1/path_z)*(getPntrToArgument(0)->get(ncols*iclose1 + i) +
                         (v4v4*dx-proj)*sp->getDerivative(ncols*iclose1 + i)/spacing -
                         dx*displace[i]) );
      zp->addDerivative( ncols*iclose2 + i, (v4v4*dx-proj)*sp->getDerivative(ncols*iclose2 + i)/(path_z*spacing) );
    }
  }
}

}
}

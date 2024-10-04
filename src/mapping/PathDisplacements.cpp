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

\par Examples

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
  keys.use("ARG");
  PathProjectionCalculator::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the average displacements should be collected and added to the average displacements");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50 percent in the average. This option may increase convergence by allowing to forget the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.setValueDescription("vector containing the average displacement between the trajectory and each of the landmarks that makes up the path");
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
  std::vector<unsigned> shape(2);
  shape[0]=nrows;
  shape[1]=ncols;
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->buildDataStore();
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
  if( iclose3<0 || iclose3>=nrows ) {
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

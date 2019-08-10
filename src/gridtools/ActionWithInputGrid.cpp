/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "ActionWithInputGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void ActionWithInputGrid::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true),
  set_zero_outside_range(false)
{
  if( getNumberOfArguments()!=1 ) error("should be exactly one argument to this action");
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) error("input to this action should be a grid");

  unsigned dimension = getPntrToArgument(0)->getRank();
  std::vector<std::string> argn( dimension ), min( dimension ), max( dimension ); std::string gtype;
  std::vector<unsigned> nbin( dimension ); std::vector<double> spacing( dimension ); std::vector<bool> ipbc( dimension );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
  if( gtype=="flat" ) gridobject.setup( "flat", ipbc, 0, 0.0 );
  else if( gtype=="fibonacci" ) gridobject.setup( "fibonacci", ipbc, nbin[0], spacing[0] );
  else plumed_merror("unknown grid type");
}

double ActionWithInputGrid::getFunctionValueAndDerivatives( const std::vector<double>& x, std::vector<double>& der ) const {
  plumed_dbg_assert( gridobject.getGridType()=="flat" ); unsigned dimension = gridobject.getDimension();
  if( set_zero_outside_range && !gridobject.inbounds( x ) ) return 0.0;

  double X,X2,X3,value=0; der.assign(der.size(),0.0);
  std::vector<double> fd(dimension);
  std::vector<double> C(dimension);
  std::vector<double> D(dimension);
  std::vector<double> dder(dimension);

  std::vector<unsigned> nindices(dimension);
  std::vector<unsigned> indices(dimension); gridobject.getIndices( x, indices );
  std::vector<unsigned> neigh; gridobject.getSplineNeighbors( gridobject.getIndex(indices), neigh );
  std::vector<double> xfloor(dimension); gridobject.getGridPointCoordinates( gridobject.getIndex(x), nindices, xfloor );

  // loop over neighbors
  for(unsigned int ipoint=0; ipoint<neigh.size(); ++ipoint) {
    double grid=getFunctionValue( neigh[ipoint] );
    for(unsigned j=0; j<dimension; ++j) dder[j] = getPntrToArgument(0)->getGridDerivative( neigh[ipoint], j );

    gridobject.getIndices( neigh[ipoint], nindices );
    double ff=1.0;
    for(unsigned j=0; j<dimension; ++j) {
      int x0=1;
      if(nindices[j]==indices[j]) x0=0;
      double ddx=gridobject.getGridSpacing()[j];
      X=fabs((x[j]-xfloor[j])/ddx-(double)x0);
      X2=X*X;
      X3=X2*X;
      double yy;
      if(fabs(grid)<0.0000001) yy=0.0;
      else yy=-dder[j]/grid;
      C[j]=(1.0-3.0*X2+2.0*X3) - (x0?-1.0:1.0)*yy*(X-2.0*X2+X3)*ddx;
      D[j]=( -6.0*X +6.0*X2) - (x0?-1.0:1.0)*yy*(1.0-4.0*X +3.0*X2)*ddx;
      D[j]*=(x0?-1.0:1.0)/ddx;
      ff*=C[j];
    }
    for(unsigned j=0; j<dimension; ++j) {
      fd[j]=D[j];
      for(unsigned i=0; i<dimension; ++i) if(i!=j) fd[j]*=C[i];
    }
    value+=grid*ff;
    for(unsigned j=0; j<dimension; ++j) der[j]+=grid*fd[j];
  }
  return value;
}

void ActionWithInputGrid::setupGridObject() {
  plumed_assert( firststep ); 
  if( gridobject.getGridType()=="flat" ) {
    unsigned dimension = getPntrToArgument(0)->getRank();
    std::vector<std::string> argn( dimension ), min( dimension ), max( dimension ); std::string gtype;
    std::vector<unsigned> nbin( dimension ); std::vector<double> spacing( dimension ); std::vector<bool> ipbc( dimension );
    (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
    gridobject.setBounds( min, max, nbin, spacing );
  }
}

void ActionWithInputGrid::doTheCalculation() {
  if( firststep ) { setupGridObject(); firststep=false; finishOutputSetup(); }
  runAllTasks();
  jobsAfterLoop();
}

void ActionWithInputGrid::calculate() {
  if( actionInChain() ) return ;
  doTheCalculation();
}

void ActionWithInputGrid::update() {
  if( skipUpdate() ) return;
  doTheCalculation();
}

void ActionWithInputGrid::runFinalJobs() {
  if( skipUpdate() ) return;
  doTheCalculation();
}

}
}


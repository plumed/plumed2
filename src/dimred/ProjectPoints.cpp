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

\par Examples

*/
//+ENDPLUMEDOC

class ProjectPoints : public ActionWithVector {
private:
  double cgtol;
  unsigned dimout;
  mutable std::vector<unsigned> rowstart;
  std::vector<SwitchingFunction> switchingFunction;
  ConjugateGradient<ProjectPoints> myminimiser;
  void getProjection( const unsigned& current, std::vector<double>& point ) const ;
public:
  static void registerKeywords( Keywords& keys );
  ProjectPoints( const ActionOptions& );
  unsigned getNumberOfDerivatives() {
    return 0;
  }
  void prepare() override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(ProjectPoints,"PROJECT_POINTS")

void ProjectPoints::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.use("ARG");
  keys.add("numbered","TARGET","the matrix of target quantities that you would like to match");
  keys.add("numbered","FUNC","a function that is applied on the distances between the points in the low dimensional space");
  keys.add("numbered","WEIGHTS","the matrix with the weights of the target quantities");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
  keys.addOutputComponent("coord","default","the coordinates of the points in the low dimensional space");
}


ProjectPoints::ProjectPoints( const ActionOptions& ao ) :
  Action(ao),
  ActionWithVector(ao),
  rowstart(OpenMP::getNumThreads()),
  myminimiser( this ) {
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
    } else if( ntoproj!=target[0]->getShape()[1] ) {
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
    target[0]->buildDataStore();
    weights[0]->buildDataStore();
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
  std::vector<unsigned> shape(1);
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
  parse("CGTOL",cgtol);
  log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  requestArguments( args );
  checkRead();
}

void ProjectPoints::prepare() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    return;
  }

  std::vector<unsigned> shape(1);
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

void ProjectPoints::getProjection( const unsigned& current, std::vector<double>& point ) const {
  Value* targ = getPntrToArgument( dimout );
  unsigned nland = getPntrToArgument(0)->getShape()[0];
  unsigned base = current;
  if( targ->getRank()==2 ) {
    base = current*targ->getShape()[1];
  }
  unsigned closest=0;
  double mindist = targ->get( base );
  for(unsigned i=1; i<nland; ++i) {
    double dist = targ->get( base + i );
    if( dist<mindist ) {
      mindist=dist;
      closest=i;
    }
  }
  // Put the initial guess near to the closest landmark  -- may wish to use grid here again Sandip??
  Random random;
  random.setSeed(-1234);
  for(unsigned j=0; j<dimout; ++j) {
    point[j] = getPntrToArgument(j)->get(closest) + (random.RandU01() - 0.5)*0.01;
  }
  // And do the optimisation
  rowstart[OpenMP::getThreadNum()]=current;
  if( targ->getRank()==2 ) {
    rowstart[OpenMP::getThreadNum()] = current*targ->getShape()[1];
  }
  myminimiser.minimise( cgtol, point, &ProjectPoints::calculateStress );
}

void ProjectPoints::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> point( dimout );
  getProjection( current, point );
  for(unsigned j=0; j<dimout; ++j) {
    myvals.setValue( getConstPntrToComponent(j)->getPositionInStream(), point[j] );
  }
}

void ProjectPoints::calculate() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    std::vector<double> point( dimout );
    getProjection( 0, point );
    for(unsigned i=0; i<dimout; ++i) {
      getPntrToComponent(i)->set(point[i]);
    }
  } else {
    runAllTasks();
  }
}

}
}

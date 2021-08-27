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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/ConjugateGradient.h"
#include "tools/SwitchingFunction.h"
#include "tools/OpenMP.h"
#include "tools/Random.h" 

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED PROJECT_POINTS
/*

*/
//+ENDPLUMEDOC

class ProjectPoints : 
public ActionWithValue,
public ActionWithArguments {
private: 
  double cgtol;
  bool updateWasRun;
  mutable std::vector<unsigned> rowstart;
  std::vector<SwitchingFunction> switchingFunction;
  ConjugateGradient<ProjectPoints> myminimiser;
public:
  static void registerKeywords( Keywords& keys );
  ProjectPoints( const ActionOptions& );
  unsigned getNumberOfDerivatives() const { return 0; }
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
  void calculate() {}
  void apply() {}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(ProjectPoints,"PROJECT_POINTS")

void ProjectPoints::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); 
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG"); 
  keys.add("numbered","TARGET","the matrix of target quantities that you would like to match");
  keys.add("numbered","FUNC","a function that is applied on the distances between the points in the low dimensional space");
  keys.add("numbered","WEIGHTS","the matrix with the weights of the target quantities");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
  keys.addOutputComponent("coord","default","the coordinates of the points in the low dimensional space");
}


ProjectPoints::ProjectPoints( const ActionOptions& ao ) :
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  updateWasRun(false),
  rowstart(OpenMP::getNumThreads()),
  myminimiser( this )
{
  unsigned nvals=0; for(unsigned i=arg_ends[0];i<arg_ends[1];++i) nvals += getPntrToArgument(i)->getNumberOfValues();
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      unsigned tvals=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tvals += getPntrToArgument(j)->getNumberOfValues();
      if( nvals!=tvals ) error("mismatch between numbers of projections");
  }
  std::vector<Value*> args( getArguments() ), target, weights; std::string sfd, errors; unsigned ntoproj=0;
  // Read in target "distances" and target weights
  for(unsigned i=1;;++i) {
      target.resize(0); if( !parseArgumentList("TARGET",i,target) ) break; 
      std::string inum; Tools::convert( i, inum ); 
      if( target.size()!=1 ) error("should only be one value in input to TARGET" + inum );
      if( (target[0]->getRank()!=1 && target[0]->getRank()!=2) || target[0]->hasDerivatives() ) error("input to TARGET" + inum + " keyword should be a vector/matrix");
      if( target[0]->getShape()[0]!=nvals ) error("number of rows in target matrix should match number of input coordinates");  
      if( i==1 && target[0]->getRank()==1 ) { ntoproj = 1; }
      else if( ntoproj==1 && target[0]->getRank()!=1 ) error("mismatch between numbers of target distances"); 
      else if( i==1 ) ntoproj = target[0]->getShape()[1]; 
      else if( ntoproj!=target[0]->getShape()[1] ) error("mismatch between numbers of target distances");
      if( !parseArgumentList("WEIGHTS",i,weights) ) error("missing WEIGHTS" + inum + " keyword in input");
      if( weights.size()!=1 ) error("should only be one value in input to WEIGHTS" + inum );
      if( weights[0]->getRank()!=1 || weights[0]->hasDerivatives() ) error("input to WEIGHTS" + inum + " keyword should be a vector");
      if( weights[0]->getShape()[0]!=nvals ) error("number of weights should match number of input coordinates");
      args.push_back( target[0] ); args.push_back( weights[0] );
      bool has_sf = parseNumbered("FUNC",i,sfd); switchingFunction.push_back( SwitchingFunction() );
      if( !has_sf ) {  
         switchingFunction[i-1].set( "CUSTOM FUNC=1-sqrt(x2) R_0=1.0", errors );
      } else {
         switchingFunction[i-1].set( sfd, errors );
         if( errors.length()!=0 ) error("problem reading switching function description " + errors); 
      }
      log.printf("  %sth term seeks to match tranformed distances with those in matrix %s \n", inum.c_str(), target[0]->getName().c_str() );
      log.printf("  in %sth term distances are transformed by 1-switching function with r_0=%s \n", inum.c_str(), switchingFunction[i-1].description().c_str() );
      log.printf("  in %sth term weights of matrix elements in stress function are given by %s \n", inum.c_str(), weights[0]->getName().c_str() );
  }
  std::vector<unsigned> shape(1); shape[0]=ntoproj; 
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      std::string num; Tools::convert( i+1, num ); addComponent( "coord-" + num, shape ); 
      componentIsNotPeriodic( "coord-" + num ); getPntrToOutput( i )->alwaysStoreValues();
      if( getPntrToArgument(arg_ends[i])->isTimeSeries() ) getPntrToOutput(i)->makeTimeSeries(); 
  }
  // Create a list of tasks to perform
  for(unsigned i=0;i<ntoproj;++i) addTaskToList(i);
  parse("CGTOL",cgtol); log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  requestArguments( args, false ); checkRead();
}

double ProjectPoints::calculateStress( const std::vector<double>& pp, std::vector<double>& der ) {
  unsigned nmatrices = ( getNumberOfArguments() - arg_ends[arg_ends.size()-1] ) / 2; double stress=0;

  unsigned t=OpenMP::getThreadNum();
  std::vector<double> dtmp( pp.size() ); unsigned nv = 1, nland = getPntrToArgument( arg_ends[arg_ends.size()-1] )->getShape()[0];
  if( getPntrToArgument( arg_ends[arg_ends.size()-1] )->getRank()==2 ) nv = getPntrToArgument( arg_ends[arg_ends.size()-1] )->getShape()[1]; 
  for(unsigned i=0;i<nland;++i) {
      // Calculate distance in low dimensional space
      double dd2 = 0; for(unsigned k=0; k<pp.size(); ++k) { dtmp[k] = pp[k] - retrieveRequiredArgument( k, i ); dd2 += dtmp[k]*dtmp[k]; }

      for(unsigned k=0; k<nmatrices; ++k ) {
          // Now do transformations and calculate differences
          double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
          // Get the weight for this connection 
          double weight = getPntrToArgument( arg_ends[arg_ends.size()-1] + 2*k + 1 )->get( i );
          // Get the difference for the connection
          double fdiff = fd - getPntrToArgument( arg_ends[arg_ends.size()-1] + 2*k )->get( rowstart[t]+i*nv );
          // Calculate derivatives
          double pref = -2.*weight*fdiff*df; for(unsigned n=0; n<pp.size(); ++n) der[n]+=pref*dtmp[n]; 
          // Accumulate the total stress
          stress += weight*fdiff*fdiff;
      }
  }
  return stress;
}

void ProjectPoints::performTask( const unsigned& current, MultiValue& myvals ) const {
  Value* targ = getPntrToArgument( arg_ends[arg_ends.size()-1] ); 
  unsigned nland = targ->getShape()[0], nv=1; 
  if( targ->getRank()==2 ) nv = targ->getShape()[1];
  unsigned closest=0; double mindist = targ->get( current ); 
  for(unsigned i=1; i<nland; ++i) {
    double dist = targ->get( current + i*nv );
    if( dist<mindist ) { mindist=dist; closest=i; }
  }
  // Put the initial guess near to the closest landmark  -- may wish to use grid here again Sandip??
  Random random; random.setSeed(-1234); unsigned nlow = arg_ends.size()-1; std::vector<double> point( nlow );
  for(unsigned j=0; j<nlow; ++j) point[j] = retrieveRequiredArgument( j, closest ) + (random.RandU01() - 0.5)*0.01;
  // And do the optimisation
  rowstart[OpenMP::getThreadNum()] = current; myminimiser.minimise( cgtol, point, &ProjectPoints::calculateStress );
  for(unsigned j=0;j<nlow;++j) myvals.setValue( getPntrToOutput(j)->getPositionInStream(), point[j] );
}

void ProjectPoints::update() {
  updateWasRun=true; runAllTasks();
}

void ProjectPoints::runFinalJobs() {
  if( updateWasRun ) return ;
  // Resize all the output stuff
  std::vector<unsigned> shape(1); shape[0]=1;
  if( getPntrToArgument( arg_ends[arg_ends.size()-1] )->getRank()==2 ) shape[0]=getPntrToArgument( arg_ends[arg_ends.size()-1] )->getShape()[1]; 
  for(unsigned i=0;i<arg_ends.size()-1;++i) getPntrToOutput(i)->setShape( shape ); 
  for(unsigned i=0;i<shape[0];++i) addTaskToList(i);
  runAllTasks();
}

}
}

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

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED ARRANGE_POINTS
/*

*/
//+ENDPLUMEDOC

class ArrangePoints : 
public ActionWithValue,
public ActionWithArguments {
private: 
  double cgtol;
  bool updateWasRun;
  std::vector<SwitchingFunction> switchingFunction;
  void checkInputMatrix( const std::string& key, const unsigned& nvals, const std::vector<Value*>& mat ) const ;
protected:
  double calculateFullStress( const std::vector<double>& p, std::vector<double>& d );
public:
  static void registerKeywords( Keywords& keys );
  ArrangePoints( const ActionOptions& );
  unsigned getNumberOfDerivatives() const { return 0; }
  void calculate() {}
  virtual void optimize( std::vector<double>& pos );
  void apply() {}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(ArrangePoints,"ARRANGE_POINTS")

void ArrangePoints::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); 
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("numbered","TARGET","the matrix of target quantities that you would like to match");
  keys.add("numbered","FUNC","a function that is applied on the distances between the points in the low dimensional space");
  keys.add("numbered","WEIGHTS","the matrix with the weights of the target quantities");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
  keys.addOutputComponent("coord","default","the coordinates of the points in the low dimensional space");
}


ArrangePoints::ArrangePoints( const ActionOptions& ao ) :
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  updateWasRun(false)
{
  std::vector<unsigned> shape(1); shape[0]=0; 
  for(unsigned i=arg_ends[0];i<arg_ends[1];++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      unsigned tvals=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
      if( shape[0]!=tvals ) error("mismatch between sizes of input coordinates");
      std::string num; Tools::convert( i+1, num ); addComponent( "coord-" + num, shape ); 
      componentIsNotPeriodic( "coord-" + num ); getPntrToOutput( i )->alwaysStoreValues();
      if( getPntrToArgument(arg_ends[i])->isTimeSeries() ) getPntrToOutput(i)->makeTimeSeries();
  }
  std::vector<Value*> args( getArguments() ), target, weights; std::string sfd, errors;
  // Read in target "distances" and target weights
  for(unsigned i=1;;++i) {
      target.resize(0); if( !parseArgumentList("TARGET",i,target) ) break; 
      std::string inum; Tools::convert( i, inum ); checkInputMatrix( "TARGET" + inum, shape[0], target );
      if( !parseArgumentList("WEIGHTS",i,weights) ) error("missing WEIGHTS" + inum + " keyword in input");
      checkInputMatrix( "WEIGHTS" + inum, shape[0], target ); 
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
  parse("CGTOL",cgtol); log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  requestArguments( args, false ); checkRead();
}

void ArrangePoints::checkInputMatrix( const std::string& key, const unsigned& nvals, const std::vector<Value*>& mat ) const {
  if( mat.size()!=1 ) error("should only be one value in input to " + key );
  if( mat[0]->getRank()!=2 || mat[0]->hasDerivatives() ) error("input to " + key + " keyword should be a matrix");
  if( mat[0]->getShape()[0]!=nvals || mat[0]->getShape()[1]!=nvals ) error("input to " + key + " keyword has the wrong size");
}

double ArrangePoints::calculateFullStress( const std::vector<double>& p, std::vector<double>& d ) {
  // Zero derivative and stress accumulators
  for(unsigned i=0; i<p.size(); ++i) d[i]=0.0; 
  double stress=0; std::vector<double> dtmp( arg_ends.size()-1 );

  unsigned nlow = arg_ends.size()-1;
  unsigned nmatrices = ( getNumberOfArguments() - arg_ends[arg_ends.size()-1] ) / 2;
  std::vector<unsigned> shape( getPntrToArgument( arg_ends[arg_ends.size()-1] )->getShape() );
  for(unsigned i=1; i<shape[0]; ++i) {
    for(unsigned j=0; j<i; ++j) {
      // Calculate distance in low dimensional space
      double dd2=0; for(unsigned k=0; k<nlow; ++k) { dtmp[k]=p[nlow*i+k] - p[nlow*j+k]; dd2+=dtmp[k]*dtmp[k]; }

      for(unsigned k=0; k<nmatrices; ++k ) {
          // Now do transformations and calculate differences
          double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
          // Get the weight for this connection 
          double weight = getPntrToArgument( arg_ends[arg_ends.size()-1] + 2*k + 1 )->get( shape[0]*i+j ); 
          // Get the difference for the connection
          double fdiff = fd - getPntrToArgument( arg_ends[arg_ends.size()-1] + 2*k )->get( shape[0]*i+j );
          // Calculate derivatives
          double pref = -2.*weight*fdiff*df;
          for(unsigned n=0; n<nlow; ++n) { double dterm=pref*dtmp[n]; d[nlow*i+n]+=dterm; d[nlow*j+n]-=dterm; }
          // Accumulate the total stress
          stress += weight*fdiff*fdiff;
      }
    }
  }
  return stress;
}

void ArrangePoints::optimize( std::vector<double>& pos ) {
  ConjugateGradient<ArrangePoints> mycgminimise( this );
  mycgminimise.minimise( cgtol, pos, &ArrangePoints::calculateFullStress );
}

void ArrangePoints::update() {
  updateWasRun=true;
  // Retrive the initial value
  unsigned nvals = getPntrToArgument( arg_ends[arg_ends.size()-1] )->getShape()[0];
  std::vector<double> pos( (arg_ends.size()-1)*nvals ); unsigned nlow = arg_ends.size()-1;
  for(unsigned i=0;i<nvals;++i) {
      for(unsigned j=0;j<nlow;++j) pos[ nlow*i + j ] = retrieveRequiredArgument( j, i ); 
  }
  // Do the optimization
  optimize( pos );
  // And set the final values
  for(unsigned i=0;i<nvals;++i) {
      for(unsigned j=0;j<nlow;++j) getPntrToOutput(j)->set( i, pos[nlow*i+j] );
  }
}

void ArrangePoints::runFinalJobs() {
   if( updateWasRun ) return ; 
   // Resize all the output stuff
   std::vector<unsigned> shape(1); shape[0]=0; 
   for(unsigned i=arg_ends[0];i<arg_ends[1];++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
   for(unsigned i=0;i<arg_ends.size()-1;++i) { 
       unsigned tvals=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
       if( shape[0]!=tvals ) error("mismatch between sizes of input coordinates"); 
       getPntrToOutput(i)->setShape( shape ); 
   }
   update();
}

}
}

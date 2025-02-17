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
#include "gridtools/GridSearch.h"
#include "SMACOF.h"

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED ARRANGE_POINTS
/*
Arrange points in a low dimensional space so that the (transformed) distances between points in the low dimensional space match the dissimilarities provided in an input matrix.

\par Examples

*/
//+ENDPLUMEDOC

class ArrangePoints :
  public ActionWithValue,
  public ActionWithArguments {
private:
  unsigned dimout, maxiter, ncycles, current_index;
  double cgtol, gbuf;
  std::vector<unsigned> npoints, nfgrid;
  std::vector<double> mypos;
  double smacof_tol, smacof_reg;
  int dist_target;
  enum {conjgrad,pointwise,smacof} mintype;
  std::vector<SwitchingFunction> switchingFunction;
  void checkInputMatrix( const std::string& key, const unsigned& nvals, const std::vector<Value*>& mat ) const ;
  double recalculateSmacofWeights( const std::vector<double>& p, SMACOF& mysmacof ) const ;
protected:
  double calculateStress( const std::vector<double>& p, std::vector<double>& d );
  double calculateFullStress( const std::vector<double>& p, std::vector<double>& d );
public:
  static void registerKeywords( Keywords& keys );
  ArrangePoints( const ActionOptions& );
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void prepare() override ;
  void calculate() override ;
  virtual void optimize( std::vector<double>& pos );
  void apply() override {}
};

PLUMED_REGISTER_ACTION(ArrangePoints,"ARRANGE_POINTS")

void ArrangePoints::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","vector","the initial positions for the projections");
  keys.addInputKeyword("numbered","TARGET","matrix","the matrix of target quantities that you would like to match");
  keys.add("numbered","FUNC","a function that is applied on the distances between the points in the low dimensional space");
  keys.addInputKeyword("numbered","WEIGHTS","matrix","the matrix with the weights of the target quantities");
  keys.add("compulsory","MINTYPE","conjgrad","the method to use for the minimisation");
  keys.add("compulsory","MAXITER","1000","maximum number of optimization cycles for optimisation algorithms");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
  keys.add("compulsory","NCYCLES","5","the number of cycles of global optimization to attempt");
  keys.add("compulsory","BUFFER","1.1","grid extent for search is (max projection - minimum projection) multiplied by this value");
  keys.add("compulsory","CGRID_SIZE","10","number of points to use in each grid direction");
  keys.add("compulsory","FGRID_SIZE","0","interpolate the grid onto this number of points -- only works in 2D");
  keys.add("compulsory","SMACTOL","1E-4","the tolerance for the smacof algorithm");
  keys.add("compulsory","SMACREG","0.001","this is used to ensure that we don't divide by zero when updating weights for SMACOF algorithm");
  keys.addOutputComponent("coord","default","vector","the coordinates of the points in the low dimensional space");
}


ArrangePoints::ArrangePoints( const ActionOptions& ao ) :
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  current_index(0),
  dist_target(-1) {
  dimout = getNumberOfArguments();
  std::vector<unsigned> shape(1);
  shape[0]=getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( shape[0]!=getPntrToArgument(i)->getNumberOfValues() ) {
      error("mismatch between sizes of input coordinates");
    }
    std::string num;
    Tools::convert( i+1, num );
    addComponent( "coord-" + num, shape );
    componentIsNotPeriodic( "coord-" + num );
    getPntrToArgument(i)->buildDataStore();
  }
  std::vector<Value*> args( getArguments() ), target, weights;
  std::string sfd, errors;
  // Read in target "distances" and target weights
  for(unsigned i=1;; ++i) {
    target.resize(0);
    if( !parseArgumentList("TARGET",i,target) ) {
      break;
    }
    std::string inum;
    Tools::convert( i, inum );
    checkInputMatrix( "TARGET" + inum, shape[0], target );
    if( !parseArgumentList("WEIGHTS",i,weights) ) {
      error("missing WEIGHTS" + inum + " keyword in input");
    }
    checkInputMatrix( "WEIGHTS" + inum, shape[0], weights );
    args.push_back( target[0] );
    args.push_back( weights[0] );
    bool has_sf = parseNumbered("FUNC",i,sfd);
    switchingFunction.push_back( SwitchingFunction() );
    if( !has_sf ) {
      switchingFunction[i-1].set( "CUSTOM FUNC=1-sqrt(x2) R_0=1.0", errors );
      dist_target=i-1;
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
  std::string mtype;
  parse("MINTYPE",mtype);
  if( mtype=="conjgrad" ) {
    mintype=conjgrad;
    log.printf("  minimimising stress function using conjugate gradients\n");
  } else if( mtype=="pointwise") {
    mintype=pointwise;
    log.printf("  minimimising stress function using pointwise global optimisation\n");
    npoints.resize(dimout);
    nfgrid.resize(dimout);
    parseVector("CGRID_SIZE",npoints);
    parse("BUFFER",gbuf);
    parse("NCYCLES",ncycles);
    parseVector("FGRID_SIZE",nfgrid);
    if( nfgrid[0]!=0 && dimout!=2 ) {
      error("interpolation only works in two dimensions");
    }
    log.printf("  doing %u cycles of global optimization sweeps\n",ncycles);
    log.printf("  using coarse grid of points that is %u",npoints[0]);
    for(unsigned j=1; j<npoints.size(); ++j) {
      log.printf(" by %u",npoints[j]);
    }
    log.printf("\n  grid is %f times larger than the difference between the position of the minimum and maximum projection \n",gbuf);
    if( nfgrid[0]>0 ) {
      log.printf("  interpolating stress onto grid of points that is %u",nfgrid[0]);
      for(unsigned j=1; j<nfgrid.size(); ++j) {
        log.printf(" by %u",nfgrid[j]);
      }
      log.printf("\n");
    }
  } else if( mtype=="smacof" ) {
    mintype=smacof;
    if( dist_target<0 ) {
      error("one of targets must be distances in order to use smacof");
    }
    log.printf("  minimising stress fucntion using smacof\n");
    parse("SMACTOL",smacof_tol);
    parse("SMACREG",smacof_reg);
    log.printf("  tolerance for smacof algorithms equals %f \n", smacof_tol);
    log.printf("  using %f as regularisation parameter for weights in smacof algorithm\n", smacof_reg);
  } else {
    error("invalid MINTYPE");
  }
  if( mintype!=smacof) {
    parse("CGTOL",cgtol);
    log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  }
  parse("MAXITER",maxiter);
  log.printf("  maximum number of iterations for minimimization algorithms equals %d \n", maxiter );
  requestArguments( args );
  checkRead();
}

void ArrangePoints::checkInputMatrix( const std::string& key, const unsigned& nvals, const std::vector<Value*>& mat ) const {
  mat[0]->buildDataStore();
  if( mat.size()!=1 ) {
    error("should only be one value in input to " + key );
  }
  if( mat[0]->getRank()!=2 || mat[0]->hasDerivatives() ) {
    error("input to " + key + " keyword should be a matrix");
  }
  if( mat[0]->getShape()[0]!=nvals || mat[0]->getShape()[1]!=nvals ) {
    error("input to " + key + " keyword has the wrong size");
  }
}

double ArrangePoints::calculateStress( const std::vector<double>& p, std::vector<double>& d ) {
  double stress=0;
  for(unsigned i=0; i<p.size(); ++i) {
    d[i]=0.0;
  }
  std::vector<double> dtmp(dimout);
  std::vector<unsigned> shape( getPntrToArgument( dimout )->getShape() );
  unsigned targi=shape[0]*current_index;
  unsigned nmatrices = ( getNumberOfArguments() - dimout ) / 2;
  for(unsigned i=0; i<shape[0]; ++i) {
    if( i==current_index ) {
      continue ;
    }
    // Calculate distance in low dimensional space
    double dd2=0;
    for(unsigned k=0; k<dimout; ++k) {
      dtmp[k]=p[k] - mypos[dimout*i+k];
      dd2+=dtmp[k]*dtmp[k];
    }

    for(unsigned k=0; k<nmatrices; ++k ) {
      // Now do transformations and calculate differences
      double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
      // Get the weight for this connection
      double weight = 0;
      for(unsigned j=0; j<shape[0]; ++j) {
        weight += getPntrToArgument( dimout + 2*k + 1 )->get( shape[0]*i+j );
      }
      // Get the difference for the connection
      double fdiff = fd - getPntrToArgument( dimout + 2*k )->get( targi+i );
      // Calculate derivatives
      double pref = -2.*weight*fdiff*df;
      for(unsigned n=0; n<dimout; ++n) {
        d[n] += pref*dtmp[n];
      }
      // Accumulate the total stress
      stress += weight*fdiff*fdiff;
    }
  }
  return stress;
}

double ArrangePoints::calculateFullStress( const std::vector<double>& p, std::vector<double>& d ) {
  // Zero derivative and stress accumulators
  for(unsigned i=0; i<p.size(); ++i) {
    d[i]=0.0;
  }
  double stress=0;
  std::vector<double> dtmp( dimout );

  unsigned nmatrices = ( getNumberOfArguments() - dimout ) / 2;
  std::vector<unsigned> shape( getPntrToArgument( dimout )->getShape() );
  for(unsigned i=1; i<shape[0]; ++i) {
    for(unsigned j=0; j<i; ++j) {
      // Calculate distance in low dimensional space
      double dd2=0;
      for(unsigned k=0; k<dimout; ++k) {
        dtmp[k]=p[dimout*i+k] - p[dimout*j+k];
        dd2+=dtmp[k]*dtmp[k];
      }

      for(unsigned k=0; k<nmatrices; ++k ) {
        // Now do transformations and calculate differences
        double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
        // Get the weight for this connection
        double weight = getPntrToArgument( dimout + 2*k + 1 )->get( shape[0]*i+j );
        // Get the difference for the connection
        double fdiff = fd - getPntrToArgument( dimout + 2*k )->get( shape[0]*i+j );
        // Calculate derivatives
        double pref = -2.*weight*fdiff*df;
        for(unsigned n=0; n<dimout; ++n) {
          double dterm=pref*dtmp[n];
          d[dimout*i+n]+=dterm;
          d[dimout*j+n]-=dterm;
        }
        // Accumulate the total stress
        stress += weight*fdiff*fdiff;
      }
    }
  }
  return stress;
}

double ArrangePoints::recalculateSmacofWeights( const std::vector<double>& p, SMACOF& mysmacof ) const {
  double stress=0, totalWeight=0;
  unsigned nmatrices = ( getNumberOfArguments() - dimout ) / 2;
  std::vector<unsigned> shape( getPntrToArgument( dimout )->getShape() );
  for(unsigned i=1; i<shape[0]; ++i) {
    for(unsigned j=0; j<i; ++j) {
      // Calculate distance in low dimensional space
      double dd2=0;
      for(unsigned k=0; k<dimout; ++k) {
        double dtmp=p[dimout*i+k] - p[dimout*j+k];
        dd2+=dtmp*dtmp;
      }
      // Calculate difference between target difference and true difference
      double wval=0, dd1 = sqrt(dd2);
      double diff = mysmacof.getDistance(i,j) - dd1;

      for(unsigned k=0; k<nmatrices; ++k ) {
        // Don't need to do anything for distances we are matching
        if( k==dist_target ) {
          continue;
        }
        // Now do transformations and calculate differences
        double df, fd = 1. - switchingFunction[k].calculateSqr( dd2, df );
        // Get the weight for this connection
        double weight = getPntrToArgument( dimout + 2*k + 1 )->get( shape[0]*i+j );
        // Get the difference for the connection
        double fdiff = getPntrToArgument( dimout + 2*k )->get( shape[0]*i+j ) - fd;
        // Now set the weight if difference in distance is larger than regularisation parameter
        if( fabs(diff)>smacof_reg  ) {
          wval -= weight*fdiff*df*dd1 / diff;
        }
        // And the total stress and weights
        stress += weight*fdiff*fdiff;
        totalWeight += weight;
      }
      mysmacof.setWeight( j, i, wval );
      mysmacof.setWeight( i, j, wval );
    }
  }
  return stress / totalWeight;
}

void ArrangePoints::optimize( std::vector<double>& pos ) {
  ConjugateGradient<ArrangePoints> mycgminimise( this );
  if( mintype==conjgrad ) {
    mycgminimise.minimise( cgtol, pos, &ArrangePoints::calculateFullStress );
  } else if( mintype==pointwise ) {
    unsigned nvals=getPntrToArgument( dimout )->getShape()[0];
    std::vector<double> gmin( dimout ), gmax( dimout ), mypoint( dimout );
    // Find the extent of the grid
    for(unsigned j=0; j<dimout; ++j) {
      gmin[j]=gmax[j]=pos[j];
    }
    for(unsigned i=1; i<nvals; ++i) {
      for(unsigned j=0; j<dimout; ++j) {
        if( pos[dimout*i+j] < gmin[j] ) {
          gmin[j] = pos[dimout*i+j];
        }
        if( pos[dimout*i+j] > gmax[j] ) {
          gmax[j] = pos[dimout*i+j];
        }
      }
    }
    for(unsigned j=0; j<dimout; ++j) {
      double gbuffer = 0.5*gbuf*( gmax[j]-gmin[j] ) - 0.5*( gmax[j]- gmin[j] );
      gmin[j]-=gbuffer;
      gmax[j]+=gbuffer;
    }
    mypos.resize( pos.size() );
    for(unsigned i=0; i<mypos.size(); ++i) {
      mypos[i] = pos[i];
    }
    gridtools::GridSearch<ArrangePoints> mygridsearch( gmin, gmax, npoints, nfgrid, this );
    // Run multiple loops over all projections
    for(unsigned i=0; i<ncycles; ++i) {
      for(unsigned j=0; j<nvals; ++j) {
        // Setup target distances and target functions for calculate stress
        current_index=j;

        // Find current projection of jth point
        for(unsigned k=0; k<dimout; ++k) {
          mypoint[k]=mypos[j*dimout+k];
        }
        // Minimise using grid search
        bool moved=mygridsearch.minimise( mypoint, &ArrangePoints::calculateStress );
        if( moved ) {
          // Reassign the new projection
          for(unsigned k=0; k<dimout; ++k) {
            mypos[dimout*j+k]=mypoint[k];
          }
          // Minimise output using conjugate gradient
          mycgminimise.minimise( cgtol, mypos, &ArrangePoints::calculateFullStress );
        }
      }
      for(unsigned i=0; i<mypos.size(); ++i) {
        pos[i] = mypos[i];
      }
    }
  } else if( mintype==smacof ) {
    SMACOF mysmacof( getPntrToArgument( dimout + 2*dist_target) );
    double stress = recalculateSmacofWeights( pos, mysmacof );

    for(unsigned i=0; i<maxiter; ++i) {
      // Optimise using smacof and current weights
      mysmacof.optimize( smacof_tol, maxiter, pos );
      // Recalculate weights matrix and sigma
      double newsig = recalculateSmacofWeights( pos, mysmacof );
      // Test whether or not the algorithm has converged
      if( fabs( newsig - stress )<smacof_tol ) {
        break;
      }
      // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
      stress=newsig;
    }
  }
}

void ArrangePoints::prepare() {
  // Make sure all the components are the right size
  std::vector<unsigned> shape(1,getPntrToArgument( dimout )->getShape()[0]);
  for(unsigned j=0; j<dimout; ++j) {
    if( getPntrToComponent(j)->getShape()[0]!=shape[0] ) {
      getPntrToComponent(j)->setShape( shape );
    }
  }
}

void ArrangePoints::calculate() {
  // Retrive the initial value
  unsigned nvals = getPntrToArgument( dimout )->getShape()[0];
  std::vector<double> pos( dimout*nvals );
  for(unsigned i=0; i<nvals; ++i) {
    for(unsigned j=0; j<dimout; ++j) {
      pos[ dimout*i + j ] = getPntrToArgument(j)->get(i);
    }
  }
  // Do the optimization
  optimize( pos );
  // And set the final values
  for(unsigned i=0; i<nvals; ++i) {
    for(unsigned j=0; j<dimout; ++j) {
      getPntrToComponent(j)->set( i, pos[dimout*i+j] );
    }
  }
}

}
}

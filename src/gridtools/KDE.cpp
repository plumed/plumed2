/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "ActionWithGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "core/PbcAction.h"
#include "tools/HistogramBead.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

//+PLUMEDOC ANALYSIS KDE
/*
Create a histogram from the input scalar/vector/matrix using KDE

\par Examples


*/
//+ENDPLUMEDOC

//+PLUMEDOC ANALYSIS SPHERICAL_KDE
/*
Create a histogram from the input scalar/vector/matrix using SPHERICAL_KDE

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class KDE : public ActionWithGrid {
private:
  double hh;
  bool hasheight;
  bool ignore_out_of_bounds, fixed_width;
  double dp2cutoff;
  std::string kerneltype;
  GridCoordinatesObject gridobject;
  std::vector<std::string> gmin, gmax;
  std::vector<double> center;
  std::vector<double> gspacing;
  unsigned num_neigh, bwargno;
  std::vector<Value> grid_diff_value;
  std::vector<unsigned> nbin, nneigh, neighbors;
  unsigned numberOfKernels, nbins;
  SwitchingFunction switchingFunction;
  double von_misses_concentration, von_misses_norm;
  void setupNeighborsVector();
  void retrieveArgumentsAndHeight( const MultiValue& myvals, std::vector<double>& args, double& height ) const ;
  double evaluateKernel( const std::vector<double>& gpoint, const std::vector<double>& args, const double& height, std::vector<double>& der ) const ;
  void setupHistogramBeads( std::vector<HistogramBead>& bead ) const ;
  double evaluateBeadValue( std::vector<HistogramBead>& bead, const std::vector<double>& gpoint, const std::vector<double>& args, const double& height, std::vector<double>& der ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit KDE(const ActionOptions&ao);
  std::vector<std::string> getGridCoordinateNames() const override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  unsigned getNumberOfDerivatives() override;
  void setupOnFirstStep( const bool incalc ) override ;
  void getNumberOfTasks( unsigned& ntasks ) override ;
  void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) override ;
  int checkTaskStatus( const unsigned& taskno, int& flag ) const override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
  void updateForceTasksFromValue( const Value* myval, std::vector<unsigned>& force_tasks ) const override ;
  void gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const override ;
};

PLUMED_REGISTER_ACTION(KDE,"KDE")
PLUMED_REGISTER_ACTION(KDE,"SPHERICAL_KDE")

void KDE::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.use("ARG");
  keys.add("optional","HEIGHTS","this keyword takes the label of an action that calculates a vector of values.  The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.add("optional","VOLUMES","this keyword take the label of an action that calculates a vector of values.  The elements of this vector "
           "divided by the volume of the Gaussian are used as weights for the Gaussians");
  // Keywords for KDE
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("optional","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","METRIC","the inverse covariance to use for the kernels that are added to the grid");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.addFlag("IGNORE_IF_OUT_OF_RANGE",false,"if a kernel is outside of the range of the grid it is safe to ignore");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  // Keywords for spherical KDE
  keys.add("compulsory","CONCENTRATION","the concentration parameter for Von Mises-Fisher distributions (only required for SPHERICAL_KDE)");
  keys.setValueDescription("a function on a grid that was obtained by doing a Kernel Density Estimation using the input arguments");
}

KDE::KDE(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  hasheight(false),
  fixed_width(false) {
  std::vector<unsigned> shape( getNumberOfArguments() );
  center.resize( getNumberOfArguments() );
  numberOfKernels=getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=1; i<shape.size(); ++i) {
    if( numberOfKernels!=getPntrToArgument(i)->getNumberOfValues() ) {
      error("mismatch between numbers of values in input arguments");
    }
  }

  bool weights_are_volumes=true;
  std::vector<std::string> weight_str;
  parseVector("VOLUMES",weight_str);
  if( weight_str.size()==0 ) {
    parseVector("HEIGHTS",weight_str);
    if( weight_str.size()>0 ) {
      weights_are_volumes=false;
    }
  }
  hasheight=(weight_str.size()==1);
  if( weight_str.size()>1 ) {
    error("only one scalar/vector/matrix should be input to HEIGHTS");
  }

  if( getName()=="KDE" ) {
    parse("KERNEL",kerneltype);
    if( kerneltype!="DISCRETE" ) {
      std::string bandwidth;
      std::vector<std::string> bwidths;
      parseVector("BANDWIDTH",bwidths);
      if( bwidths.size()>0 ) {
        std::string band="VALUES=" + bwidths[0];
        for(unsigned i=0; i<bwidths.size(); ++i) {
          if( i>0 ) {
            band += "," + bwidths[i];
          }
        }
        plumed.readInputLine( getLabel() + "_sigma: CONSTANT " + band );
        plumed.readInputLine( getLabel() + "_cov: CUSTOM ARG=" + getLabel() + "_sigma FUNC=x*x PERIODIC=NO" );
        plumed.readInputLine( getLabel() + "_icov: CUSTOM ARG=" + getLabel() + "_cov FUNC=1/x PERIODIC=NO" );
        bandwidth = getLabel() + "_icov";

        if( (kerneltype=="gaussian" || kerneltype=="GAUSSIAN") && weights_are_volumes ) {
          std::string pstr;
          Tools::convert( sqrt(pow(2*pi,bwidths.size())), pstr );
          plumed.readInputLine( getLabel() + "_bwprod: PRODUCT ARG=" + getLabel() + "_cov");
          plumed.readInputLine( getLabel() + "_vol: CUSTOM ARG=" + getLabel() + "_bwprod FUNC=(sqrt(x)*" + pstr + ") PERIODIC=NO");
          if( hasheight ) {
            plumed.readInputLine( getLabel() + "_height: CUSTOM ARG=" + weight_str[0] + "," + getLabel() + "_vol FUNC=x/y PERIODIC=NO");
          } else {
            plumed.readInputLine( getLabel() + "_height: CUSTOM ARG=" + getLabel() + "_vol FUNC=1/x PERIODIC=NO");
          }
          hasheight=true;
          weight_str.resize(1);
          weight_str[0] = getLabel() + "_height";
        }
      } else {
        parse("METRIC",bandwidth);
      }
      weight_str.push_back( bandwidth );
    }
  }
  if( weight_str.size()>0 ) {
    std::vector<Value*> weight_args;
    ActionWithArguments::interpretArgumentList( weight_str, plumed.getActionSet(), this, weight_args );
    std::vector<Value*> args( getArguments() );
    args.push_back( weight_args[0] );
    if( hasheight && weight_args[0]->getNumberOfValues()>1 && numberOfKernels!=weight_args[0]->getNumberOfValues() ) {
      error("mismatch between numbers of values in input arguments and HEIGHTS");
    }

    if( weight_str.size()==2 ) {
      log.printf("  quantities used for weights are : %s \n", weight_str[0].c_str() );
      args.push_back( weight_args[1] );
      if( weight_args[1]->getRank()==1 && weight_args[1]->getNumberOfValues()!=shape.size() ) {
        error("size of bandwidth vector is incorrect");
      }
      if( weight_args[1]->getRank()>2 ) {
        error("bandwidths cannot have rank greater than 2");
      }
      bwargno=args.size()-1;
      log.printf("  bandwidths are taken from : %s \n", weight_str[1].c_str() );
    } else if( !hasheight ) {
      if( weight_args[0]->getRank()==1 && weight_args[0]->getNumberOfValues()!=shape.size() ) {
        error("size of bandwidth vector is incorrect");
      }
      if( weight_args[0]->getRank()>2 ) {
        error("bandwidths cannot have rank greater than 2");
      }
      bwargno=args.size()-1;
      log.printf("  bandwidths are taken from : %s \n", weight_str[0].c_str() );
    } else if ( weight_str.size()==1 ) {
      log.printf("  quantities used for weights are : %s \n", weight_str[0].c_str() );
    } else {
      error("only one scalar/vector/matrix should be input to HEIGHTS");
    }
    requestArguments( args );
  }

  if( getName()=="KDE" ) {
    bool hasauto=false;
    gmin.resize( shape.size() );
    gmax.resize( shape.size() );
    parseVector("GRID_MIN",gmin);
    parseVector("GRID_MAX",gmax);
    for(unsigned i=0; i<gmin.size(); ++i) {
      if( gmin[i]=="auto" ) {
        log.printf("  for %dth coordinate min and max are set from cell directions \n", (i+1) );
        hasauto=true;  // We need to do a preparation step to set the grid from the box size
        if( gmax[i]!="auto" ) {
          error("if gmin is set from box vectors gmax must also be set in the same way");
        }
        if( getPntrToArgument(i)->isPeriodic() ) {
          if( gmin[i]=="auto" ) {
            getPntrToArgument(i)->getDomain( gmin[i], gmax[i] );
          } else {
            std::string str_min, str_max;
            getPntrToArgument(i)->getDomain( str_min, str_max );
            if( str_min!=gmin[i] || str_max!=gmax[i] ) {
              error("all periodic arguments should have the same domain");
            }
          }
        } else if( getPntrToArgument(i)->getName().find(".")!=std::string::npos ) {
          std::size_t dot = getPntrToArgument(i)->getName().find_first_of(".");
          std::string name = getPntrToArgument(i)->getName().substr(dot+1);
          if( name!="x" && name!="y" && name!="z" ) {
            error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
          }
        } else {
          error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
        }
      } else {
        log.printf("  for %dth coordinate min is set to %s and max is set to %s \n", (i+1), gmin[i].c_str(), gmax[i].c_str() );
      }
    }
    if( hasauto && gmin.size()>3 ) {
      error("can only set GRID_MIN and GRID_MAX automatically if components of distance are used in input");
    }

    parseVector("GRID_BIN",nbin);
    parseVector("GRID_SPACING",gspacing);
    parse("CUTOFF",dp2cutoff);
    if( kerneltype.find("bin")==std::string::npos && kerneltype!="DISCRETE" ) {
      std::string errors;
      switchingFunction.set( kerneltype + " R_0=1.0 NOSTRETCH", errors );
      if( errors.length()!=0 ) {
        error("problem reading switching function description " + errors);
      }
    }

    if( nbin.size()!=shape.size() && gspacing.size()!=shape.size() ) {
      error("GRID_BIN or GRID_SPACING must be set");
    }
    // Create a value
    std::vector<bool> ipbc( shape.size() );
    for(unsigned i=0; i<shape.size(); ++i) {
      if( getPntrToArgument( i )->isPeriodic() || gmin[i]=="auto" ) {
        ipbc[i]=true;
      } else {
        ipbc[i]=false;
      }
    }
    gridobject.setup( "flat", ipbc, 0, 0.0 );
  } else {
    if( shape.size()!=3 ) {
      error("should have three coordinates in input to this action");
    }

    parse("GRID_BIN",nbins);
    log.printf("  setting number of bins to %d \n", nbins );
    parse("CONCENTRATION",von_misses_concentration);
    fixed_width=true;
    von_misses_norm = von_misses_concentration / ( 4*pi*sinh( von_misses_concentration ) );
    log.printf("  setting concentration parameter to %f \n", von_misses_concentration );

    // Create a value
    std::vector<bool> ipbc( shape.size(), false );
    double fib_cutoff = std::log( epsilon / von_misses_norm ) / von_misses_concentration;
    gridobject.setup( "fibonacci", ipbc, nbins, fib_cutoff );
    checkRead();

    // Setup the grid
    shape[0]=nbins;
    shape[1]=shape[2]=1;
  }
  parseFlag("IGNORE_IF_OUT_OF_RANGE",ignore_out_of_bounds);
  if( ignore_out_of_bounds ) {
    log.printf("  ignoring kernels that are outside of grid \n");
  }
  addValueWithDerivatives( shape );
  setNotPeriodic();
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  // Make sure we store all the arguments
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    getPntrToArgument(i)->buildDataStore();
  }
  // Check for task reduction
  updateTaskListReductionStatus();
  setupOnFirstStep( false );
}

void KDE::setupOnFirstStep( const bool incalc ) {
  if( getName()=="SPHERICAL_KDE" ) {
    return ;
  }

  for(unsigned i=0; i<getNumberOfDerivatives(); ++i) {
    if( gmin[i]=="auto" && incalc ) {
      double lcoord, ucoord;
      PbcAction* bv = plumed.getActionSet().selectWithLabel<PbcAction*>("Box");
      Tensor box( bv->getPbc().getBox() );
      std::size_t dot = getPntrToArgument(i)->getName().find_first_of(".");
      std::string name = getPntrToArgument(i)->getName().substr(dot+1);
      if( name=="x" ) {
        lcoord=-0.5*box(0,0);
        ucoord=0.5*box(0,0);
      } else if( name=="y" ) {
        lcoord=-0.5*box(1,1);
        ucoord=0.5*box(1,1);
      } else if( name=="z" ) {
        lcoord=-0.5*box(2,2);
        ucoord=0.5*box(2,2);
      } else {
        plumed_error();
      }
      // And convert to strings for bin and bmax
      Tools::convert( lcoord, gmin[i] );
      Tools::convert( ucoord, gmax[i] );
    }
    if( incalc ) {
      grid_diff_value.push_back( Value() );
      if( gridobject.isPeriodic(i) ) {
        grid_diff_value[i].setDomain( gmin[i], gmax[i] );
      } else {
        grid_diff_value[i].setNotPeriodic();
      }
    }
  }
  // And setup the grid object
  gridobject.setBounds( gmin, gmax, nbin, gspacing );
  std::vector<unsigned> shape( gridobject.getNbin(true) );
  getPntrToComponent(0)->setShape( shape );
  bool hasauto=false;
  for(unsigned i=0; i<gmin.size(); ++i) {
    if(gmin[i]=="auto" || gmax[i]=="auto" ) {
      hasauto=true;
      break;
    }
  }
  // And setup the neighbors
  if( !hasauto && kerneltype!="DISCRETE" && getPntrToArgument(bwargno)->isConstant() ) {
    fixed_width=true;
    setupNeighborsVector();
  }
}

void KDE::setupNeighborsVector() {
  if( kerneltype!="DISCRETE" ) {
    std::vector<double> support(gmin.size(),0);
    nneigh.resize( gmin.size() );
    if( kerneltype.find("bin")!=std::string::npos ) {
      std::size_t dd = kerneltype.find("-bin");
      HistogramBead bead;
      bead.setKernelType( kerneltype.substr(0,dd) );
      Value* bw_arg=getPntrToArgument(bwargno);
      if( bw_arg->getRank()<2 ) {
        for(unsigned i=0; i<support.size(); ++i) {
          bead.set( 0, gridobject.getGridSpacing()[i], 1./sqrt(bw_arg->get(i)) );
          support[i] = bead.getCutoff();
          nneigh[i] = static_cast<unsigned>( ceil( support[i]/gridobject.getGridSpacing()[i] ));
        }
      } else {
        plumed_error();
      }
    } else {
      Value* bw_arg=getPntrToArgument(bwargno);
      if( bw_arg->getRank()<2 ) {
        for(unsigned i=0; i<support.size(); ++i) {
          support[i] = sqrt(2.0*dp2cutoff)*(1.0/sqrt(bw_arg->get(i)));
          nneigh[i] = static_cast<unsigned>( ceil( support[i] / gridobject.getGridSpacing()[i] ) );
        }
      } else if( bw_arg->getRank()==2 ) {
        Matrix<double> metric(support.size(),support.size());
        unsigned k=0;
        for(unsigned i=0; i<support.size(); ++i) {
          for(unsigned j=0; j<support.size(); ++j) {
            metric(i,j)=bw_arg->get(k);
            k++;
          }
        }
        Matrix<double> myautovec(support.size(),support.size());
        std::vector<double> myautoval(support.size());
        diagMat(metric,myautoval,myautovec);
        double maxautoval=1/myautoval[0];
        unsigned ind_maxautoval=0;
        for(unsigned i=1; i<support.size(); i++) {
          double neweig=1/myautoval[i];
          if(neweig>maxautoval) {
            maxautoval=neweig;
            ind_maxautoval=i;
          }
        }
        for(unsigned i=0; i<support.size(); i++) {
          support[i] = sqrt(2.0*dp2cutoff)*fabs(sqrt(maxautoval)*myautovec(i,ind_maxautoval));
          nneigh[i] = static_cast<unsigned>( ceil( support[i] / gridobject.getGridSpacing()[i] ) );
        }
      } else {
        plumed_error();
      }
    }
    for(unsigned i=0; i<gridobject.getDimension(); ++i) {
      double fmax, fmin;
      Tools::convert( gridobject.getMin()[i], fmin );
      Tools::convert( gridobject.getMax()[i], fmax );
      if( gridobject.isPeriodic(i) && 2*support[i]>(fmax-fmin) ) {
        error("bandwidth is too large for periodic grid");
      }
    }
  }
}

unsigned KDE::getNumberOfDerivatives() {
  return gridobject.getDimension();
}

std::vector<std::string> KDE::getGridCoordinateNames() const {
  std::vector<std::string> names( gridobject.getDimension() );
  for(unsigned i=0; i<names.size(); ++i) {
    names[i] = getPntrToArgument(i)->getName();
  }
  return names;
}

const GridCoordinatesObject& KDE::getGridCoordinatesObject() const {
  return gridobject;
}

void KDE::areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {
  if( numberOfKernels==1 || (hasheight && getPntrToArgument(gridobject.getDimension())->getRank()>0) ) {
    task_reducing_actions.push_back(this);
  }
}

void KDE::getNumberOfTasks( unsigned& ntasks ) {
  if( !fixed_width ) {
    setupNeighborsVector();
  }
  ntasks = numberOfKernels = getPntrToArgument(0)->getNumberOfValues();
  if( numberOfKernels>1 ) {
    return;
  }

  hh = 1.0;
  if( hasheight ) {
    hh = getPntrToArgument(gridobject.getDimension())->get();
  }
  for(unsigned i=0; i<center.size(); ++i) {
    center[i]=getPntrToArgument(i)->get();
  }
  if( !ignore_out_of_bounds && !gridobject.inbounds( center ) ) {
    //if( fabs(height)>epsilon ) warning("bounds are possibly set too small as hills with substantial heights are being ignored");
    return;
  }
  if( kerneltype=="DISCRETE" ) {
    num_neigh=1;
    neighbors.resize(1);
    for(unsigned i=0; i<center.size(); ++i) {
      center[i] += 0.5*gridobject.getGridSpacing()[i];
    }
    neighbors[0]=gridobject.getIndex( center );
  } else {
    gridobject.getNeighbors( center, nneigh, num_neigh, neighbors );
  }
  ntasks = getPntrToComponent(0)->getNumberOfValues();
  return;
}

int KDE::checkTaskStatus( const unsigned& taskno, int& flag ) const {
  if( numberOfKernels>1 ) {
    if( hasheight && getPntrToArgument(gridobject.getDimension())->getRank()>0
        && fabs(getPntrToArgument(gridobject.getDimension())->get(taskno))<epsilon ) {
      return 0;
    }
    return 1;
  }
  for(unsigned i=0; i<num_neigh; ++i) {
    if( taskno==neighbors[i] ) {
      return 1;
    }
  }
  return 0;
}

void KDE::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( numberOfKernels==1 ) {
    double newval;
    std::vector<double> args( gridobject.getDimension() ), der( gridobject.getDimension() );
    unsigned valout = getConstPntrToComponent(0)->getPositionInStream();
    gridobject.getGridPointCoordinates( current, args );
    if( getName()=="KDE" ) {
      if( kerneltype=="DISCRETE" ) {
        newval = 1.0;
      } else if( kerneltype.find("bin")!=std::string::npos ) {
        double val=hh;
        std::size_t dd = kerneltype.find("-bin");
        HistogramBead bead;
        bead.setKernelType( kerneltype.substr(0,dd) );
        Value* bw_arg=getPntrToArgument(bwargno);
        for(unsigned j=0; j<args.size(); ++j) {
          if( gridobject.isPeriodic(j) ) {
            double lcoord,  ucoord;
            Tools::convert( gmin[j], lcoord );
            Tools::convert( gmax[j], ucoord );
            bead.isPeriodic( lcoord, ucoord );
          } else {
            bead.isNotPeriodic();
          }
          if( bw_arg->getRank()<2 ) {
            bead.set( args[j], args[j]+gridobject.getGridSpacing()[j], 1/sqrt(bw_arg->get(j)) );
          } else if( bw_arg->getRank()==2 ) {
            plumed_error();
          }
          double contr = bead.calculateWithCutoff( args[j], der[j] );
          val = val*contr;
          der[j] = der[j] / contr;
        }
        for(unsigned j=0; j<args.size(); ++j) {
          der[j] *= val;
        }
        newval=val;
      } else {
        newval = evaluateKernel( args, center, hh, der );
      }
    } else {
      double dot=0;
      for(unsigned i=0; i<der.size(); ++i) {
        dot += args[i]*center[i];
      }
      newval = hh*von_misses_norm*exp( von_misses_concentration*dot );
      for(unsigned i=0; i<der.size(); ++i) {
        der[i] = von_misses_concentration*newval*args[i];
      }
    }
    myvals.setValue( valout, newval );
    for(unsigned i=0; i<der.size(); ++i) {
      myvals.addDerivative( valout, i, der[i] );
      myvals.updateIndex( valout, i );
    }
  }
}

void KDE::retrieveArgumentsAndHeight( const MultiValue& myvals, std::vector<double>& args, double& height ) const {
  height=1.0;
  for(unsigned i=0; i<args.size(); ++i) {
    args[i]=getPntrToArgument(i)->get( myvals.getTaskIndex() );
  }
  if( hasheight && getPntrToArgument(args.size())->getRank()==0 ) {
    height = getPntrToArgument( args.size() )->get();
  } else if( hasheight ) {
    height = getPntrToArgument( args.size() )->get( myvals.getTaskIndex() );
  }
}

double KDE::evaluateKernel( const std::vector<double>& gpoint, const std::vector<double>& args, const double& height, std::vector<double>& der ) const {
  double r2=0, hval = height;
  Value* bw_arg=getPntrToArgument(bwargno);
  if( bw_arg->getRank()<2 ) {
    for(unsigned j=0; j<der.size(); ++j) {
      double tmp = -grid_diff_value[j].difference( gpoint[j], args[j] );
      der[j] = tmp*bw_arg->get(j);
      r2 += tmp*der[j];
    }
  } else if( bw_arg->getRank()==2 ) {
    for(unsigned j=0; j<der.size(); ++j) {
      der[j]=0;
      double dp_j, dp_k;
      dp_j = -grid_diff_value[j].difference( gpoint[j], args[j] );
      for(unsigned k=0; k<der.size(); ++k ) {
        if(j==k) {
          dp_k = dp_j;
        } else {
          dp_k = -grid_diff_value[k].difference( gpoint[k], args[k] );
        }
        der[j] += bw_arg->get(j*der.size()+k)*dp_k;
        r2 += dp_j*dp_k*bw_arg->get(j*der.size()+k);
      }
    }
  } else {
    plumed_error();
  }
  double dval, val=hval*switchingFunction.calculateSqr( r2, dval );
  dval *= hval;
  for(unsigned j=0; j<der.size(); ++j) {
    der[j] *= dval;
  }
  return val;
}

void KDE::setupHistogramBeads( std::vector<HistogramBead>& bead ) const {
  std::size_t dd = kerneltype.find("-bin");
  std::string ktype=kerneltype.substr(0,dd);
  for(unsigned j=0; j<bead.size(); ++j) {
    bead[j].setKernelType( ktype );
    if( gridobject.isPeriodic(j) ) {
      double lcoord,  ucoord;
      Tools::convert( gmin[j], lcoord );
      Tools::convert( gmax[j], ucoord );
      bead[j].isPeriodic( lcoord, ucoord );
    } else {
      bead[j].isNotPeriodic();
    }
  }
}

double KDE::evaluateBeadValue( std::vector<HistogramBead>& bead, const std::vector<double>& gpoint, const std::vector<double>& args,
                               const double& height, std::vector<double>& der ) const {
  double val=height;
  std::vector<double> contr( args.size() );
  Value* bw_arg=getPntrToArgument(bwargno);
  if( bw_arg->getRank()<2 ) {
    for(unsigned j=0; j<args.size(); ++j) {
      bead[j].set( gpoint[j], gpoint[j]+gridobject.getGridSpacing()[j], 1/sqrt(bw_arg->get(j)) );
      contr[j] = bead[j].calculateWithCutoff( args[j], der[j] );
      val = val*contr[j];
    }
  } else {
    plumed_error();
  }
  for(unsigned j=0; j<args.size(); ++j) {
    if( fabs(contr[j])>epsilon ) {
      der[j] *= val / contr[j];
    }
  }
  return val;
}

void KDE::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                             const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 );
  if( numberOfKernels==1 ) {
    unsigned istart = bufstart + (1+gridobject.getDimension())*code;
    unsigned valout = getConstPntrToComponent(0)->getPositionInStream();
    buffer[istart] += myvals.get( valout );
    for(unsigned i=0; i<gridobject.getDimension(); ++i) {
      buffer[istart+1+i] += myvals.getDerivative( valout, i );
    }
    return;
  }
  std::vector<double> args( gridobject.getDimension() );
  double height;
  retrieveArgumentsAndHeight( myvals, args, height );
  if( !ignore_out_of_bounds && !gridobject.inbounds( args ) ) {
    // if( fabs(height)>epsilon ) warning("bounds are possibly set too small as hills with substantial heights are being ignored");
    return ;
  }
  // Add the kernel to the grid
  unsigned num_neigh;
  std::vector<unsigned> neighbors;
  if( kerneltype!="DISCRETE" ) {
    gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  }
  std::vector<double> der( args.size() ), gpoint( args.size() );
  if( fabs(height)>epsilon ) {
    if( getName()=="KDE" ) {
      if( kerneltype=="DISCRETE" ) {
        std::vector<double> newargs( args.size() );
        for(unsigned i=0; i<args.size(); ++i) {
          newargs[i] = args[i] + 0.5*gridobject.getGridSpacing()[i];
        }
        plumed_assert( bufstart + gridobject.getIndex( newargs )*(1+args.size())<buffer.size() );
        buffer[ bufstart + gridobject.getIndex( newargs )*(1+args.size()) ] += height;
      } else if( kerneltype.find("bin")!=std::string::npos ) {
        std::vector<HistogramBead> bead( args.size() );
        setupHistogramBeads( bead );
        for(unsigned i=0; i<num_neigh; ++i) {
          gridobject.getGridPointCoordinates( neighbors[i], gpoint );
          double val = evaluateBeadValue( bead, gpoint, args, height, der );
          buffer[ bufstart + neighbors[i]*(1+der.size()) ] += val;
          for(unsigned j=0; j<der.size(); ++j) {
            buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += val*der[j];
          }
        }
      } else {
        for(unsigned i=0; i<num_neigh; ++i) {
          gridobject.getGridPointCoordinates( neighbors[i], gpoint );
          buffer[ bufstart + neighbors[i]*(1+der.size()) ] += evaluateKernel( gpoint, args, height, der );
          for(unsigned j=0; j<der.size(); ++j) {
            buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += der[j];
          }
        }
      }
    } else {
      for(unsigned i=0; i<num_neigh; ++i) {
        gridobject.getGridPointCoordinates( neighbors[i], gpoint );
        double dot=0;
        for(unsigned j=0; j<gpoint.size(); ++j) {
          dot += args[j]*gpoint[j];
        }
        double newval = height*von_misses_norm*exp( von_misses_concentration*dot );
        buffer[ bufstart + neighbors[i]*(1+gpoint.size()) ] += newval;
        for(unsigned j=0; j<gpoint.size(); ++j) {
          buffer[ bufstart + neighbors[i]*(1+gpoint.size()) + 1 + j ] += von_misses_concentration*newval*gpoint[j];
        }
      }
    }
  }
}

void KDE::updateForceTasksFromValue( const Value* myval, std::vector<unsigned>& force_tasks ) const {
  if( !myval->forcesWereAdded() ) {
    return ;
  }
  if( numberOfKernels==1 ) {
    plumed_error();
  }

  int flag=1;
  for(unsigned i=0; i<numberOfKernels; ++i) {
    if( checkTaskStatus( i, flag ) ) {
      force_tasks.push_back(i);
    }
  }
}

void KDE::gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( numberOfKernels==1 ) {
    plumed_error();
    return;
  }
  double height;
  std::vector<double> args( gridobject.getDimension() );
  retrieveArgumentsAndHeight( myvals, args, height );
  unsigned num_neigh;
  std::vector<unsigned> neighbors;
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  std::vector<double> der( args.size() ), gpoint( args.size() );
  unsigned hforce_start = 0;
  for(unsigned j=0; j<der.size(); ++j) {
    hforce_start += getPntrToArgument(j)->getNumberOfStoredValues();
  }
  if( fabs(height)>epsilon ) {
    if( getName()=="KDE" ) {
      if( kerneltype.find("bin")!=std::string::npos ) {
        std::vector<HistogramBead> bead( args.size() );
        setupHistogramBeads( bead );
        for(unsigned i=0; i<num_neigh; ++i) {
          gridobject.getGridPointCoordinates( neighbors[i], gpoint );
          double val = evaluateBeadValue( bead, gpoint, args, height, der );
          double fforce = getConstPntrToComponent(0)->getForce( neighbors[i] );
          if( hasheight && getPntrToArgument(args.size())->getRank()==0 ) {
            forces[ hforce_start ] += val*fforce / height;
          } else if( hasheight ) {
            forces[ hforce_start + getPntrToArgument(args.size())->getIndexInStore(itask) ] += val*fforce / height;
          }
          unsigned n=0;
          for(unsigned j=0; j<der.size(); ++j) {
            forces[n + getPntrToArgument(j)->getIndexInStore(itask)] += der[j]*fforce;
            n += getPntrToArgument(j)->getNumberOfStoredValues();
          }
        }
      } else {
        for(unsigned i=0; i<num_neigh; ++i) {
          gridobject.getGridPointCoordinates( neighbors[i], gpoint );
          double val = evaluateKernel( gpoint, args, height, der ), fforce = getConstPntrToComponent(0)->getForce( neighbors[i] );
          if( hasheight && getPntrToArgument(args.size())->getRank()==0 ) {
            forces[ hforce_start ] += val*fforce / height;
          } else if( hasheight ) {
            forces[ hforce_start + getPntrToArgument(args.size())->getIndexInStore(itask) ] += val*fforce / height;
          }
          unsigned n=0;
          for(unsigned j=0; j<der.size(); ++j) {
            forces[n + getPntrToArgument(j)->getIndexInStore(itask)] += -der[j]*fforce;
            n += getPntrToArgument(j)->getNumberOfStoredValues();
          }
        }
      }
    } else {
      for(unsigned i=0; i<num_neigh; ++i) {
        gridobject.getGridPointCoordinates( neighbors[i], gpoint );
        double dot=0;
        for(unsigned j=0; j<gpoint.size(); ++j) {
          dot += args[j]*gpoint[j];
        }
        double fforce = myval->getForce( neighbors[i] );
        double newval = height*von_misses_norm*exp( von_misses_concentration*dot );
        if( hasheight && getPntrToArgument(args.size())->getRank()==0 ) {
          forces[ hforce_start ] += newval*fforce / height;
        } else if( hasheight ) {
          forces[ hforce_start + getPntrToArgument(args.size())->getIndexInStore(itask) ] += newval*fforce / height;
        }
        unsigned n=0;
        for(unsigned j=0; j<gpoint.size(); ++j) {
          forces[n + getPntrToArgument(j)->getIndexInStore(itask)] += von_misses_concentration*newval*gpoint[j]*fforce;
          n += getPntrToArgument(j)->getNumberOfStoredValues();
        }
      }
    }
  }
}

}
}

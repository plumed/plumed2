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
#include "KDE.h"
#include "tools/OpenMP.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Pbc.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDCALC HISTOGRAM
/*
Accumulate the average probability density along a few CVs from a trajectory.

When using this method it is supposed that you have some collective variable \f$\zeta\f$ that
gives a reasonable description of some physical or chemical phenomenon.  As an example of what we
mean by this suppose you wish to examine the following SN2 reaction:

\f[
 \textrm{OH}^- + \textrm{CH}_3Cl  \rightarrow \textrm{CH}_3OH + \textrm{Cl}^-
\f]

The distance between the chlorine atom and the carbon is an excellent collective variable, \f$\zeta\f$,
in this case because this distance is short for the reactant, \f$\textrm{CH}_3Cl\f$, because the carbon
and chlorine are chemically bonded, and because it is long for the product state when these two atoms are
not chemically bonded.  We thus might want to accumulate the probability density, \f$P(\zeta)\f$, as a function of this distance
as this will provide us with information about the overall likelihood of the reaction.   Furthermore, the
free energy, \f$F(\zeta)\f$, is related to this probability density via:

\f[
F(\zeta) = - k_B T \ln P(\zeta)
\f]

Accumulating these probability densities is precisely what this Action can be used to do.  Furthermore, the conversion
of the histogram to the free energy can be achieved by using the method \ref CONVERT_TO_FES.

We calculate histograms within PLUMED using a method known as kernel density estimation, which you can read more about here:

https://en.wikipedia.org/wiki/Kernel_density_estimation

In PLUMED the value of \f$\zeta\f$ at each discrete instant in time in the trajectory is accumulated.  A kernel, \f$K(\zeta-\zeta(t'),\sigma)\f$,
centered at the current value, \f$\zeta(t)\f$, of this quantity is generated with a bandwith \f$\sigma\f$, which
is set by the user.  These kernels are then used to accumulate the ensemble average for the probability density:

\f[
\langle P(\zeta) \rangle = \frac{ \sum_{t'=0}^t w(t') K(\zeta-\zeta(t'),\sigma) }{ \sum_{t'=0}^t w(t') }
\f]

Here the sums run over a portion of the trajectory specified by the user.  The final quantity evalulated is a weighted
average as the weights, \f$w(t')\f$, allow us to negate the effect any bias might have on the region of phase space
sampled by the system.  This is discussed in the section of the manual on \ref Analysis.

A discrete analogue of kernel density estimation can also be used.  In this analogue the kernels in the above formula
are replaced by dirac delta functions.   When this method is used the final function calculated is no longer a probability
density - it is instead a probability mass function as each element of the function tells you the value of an integral
between two points on your grid rather than the value of a (continuous) function on a grid.

Additional material and examples can be also found in the tutorials \ref belfast-1 and \ref lugano-1.

\par A note on block averaging and errors

Some particularly important
issues related to the convergence of histograms and the estimation of error bars around the ensemble averages you calculate are covered in \ref trieste-2.
The technique for estimating error bars that is known as block averaging is introduced in this tutorial.  The essence of this technique is that
the trajectory is split into a set of blocks and separate ensemble averages are calculated from each separate block of data.  If \f$\{A_i\}\f$ is
the set of \f$N\f$ block averages that are obtained from this technique then the final error bar is calculated as:

\f[
\textrm{error} = \sqrt{ \frac{1}{N} \frac{1}{N-1} \sum_{i=1}^N (A_i^2 - \langle A \rangle )^2 } \qquad \textrm{where} \qquad \langle A \rangle = \frac{1}{N} \sum_{i=1}^N A_i
\f]

If the simulation is biased and reweighting is performed then life is a little more complex as each of the block averages should be calculated as a
weighted average.  Furthermore, the weights should be taken into account when the final ensemble and error bars are calculated.  As such the error should be:

\f[
\textrm{error} = \sqrt{ \frac{1}{N} \frac{\sum_{i=1}^N W_i }{\sum_{i=1}^N W_i - \sum_{i=1}^N W_i^2 / \sum_{i=1}^N W_i} \sum_{i=1}^N W_i (A_i^2 - \langle A \rangle )^2 }
\f]

where \f$W_i\f$ is the sum of all the weights for the \f$i\f$th block of data.

If we wish to caclulate a normalized histogram we must calculate ensemble averages from our biased simulation using:
\f[
 \langle H(x) \rangle = \frac{\sum_{t=1}^M w_t K( x - x_t,\sigma) }{\sum_{t=1}^M w_t}
\f]
where the sums runs over the trajectory, \f$w_t\f$ is the weight of the \f$t\f$th trajectory frame, \f$x_t\f$ is the value of the cv for the \f$t\f$th
trajectory frame and \f$K\f$ is a kernel function centered on \f$x_t\f$ with bandwidth \f$\sigma\f$.  The quantity that is evaluated is the value of the
normalized histogram at point \f$x\f$.  The following ensemble average will be calculated if you use the NORMALIZATION=true option in HISTOGRAM.
If the ensemble average is calculated in this way we must calculate the associated error bars from our block averages using the second of the expressions
above.

A number of works have shown that when biased simulations are performed it is often better to calculate an unormalized estimate of the histogram using:
\f[
\langle H(x) \rangle = \frac{1}{M} \sum_{t=1}^M w_t K( x - x_t,\sigma)
\f]
instead of the expression above.  As such this is what is done by default in HISTOGRAM or if the NORMALIZATION=ndata option is used.
When the histogram is calculated in this second way the first of the two formula above can be used when calculating error bars from
block averages.

\par Examples

The following input monitors two torsional angles during a simulation
and outputs a continuos histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  USE_ALL_DATA
  KERNEL=DISCRETE
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
Notice how the CLEAR keyword is used here and how it is not used in the
previous example.

\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 CLEAR=100000
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  GRID_WFILE=histo
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(KDE,"KDE_CALC")

void KDE::registerKeywords( Keywords& keys ) {
  HistogramBase::registerKeywords( keys );
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","METRIC","the inverse covariance to use for the kernels that are added to the grid");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.addFlag("IGNORE_IF_OUT_OF_RANGE",false,"if a kernel is outside of the range of the grid it is safe to ignore");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

KDE::KDE(const ActionOptions&ao):
  Action(ao),
  HistogramBase(ao),
  firststep(false),
  ignore_out_of_bounds(false),
  gmin( getNumberOfDerivatives() ),
  gmax( getNumberOfDerivatives() )
{
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
  for(unsigned i=0; i<gmin.size(); ++i) {
    if( gmin[i]=="auto" ) {
      log.printf("  for %dth coordinate min and max are set from cell directions \n", (i+1) );
      firststep=true;  // We need to do a preparation step to set the grid from the box size
      if( gmax[i]!="auto" ) error("if gmin is set from box vectors gmax must also be set in the same way");
      if( arg_ends.size()==0 ) {
        if( getPntrToArgument(i)->isPeriodic() ) {
          if( gmin[i]=="auto" ) getPntrToArgument(i)->getDomain( gmin[i], gmax[i] );
          else {
            std::string str_min, str_max; getPntrToArgument(i)->getDomain( str_min, str_max );
            if( str_min!=gmin[i] || str_max!=gmax[i] ) error("all periodic arguments should have the same domain");
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
        for(unsigned j=arg_ends[i]; j<arg_ends[i+1]; ++j) {
          if ( getPntrToArgument(j)->isPeriodic() ) {
            if( gmin[i]=="auto" ) getPntrToArgument(j)->getDomain( gmin[i], gmax[i] );
            else {
              std::string str_min, str_max; getPntrToArgument(j)->getDomain( str_min, str_max );
              if( str_min!=gmin[i] || str_max!=gmax[i] ) error("all periodic arguments should have the same domain");
            }
          } else if( getPntrToArgument(j)->getName().find(".")!=std::string::npos ) {
            std::size_t dot = getPntrToArgument(j)->getName().find_first_of(".");
            std::string name = getPntrToArgument(j)->getName().substr(dot+1);
            if( name!="x" && name!="y" && name!="z" ) {
              error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
            }
          } else {
            error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
          }
        }
      }
    } else {
      log.printf("  for %dth coordinate min is set to %s and max is set to %s \n", (i+1), gmin[i].c_str(), gmax[i].c_str() );
    }
  }
  if( firststep && gmin.size()>3 ) error("can only set GRID_MIN and GRID_MAX automatically if components of distance are used in input");

  parseVector("GRID_BIN",nbin); parseVector("GRID_SPACING",gspacing); parse("CUTOFF",dp2cutoff);
  parse("KERNEL",kerneltype); 
  if( kerneltype!="DISCRETE" ) {
      std::vector<std::string> bandwidth(1); parseVector("METRIC",bandwidth);
      std::vector<Value*> bw_args; interpretArgumentList( bandwidth, bw_args );
      if( bw_args[0]->hasDerivatives() ) error("bandwidth should not have derivatives");
      if( bw_args[0]->getRank()==1 && bw_args[0]->getNumberOfValues( getLabel() )!=getNumberOfDerivatives() ) error("size of bandwidth vector is incorrect");
      if( bw_args[0]->getRank()==2 ) error("not implemented KDE with matrices for bandwidth yet");
      log.printf("  bandwidths are taken from : %s \n", bandwidth[0].c_str() );
      std::vector<Value*> args( getArguments() ); args.push_back( bw_args[0] );
      requestArguments( args, true );
  }
  createTaskList();

  if( kerneltype.find("bin")==std::string::npos && kerneltype!="DISCRETE" ) {
     std::string errors; switchingFunction.set( kerneltype + " R_0=1.0 NOSTRETCH RETURN_DERIV", errors ); 
     if( errors.length()!=0 ) error("problem reading switching function description " + errors);
  }

  if( nbin.size()!=getNumberOfDerivatives() && gspacing.size()!=getNumberOfDerivatives() ) error("GRID_BIN or GRID_SPACING must be set");
  cval.resize( gmin.size() );

  parseFlag("IGNORE_IF_OUT_OF_RANGE",ignore_out_of_bounds);
  if( ignore_out_of_bounds ) log.printf("  ignoring kernels that are outside of grid \n");
  // Create a value
  std::vector<bool> ipbc( getNumberOfDerivatives() );
  for(unsigned i=0; i<getNumberOfDerivatives(); ++i) {
    unsigned k=i; if( arg_ends.size()>0 ) k=arg_ends[i];
    if( getPntrToArgument( k )->isPeriodic() || gmin[i]=="auto" ) ipbc[i]=true;
    else ipbc[i]=false;
  }
  gridobject.setup( "flat", ipbc, 0, 0.0 ); checkRead();

  // Setup the grid if we are not using automatic bounds
  if( !firststep ) {
    gridobject.setBounds( gmin, gmax, nbin, gspacing );
    std::vector<unsigned> shape( gridobject.getNbin(true) );
    for(unsigned i=0; i<gmin.size(); ++i) {
        grid_diff_value.push_back( Value() );
        if( gridobject.isPeriodic(i) ) grid_diff_value[i].setDomain( gmin[i], gmax[i] );
        else grid_diff_value[i].setNotPeriodic();
    }
    addValueWithDerivatives( shape ); setupNeighborsVector();
  } else {
    std::vector<unsigned> shape( getNumberOfDerivatives(), 1 );
    addValueWithDerivatives( shape );
  }
}

void KDE::setupNeighborsVector() {
  if( kerneltype!="DISCRETE" ) {
    std::vector<double> point(gmin.size(), 0), support(gmin.size(),0); nneigh.resize( gmin.size() );
    if( kerneltype.find("bin")!=std::string::npos ) {
      std::size_t dd = kerneltype.find("-bin"); 
      HistogramBead bead; bead.setKernelType( kerneltype.substr(0,dd) );
      Value* bw_arg=getPntrToArgument(arg_ends[arg_ends.size()-1]);
      if( bw_arg->getRank()<2 ) {
          for(unsigned i=0; i<point.size(); ++i) {
            bead.set( 0, gridobject.getGridSpacing()[i], 1./sqrt(bw_arg->get(i)) );
            support[i] = bead.getCutoff(); nneigh[i] = static_cast<unsigned>( ceil( support[i]/gridobject.getGridSpacing()[i] ));
          } 
      } else plumed_error();
    } else {
      Value* bw_arg=getPntrToArgument(arg_ends[arg_ends.size()-1]);
      if( bw_arg->getRank()<2 ) {
          for(unsigned i=0; i<support.size(); ++i) {
             support[i] = sqrt(2.0*dp2cutoff)*(1.0/sqrt(bw_arg->get(i)));
             nneigh[i] = static_cast<unsigned>( ceil( support[i] / gridobject.getGridSpacing()[i] ) );
          }
      } else plumed_error();
    }
    for(unsigned i=0; i<gridobject.getDimension(); ++i) {
      double fmax, fmin; Tools::convert( gridobject.getMin()[i], fmin ); Tools::convert( gridobject.getMax()[i], fmax );
      if( gridobject.isPeriodic(i) && 2*support[i]>(fmax-fmin) ) error("bandwidth is too large for periodic grid");
    }
  }
}

void KDE::completeGridObjectSetup() {
  if( firststep ) {
    for(unsigned i=0; i<getNumberOfDerivatives(); ++i) {
      if( gmin[i]=="auto" ) {
        double lcoord, ucoord; Tensor box( plumed.getAtoms().getPbc().getBox() );
        std::size_t dot = getPntrToArgument(i)->getName().find_first_of(".");
        std::string name = getPntrToArgument(i)->getName().substr(dot+1);
        if( name=="x" ) { lcoord=-0.5*box(0,0); ucoord=0.5*box(0,0); }
        else if( name=="y" ) { lcoord=-0.5*box(1,1); ucoord=0.5*box(1,1); }
        else if( name=="z" ) { lcoord=-0.5*box(2,2); ucoord=0.5*box(2,2); }
        else plumed_error();
        // And convert to strings for bin and bmax
        Tools::convert( lcoord, gmin[i] ); Tools::convert( ucoord, gmax[i] );
      }
      grid_diff_value.push_back( Value() );
      if( gridobject.isPeriodic(i) ) grid_diff_value[i].setDomain( gmin[i], gmax[i] );
      else grid_diff_value[i].setNotPeriodic();
    }
    // And setup the grid object
    gridobject.setBounds( gmin, gmax, nbin, gspacing );
    std::vector<unsigned> shape( gridobject.getNbin(true) );
    getPntrToComponent(0)->setShape( shape );
    // And create the tasks
    if( one_kernel_at_a_time ) {
      for(unsigned i=0; i<gridobject.getNumberOfPoints(); ++i) addTaskToList(i);
    }
    // And setup the neighbors
    setupNeighborsVector();
    // And never do this again
    firststep=false;
  }
}

void KDE::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                      std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                                      std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  bool isdists=dumpcube; double units=1.0; gtype="flat";
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    unsigned k=i; if( arg_ends.size()>0 ) k = arg_ends[i];
    if( getPntrToArgument( k )->getName().find(".")==std::string::npos ) { isdists=false; break; }
    std::size_t dot = getPntrToArgument( k )->getName().find(".");
    std::string name = getPntrToArgument( k )->getName().substr(dot+1);
    if( name!="x" && name!="y" && name!="z" ) { isdists=false; break; }
  }
  if( isdists ) {
    if( plumed.getAtoms().usingNaturalUnits() ) units = 1.0/0.5292;
    else units = plumed.getAtoms().getUnits().getLength()/.05929;
  }
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    unsigned k=i; if( arg_ends.size()>0 ) k = arg_ends[i];
    argn[i] = getPntrToArgument( k )->getName(); double gmin, gmax;
    if( isdists && gridobject.getMin().size()>0 ) {
      Tools::convert( gridobject.getMin()[i], gmin ); Tools::convert( gmin*units, min[i] );
      Tools::convert( gridobject.getMax()[i], gmax ); Tools::convert( gmax*units, max[i] );
    } else if( gridobject.getMin().size()>0 ) { min[i] = gridobject.getMin()[i]; max[i] = gridobject.getMax()[i]; }
    if( nbin.size()>0 ) out_nbin[i]=nbin[i];
    if( gspacing.size()>0 ) spacing[i]=units*gspacing[i];
    pbc[i]=gridobject.isPeriodic(i);
  }
}

void KDE::buildSingleKernel( std::vector<unsigned>& tflags, const double& height, std::vector<double>& args ) {
  if( kerneltype=="DISCRETE" ) {
    for(unsigned i=0; i<args.size(); ++i) args[i] += 0.5*gridobject.getGridSpacing()[i];
    tflags[ gridobject.getIndex( args ) ] = 1; return;
  } else { 
    cheight = height; for(unsigned i=0; i<args.size(); ++i) cval[i] = args[i];
  } 
  unsigned num_neigh; std::vector<unsigned> neighbors;
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0; i<num_neigh; ++i) tflags[ neighbors[i] ] = 1;
}

double KDE::calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const {
  if( kerneltype=="DISCRETE" ) return 1.0;

  if( kerneltype.find("bin")!=std::string::npos ) {
    double val=cheight; std::size_t dd = kerneltype.find("-bin");
    HistogramBead bead; bead.setKernelType( kerneltype.substr(0,dd) ); Value* bw_arg=getPntrToArgument(arg_ends[arg_ends.size()-1]);
    for(unsigned j=0; j<args.size(); ++j) {
      if( gridobject.isPeriodic(j) ) {
        double lcoord,  ucoord; Tools::convert( gmin[j], lcoord );
        Tools::convert( gmax[j], ucoord ); bead.isPeriodic( lcoord, ucoord );
      } else bead.isNotPeriodic();
      if( bw_arg->getRank()<2 ) bead.set( args[j], args[j]+gridobject.getGridSpacing()[j], 1/sqrt(bw_arg->get(j)) );
      else if( bw_arg->getRank()==2 ) plumed_error();
      double contr = bead.calculateWithCutoff( cval[j], der[j] );
      val = val*contr; der[j] = der[j] / contr;
    }
    for(unsigned j=0; j<args.size(); ++j) der[j] *= val; return val;
  } else {
    return evaluateKernel( args, cval, cheight, der );
  }
}

double KDE::evaluateKernel( const std::vector<double>& gpoint, const std::vector<double>& args, const double& height, std::vector<double>& der ) const {
  double r2=0, hval = height; Value* bw_arg=getPntrToArgument(arg_ends[arg_ends.size()-1]);
  if( bw_arg->getRank()<2 ) {
      for(unsigned j=0; j<der.size(); ++j) {
         double tmp = -grid_diff_value[j].difference( gpoint[j], args[j] );
         der[j] = tmp*bw_arg->get(j); r2 += tmp*der[j]; 
      }
  } else plumed_error();
  double dval, val=hval*switchingFunction.calculateSqr( r2, dval ); 
  dval *= hval; for(unsigned j=0; j<der.size(); ++j) der[j] *= dval;
  return val;
}

void KDE::setupHistogramBeads( std::vector<HistogramBead>& bead ) const {
  std::size_t dd = kerneltype.find("-bin");
  std::string ktype=kerneltype.substr(0,dd);
  for(unsigned j=0; j<bead.size(); ++j) {
    bead[j].setKernelType( ktype );
    if( gridobject.isPeriodic(j) ) {
      double lcoord,  ucoord; Tools::convert( gmin[j], lcoord );
      Tools::convert( gmax[j], ucoord ); bead[j].isPeriodic( lcoord, ucoord );
    } else bead[j].isNotPeriodic();
  }
}

double KDE::evaluateBeadValue( std::vector<HistogramBead>& bead, const std::vector<double>& gpoint, const std::vector<double>& args,
                                     const double& height, std::vector<double>& der ) const {
  double val=height; std::vector<double> contr( args.size() ); Value* bw_arg=getPntrToArgument(arg_ends[arg_ends.size()-1]);
  if( bw_arg->getRank()<2 ) {
      for(unsigned j=0; j<args.size(); ++j) {
        bead[j].set( gpoint[j], gpoint[j]+gridobject.getGridSpacing()[j], 1/sqrt(bw_arg->get(j)) );
        contr[j] = bead[j].calculateWithCutoff( args[j], der[j] );
        val = val*contr[j];
      }
  } else plumed_error();
  for(unsigned j=0; j<args.size(); ++j) {
    if( fabs(contr[j])>epsilon ) der[j] *= val / contr[j];
  }
  return val;
}

void KDE::addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( kerneltype=="DISCRETE" ) {
      std::vector<double> newargs( args.size() );
      for(unsigned i=0; i<args.size(); ++i) newargs[i] = args[i] + 0.5*gridobject.getGridSpacing()[i];
      buffer[ bufstart + gridobject.getIndex( newargs )*(1+args.size()) ] += height;
      return;
  }
  if( ignore_out_of_bounds && !gridobject.inbounds( args ) ) return ;

  unsigned num_neigh; std::vector<unsigned> neighbors;
  std::vector<double> gpoint( args.size() ), der( args.size() );
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  if( kerneltype.find("bin")!=std::string::npos ) {
    std::vector<HistogramBead> bead( args.size() ); setupHistogramBeads( bead );
    for(unsigned i=0; i<num_neigh; ++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint );
      double val = evaluateBeadValue( bead, gpoint, args, height, der );
      buffer[ bufstart + neighbors[i]*(1+der.size()) ] += val;
      for(unsigned j=0; j<der.size(); ++j) buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += val*der[j];
    }
  } else {
    for(unsigned i=0; i<num_neigh; ++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint ); 
      buffer[ bufstart + neighbors[i]*(1+der.size()) ] += evaluateKernel( gpoint, args, height, der );               
      for(unsigned j=0; j<der.size(); ++j) buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += der[j]; 
    }
  }
}

void KDE::addKernelForces( const unsigned& heights_index, const unsigned& itask, const std::vector<double>& args, const unsigned& htask, const double& height, std::vector<double>& forces ) const {
  plumed_assert( kerneltype!="DISCRETE" );
  unsigned num_neigh; std::vector<unsigned> neighbors;
  std::vector<double> gpoint( args.size() ), der( args.size() );
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  if( kerneltype.find("bin")!=std::string::npos ) {
    std::vector<HistogramBead> bead( args.size() ); setupHistogramBeads( bead );
    for(unsigned i=0; i<num_neigh; ++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint );
      double val = evaluateBeadValue( bead, gpoint, args, height, der ); double fforce = getPntrToOutput(0)->getForce( neighbors[i] );
      if( heights_index==2 ) forces[ args.size()*numberOfKernels + htask ] += val*fforce / height;
      unsigned n=itask; for(unsigned j=0; j<der.size(); ++j) { forces[n] += der[j]*fforce; n += numberOfKernels; }
    }
  } else {
    for(unsigned i=0; i<num_neigh; ++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint ); 
      double val = evaluateKernel( gpoint, args, height, der ), fforce = getPntrToOutput(0)->getForce( neighbors[i] );
      if( heights_index==2 ) forces[ args.size()*numberOfKernels + htask ] += val*fforce / height;
      unsigned n=itask; for(unsigned j=0; j<der.size(); ++j) { forces[n] += -der[j]*fforce; n += numberOfKernels; }
    }
  }
}

}
}

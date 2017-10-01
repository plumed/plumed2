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
#include "HistogramBase.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"

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

class Histogram : public HistogramBase {
private:
  KernelFunctions* kernel;
  std::string kerneltype;
  std::vector<unsigned> nbin, nneigh;
  std::vector<double> min, max, bandwidth;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit Histogram(const ActionOptions&ao);
  ~Histogram();
  void getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin, std::vector<bool>& pbc ) const ;
  void buildSingleKernel( std::vector<unsigned>& tflags, std::vector<double>& args );
  double calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const ;
  void addKernelToGrid( const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const ;
};

PLUMED_REGISTER_ACTION(Histogram,"KDE")
PLUMED_REGISTER_SHORTCUT(Histogram,"HISTOGRAM")

void Histogram::shortcutKeywords( Keywords& keys ){
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
                                    "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
  keys.add("compulsory","NORMALIZATION","true","This controls how the data is normalized it can be set equal to true, false or ndata.  The differences between "
                                               "these options are explained in the manual page for \\ref HISTOGRAM");
  keys.add("optional","UPDATE_FROM","Only update this action from this time");
  keys.add("optional","UPDATE_UNTIL","Only update this action until this time");
}

void Histogram::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ){
  if( words[0]=="HISTOGRAM" ) {
      // Make the kde object
      std::vector<std::string> hist_words; hist_words.push_back( lab + "_kde:");
      hist_words.push_back("KDE"); for(unsigned i=1;i<words.size();++i) hist_words.push_back( words[i] );
      actions.push_back( hist_words );
      // And the averaging object
      std::vector<std::string> av_words; av_words.push_back( lab + ":");
      av_words.push_back("AVERAGE"); av_words.push_back("ARG=" + lab + "_kde" );
      av_words.push_back("STRIDE=" + keys.find("STRIDE")->second ); 
      av_words.push_back("CLEAR=" + keys.find("CLEAR")->second );
      av_words.push_back("NORMALIZATION=" + keys.find("NORMALIZATION")->second );
      if( keys.count("LOGWEIGHTS") ) av_words.push_back("LOGWEIGHTS=" + keys.find("LOGWEIGHTS")->second );
      if( keys.count("UPDATE_FROM") ) av_words.push_back("UPDATE_FROM=" + keys.find("UPDATE_FROM")->second );
      if( keys.count("UPDATE_UNTIL") ) av_words.push_back("UPDATE_UNTIL=" + keys.find("UPDATE_UNTIL")->second );
      actions.push_back( av_words );
  }
}

void Histogram::registerKeywords( Keywords& keys ) {
  HistogramBase::registerKeywords( keys );
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

Histogram::Histogram(const ActionOptions&ao):
  Action(ao),
  HistogramBase(ao),
  kernel(NULL),
  bandwidth( arg_ends.size()-1 )
{
  std::vector<std::string> gmin( arg_ends.size()-1 ), gmax( arg_ends.size()-1 );
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  parse("KERNEL",kerneltype); if( kerneltype!="DISCRETE" ) parseVector("BANDWIDTH",bandwidth); 
  if( nbin.size()!=(arg_ends.size()-1) && gspacing.size()!=(arg_ends.size()-1) ) error("GRID_BIN or GRID_SPACING must be set");

  // Create a value
  std::vector<bool> ipbc( gmin.size() ); 
  for(unsigned i=0;i<arg_ends.size()-1;++i) ipbc[i] = getPntrToArgument( arg_ends[i] )->isPeriodic();
  gridobject.setup( "flat", ipbc, 0, 0.0 ); gridobject.setBounds( gmin, gmax, nbin, gspacing ); 
  std::vector<unsigned> shape( gridobject.getNbin(true) ); 
  addValueWithDerivatives( shape ); checkRead();

  // Now setup everything for the kernel
  std::vector<double> point(gmin.size(), 0);
  if( kerneltype!="DISCRETE" ) {
      KernelFunctions kernel( point, bandwidth, kerneltype, false, 1.0, true ); 
      nneigh=kernel.getSupport( gridobject.getGridSpacing() ); 
      std::vector<double> support( kernel.getContinuousSupport() );
      for(unsigned i=0; i<gmin.size(); ++i) {
        double fmax, fmin; Tools::convert( gmin[i], fmin ); Tools::convert( gmax[i], fmax );
        if( ipbc[i] && 2*support[i]>(fmax-fmin) ) error("bandwidth is too large for periodic grid");
      }
  } 
}

Histogram::~Histogram() {
  if( kernel ) delete kernel;
}

void Histogram::getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min,
                                      std::vector<std::string>& max, std::vector<unsigned>& nbin, std::vector<bool>& pbc ) const {
  std::vector<unsigned> nn( gridobject.getNbin( false ) );
  for(unsigned i=0;i<getPntrToOutput(0)->getRank();++i) {
      argn[i] = getPntrToArgument( arg_ends[i] )->getName();
      min[i]=gridobject.getMin()[i]; max[i]=gridobject.getMax()[i]; nbin[i]=nn[i]; pbc[i]=gridobject.isPeriodic(i); 
  } 
} 

void Histogram::buildSingleKernel( std::vector<unsigned>& tflags, std::vector<double>& args ) {
  if( kerneltype=="DISCRETE" ) {
      for(unsigned i=0; i<args.size(); ++i) args[i] += 0.5*gridobject.getGridSpacing()[i];
      tflags[ gridobject.getIndex( args ) ] = 1; 
  } else {
      kernel = new KernelFunctions( args, bandwidth, kerneltype, false, 1.0, true );
      unsigned num_neigh; std::vector<unsigned> neighbors;
      gridobject.getNeighbors( args, nneigh, num_neigh, neighbors ); 
      for(unsigned i=0;i<num_neigh;++i) tflags[ neighbors[i] ] = 1;
  }
}

double Histogram::calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const {
  if( kerneltype=="DISCRETE" ) return 1.0;

  std::vector<Value*> vv; std::string str_min, str_max; 
  for(unsigned i=0; i<der.size(); ++i) {
      vv.push_back( new Value() );
      if( getPntrToArgument(i)->isPeriodic() ) {
          getPntrToArgument(i)->getDomain( str_min, str_max );
          vv[i]->setDomain( str_min, str_max );
      } else vv[i]->setNotPeriodic();
      vv[i]->set( args[i] );
  }
  double val = kernel->evaluate( vv, der, true );
  for(unsigned i=0; i<der.size(); ++i) delete vv[i]; 
  return val;
}

void Histogram::addKernelToGrid( const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const {
  std::vector<Value*> vv; std::string str_min, str_max; 
  for(unsigned i=0; i<args.size(); ++i) {
      vv.push_back( new Value() );
      if( getPntrToArgument(i)->isPeriodic() ) {
          getPntrToArgument(i)->getDomain( str_min, str_max );
          vv[i]->setDomain( str_min, str_max );
      } else vv[i]->setNotPeriodic();
  }
  KernelFunctions* kk = new KernelFunctions( args, bandwidth, kerneltype, false, 1.0, true );
  std::vector<double> gpoint( args.size() ), der( args.size() );
  unsigned num_neigh; std::vector<unsigned> neighbors; 
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0;i<num_neigh;++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint );
      for(unsigned j=0; j<der.size(); ++j) vv[j]->set( gpoint[j] );
      buffer[ bufstart + neighbors[i]*(1+der.size()) ] += kk->evaluate( vv, der, true );
      for(unsigned j=0; j<der.size(); ++j) buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += der[j];
  } 
  delete kk; for(unsigned i=0; i<der.size(); ++i) delete vv[i];
}

}
}

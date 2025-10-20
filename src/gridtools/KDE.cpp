/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS KDE
/*
Create a histogram from the input scalar/vector/matrix using KDE

This action can be used to construct instantaneous distributions for quantities by using [kernel density esstimation](https://en.wikipedia.org/wiki/Kernel_density_estimation).
The input arguments must all have the same rank and size but you can use a scalar, vector or matrix in input.  The distribution
of this quantity on a grid is then computed using kernel density estimation.

The following example demonstrates how this action can be used with a scalar as input:

```plumed
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_BIN=100 BANDWIDTH=0.2
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

This input outputs a different file on every time step. These files contain a function stored on a grid.  The function output in this case
consists of a single Gaussian with $\sigma=0.2$ that is centered on the instantaneous value of the distance between atoms 1 and 2.  Obviously,
you are unlikely to use an input like the one above. The more usual thing to do would be to accumulate the histogram over the course of a
few trajectory frames using the [ACCUMULATE](ACCUMULATE.md) command as has been done in the input below, which estimates a histogram as a function
of two collective variables:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
kde: KDE ARG=d1,d2 GRID_MIN=0.0,0.0 GRID_MAX=1.0,1.0 GRID_BIN=100,100 BANDWIDTH=0.2,0.2
histo: ACCUMULATE ARG=kde STRIDE=1
DUMPGRID ARG=histo FILE=histo.grid STRIDE=10000
```

Notice, that you can also achieve something similar by using the [HISTOGRAM](HISTOGRAM.md) shortcut.

## Controlloing the grid

If you prefer to specify the grid spacing rather than the number of bins you can do so using the GRID_SPACING keyword as shown below:

```plumed
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_SPACING=0.01 BANDWIDTH=0.2
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

If $x$ is one of the input arguments to the KDE action and $x<g_{min}$ or $x>g_{max}$, where $g_{min}$ and $g_{max}$ are the minimum
and maximum values on the grid for that argument that were specified using GRID_MIN and GRID_MAX, then by PLUMED will crash.

Notice also that when you use Gaussian kernels to accumulate a denisty as in the input above you need to define a cutoff beyond, which the
Gaussian (which is a function with infinite support) is assumed not to contribute to the accumulated density.  When setting this cutoff you
set the value of $x$ in the following expression $\sigma \sqrt{2*x}$, where $\sigma$ is the bandwidth.  By default $x$ is set equal to 6.25 but
you can change this value by using the CUTOFF keyword as shown below:

```plumed
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_SPACING=0.01 BANDWIDTH=0.2 CUTOFF=6.25
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

## Constructing the density

If you are performing a simulation in the NVT ensemble and wish to look at the density as a function of position in the cell you can use an input like the one shown below:

```plumed
a: FIXEDATOM AT=0,0,0
dens: DISTANCES ATOMS=1-100 ORIGIN=a COMPONENTS
kde: KDE ARG=dens.x,dens.y,dens.z GRID_BIN=100,100,100 BANDWIDTH=0.05,0.05,0.05
DUMPGRID ARG=kde STRIDE=1 FILE=density
```

Notice that you do not need to specify GRID_MIN and GRID_MAX values with this input. In this case PLUMED gets the extent of the grid from the cell vectors during the first
step of the simulation.

## Specifying a non diagonal bandwidth

If for any reason you want to use a bandwidth that is not diagonal when doing kensity density estimation you can do by using an input similar to the one shown below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
kde: KDE ...
  ARG=d1,d2 GRID_MIN=0.0,0.0
  GRID_MAX=1.0,1.0 GRID_BIN=100,100
  BANDWIDTH=0.2,0.1,0.1,0.2 HEIGHTS=1
...
histo: ACCUMULATE ARG=kde STRIDE=1
DUMPGRID ARG=histo FILE=histo.grid STRIDE=10000
```

As there are two arguments for this KDE action the four numbers passed in the bandwdith parameter are interepretted as a $2\times 2$ matrix.
Notice that you can also pass the information for the bandwidth in from another argument as has been done here:

```plumed
m: CONSTANT VALUES=0.2,0.1,0.1,0.2 NROWS=2 NCOLS=2

d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
kde: KDE ...
  ARG=d1,d2 GRID_MIN=0.0,0.0
  GRID_MAX=1.0,1.0 GRID_BIN=100,100
  BANDWIDTH=m HEIGHTS=1
...
histo: ACCUMULATE ARG=kde STRIDE=1
DUMPGRID ARG=histo FILE=histo.grid STRIDE=10000
```

In this case the input is equivalent to the first input above and the bandwidth is a constant.  You could, however, also use a non-constant value as input to the BANDWIDTH keyword.

## Working with vectors and scalars

If the input to your KDE action is a set of scalars it appears odd to separate the process of computing the KDE from the process of accumulating the histogram. However, if
you are using vectors as in the example below, this division can be helpful.

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_BIN=100 BANDWIDTH=0.2
```

In the papea cited in the bibliography below, the [KL_ENTROPY](KL_ENTROPY.md) between the instantaneous distribution of CVs and a reference distribution was introduced
as a collective variable. As is detailed in the documentation for that action, the ability to calculate the instaneous histogram from an input vector is essential to
reproducing these calculations.

Notice that you can also use a one or multiple matrices in the input for a KDE object.  The example below uses the angles between the z axis and set of bonds aroud two
atoms:

```plumed
d1: DISTANCE_MATRIX GROUPA=1,2 GROUPB=3-10 COMPONENTS
phi: CUSTOM ARG=d1.z,d1.w FUNC=acos(x/y) PERIODIC=NO
kde: KDE ARG=phi GRID_MIN=0 GRID_MAX=pi GRID_BIN=200 BANDWIDTH=0.1
```

## Using different weights

In all the inputs above the kernels that are added to the grid on each step are Gaussians with that are normalised so that their integral over all space is one. If you want your
Gaussians to have a particular height you can use the HEIGHT keyword as illustrated below:

```plumed
d1: CONTACT_MATRIX GROUPA=1,2 GROUPB=3-10 SWITCH={RATIONAL R_0=0.1} COMPONENTS
mag: CUSTOM ARG=d1.x,d1.y,d1.z FUNC=x*x+y*y+z*z PERIODIC=NO
phi: CUSTOM ARG=d1.z,mag FUNC=acos(x/sqrt(y)) PERIODIC=NO
kde: KDE ARG=phi GRID_MIN=0 GRID_MAX=pi HEIGHTS=d1.w GRID_BIN=200 BANDWIDTH=0.1
```

As indicated above, the HEIGHTS keyword should be passed a Value that has the same rank and size as the arguments that are passed using the ARG keyword. Each of the Gaussian kernels
that are added to the grid in this case have a value equal to the weight at the maximum of the function.

Notice that you can also use the VOLUMES keyword in a similar way as shown below:

```plumed
d1: CONTACT_MATRIX GROUPA=1,2 GROUPB=3-10 SWITCH={RATIONAL R_0=0.1} COMPONENTS
mag: CUSTOM ARG=d1.x,d1.y,d1.z FUNC=x*x+y*y+z*z PERIODIC=NO
phi: CUSTOM ARG=d1.z,mag FUNC=acos(x/sqrt(y)) PERIODIC=NO
kde: KDE ARG=phi GRID_MIN=0 GRID_MAX=pi VOLUMES=d1.w GRID_BIN=200 BANDWIDTH=0.1
```

Now, however, the integral of the Gaussians over all space are equal to the elements of d1.w.

*/
//+ENDPLUMEDOC

//+PLUMEDOC ANALYSIS SPHERICAL_KDE
/*
Create a histogram from the input scalar/vector/matrix using SPHERICAL_KDE

This action operates similarly to [KDE](KDE.md) but it is designed to be used for investigating [directional statistics]().
It is particularly useful if you are looking at the distribution of bond vectors as illustrated in the input below:

```plumed
# Calculate all the bond vectors
d1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1} COMPONENTS
# Normalise the bond vectors
mag: CUSTOM ARG=d1.x,d1.y,d1.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
d1x: CUSTOM ARG=d1.x,mag FUNC=x/y PERIODIC=NO
d1y: CUSTOM ARG=d1.y,mag FUNC=x/y PERIODIC=NO
d1z: CUSTOM ARG=d1.z,mag FUNC=x/y PERIODIC=NO
# And construct the KDE
kde: SPHERICAL_KDE ARG=d1x,d1y,d1z HEIGHTS=d1.w CONCENTRATION=100 GRID_BIN=144
```

Each bond vector here contributes a [Fisher von-Mises kernel](https://en.wikipedia.org/wiki/Von_Mises–Fisher_distribution) to the spherical grid.  This spherical grid is constructed
using a [Fibonnacci sphere algorithm](https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere) so the number of specified using the GRID_BIN keyword must be a Fibonacci number.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

template <class K, class P>
class KDEGridTools {
public:
  double dp2cutoff;
  std::vector<double> gspacing;
  std::vector<std::size_t> nbin;
  std::vector<std::string> gmin, gmax;
  static void registerKeywords( Keywords& keys );
  static void readHeightKeyword( bool canusevol, std::size_t nargs, const std::vector<std::string>& bw, ActionWithArguments* action );
  static void readBandwidth( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw );
  static void readBandwidthKeyword( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw, std::vector<Value*>& bwargs );
  static void readBandwidthAndHeight( const P& params, ActionWithArguments* action );
  static void convertHeightsToVolumes( const std::size_t& nargs, const std::vector<std::string>& bw, const std::string& volstr, ActionWithArguments* action );
  static void readGridParameters( KDEGridTools<K,P>& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape );
  static void setupGridBounds( KDEGridTools<K,P>& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval );
  static void getDiscreteSupport( const KDEGridTools<K,P>& g, P& p, const K& kp, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject );
  static void getNeighbors( const P& p, K& kp, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors );
};

template <class K, class P>
void KDEGridTools<K,P>::registerKeywords( Keywords& keys ) {
  keys.add("optional","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("optional","VOLUMES","this keyword take the label of an action that calculates a vector of values.  The elements of this vector "
           "divided by the volume of the Gaussian are used as weights for the Gaussians");
  keys.add("optional","HEIGHTS","this keyword takes the label of an action that calculates a vector of values. The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
}

template <class K, class P>
void KDEGridTools<K,P>::readHeightKeyword( bool canusevol, std::size_t nargs, const std::vector<std::string>& bw, ActionWithArguments* action ) {
  std::string weight_str;
  action->parse("HEIGHTS",weight_str);
  std::string str_nvals;
  if( (action->getPntrToArgument(0))->getRank()==2 ) {
    std::string nr, nc;
    Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
    Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
    str_nvals = nr + "," + nc;
  } else {
    Tools::convert( (action->getPntrToArgument(0))->getNumberOfValues(), str_nvals );
  }
  if( weight_str.length()>0 ) {
    KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( weight_str, action, "_heights", true );
  } else if( canusevol ) {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_volumes: ONES SIZE=" + str_nvals ), false );
    KDEGridTools<K,P>::convertHeightsToVolumes(nargs,bw,action->getLabel() + "_volumes",action);
  } else {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: ONES SIZE=" + str_nvals ), false );
    KDEHelper<K,P,KDEGridTools<K,P>>::addArgument( action->getLabel() + "_heights", action );
  }
}

template <>
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::readBandwidthAndHeight( const DiscreteKernel& params, ActionWithArguments* action ) {
  std::size_t nargs = action->getNumberOfArguments();
  KDEGridTools<DiagonalKernelParams,DiscreteKernel>::readHeightKeyword( false, nargs, std::vector<std::string>(), action );
}

template<class K, class P>
void KDEGridTools<K, P>::readBandwidthKeyword( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw, std::vector<Value*>& bwargs ) {
  action->parseVector("BANDWIDTH",bw);
  if( nargs>1 && bw.size()==1 ) {
    ActionWithArguments::interpretArgumentList( bw, action->plumed.getActionSet(), action, bwargs );
    if( bwargs.size()!=1 ) {
      action->error("invalid bandwidth found");
    }
    // Create a vector of ones with the right size
    std::string nvals;
    if( (action->getPntrToArgument(0))->getRank()==2 ) {
      std::string nr, nc;
      Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
      Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
      nvals = nr + "," + nc;
    } else {
      Tools::convert( (action->getPntrToArgument(0))->getNumberOfValues(), nvals );
    }
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_bwones: ONES SIZE=" + nvals ), false );
  }
}

template <class K, class P>
void KDEGridTools<K,P>::readBandwidth( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw ) {
  plumed_assert( typeid(K)==typeid(DiagonalKernelParams) );
  std::vector<Value*> bwargs;
  readBandwidthKeyword( nargs, action, bw, bwargs );
  if( nargs>1 && bw.size()==1 ) {
    if( bwargs[0]->getRank()!=1 || bwargs[0]->getNumberOfValues()!=nargs ) {
      action->error("invalid input for bandwidth parameter");
    }
    std::string str_i;
    bw.resize( nargs );
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_scalar_bw" + str_i + ": SELECT_COMPONENTS ARG=" + bwargs[0]->getName() + " COMPONENTS=" + str_i ), false );
      bw[i] = action->getLabel() + "_bw" + str_i;
      action->plumed.readInputWords( Tools::getWords( bw[i] + ": CUSTOM ARG=" + action->getLabel() + "_bwones," + action->getLabel() + "_scalar_bw" + str_i + " FUNC=x*y PERIODIC=NO"), false );
    }
  } else if( bw.size()==nargs ) {
    double bwval;
    std::string str_i;
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      if( Tools::convertNoexcept( bw[i], bwval ) && fabs(bwval)<epsilon ) {
        KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( bw[i], action, "_bwz" + str_i, true );
      } else {
        KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( bw[i], action, "_bw" + str_i, true );
      }
    }
  } else {
    action->error("wrong number of arguments specified in input to bandwidth parameter");
  }
}

template <>
void KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>::readBandwidthAndHeight( const HistogramBeadKernel& params, ActionWithArguments* action ) {
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidth( nargs, action, bw );
  KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>::readHeightKeyword( false, nargs, bw, action );
}

template <>
void KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>::readBandwidthAndHeight( const RegularKernel<DiagonalKernelParams>& params, ActionWithArguments* action ) {
  std::vector<Value*> bwargs;
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidth( nargs, action, bw );
  std::string volstr;
  action->parse("VOLUMES",volstr);
  if( volstr.length()>0 ) {
    if( !params.canusevol ) {
      action->error("cannot use normalized kernels with selected kernel type");
    }
    // Check if we are using Gaussian kernels
    KDEHelper<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>>::readKernelParameters( volstr, action, "_volumes", false );
    convertHeightsToVolumes(nargs, bw, volstr, action);
  } else {
    KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>::readHeightKeyword( params.canusevol, nargs, bw, action );
  }
}

template <>
void KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>::readBandwidthAndHeight( const RegularKernel<NonDiagonalKernelParams>& params, ActionWithArguments* action ) {
  std::vector<Value*> bwargs;
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidthKeyword( nargs, action, bw, bwargs );
  if( nargs>1 && bw.size()==1 ) {
    if( bwargs[0]->getRank()!=2 || bwargs[0]->getShape()[0]!=nargs || bwargs[0]->getShape()[1]!=nargs  ) {
      action->error("invalid input for bandwidth parameter");
    }
    std::string str_i, str_j;
    bw.resize( nargs*nargs );
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      for(unsigned j=0; j<nargs; ++j) {
        Tools::convert( j+1, str_j );
        action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_scalar_bw" + str_i + "_" + str_j + ": SELECT_COMPONENTS ARG=" + bwargs[0]->getName() + " COMPONENTS=" + str_i + "." + str_j ), false );
        bw[i*nargs+j] = action->getLabel() + "_bw" + str_i + "_" + str_j;
        action->plumed.readInputWords( Tools::getWords( bw[i*nargs+j] + ": CUSTOM ARG=" + action->getLabel() + "_bwones," + action->getLabel() + "_scalar_bw" + str_i + "_" + str_j + " FUNC=x*y PERIODIC=NO"), false );
      }
    }
  } else if( bw.size()==nargs*nargs ) {
    std::string str_i, str_j;
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      for(unsigned j=0; j<nargs; ++j) {
        Tools::convert( j+1, str_j );
        KDEHelper<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>>::readKernelParameters( bw[i*nargs+j], action, "_bw" + str_i + "_" + str_j, true );
      }
    }
  } else {
    action->error("wrong number of arguments specified in input to bandwidth parameter");
  }
  std::string volstr;
  action->parse("VOLUMES",volstr);
  if( volstr.length()>0 ) {
    if( !params.canusevol ) {
      action->error("cannot use normalized kernels with selected kernel type");
    }
    // Check if we are using Gaussian kernels
    KDEHelper<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>>::readKernelParameters( volstr, action, "_volumes", false );
    convertHeightsToVolumes(nargs, bw, volstr, action);
  } else {
    KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>::readHeightKeyword( params.canusevol, nargs, bw, action );
  }
}

template <class K, class P>
void KDEGridTools<K, P>::convertHeightsToVolumes( const std::size_t& nargs, const std::vector<std::string>& bw, const std::string& volstr, ActionWithArguments* action ) {
  if( bw.size()==nargs ) {
    unsigned nonzeroargs=0;
    for(unsigned i=0; i<nargs; ++i) {
      if( bw[i].find("_bwz")==std::string::npos ) {
        nonzeroargs++;
      }
    }
    std::string str_i, nargs_str;
    Tools::convert( nonzeroargs, nargs_str );
    std::string varstr = "VAR=h", funcstr = "FUNC=h/(sqrt((2*pi)^" + nargs_str + ")", argstr = "ARG=" + volstr;
    for(unsigned i=0; i<nargs; ++i) {
      if( bw[i].find("_bwz")!=std::string::npos ) {
        continue;
      }
      Tools::convert( i+1, str_i );
      varstr += ",b" + str_i;
      funcstr += "*b" + str_i;
      argstr += "," + bw[i];
    }
    funcstr += ")";
    if( (action->getPntrToArgument(0))->getNumberOfValues()==1 ) {
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO " + argstr + " " + varstr + " " + funcstr), false );
    } else {
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO " + argstr + " " + varstr + " " + funcstr + " MASK=" + volstr), false );
    }
    KDEHelper<K,P,KDEGridTools<K,P>>::addArgument( action->getLabel() + "_heights", action );
  } else if( bw.size()==nargs*nargs ) {
    action->error("have not implemented normalization parameter for non-diagonal kernels");
  }

}

template <class K, class P>
void KDEGridTools<K,P>::readGridParameters( KDEGridTools<K,P>& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape ) {
  g.gmin.resize( shape.size() );
  g.gmax.resize( shape.size() );
  action->parseVector("GRID_MIN",g.gmin);
  action->parseVector("GRID_MAX",g.gmax);
  for(unsigned i=0; i<g.gmin.size(); ++i) {
    if( g.gmin[i]=="auto" ) {
      action->log.printf("  for %dth coordinate min and max are set automatically \n", (i+1) );
      if( g.gmax[i]!="auto" ) {
        action->error("if gmin is set automatically gmax must also be set automatically");
      }
      Value* myarg = action->getPntrToArgument(i);
      if( myarg->isPeriodic() ) {
        if( g.gmin[i]=="auto" ) {
          myarg->getDomain( g.gmin[i], g.gmax[i] );
        } else {
          std::string str_min, str_max;
          myarg->getDomain( str_min, str_max );
          if( str_min!=g.gmin[i] || str_max!=g.gmax[i] ) {
            action->error("all periodic arguments should have the same domain");
          }
        }
      } else if( myarg->getName().find(".")!=std::string::npos ) {
        std::size_t dot = myarg->getName().find_first_of(".");
        std::string name = myarg->getName().substr(dot+1);
        if( name!="x" && name!="y" && name!="z" ) {
          action->error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
        }
      } else {
        action->error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
      }
    } else {
      action->log.printf("  for %dth coordinate min is set to %s and max is set to %s \n", (i+1), g.gmin[i].c_str(), g.gmax[i].c_str() );
    }
  }

  action->parseVector("GRID_BIN",g.nbin);
  action->parseVector("GRID_SPACING",g.gspacing);
  action->parse("CUTOFF",g.dp2cutoff);

  if( g.nbin.size()!=shape.size() && g.gspacing.size()!=shape.size() ) {
    action->error("GRID_BIN or GRID_SPACING must be set");
  }
  // Create a value
  std::vector<bool> ipbc( shape.size() );
  for(unsigned i=0; i<shape.size(); ++i) {
    if( (action->getPntrToArgument( i ))->isPeriodic() || g.gmin[i]=="auto" ) {
      ipbc[i]=true;
      if( g.nbin.size()==shape.size() ) {
        shape[i] = g.nbin[i];
      }
    } else {
      ipbc[i]=false;
      if( g.nbin.size()==shape.size() ) {
        shape[i] = g.nbin[i]+1;
      }
    }
  }
  gridobject.setup( "flat", ipbc, 0, 0.0 );
}

template <class K, class P>
void KDEGridTools<K,P>::setupGridBounds( KDEGridTools<K,P>& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval ) {
  for(unsigned i=0; i<gridobject.getDimension(); ++i) {
    if( g.gmin[i]=="auto" ) {
      double lcoord, ucoord;
      std::size_t dot = args[i]->getName().find_first_of(".");
      std::string name = args[i]->getName().substr(dot+1);
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
      Tools::convert( lcoord, g.gmin[i] );
      Tools::convert( ucoord, g.gmax[i] );
    }
  }
  // And setup the grid object
  gridobject.setBounds( g.gmin, g.gmax, g.nbin, g.gspacing );
  myval->setShape( gridobject.getNbin(true) );
}

template <class K, class P>
void KDEGridTools<K, P>::getDiscreteSupport( const KDEGridTools<K,P>& g, P& p, const K& kp, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  std::size_t ng = gridobject.getDimension();
  plumed_assert( nneigh.size()==ng );
  std::vector<double> support( ng );
  P::getSupport( p, kp, g.dp2cutoff, support );
  for(unsigned i=0; i<ng; ++i) {
    nneigh[i] = static_cast<unsigned>( ceil( support[i]/gridobject.getGridSpacing()[i] ));
  }
}

template <>
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::getDiscreteSupport( const KDEGridTools<DiagonalKernelParams,DiscreteKernel>& g, DiscreteKernel& p, const DiagonalKernelParams& kp, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  return;
}

template <class K, class P>
void KDEGridTools<K,P>::getNeighbors( const P& p, K& kp, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  gridobject.getNeighbors( kp.at, nneigh, num_neighbors, neighbors );
}

template <>
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::getNeighbors( const DiscreteKernel& p, DiagonalKernelParams& kp, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  num_neighbors=1;
  neighbors.resize(1);
  for(unsigned i=0; i<kp.at.size(); ++i) {
    kp.at[i] += 0.5*gridobject.getGridSpacing()[i];
  }
  neighbors[0]=gridobject.getIndex( kp.at );
}

class SphericalKDEGridTools {
public:
  std::size_t nbins;
  static void registerKeywords( Keywords& keys );
  static void readBandwidthAndHeight( const UniversalVonMisses& params, ActionWithArguments* action );
  static void readGridParameters( SphericalKDEGridTools& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape );
  static void setupGridBounds( SphericalKDEGridTools& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval ) {}
  static void getDiscreteSupport( const SphericalKDEGridTools& g, const UniversalVonMisses& p, const VonMissesKernelParams& kp, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject );
  static void getNeighbors( const UniversalVonMisses& p, const VonMissesKernelParams& kp, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors );
};

void SphericalKDEGridTools::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","CONCENTRATION","the concentration parameter for the Von Mises-Fisher distributions");
  keys.add("compulsory","HEIGHTS","1.0","this keyword takes the label of an action that calculates a vector of values. The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.add("compulsory","GRID_BIN","the number of points on the fibonacci sphere at which the density should be evaluated");
}

void SphericalKDEGridTools::readBandwidthAndHeight( const UniversalVonMisses& params, ActionWithArguments* action ) {
  // Read in the concentration parameters
  std::string von_misses_concentration;
  action->parse("CONCENTRATION",von_misses_concentration);
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::readKernelParameters( von_misses_concentration, action, "_vmconcentration", true );
  action->log.printf("  getting concentration parameters from %s \n", von_misses_concentration.c_str() );
  // Read in the heights
  std::string weight_str;
  action->parse("HEIGHTS",weight_str);
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::readKernelParameters( weight_str, action, "_volumes", false );
  if( (action->getPntrToArgument(0))->getNumberOfValues()==1 ) {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO ARG=" + weight_str + "," + action->getLabel() + "_vmconcentration FUNC=x*y/(4*pi*sinh(y))" ), false );
  } else {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO ARG=" + weight_str + "," + action->getLabel() + "_vmconcentration FUNC=x*y/(4*pi*sinh(y)) MASK=" + weight_str ), false );
  }
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::addArgument( action->getLabel() + "_heights", action );
  action->log.printf("  getting heights from %s \n", weight_str.c_str() );
}

void SphericalKDEGridTools::readGridParameters( SphericalKDEGridTools& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape ) {
  if( shape.size()!=3 ) {
    action->error("should have three coordinates in input to this action");
  }
  action->parse("GRID_BIN",g.nbins);
  action->log.printf("  setting number of bins to %zu \n", g.nbins );
  std::vector<bool> ipbc( 3, false );
  gridobject.setup( "fibonacci", ipbc, g.nbins, 0 );
  shape[0]=g.nbins;
  shape[1]=shape[2]=1;
}

void SphericalKDEGridTools::getDiscreteSupport( const SphericalKDEGridTools& g, const UniversalVonMisses& p, const VonMissesKernelParams& kp, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  plumed_assert( nneigh.size()==gridobject.getDimension() );
  std::vector<bool> ipbc( 3, false );
  double fib_cutoff = std::log( epsilon / (kp.concentration/(4*pi*sinh(kp.concentration))) ) / kp.concentration;
  gridobject.setup( "fibonacci", ipbc, gridobject.getNumberOfPoints(), fib_cutoff );
}

void SphericalKDEGridTools::getNeighbors( const UniversalVonMisses& p, const VonMissesKernelParams& kp, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  gridobject.getNeighbors( kp.at, nneigh, num_neighbors, neighbors );
}

typedef KDE<DiagonalKernelParams,DiscreteKernel,KDEGridTools<DiagonalKernelParams,DiscreteKernel>> discretekde;
PLUMED_REGISTER_ACTION(discretekde,"KDE_DISCRETE")
typedef KDE<DiagonalKernelParams,HistogramBeadKernel,KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>> beadkde;
PLUMED_REGISTER_ACTION(beadkde,"KDE_BEADS")
typedef KDE<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>> flatkde;
PLUMED_REGISTER_ACTION(flatkde,"KDE_KERNELS")
typedef KDE<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools> sphericalkde;
PLUMED_REGISTER_ACTION(sphericalkde,"SPHERICAL_KDE")
typedef KDE<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>> flatfkde;
PLUMED_REGISTER_ACTION(flatfkde,"KDE_FULLCOVAR")


class KDEShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit KDEShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(KDEShortcut,"KDE")

void KDEShortcut::registerKeywords(Keywords& keys) {
  KDE<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>>::registerKeywords( keys );
  keys.addActionNameSuffix("_DISCRETE");
  keys.addActionNameSuffix("_KERNELS");
  keys.addActionNameSuffix("_BEADS");
  keys.addActionNameSuffix("_FULLCOVAR");
}

KDEShortcut::KDEShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string kerneltype;
  parse("KERNEL",kerneltype);
  if( kerneltype=="DISCRETE" ) {
    readInputLine( getShortcutLabel() + ": KDE_DISCRETE " + convertInputLineToString() );
    return;
  }
  std::vector<std::string> args;
  parseVector("ARG", args );
  std::vector<Value*> argvals;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, argvals );
  std::string argstr = " ARG=" + argvals[0]->getName();
  for(unsigned i=1; i<argvals.size(); ++i) {
    argstr += "," + argvals[i]->getName();
  }
  std::vector<std::string> bw;
  parseVector("BANDWIDTH",bw);
  std::string bwstr = " BANDWIDTH=" + bw[0];
  for(unsigned i=1; i<bw.size(); ++i) {
    bwstr += "," + bw[i];
  }
  if( bw.size() == 1 && argvals.size()>1 ) {
    std::vector<Value*> bwargs;
    ActionWithArguments::interpretArgumentList( bw, plumed.getActionSet(), this, bwargs );
    if( bwargs.size()!=1 ) {
      error("invalid input for bandwidth parameter");
    } else if( bwargs[0]->getRank()<=1 ) {
      if( kerneltype.find("bin")==std::string::npos ) {
        readInputLine( getShortcutLabel() + ": KDE_KERNELS " + argstr + " " + bwstr + " KERNEL=" + kerneltype + " " + convertInputLineToString() );
      } else {
        std::size_t dd = kerneltype.find("-bin");
        readInputLine( getShortcutLabel() + ": KDE_BEADS " + argstr + " " + bwstr + " KERNEL=" + kerneltype.substr(0,dd) + " " + convertInputLineToString() );
      }
    } else if( bwargs[0]->getRank()==2 ) {
      readInputLine( getShortcutLabel() + ": KDE_FULLCOVAR" + argstr + " " + bwstr + " KERNEL=" + kerneltype + " " + convertInputLineToString() );
    } else {
      error("found strange rank for bandwidth parameter");
    }
  } else if( bw.size()==argvals.size() ) {
    if( kerneltype.find("bin")==std::string::npos ) {
      readInputLine( getShortcutLabel() + ": KDE_KERNELS " + argstr + " " + bwstr + " KERNEL=" + kerneltype + " " + convertInputLineToString() );
    } else {
      std::size_t dd = kerneltype.find("-bin");
      readInputLine( getShortcutLabel() + ": KDE_BEADS " + argstr + " " + bwstr + " KERNEL=" + kerneltype.substr(0,dd) + " " + convertInputLineToString() );
    }
  } else if( bw.size()==argvals.size()*argvals.size() ) {
    readInputLine( getShortcutLabel() + ": KDE_FULLCOVAR" + argstr + " " + bwstr + " KERNEL=" + kerneltype + " " + convertInputLineToString() );
  } else {
    error("invalid input for bandwidth");
  }
}

}
}

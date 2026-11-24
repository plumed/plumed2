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
#include "KDEGridTools.h"
#include "SphericalKDEGridTools.h"
#include "core/AccelerableShortcut.h"
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
as a collective variable. As is detailed in the documentation for that action, the ability to calculate the instantaneous histogram from an input vector is essential to
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

## Setting the kernel

When you run molecular dynamics and calculate $d$ collective variables using PLUMED you generate a time series of $d$ dimensional vectors $\mathbf{x}_t$. The Kernel Density Estimation
method that is implemented in this action allows you to calculate an estimate of the probabilty density for the $d$ collective variables, which is a function with $d$ arguments.
This estimate is calculated using:

$$
f(\mathbf{x}) = \sum_{t=1}^n h_t K\left[ \sqrt{(\mathbf{x}-\mathbf{x}_t)^T \mathbf{H}^{-1} (\mathbf{x}-\mathbf{x}_t)} \right]
$$

This expression assumes that you are calculating the density from $n$ trajectory frames. $\mathbf{H}$ is bandwidths matrix and $h_t$ is the height that you want to use for the $t$th of your observations.
$K$ is then a non-negative function called the kernel.  This function acts upon the [Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) between $\mathbf{x}$ and $\mathbf{x}_t$.

To provide users with flexibility, the functional form for the kernel is specified using a syntax similar to the one described in the documentation for the [LESS_THAN](LESS_THAN.md)
action. The only difference from the instructions for swtiching functions in [LESS_THAN](LESS_THAN.md) is that PLUMED automatically
sets the $r_0$ parameter of your switching function to one for you. __You should thus not include R_0 when you set the parameters of your kernel function.__

By default PLUMED uses a Gaussian kernel function as illustrated in the input below:

```plumed
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_BIN=100 BANDWIDTH=0.2 KERNEL={GAUSSIAN D_MAX=3.5355}
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

Notice that a stretching and scaling is performed on the gaussian kernel that is employed here so the function is guaranteed to be equal to zero if
$\sqrt{(\mathbf{x}-\mathbf{x}_t)^T \mathbf{H}^{-1} (\mathbf{x}-\mathbf{x}_t)}>3.5355$.  One consequence of this is that the integral of the $d$ dimensional
Gaussian will be slightly less than $(2\pi)^{d/2}\det(H)^{1/2}$. This can be a problem if you are assuming that the integral of your Gaussians is equal to this value when
normalising them.  You can partially avoid this problem by using the NOSTRETCH option as shown below:

```plumed
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_BIN=100 BANDWIDTH=0.2 KERNEL={GAUSSIAN NOSTRETCH D_MAX=3.5355}
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

This, however, introduces a discontinuity in the kernel at $d_{max}$.

__We strongly recommend setting a D_MAX parameter for your kernels__.
If you do not set a D_MAX value for your KERNEL then a cutoff of $d_{max}=\sqrt{2*6.25}\approx 3.5355$ is used. If you are calculating a 1D histogram and D_MAX is unset PLUMED
then you evaluate a kernel at centered at $a$ on the gridpoint at $b$ if $\frac{|a-b|}{s}<=\left\lceil \frac{d_{max}\sigma}{s}\right\rceil$, where $s$ is the spacing between
grid points and $\sigma$ is the bandwidth. The ceiling on the right hand side of this inequality obviously makes this s necessary but not sufficient condition for ensuring that
$\sqrt{(\mathbf{x}-\mathbf{x}_t)^T \mathbf{H}^{-1} (\mathbf{x}-\mathbf{x}_t)}<d_{max}$.  In higher dimensions this becomes even more severe.  In 2D to
determine if we need to evaluate a kernel centered at $(a,b)$ for the grid point $(c,d)$ we first determine if $\frac{|a-c|}{s}<=\left\lceil \frac{d_{max}\sigma}{s}\right\rceil$ and
if $\frac{|b-d|}{s}<=\left\lceil \frac{d_{max}\sigma}{s}\right\rceil$. With these tests there are many more points that satisfy these two criteria but that do not satisfy
$\sqrt{(\mathbf{x}-\mathbf{x}_t)^T \mathbf{H}^{-1} (\mathbf{x}-\mathbf{x}_t)}<d_{max}$. The key thing to recognise is that when D_MAX is set PLUMED sets the value of the kernel equal to zero if
$\sqrt{(\mathbf{x}-\mathbf{x}_t)^T \mathbf{H}^{-1} (\mathbf{x}-\mathbf{x}_t)}>d_{max}$ even when the NOSTRETCH option is employed. If D_MAX is not set this check is not performed and
kernels can thus contribute density to grid points that are more than $\sqrt(2*6.25)$ from the center of the kernel.

Kernel density estimation is a method for estimating a histogram but it is probably not the first method that you learned for estimating histograms. Instead you probably learned to assume that the
the density is constant between grid points. In 1D we thus express the density between grid point $x_k$ and $x_{k+1}$ using:

$$
F(x_{k+1}) - F(x_k) = \int_{x_k}^{x_{k+1}} f(x) \textrm{d}x = \int_{x_k}^{x_{k+1}} \sum_{t=1}^n h_t K\left[ \frac{|x - x_t|}{\sigma} \right] \textrm{d}x
$$

In the typical method for calculating histograms $\sigma=x_{k+1}-x_k$ and $K$ is then:

$$
K(x) = \begin{cases}
1 & \textrm{if} x<1 \\
0 & \textrm{otherwise}
\end{cases}
$$

This is the method used for estimating the histogram when you use KERNEL=DISCRETE. Other (continuous) choices for $K$ are possible and are described in the documentation for [BETWEEN](BETWEEN.md).
For example, in the following input $K$ is is a Gaussian:

```
d1: DISTANCE ATOMS=1,2
kde: KDE ARG=d1 GRID_MIN=0.0 GRID_MAX=1.0 GRID_BIN=100 BANDWIDTH=0.2 KERNEL=gaussian-bin
DUMPGRID ARG=kde STRIDE=1 FILE=kde.grid
```

Multidimensional analogues of this approach can also be used. For example, in a 2D analogue the density in the square with corners at $(x_k,y_j)$, $(x_{k+1},y_j)$, $(x_k,y_{j+1})$ and $(x_{k+1},y_{j+1})$ is computed using:

$$
[ F(x_{k+1}) - F(x_k) ][ F(y_{j+1} - F(y_j) ] = \int_{x_k}^{x_{k+1}} f(x) \textrm{d}x = \sum_{t=1}^n h_t \int_{x_k}^{x_{k+1}} K\left[ \frac{|x - x_t|}{\sigma} \right] \textrm{d}x \int_{y_j}^{y_{j+1}} K\left[ \frac{|y - y_t|}{\sigma} \right] \textrm{d}y
$$

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

typedef KDE<DiagonalKernelParams,DiscreteKernel,KDEGridTools<DiagonalKernelParams,DiscreteKernel>> discretekde;
PLUMED_REGISTER_ACTION(discretekde,"KDE_DISCRETE")
typedef KDE<DiagonalKernelParams,HistogramBeadKernel,KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>> beadkde;
PLUMED_REGISTER_ACTION(beadkde,"KDE_BEADS")
typedef KDE<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>> flatkde;
PLUMED_REGISTER_ACTION(flatkde,"KDE_KERNELS")
typedef KDE<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools> sphericalkde;
PLUMED_REGISTER_ACTION(sphericalkde,"SPHERICAL_KDE_CPU")
typedef AccelerableShortcut<sphericalkde> sphshortcut;
PLUMED_REGISTER_ACTION(sphshortcut,"SPHERICAL_KDE")
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
  keys.addActionNameSuffix("_DISCRETEACC");
  keys.addActionNameSuffix("_KERNELSACC");
  keys.addActionNameSuffix("_BEADSACC");
  keys.addActionNameSuffix("_FULLCOVARACC");
}

KDEShortcut::KDEShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  bool usegpuFLAG=false;
  parseFlag("USEGPU",usegpuFLAG);
  std::string kerneltype;
  parse("KERNEL",kerneltype);
  if( kerneltype=="DISCRETE" ) {
    readInputLine( getShortcutLabel() + ": KDE_DISCRETE" + (usegpuFLAG ? "ACC ":" ") + convertInputLineToString() );
    return;
  }
  std::vector<std::string> args;
  parseVector("ARG", args );
  std::vector<Value*> argvals;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, argvals );
  if( argvals.empty() ) {
    error("no arguments were specified with ARG");
  }
  std::string argstr = " ARG=" + argvals[0]->getName();
  for(unsigned i=1; i<argvals.size(); ++i) {
    argstr += "," + argvals[i]->getName();
  }
  std::vector<std::string> bw;
  parseVector("BANDWIDTH",bw);
  if( bw.empty() ) {
    error("BANDWIDTH is required unless KERNEL=DISCRETE");
  }
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
        readInputLine( getShortcutLabel() + ": KDE_KERNELS" + (usegpuFLAG ? "ACC ":" ") + argstr + " " + bwstr + " KERNEL={" + kerneltype + "} " + convertInputLineToString() );
      } else {
        std::size_t dd = kerneltype.find("-bin");
        readInputLine( getShortcutLabel() + ": KDE_BEADS" + (usegpuFLAG ? "ACC ":" ") + argstr + " " + bwstr + " KERNEL=" + kerneltype.substr(0,dd) + " " + convertInputLineToString() );
      }
    } else if( bwargs[0]->getRank()==2 ) {
      readInputLine( getShortcutLabel() + ": KDE_FULLCOVAR" + (usegpuFLAG ? "ACC ":" ") + argstr + " " + bwstr + " KERNEL={" + kerneltype + "} " + convertInputLineToString() );
    } else {
      error("found strange rank for bandwidth parameter");
    }
  } else if( bw.size()==argvals.size() ) {
    if( kerneltype.find("bin")==std::string::npos ) {
      readInputLine( getShortcutLabel() + ": KDE_KERNELS" + (usegpuFLAG ? "ACC ":" ") + argstr + " " + bwstr + " KERNEL={" + kerneltype + "} " + convertInputLineToString() );
    } else {
      std::size_t dd = kerneltype.find("-bin");
      readInputLine( getShortcutLabel() + ": KDE_BEADS" + (usegpuFLAG ? "ACC ":" ") + argstr + " " + bwstr + " KERNEL=" + kerneltype.substr(0,dd) + " " + convertInputLineToString() );
    }
  } else if( bw.size()==argvals.size()*argvals.size() ) {
    readInputLine( getShortcutLabel() + ": KDE_FULLCOVAR" + (usegpuFLAG ? "ACC ":" ")+ argstr + " " + bwstr + " KERNEL={" + kerneltype + "} " + convertInputLineToString() );
  } else {
    error("invalid input for bandwidth");
  }
}

}
}

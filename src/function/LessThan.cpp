/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION LESS_THAN
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

This action takes one argument, $r$ and evaluates the following function:

$$
w(r) = s(r)
$$

The function $s(r)$ here is a switching function so the output value $w$ is a number between 0 and 1.
Switching functions are typically used to determine if an input value is less than some cutoff.  The value
of $w$ the switches smoothly from one to zero as the input value $r$ crosses the threshold of interest.

The following example, shows how we can apply a switching function on the instantaneous value of the
distance between atom 1 and atom 2.

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

The output here is close to one when the distance between atoms 1 and 2 is less than 0.2 nm close to zero
for values that are greater than 0.2.

## Transforming squares of quantities

If you have computed the square of the distance you can use the flag SQUARED to indicate that the input
quantity is the square of the distance as indicated below:

```plumed
d: DISTANCE COMPONENTS ATOMS=1,2
dm: CUSTOM ARG=d.x,d.y,d.z FUNC=x*x+y*y+z*z PERIODIC=NO
s: LESS_THAN ARG=dm SQUARED SWITCH={RATIONAL R_0=0.2}
```

This option can be useful for improving performance by removing the expensive square root operations.

## Non rank zero arguments

Instead of passing a single scalar in the input to the `LESS_THAN` action you can pass a single vector as shown here:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

The input to the `LESS_THAN` action here is a vector with four elements. The output from the action `b` is similarly
a vector with four elements. In calculating the elements of this vector PLUMED applies the function described in the previous
section on each of the distances in turn. The first element of `b` thus tells you if the distance between atoms 1 and 2 is less than
0.2 nm, the second element tells you if the distance between atoms 3 and 4 is less than 0.2 nm and so on.

You can use the commands in the above example input to evaluate the number of distances that are less than a cutoff as follows:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

The final scalar that is output here is evaluated using the following summation:

$$
s = \sum_i s(d_i)
$$

where the sum over $i$ here runs over the four distances in the above expression. This scalar tells you the number of distances that are
less than 0.2 nm.

Notice that you can do something similar with a matrix as input as shown below:

```plumed
d: DISTANCE_MATRIX GROUPA=1-10 GROUPB=11-20
b: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

This input tells PLUMED to calculate the 100 distances between the atoms in the two input groups. The final value that is printed to the colvar file then
tells you how many of these distances are less than 0.2 nm.

## The MASK keyword

Consider the following input:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Compute the coordination numbers
adj: CONTACT_MATRIX GROUP=1-400 SWITCH={RATIONAL R_0=0.3 D_MAX=1.0} MASK=sphere
ones: ONES SIZE=400
coord: MATRIX_VECTOR_PRODUCT ARG=adj,ones
# Determine if coordination numbers are less than 4 or not
lt: LESS_THAN ARG=coord MASK=sphere SWITCH={RATIONAL R_0=4.0}
# And calculate how many atoms in the region of interst have a coordination number that is less than four
ltsphere: CUSTOM ARG=lt,sphere FUNC=x*y PERIODIC=NO
cv: SUM ARG=ltsphere PERIODIC=NO
PRINT ARG=cv FILE=colvar
```

This input calculates the total number of atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$ and have a
coordination number that is less than 4.  Notice that to reduce the computational expense we use the MASK keyword in the input to
[CONTACT_MATRIX](CONTACT_MATRIX.md) so PLUMED knows to not bother calculating the coordination numbers of atoms that are not within the spherical region of
interest. Further notice that we have also used this MASK keyword in the input to LESS_THAN to prevent PLUMED from transforming the coordination numbers
that have not been calculated with the switching function to get a further speed up.

Using the MASK keyword in this way is necessary if the input argument is a vector. If the input argument is a matrix then you should never need to use the
MASK keyword. PLUMED will use the sparsity pattern for the input matrix to reduce the number of transformations that are performed within LESS_THAN.

## Switching functions types

PLUMED allows you to use a range of different switching function types.  Most of these
require you to provide at least one input parameter $r_0$. The switching function is then defined so that for
$r \le d_0 \quad s(r)=1.0$ while for $r > d_0$ the function decays smoothly to 0.  By changing the switching
function you are using you change the decay that occurs in the switching function for $r$ values that are greater
that $d_0$.

You can specify the value at which the the switching function goes to zero by using the D_MAX keyword.
PLUMED will ensure that the swtiching function is definitely zero for for input values that are greater than
or equal to D_MAX by by stretching and shifting the function using the following transformations.

$
s(r) = \frac{s'(r)-s'(d_{max})}{s'(0)-s'(d_{max})}
$

The various expressions that are used for $s'(r)$ in the above expression are described in the following sections.
Scaling the switching function so that it is exactly zero at $d_{max}$ using the transformation described above
has been the default behaviour within PLUMED since version 2.2.
If you wish to use a D_MAX parameter and the $s'(r)$ functions described in the following sections directly you can use the
keyword `NOSTRETCH`. However, the advangage of performing the stretching and scaling in this way is that the function has
no discontinuities.  If you use the `NOSTRETCH` option the switching function may have a discontinuity at $d_{max}$.

### Rational switching function

The rational switching function is the most commonly used of the switching functions.  The implementation of this function within PLUMED has been carefully
optimised so calculations using this switching function should be fast.

For this function $s'(r)$ has the following functional form:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
\frac{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{n} }{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{m} } & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

The following example illustrates how you can use the rational switching function is used in practice:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.5 D_0=0.1 NN=8 MM=16 D_MAX=1.0}
```

As was discussed in earlier parts of this page, you can also specify a rational switching function using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

In this case, the paramers $d_0$, $n$ and $m$ are set to their default values of 0, 6 and 12 respectively. Meanwhile, if D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

Notice that you can choose to set a subset of the switching function parameters and to use the default values for the unset parameters by using an input like this one:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2 NN=8}
```

The input above sets $n=8$ and $m=2n=16$. Meanwhile, $d_0$ is set to its default value of 0.

Lastly, note that you can also specify the parameters for the rational switching function without using the SWITCH keyword as indicated below:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d D_0=0.1 R_0=0.2 NN=6 MM=12
```

The newer syntax using the SWITCH keyword is better than this option as there is no way to set the D_MAX parameter this way and because if you use this syntax you are forced to
use the RATIONAL switching function type.

### Exponential switching function

The way that the exponential switching function can be used is illustrated in the following example input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={EXP D_0=0.1 R_0=0.2 D_MAX=1.0}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
\exp\left(-\frac{ r - d_0 }{ r_0 }\right) & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

You can also specify that an exponential switching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={EXP R_0=0.2}
```

For this input $d_0$ is set to its default value of 0.  Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

### Gaussian switching function

The way that the gaussian switching function can be used is illustrated in the following example input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={GAUSSIAN D_0=0.1 R_0=0.2 D_MAX=1.0}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
\exp\left(-\frac{ (r - d_0)^2 }{ 2r_0^2 }\right) & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

You can also specify that an gaussian switching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={GAUSSIAN R_0=0.2}
```

For this input $d_0$ is set to its default value of 0.  Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

### Sketch-map switching function

You can specify that the sketch-map switching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={SMAP D_0=1 R_0=4 A=3 B=2 D_MAX=10}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
\left[ 1 + ( 2^{a/b} -1 )\left( \frac{r-d_0}{r_0} \right)^a \right]^{-b/a} & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

This type of swtiching function would be used with PLUMED's implementation of [SKETCHMAP](SKETCHMAP.md). If you are
performing a sketch-map projection an input like the one below would be more common:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={SMAP R_0=4 A=3 B=2}
```

With this input $d_0$ is set to its default value of 0.  Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

### Q switching function

You can specify that the Q switching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={Q D_0=0 REF=0.498366 BETA=50.0 LAMBDA=1.8 R_0=0.01 D_MAX=20}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
s(r) = \frac{1}{1 + \exp(\beta(r_{ij} - \lambda r_{ij}^0))} & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

The value of $r_{ij}^0$ is specified in the above input by the `REF` keyword.  Notice that you can also specify this type of switching function using
the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={Q REF=0.498366 R_0=0.01 }
```

With this input $_0$, $\lambda$ and $\beta$ are set equal to their default values of 0, 1.8 and 50 nm$^{-1}$ respectively. Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed. These default values for $\lambda$ and $\beta$ are suitable if you are simulating all atom models.  However, if you are using
a coarse grainded models values of $\lambda=1.5$ and $\beta=50 nm^-1$ are more appropriate.

### Cubic switching function

You can specify that the cubic swtiching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={CUBIC D_0=0.1 D_MAX=0.5}
```

If this type of expression is used then the $s(r)$ is calculated as:

$$
s(r) = (y-1)^2(1+2y) \qquad \textrm{where} \quad y = \frac{r - d_0}{d_{max}-d_{0}}
$$

No stretching is required for this type of switching function as its functional form ensures that it is zero at $d_{max}$.

### Tanh switching function

You can specify that the cubic swtiching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={TANH D_0=0.1 R_0=0.1 D_MAX=1.0}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s'(r) = \begin{cases}
1 & \textrm{if} \quad r<d_0 \\
1 - \tanh\left( \frac{ r - d_0 }{ r_0 } \right) & \textrm{if} \quad d_0 < r < d_{max} \\
0 & \textrm{otherwise}
\end{cases}
$$

You can also specify that a tanh switching function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={TANH R_0=0.2}
```

For this input $d_0$ is set to its default value of 0.  Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

### Cosinus switching function

You can specify that the cosinum function is to be used by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={COSINUS D_0=0.1 R_0=0.1}
```

The $s'(r)$ that is used in this input is given by the following expression:

$$
s(r) =\left\{\begin{array}{ll}
   1                                                           & \mathrm{if } r \leq d_0 \\
   0.5 \left( \cos ( \frac{ r - d_0 }{ r_0 } \pi ) + 1 \right) & \mathrm{if } d_0 < r\leq d_0 + r_0 \\
   0
\end{array} \right\}
$$

No stretching is required for this type of switching function as its functional form ensures that it is zero at $d_0+r_0$

### Custom switching function

If you want to use a switching function that is different to all of the functions listed above you can use the cusom option that is illustrated below:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={CUSTOM FUNC=1/(1+x^6) D_0=0.1 R_0=0.1 D_MAX=10}
```

This option allows you to use the lepton library that is used to implement the [CUSTOM](CUSTOM.md) action to evaluate your switching function. With this
option you specify the functional form for $s'(r)$ by providing a function of `x` to the `FUNC` keyword.  The `x` that enters the switching function
definition that you provide is given by:

$$
x = \frac{r - d_0}{r_0}
$$

Notice, that you can also express your switching function in terms of $x^2$ as illustrated in the input below:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={CUSTOM FUNC=1/(1+x2^3) D_0=0.1 R_0=0.1 D_MAX=10}
```

which is equivalent to the earlier input with the CUSTOM switching function.  However, using `x2` in place of `x` is often more computationally efficient
as you can avoid an expensive square root operation.  Notice, lastly, that just as there is a [MATHEVAL](MATHEVAL.md) action that is equivalent to [CUSTOM](CUSTOM.md),
there is also a matheval equivalent that you can use here:

```plumed
d: DISTANCE ATOMS=1,2
s: LESS_THAN ARG=d SWITCH={MATHEVAL FUNC=1/(1+x2^3) R_0=0.1}
```

For this input $d_0$ is set to its default value of 0. Furthermore, as D_MAX is unset the stretching and scaling that
was described in the previous section is not performed.

!!! caution "performance of CUSTOM"

    With the default implementation CUSTOM is slower than other functions
    (e.g., it is slower than an equivalent RATIONAL function by approximately a factor 2).
    You can find information on how to improve its performance in the documenation for [CUSTOM](CUSTOM.md)

*/
//+ENDPLUMEDOC

class LessThan {
public:
  bool squared;
#ifdef __PLUMED_USE_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_USE_OPENACC
  static void registerKeywords( Keywords& keys );
  static void read( LessThan& func,
                    ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const LessThan& func,
                    bool noderiv,
                    const View<const double>& args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],squared)
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(squared,this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

typedef FunctionShortcut<LessThan> LessThanShortcut;
PLUMED_REGISTER_ACTION(LessThanShortcut,"LESS_THAN")
typedef FunctionOfScalar<LessThan> ScalarLessThan;
PLUMED_REGISTER_ACTION(ScalarLessThan,"LESS_THAN_SCALAR")
typedef FunctionOfVector<LessThan> VectorLessThan;
PLUMED_REGISTER_ACTION(VectorLessThan,"LESS_THAN_VECTOR")
typedef FunctionOfMatrix<LessThan> MatrixLessThan;
PLUMED_REGISTER_ACTION(MatrixLessThan,"LESS_THAN_MATRIX")

void LessThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the input is less than a threshold");
}

void LessThan::read( LessThan& func, ActionWithArguments* action, FunctionOptions& options ) {
  options.derivativeZeroIfValueIsZero = true;
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }


  std::string errors;
  std::string sfinput;
  action->parse("SWITCH",sfinput);
  if(sfinput.length()>0) {
    func.switchingFunction.set(sfinput,errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    action->parse("R_0",r0);
    if(r0<=0.0) {
      action->error("R_0 should be explicitly specified and positive");
    }
    action->parse("D_0",d0);
    action->parse("NN",nn);
    action->parse("MM",mm);
    func.switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<func.switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",func.squared);
  if( func.squared ) {
    action->log<<"  input quantity is square of quantity that switching function acts upon\n";
  }
}

void LessThan::calc( const LessThan& func,
                     bool noderiv,
                     const View<const double>& args,
                     FunctionOutput& funcout ) {
  double d;
  if( func.squared ) {
    funcout.values[0] = func.switchingFunction.calculateSqr( args[0], d );
  } else {
    funcout.values[0] = func.switchingFunction.calculate( args[0], d );
  }
  if( !noderiv ) {
    funcout.derivs[0][0] = args[0]*d;
  }
}

}
}



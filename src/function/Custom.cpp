/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "Custom.h"
#include "core/ActionRegister.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "tools/OpenMP.h"
#include "tools/LeptonCall.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION CUSTOM
/*
Calculate a combination of variables using a custom expression.

CUSTOM is one of the most useful actions in PLUMED.  This action takes in a list of arguments and then uses
the [lepton mathematical expression parser](https://simtk.org/projects/lepton) to evaluate a user defined function
of these input arguments. We can thus use this action in the input below to perform a metadynamics
simulation using the difference between two distances as a CV.

```plumed
dAB: DISTANCE ATOMS=10,12
dAC: DISTANCE ATOMS=10,15
diff: CUSTOM ARG=dAB,dAC FUNC=y-x PERIODIC=NO
# notice: the previous line could be replaced with the following
# diff: COMBINE ARG=dAB,dAC COEFFICIENTS=-1,1
METAD ARG=diff SIGMA=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100
```

The particular function that should be evaluated from the input arguments is specified using the `FUNC` keyword.
The function provided to the `FUNC` keyword is written in terms of `x` and `y` in the input above.  `x` here deonetes the
first argument provided to the ARG keyword, `dAB`, while `y` is the second argument provided to the `ARG` keyword, `dAC`.

## The VAR keyword

If you wish, you can rewrite the example input above more transparantly as:

```plumed
dAB: DISTANCE ATOMS=10,12
dAC: DISTANCE ATOMS=10,15
diff: CUSTOM ARG=dAB,dAC VAR=dAB,dAC FUNC=dAC-dAB PERIODIC=NO
METAD ARG=diff SIGMA=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100
```

By using the `VAR` keyword here we ensure that we can use the labels of the arguments in the mathematical expression that
we provide to the `FUNC` argument.  Notice that, if you have four or more arguments you must use the `VAR` keyword.  With
three or less arguments the `VAR` keyword can be ommitted and the symbols `x`, `y` and `z` can be used to denote the first,
second and third arguments respectively.

The following input illustrates a case where using the `VAR` keyword is essential.  This input tells PLUMED to
print the angle between the vector connecting atoms 1,2 and the vector connecting atoms 2,3.

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
d2: DISTANCE ATOMS=2,3 COMPONENTS
theta: CUSTOM ...
  ARG=d1.x,d1.y,d1.z,d2.x,d2.y,d2.z
  VAR=ax,ay,az,bx,by,bz
  FUNC=acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz)))
  PERIODIC=NO
...

PRINT ARG=theta
```

## Functions and constants

The previous examples demonstrates how CUSTOM can be used to evaluate mathematical expressions that involve taking powers (^), adding (+),
multiplying (*), dividing (/) and subtracting (-) variables and can control the order in which operations are performed by using parenthesis.
In addition to these basic binary operations you can also use a range of mathematical functions in your expressions.  The following table lists all the mathematical functions
that you can use in in the input to the `FUNC` keyword for CUSTOM.

| Function | Description |
|:--------:|:------------|
| sqrt(x)       | The square root of x |
| exp(x)        | The exponential of x i.e. $e^{x}$ |
| log(x)        | The natural logarithm of x |
| sin(x)        | The sine of x |
| cos(x)        | The cosine of x |
| sec(x)        | The reciprocal of the cosine of $x$. $\frac{1}{\cos(x)}$ |
| csc(x)        | The reciprocal of the sine of $x$. $\frac{1}{\sin(x)}$ |
| tan(x)        | The tangent of x i.e. $\frac{\sin(x)}{\cos(x)}$ |
| cot(x)        | The reciprocal of the tangent of $x$. $\frac{1}{\tan(x)}$ |
| asin(x)       | The principal arc sine of x. Returns $-\frac{\pi}{2} \le y \le \frac{\pi}{2}$ which gives $x = \sin(y)$ |
| acos(x)       | The principal arc cosine of x. Returns $0 \le y \le \pi$ which gives $x=\cos(y)$ |
| atan(x)       | The principal arc tangent of x.  Returns $-\frac{\pi}{2} \le y \le \frac{\pi}{2}$ which gives $x = \tan(y)$ |
| atan2(x,y)    | The principal arg tangent of $\frac{x}{y}$.  Returns $-\pi \le z \le \pi$ which gives $\frac{x}{y} = \tan(z)$ |
| sinh(x)       | The hyperbolic sine of $x$ |
| cosh(x)       | the hyperbolic cosine of $x$ |
| tanh(x)       | The hyperbolic tangent of $x$ |
| erf(x)        | The [error function](https://en.wikipedia.org/wiki/Error_function) $\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2} \textrm{d}t$ |
| erfc(x)       | The complementary error function $1-\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2} \textrm{d}t$ |
| step(x)       | 1 if $x \ge 0$ and 0 otherwise |
| delta(x)      | inf if $x=0$ and 0 otherwise |
| nandelta(x)   | nan if $x=0$ and 0 otherwise |
| square(x)     | The square of x i.e. $x^2$ |
| cube(x)       | The cube of x i.e. $x^3$ |
| recip(x)      | The reciprocal of x i.e. $\frac{1}{x}$ |
| min(x,y)      | If $x<y$ this function returns $x$.  If $y\le x$ this function returns $y$ |
| max(x,y)      | If $x>y$ this function returns $x$.  If $y\ge x$ this function returns $y$ |
| abs(x)        | The absolute value of x |
| floor(x)      | The largest integer that is less than $x$ |
| ceil(x)       | The smallest integer that is greater than $x$ |
| select(x,y,z) | If $x==0$ this returns $z$ otherwise this returns $y$ |
| acot(x)       | Returns the value of $-\frac{\pi}{2} \le y \le \frac{\pi}{2}$ which gives $\frac{1}{x}=\tan(y)$ |
| asec(x)       | Returns the value of $0 \le y \le \pi$ which gives $\frac{1}{x}=\cos(y)$ |
| acsc(x)       | Returns the value of $-\frac{\pi}{2} \le y \le \frac{\pi}{2}$ which gives $\frac{1}{x}=\sin(y)$ |
| coth(x)       | The recipricoal of the hyperbolic tangent of $x$. $\frac{1}{\tanh(x)}$ |
| sech(x)       | The recipricoal of the hyperbolic cosine of $x$. $\frac{1}{\cosh(x)}$ |
| csch(x)       | The recipricoal of the hyperbolic sine of $x$. $\frac{1}{\sinh(x)}$ |
| asinh(x)      | The nonnegative area hyperbolic sine of $x$. Returns $y$ which gives $x=\sinh(y)$ |
| acosh(x)      | The nonnegative area hyperbolic cosine of $x$. Returns $0\le y \le \infty$ which gives $x=\cosh(y)$ |
| atanh(x)      | The nonnegative area hyperbolic tangent of $-1 \le x \le 1$.  Returns $y$ which gives $x=\tanh(y)$. |
| acoth(x)      | The inverse hyperbolic tangent of $x$ is calculated as $\frac{1}{2}\ln\left( \frac{x+1}{x-1} \right)$ |
| asech(x)      | The inverse hyperbolic secant of $x$ is calculated as $\log\left( \sqrt{\frac{1}{x}-1}\sqrt{\frac{1}{x}+1} + \frac{1}{x}\right)$ |
| acsch(x)      | The inverse hyperbolic cosecant of $x$ is calculated as $\log\left(\frac{1}{x}+\sqrt{\frac{1}{x^2}+1}\right)$ |

Notice, that you can also incorporate the following constants in the expressions that are used in the input to FUNC:

| Symbol | Description |
| :----- |:------------|
| `e` | Euler's number - the base for the natural logarithm |
| `log2e` | $1 / \log(2)$ |
| `log10e` | $1 / \log(10)$ |
| `ln2` | $\log(2)$ |
| `ln10` | $\log(10)$ |
| `pi`   | the circle constant $\pi$ |
| `pi_2` | $\pi / 2$ |
| `pi_4` | $\pi / 4$ |
| `sqrt2 | $\sqrt(2)$ |
| `sqrt1_2 ` | $\sqrt(0.5)$ |

The way one of these constants can be used in an expression is illustrated in the following example:

```plumed
d: DISTANCE ATOMS=1,2
# This function is evaluating the area of the circle whose radius is
# given by the distance between atoms 1 and 2.
c: CUSTOM ARG=d FUNC=pi*x*x PERIODIC=NO
PRINT ARG=c FILE=colvar
```

## The PERIODIC keyword

If the function you have evaluated has a periodic codomain then you have to use the PERIODIC keyword to specify the codomain of the function.
To see how to use this keyword consider the following input:

```plumed
d: DISTANCE ATOMS=1,2 COMPONENTS
theta: CUSTOM ARG=d.x,d.y FUNC=atan2(y,x) PERIODIC=-pi,+pi
PRINT ARG=theta FILE=colvar
```

The above input calculates the torsional angle between the vector connecting atoms 1 and 2 and positive $x$ direction around the $z$ axis.
This quantity must lie between $-\pi$ and $\pi$ so we use the PERIODIC keyword to tell PLUMED this fact about the output quantity.

Notice that PLUMED cannot automatically detect if the codomain of the function is periodic.  You thus always have to specifiy the
periodicity of the codomain whenever you use the CUSTOM command.  However, the codomain of most functions that people use are not periodic
so you can normally just write `PERIODIC=NO`

## The step function

The `step` operation (that is the Heaviside function) from the table above allow you to use if clauses in CUSTOM actions.
As discussed in the table `step(x)` is 1 when `x` is positive and `0` when `x` is negative.  So the function `step(1-x)` is
1 if x is less than one and zero otherwise.

Using the `step` operation multiple times in a function allows you to perform logical operations.  For example, the equivalent of the AND operator
is the product so, for example, `step(1.0-x)*step(x-0.5)` is only equal to 1 when x is greater than 0.5 AND less than 1.0.
By a similar logic we can use the function `1-step(1.0-x)*step(x-0.5)` to create a value that is 1 if x is less than 0.5 OR
greater than 1.0.

The example below illustrtes how we can put these ideas of using the `step` operation into practise.

```plumed
d: DISTANCE ATOMS=10,15
m: CUSTOM ARG=d FUNC=0.5*step(0.5-x)+x*step(x-0.5) PERIODIC=NO
# check the function you are applying:
PRINT ARG=d,m FILE=checkme
RESTRAINT ARG=d AT=0.5 KAPPA=10.0
```

The meaning of the function `0.5*step(0.5-x)+x*step(x-0.5)` in this example is:
- If x<0.5 (step(0.5-x)!=0) use 0.5
- If x>0.5 (step(x-0.5)!=0) use x
Notice that the same result can be achieved by using [UPPER_WALLS](UPPER_WALLS.md)
However, with CUSTOM you can create much more complex definitions.

Notice that we can apply a force on the value `m` by using the [RESTRAINT](RESTRAINT.md) command as the function we defined in the expression that was passed to the `FUNC` keyword is continuous.
In general, however, you must be careful when using the `step`, `delta`, `nandelta` and `select` functions as you can easily write expression for
discontinuous functions by using these operations. If you want to check if a function you have created using `step`
is continuous you can easily plot the function in gnuplot by using a commands like those shown below

````
# this allow to step function to be used in gnuplot:
gnuplot> step(x)=0.5*(erf(x*10000000)+1)
# here you can test your function
gnuplot> p 0.5*step(0.5-x)+x*step(x-0.5)
````

## Using TIME as a variable

Notice that you can use CUSTOM to implement a [MOVINGRESTRAINT](MOVINGRESTRAINT.md) as shown below.

```plumed
t: TIME
d: DISTANCE ATOMS=1,2
f: CUSTOM ARG=d,t FUNC=100*(x-((0.2-0.1)*y/100))^2 PERIODIC=NO
BIASVALUE ARG=f
```

In a 100~ps simulation that was run with this input the distance beetween atom 1 and 2 would be forced to increase from 0.1 to 0.2 nm.

## Using CUSTOM in shortcuts

The CUSTOM action is used in many of the [shortcut actions](shortcuts.md) that are implemented in PLUMED. We think that using this action in these places is beneficial as it ensures that the mathematical
expressions that are used in the method are visible in the log. We have found that there are many actions that
are used in relatively few papers. When implementing these actions we think that sharing implementations of these methods that are comprehensible is more important than sharing methods that are
fast.

As an example of why this is useful consider some of the variants of the DISTANCE keyword that were present in PLUMED 1.3. These variants allowed you to compute the distance between a point and a line defined by
two other points or the progression along that line.  In PLUMED 2.10 we can implement these variables using the following input file:

```plumed
# take center of atoms 1 to 10 as reference point 1
p1: CENTER ATOMS=1-10
# take center of atoms 11 to 20 as reference point 2
p2: CENTER ATOMS=11-20
# take center of atoms 21 to 30 as reference point 3
p3: CENTER ATOMS=21-30

# compute distances
d12: DISTANCE ATOMS=p1,p2
d13: DISTANCE ATOMS=p1,p3
d23: DISTANCE ATOMS=p2,p3

# compute progress variable of the projection of point p3
# along the vector joining p1 and p2
# notice that progress is measured from the middle point
onaxis: CUSTOM ARG=d13,d23,d12 FUNC=(0.5*(y^2-x^2)/z) PERIODIC=NO

# compute between point p3 and the vector joining p1 and p2
fromaxis: CUSTOM ARG=d13,d23,d12,onaxis VAR=x,y,z,o FUNC=(0.5*(y^2+x^2)-o^2-0.25*z^2) PERIODIC=NO

PRINT ARG=onaxis,fromaxis
```

The equations in this input were also used to combine [RMSD](RMSD.md) values from different snapshots of a protein so as to define
progression (S) and distance (Z) variables in the paper that is cited in the bibliography.  We can understand how these expressions
are derived by noting that $x$, $y$ and $z$ are the distances between atoms 1 and 3, 2 and 3 and 1 and 2 respectively.  The projection
of the vector connecting atom 1 to atom 3 onto the vector connecting atom 1 to atom 2 is thus $x\cos(\theta)$, where theta is the angle
between the vector connecting atoms 1 and 3 and the vector connecting atoms 1 and 2.  We can arrive at the following expression for $x\cos(\theta)$
by rearranging the cosine rule:

$$
x\cos(\theta) = \frac{y^2 - x^2}{z} - \frac{z}{2}
$$

Notice that the value called `onaxis` in the above input is thus $o=x\cos(\theta) + \frac{z}{2}$.  Adding the factor of $\frac{z}{2}$ ensures that the origin
is at the center of the bond connecting atom 1 to atom 2.

The value `fromaxis` measures the square of the distance from the the line.  It is calculated using pythagoras theorem as follows:

$$
f^2 = y^2 - x^2\cos^2(\theta)
$$

Inserting $x\cos(\theta) = o - \frac{z}{2}$ into this expression gives:

$$
f^2 = y^2 - (o -\frac{z}{2})^2 = y^2 - o^2 + oz - \frac{z^2}{4}
$$

Inserting the fact that $oz = \frac{y^2 - x^2}{2}$, which comes from the expression for $o$ that was used to calculate `onaxis`, gets us to the expression that is used to calculate `fromaxis`.

## CUSTOM with vector arguments

The examples above have shown how CUSTOM operates when the input arguments are scalars. You can also pass vectors in the argument to a CUSTOM action. The function you define will then
be applied to each element of the vector in turn.  So, for example in the input below a vector that contains three angles is passed to the CUSTOM action. The CUSTOM action calculates the
cosine of these three angles and outputs them in a three dimensional vector called `c` that is printed to the colvar file.

```plumed
a: ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
c: CUSTOM ARG=a FUNC=cos(x) PERIODIC=NO
PRINT ARG=c FILE=colvar
```

You are not confined to passing a single vector to CUSTOM. The input below shows what you can do by pass two vectors with the same numbers of elements. The first element of the vector
output by the CUSTOM action here contains the projection of the vector connecting atom 1 and 2 on the vector connecting atom 2 and 3, the second element contains the projection of the
vector connecting atoms 4 and 5 on the vector connecting atoms 5 and 6 and the final element contains the projection of the vector connecting atoms 7 and 8 on the vector connecting atoms
8 and 9.

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=4,5 ATOMS3=7,8
a: ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
c: CUSTOM ARG=d,a FUNC=x*cos(y) PERIODIC=NO
PRINT ARG=c FILE=colvar
```

Notice, when we multiply two vectors in CUSTOM the output is a vector.  This product that emerges from using a CUSTOM action is __not__ the scalar or cross product of the input vectors.

Lastly, notice that you can pass a mixture of vectors and scalars in the input to a CUSTOM action. However, the labels of the scalar values that appear in the input for the ARG keyword
must appear after the labels of vector vectors as shown below:

```plumed
d: DISTANCE ATOMS=1,2
a: ANGLE ATOMS1=1,2,3 ATOMS2=1,2,4 ATOMS3=1,2,5 ATOMS4=1,2,6
c: CUSTOM ARG=a,d FUNC=y*cos(x) PERIODIC=NO
PRINT ARG=c FILE=colvar
```

The multiplication of a scalar by a vector in the above input is done in [the usual way](https://en.wikipedia.org/wiki/Scalar_multiplication).
Similarly, dividing a vector by a scalar is equivalent to multiplying the vector by the reciprocal of the scalar.  If you write an expression that adds or subtract a
scalar from a vector addition is performed.  To understand why consider the following example that adds the constant `c` to the input vector of distances `d`. The
result of this operation is a vector `f` that contains the three input distances from `d` with 23 added to each of them.

```plumed
c: CONSTANT VALUE=23
d: DISTANCE ATOMS1=1,2 ATOMS2=4,5 ATOMS3=7,8
f: CUSTOM ARG=d,c FUNC=x+y PERIODIC=NO
PRINT ARG=f FILE=colvar
```

## CUSTOM, MASK and vector arguments

The example input below illustrates a way you can use the MASK keyword to improve PLUMED performance for inputs that use CUSTOM actions that take arguments that
are vectors:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
ow: GROUP ATOMS=1-400
hw: GROUP ATOMS=401-1200
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=ow CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate coordination with oxygen atoms
ow_ones: ONES SIZE=400
ow_mat: CONTACT_MATRIX GROUP=ow MASK=sphere SWITCH={RATIONAL R_0=0.3 D_MAX=0.8}
ow_coord: MATRIX_VECTOR_PRODUCT ARG=ow_mat,ow_ones
# And transform the coordination numbers by a switching function
ow_mt: MORE_THAN ARG=ow_coord MASK=sphere SWITCH={RATIONAL R_0=4}
# Calculate the coordination with hydrogen atoms
hw_ones: ONES SIZE=800
hw_mat: CONTACT_MATRIX GROUPA=ow GROUPB=hw MASK=sphere SWITCH={RATIONAL R_0=0.1 D_MAX=0.5}
hw_coord: MATRIX_VECTOR_PRODUCT ARG=hw_mat,hw_ones
# And transform the coordination numbers by a switching function
hw_mt: MORE_THAN ARG=hw_coord MASK=sphere SWITCH={RATIONAL R_0=2}
# And multiply the two vectors of transformed coordination numbers together
prod: CUSTOM ARG=hw_mt,ow_mt FUNC=x*y MASK=sphere PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
cv: SUM ARG=prod PERIODIC=NO
# And print out final CV to a file
PRINT ARG=cv FILE=colvar STRIDE=1
```

The variable being computed here measures the number of oxygen atoms in a spherical region centered on the point $(2.5,2.5,2.5)$ that have a coordination
that is greater than 4 with other oxygen atoms and a coordination that is greater than 2 with other hydrogen atoms.  Notice how by using the MASK keyword
in the [CONTACT_MATRIX](CONTACT_MATRIX.md), [MORE_THAN](MORE_THAN.md and CUSTOM actions we avoid calculating the order parameters that are not in the spherical
region of interest in the way that is described in more detail in the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).


## CUSTOM with matrix arguments

You can also pass matrices in the argument to a CUSTOM action. These input matrices are treated similarly to input vectors. In other words, any function you define
is applied to each element of the matrix in turn so if the input matrix is $N \times M$ the output matrix is also $N \times M$.  The following example illustrates how you
can use this functionality to calculate all the angles between a set of bond vectors.

```plumed
# Calculate the vectors connecting four atoms
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
# Calculate the norm of these four vectors
norm: CUSTOM ARG=d.x,d.y,d.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
# Now calculate the directors of the vectors
norm_x: CUSTOM ARG=d.x,norm FUNC=x/y PERIODIC=NO
norm_y: CUSTOM ARG=d.y,norm FUNC=x/y PERIODIC=NO
norm_z: CUSTOM ARG=d.z,norm FUNC=x/y PERIODIC=NO
# And combine all these directors in a matrix
stack: VSTACK ARG=norm_x,norm_y,norm_z
# Now calculate the matrix of dot products between directors (these are the cosines of the angles)
stackT: TRANSPOSE ARG=stack
cosa: MATRIX_PRODUCT ARG=stack,stackT
# And finally get the 4x4 matrix of angles and print it to a file
angles: CUSTOM ARG=cosa FUNC=acos(x) PERIODIC=NO
PRINT ARG=angles FILE=colvar
```

Notice that you can pass multiple $N\times M$ matrices in the input to a CUSTOM action as illustrated in the example below:

```plumed
c1: CONSTANT VALUES=2,3,4,5 NROWS=2 NCOLS=2
c2: CONSTANT VALUES=1,0,0,1 NROWS=2 NCOLS=2
f: CUSTOM ARG=c1,c2 FUNC=x*y PERIODIC=NO
PRINT ARG=f FILE=colvar
```

The four by four matrix that is calculated by the custom action in this input is given by:

$$
f = \left(
\begin{matrix}
2 & 0 \\
0 & 5
\end{matrix}
\right)
$$

which is the element-wise [Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)) and __not__ the matrix product.
By a similar token if you have a custom command that takes two matrices in input and use `FUNC=x/y` the Hadamard product between the matrix
x and a matrix that contains the reciprocals of each of the elements in y is computed.  If you wish to calcalculate the
[product of two matrices](https://en.wikipedia.org/wiki/Matrix_multiplication) you should use the [MATRIX_PRODUCT](MATRIX_PRODUCT.md) comamnd.
Similarly, if you want to calculate the product of a matrix and a vector you should use the [MATRIX_VECTOR_PRODUCT](MATRIX_VECTOR_PRODUCT.md)
command.

Lastly, note that you can pass a mixture of scalars and $N\times M$ matrices in the input to a CUSTOM command. As with vectors, you can think of
any scalars you pass as being converted into $N\times M$ matrix in which every element is equal to the input scalar.  Furthermore, the labels of
the input scalars must appear __after__ the labels fo the input matrices in the input for the ARG keyword.

## CUSTOM, MASK and matrix arguments

The example input below illustrates a way you can use the MASK keyword to improve PLUMED performance for inputs that use CUSTOM actions that take arguments that
are matrices:

```plumed
# Calculate the distances between atoms
d_mat: DISTANCE_MATRIX GROUP=1-64 CUTOFF=4.5 COMPONENTS
# Find the six nearest atom to each of the coordinates
nn: NEIGHBORS ARG=d_mat.w NLOWEST=6
# Evalulate the lengths of the bonds connecting each atom to its nearest neighbors
r: CUSTOM ...
   ARG=d_mat.x,d_mat.y,d_mat.z MASK=nn
   FUNC=sqrt(x*x+y*y+z*z)
   PERIODIC=NO
...
# Evaluate a function of the nearest neighbors
d_wvfunc: CUSTOM ...
   ARG=nn,d_mat.x,d_mat.y,d_mat.z,r MASK=nn
   FUNC=w*(((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3)
   VAR=w,x,y,z,r
   PERIODIC=NO
...
# Calculate the sum of fcc cubic function values for each atom
d_ones: ONES SIZE=64
d: MATRIX_VECTOR_PRODUCT ARG=d_wvfunc,d_ones
# Calculate the number of neighbours
d_denom: MATRIX_VECTOR_PRODUCT ARG=nn,d_ones
# Calculate the average value of the fcc cubic function per bonds
d_n: CUSTOM ARG=d,d_denom FUNC=x/y PERIODIC=NO
d_mean: MEAN ARG=d_n PERIODIC=NO
PRINT ARG=d_mean FILE=colv
```

Using the MASK here ensures that we only calculate the CUSTOM functions for the six nearest neighbors of each of the input atoms. However, it is important to note that
this input is rather unsual.  In general you rarely need to use the MASK keyword for CUSTOM actions that take matrices in input.  For example, if we use the following
(more conventional version) of the above CV the MASK is unecessary as the sparsity pattern of the input matrices are all identical.

```plumed
cmat: CONTACT_MATRIX GROUP=1-64 SWITCH={RATIONAL R_0=1.0 D_MAX=4.5} COMPONENTS
# Evalulate the lengths of the bonds connecting each atom to its nearest neighbors
r: CUSTOM ...
   ARG=cmat.x,cmat.y,cmat.z
   FUNC=sqrt(x*x+y*y+z*z)
   PERIODIC=NO
...
# Evaluate a function of the nearest neighbors
d_wvfunc: CUSTOM ...
   ARG=cmat.w,cmat.x,cmat.y,cmat.z,r
   FUNC=w*(((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3)
   VAR=w,x,y,z,r
   PERIODIC=NO
...
# Calculate the sum of fcc cubic function values for each atom
d_ones: ONES SIZE=64
d: MATRIX_VECTOR_PRODUCT ARG=d_wvfunc,d_ones
# Calculate the number of neighbours
d_denom: MATRIX_VECTOR_PRODUCT ARG=cmat.w,d_ones
# Calculate the average value of the fcc cubic function per bonds
d_n: CUSTOM ARG=d,d_denom FUNC=x/y PERIODIC=NO
d_mean: MEAN ARG=d_n PERIODIC=NO
PRINT ARG=d_mean FILE=colv
```

## CUSTOM with grid arguments

CUSTOM will also accept a function on a grid in input.  You might use this feature if you want to calculate a free energy from a histogram
using $F(s) = -k_BT\log[H(s)]$ as illustrated in the input below:

```plumed
x: DISTANCE ATOMS=1,2
hA: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
fE: CUSTOM ARG=hA FUNC=-2.5*log(x) PERIODIC=NO
DUMPGRID ARG=fE FILE=fes.dat
```

Another way that this functonality is routinely used in PLUMED is in calculating the density field for a symmetry function. The example below shows how you
would do this in practise. The final output here is a function on a grid that tells you the average value of the FCC order parameter at each point in the cell.

```plumed
# Evaluate the FCC order parameter for all of the atoms
fcc: FCCUBIC SPECIES=1-5184 SWITCH={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
# Calculate the distance between each atom and the origin on atom 1
dens_dist: DISTANCES ORIGIN=1 ATOMS=fcc COMPONENTS
# Do a KDE using the FCC order parameters for the weights of each gaussian and the positions of the atoms as the centers
dens_numer: KDE HEIGHTS=fcc_n ARG=dens_dist.x,dens_dist.y,dens_dist.z GRID_BIN=14,14,28 BANDWIDTH=1.0,1.0,1.0
# Estimate the density of atoms at each point in the box
ones: ONES SIZE=5184
dens_denom: KDE ARG=dens_dist.x,dens_dist.y,dens_dist.z GRID_BIN=14,14,28 HEIGHTS=ones BANDWIDTH=1.0,1.0,1.0
# Now calculate the average value of the order parameter at each point in the box
dens: CUSTOM ARG=dens_numer,dens_denom FUNC=x/y PERIODIC=NO
DUMPCUBE ARG=dens FILE=dens.cube FMT=%8.4f
```

## Making lepton library faster

If you use a lot of [CUSTOM](CUSTOM.md) actions or CUSTOM [switching functions](LESS_THAN.md) it can slow down PLUMED.
These commands use the lepton library that is included in PLUMED.
This library replaces libmatheval since PLUMED 2.5, and by itself it is significantly faster than libmatheval.
However, you can make it even faster using a [just-in-time compiler](https://github.com/asmjit/asmjit.git).
As of PLUMED 2.6, the correct version of ASMJIT is embedded in PLUMED.
As of PLUMED 2.8, ASMJIT is enabled by default on supported architectures (X86/X64).
You can disable it at runtime setting the environment variable `PLUMED_USE_ASMJIT`:

```bash
export PLUMED_USE_ASMJIT=no
```

In some case using a custom expression is almost as fast as using a hard-coded
function. For instance, with an input that contained the following lines:

```plumed
c: COORDINATION GROUPA=1-108 GROUPB=1-108 R_0=1
d_fast: COORDINATION GROUPA=1-108 GROUPB=1-108 SWITCH={CUSTOM FUNC=1/(1+x2^3) R_0=1}
```

I (GB) obtained the following timings (on a Macbook laptop):

````
...
PLUMED: 4A  1 c                                          108     0.126592     0.001172     0.000701     0.002532
PLUMED: 4A  2 d_fast                                      108     0.135210     0.001252     0.000755     0.002623
...
````

Notice the usage of `x2` as a variable for the switching function, which
avoids an unnecessary square root calculation (this is done automatically by the hard-coded switching functions
when you use only even powers). The asmjit calculation (`d_fast`) takes less than 10% more than the hard-coded
one (`c`).


*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Custom> CustomShortcut;
PLUMED_REGISTER_ACTION(CustomShortcut,"CUSTOM")
PLUMED_REGISTER_ACTION(CustomShortcut,"MATHEVAL")
typedef FunctionOfScalar<Custom> ScalarCustom;
PLUMED_REGISTER_ACTION(ScalarCustom,"CUSTOM_SCALAR")
PLUMED_REGISTER_ACTION(ScalarCustom,"MATHEVAL_SCALAR")
typedef FunctionOfVector<Custom> VectorCustom;
PLUMED_REGISTER_ACTION(VectorCustom,"CUSTOM_VECTOR")
PLUMED_REGISTER_ACTION(VectorCustom,"MATHEVAL_VECTOR")
typedef FunctionOfMatrix<Custom> MatrixCustom;
PLUMED_REGISTER_ACTION(MatrixCustom,"CUSTOM_MATRIX")
PLUMED_REGISTER_ACTION(MatrixCustom,"MATHEVAL_MATRIX")

//+PLUMEDOC FUNCTION MATHEVAL
/*
An alias to the CUSTOM function that can also be used to calaculate combinations of variables using a custom expression.

The documentation for this action is identical to that for [CUSTOM](CUSTOM.md).  You can thus use it to evaluate a evaluate an arbitrary function as in the following example input:

```plumed
dAB: DISTANCE ATOMS=10,12
dAC: DISTANCE ATOMS=10,15
diff: MATHEVAL ARG=dAB,dAC FUNC=y-x PERIODIC=NO
# notice: the previous line could be replaced with the following
# diff: COMBINE ARG=dAB,dAC COEFFICIENTS=-1,1
METAD ARG=diff SIGMA=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100
```

This alias is kept in order to maintain compatibility with previous PLUMED versions.
However, notice that as of PLUMED 2.5 the libmatheval library is not linked anymore,
and that the MATHEVAL action evaluates functions [the Lepton library](https://simtk.org/projects/lepton).

*/
//+ENDPLUMEDOC

void Custom::registerKeywords(Keywords& keys) {
  if( keys.getDisplayName()=="MATHEVAL") {
    keys.setDeprecated("CUSTOM");
  }
  keys.use("PERIODIC");
  keys.add("compulsory","FUNC","the function you wish to evaluate");
  keys.add("optional","VAR","the names to give each of the arguments in the function.  If you have up to three arguments in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("scalar/vector/matrix/grid","an arbitrary function");
  keys.addDOI("10.1093/nar/gkv872");
}

void Custom::read( Custom& f, ActionWithArguments* action, FunctionOptions& options ) {
  // Read in the variables
  unsigned nargs = action->getNumberOfArguments();
  ActionWithVector* av=dynamic_cast<ActionWithVector*>(action);
  if( av && av->getNumberOfMasks()>0 ) {
    nargs = nargs - av->getNumberOfMasks();
  }
  action->parseVector("VAR",f.var);
  action->parse("FUNC",f.func);
  if(f.var.size()==0) {
    f.var.resize(nargs);
    if(f.var.size()>3) {
      action->error("Using more than 3 arguments you should explicitly write their names with VAR");
    }
    if(f.var.size()>0) {
      f.var[0]="x";
    }
    if(f.var.size()>1) {
      f.var[1]="y";
    }
    if(f.var.size()>2) {
      f.var[2]="z";
    }
  }
  if(f.var.size()!=nargs) {
    action->error("Size of VAR array should be the same as number of arguments");
  }
  // Check for operations that are not multiplication (this can probably be done much more cleverly)
  bool onlymultiplication = f.func.find("*")!=std::string::npos;
  // Find first bracket in expression
  if( f.func.find("(")!=std::string::npos ) {
    std::size_t br = f.func.find_first_of("(");
    std::string subexpr=f.func.substr(0,br);
    onlymultiplication = f.func.find("*")!=std::string::npos;
    if( subexpr.find("/")!=std::string::npos ) {
      std::size_t sl = f.func.find_first_of("/");
      std::string aa = subexpr.substr(0,sl);
      subexpr=aa;
    }
    if( subexpr.find("+")!=std::string::npos || subexpr.find("-")!=std::string::npos ) {
      onlymultiplication=false;
    }
    // Now work out which vars are in multiplication
    if( onlymultiplication ) {
      for(unsigned i=0; i<f.var.size(); ++i) {
        if( subexpr.find(f.var[i])!=std::string::npos &&
            action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
          f.check_multiplication_vars.push_back(i);
        }
      }
    }
  } else if( f.func.find("/")!=std::string::npos ) {
    onlymultiplication=true;
    if( f.func.find("+")!=std::string::npos || f.func.find("-")!=std::string::npos ) {
      onlymultiplication=false;
    }
    if( onlymultiplication ) {
      std::size_t br = f.func.find_first_of("/");
      std::string subexpr=f.func.substr(0,br);
      for(unsigned i=0; i<f.var.size(); ++i) {
        if( subexpr.find(f.var[i])!=std::string::npos &&
            action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
          f.check_multiplication_vars.push_back(i);
        }
      }
    }
  } else if( f.func.find("+")!=std::string::npos || f.func.find("-")!=std::string::npos ) {
    onlymultiplication=false;
  } else {
    for(unsigned i=0; i<f.var.size(); ++i) {
      if( action->getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
        f.check_multiplication_vars.push_back(i);
      }
    }
  }
  if( f.check_multiplication_vars.size()>0 ) {
    action->log.printf("  optimizing implementation as function only involves multiplication \n");
  }

  action->log.printf("  with function : %s\n",f.func.c_str());
  action->log.printf("  with variables :");
  for(unsigned i=0; i<f.var.size(); i++) {
    action->log.printf(" %s",f.var[i].c_str());
  }
  action->log.printf("\n");
  f.function.set( f.func, f.var, action );
  std::vector<double> zeros( nargs, 0 );
  double fval = fabs(f.function.evaluate(zeros));
  f.zerowhenallzero=(fval<epsilon );
  if( f.zerowhenallzero ) {
    action->log.printf("  not calculating when all arguments are zero \n");
  }
  options.derivativeZeroIfValueIsZero = f.check_multiplication_vars.size()>0;
}

std::string Custom::getFunctionString( const std::string& func ) {
  std::string outstr="", mult="*";
  for(unsigned i=0; i<func.length(); ++i) {
    if( func[i]==mult[0] ) {
      outstr += "\\*";
    } else {
      outstr += func[i];
    }
  }
  return outstr;
}

void Custom::calc( const Custom& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout ) {
  if( args.size()>1 ) {
    bool allzero=false;
    if( func.check_multiplication_vars.size()>0 ) {
      for(unsigned i=0; i<func.check_multiplication_vars.size(); ++i) {
        if( fabs(args[func.check_multiplication_vars[i]])<epsilon ) {
          allzero=true;
          break;
        }
      }
    } else if( func.zerowhenallzero ) {
      allzero=(fabs(args[0])<epsilon);
      for(unsigned i=1; i<args.size(); ++i) {
        if( fabs(args[i])>epsilon ) {
          allzero=false;
          break;
        }
      }
    }
    if( allzero ) {
      funcout.values[0]=0;
      if( noderiv ) {
        return;
      }
      for(unsigned i=0; i<args.size(); i++) {
        funcout.derivs[0][i] = 0.0;
      }
      return;
    }
  }
  funcout.values[0] = func.function.evaluate( args );
  if( !noderiv ) {
    for(unsigned i=0; i<args.size(); i++) {
      funcout.derivs[0][i] = func.function.evaluateDeriv( i, args );
    }
  }
}

}
}



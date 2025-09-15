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
#include "SphericalHarmonic.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfMatrix.h"
#include "core/ActionRegister.h"


namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR SPHERICAL_HARMONIC
/*
Calculate the values of all the spherical harmonic funtions for a particular value of l.

This action allows you to the all the [spherical harmonic](https://en.wikipedia.org/wiki/Spherical_harmonics) functions for a particular input
$L$ value.  As discussed in more detail in the article provided above the spherical harmonics this action calculates have the following general form:

$$
Y_l^m(\theta,\phi) = wNe^{im\phi} P_l^m(\cos\theta)
$$

where $N$ is a normalisation constant, $P_l^m$ is tn associated Legendre Polynomial and $w$ is an optional weight that can be passed in an input argument or simply set equal to one.
$e^{i\phi}$ and $\cos(\theta) are computed from the other input arguments, $x, y$ and $z$ as follows:

$$
e^{i\phi} = \frac{x}{r} + i\frac{y}{r} \qquad \textrm{and} \qquad \cos(\theta) = \frac{z}{r} \qquad \textrm{where} \qquad r = \sqrt{x^2 + y^2 + z^2}
$$

At present, the arguments for this action must be matrices. However, it would be easy to add functionality that would allow you to compute this function for scalar or vector input.
The arguments must all have the same shape as the two output components will also be matrices that are
calculated by applying the function above to each of the elements of the input matrix in turn.  The number of components output will be equal to $2(2L+1)$ and will contain
the real and imaginary parts of the $Y_l^m$ functions with the the $2l+1$ possible $m$ values.

The following intput provides an example that demonstrates how this function is used:

```plumed
d: DISTANCE_MATRIX GROUP=1-10 COMPONENTS
c: SPHERICAL_HARMONIC L=1 ARG=d.x,d.y,d.z
PRINT ARG=c.rm-n1 FILE=real_part_m-1
PRINT ARG=c.im-n1 FILE=imaginary_part_m-1
PRINT ARG=c.rm-0 FILE=real_part_m0
PRINT ARG=c.im-0 FILE=imaginary_part_m0
PRINT ARG=c.rm-p1 FILE=real_part_m+1
PRINT ARG=c.im-p1 FILE=imaginary_part_m+1
```

The DISTANCE_MATRIX command in the above input computes 3 $10\times10$ matrices.  These 3 $10\times10$ matrices are used in the input to the sphierical harmonic command,
which in turn outputs 6 $10\times10$ matrices that contain the real and imaginary parts when the three spherical harmonic functions with $l=1$ are applied element-wise to the above input. These six $10\times10$
matrices are then output to six separate files.

In the above example the weights for every distance is set equal to one.  The following example shows how an argument can be used to set the $w$ values to use when computing the function
above.

```plumed
s: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0}
sc: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0} COMPONENTS
c: SPHERICAL_HARMONIC L=1 ARG=sc.x,sc.y,sc.z,s
PRINT ARG=c.rm-n1 FILE=real_part_m-1
PRINT ARG=c.im-n1 FILE=imaginary_part_m-1
PRINT ARG=c.rm-0 FILE=real_part_m0
PRINT ARG=c.im-0 FILE=imaginary_part_m0
PRINT ARG=c.rm-p1 FILE=real_part_m+1
PRINT ARG=c.im-p1 FILE=imaginary_part_m+1
```

This function is used in the calculation of the Steinhardt order parameters, which are described in detail [here](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html).

*/
//+ENDPLUMEDOC

typedef function::FunctionShortcut<SphericalHarmonic> SpHarmShortcut;
PLUMED_REGISTER_ACTION(SpHarmShortcut,"SPHERICAL_HARMONIC")
typedef function::FunctionOfMatrix<SphericalHarmonic> MatrixSpHarm;
PLUMED_REGISTER_ACTION(MatrixSpHarm,"SPHERICAL_HARMONIC_MATRIX")

}
}


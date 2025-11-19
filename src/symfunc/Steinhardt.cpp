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
#include "CoordinationNumbers.h"
#include "core/ActionShortcut.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR Q1
/*
Calculate 1st order Steinhardt parameters

The 1st order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered with the atoms aranged on a line.  The Steinhardt parameter for atom, $i$ is complex vector whose components are
calculated using the following formula:

$$
q_{1m}(i) = \sum_j \sigma( r_{ij} ) Y_{1m}(\mathbf{r}_{ij})
$$

where $Y_{1m}$ is one of the 1st order spherical harmonics so $m$ is a number that runs from $-1$ to
$+1$.  The function $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  The parameters of this function should be set so that it the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html), the Steinhardt parameters can
be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

$$
Q_1(i) = \frac{1}{\sum_j \sigma(r_{ij}) } \sqrt{ \sum_{m=-1}^1 q_{1m}(i)^{*} q_{1m}(i) }
$$

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. The following input illustrates
how by averaging the value of this norm over all the atoms in the system you can measure the global degree of order in the system:

```plumed
q1: Q1 SPECIES=1-64 D_0=1.3 R_0=0.2
q1_mean: MEAN ARG=q1 PERIODIC=NO
PRINT ARG=q1_mean FILE=colvar
```

In the above input the rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
q1: Q1 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
q1_mean: MEAN ARG=q1 PERIODIC=NO
PRINT ARG=q1_mean FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atoms

The command below could be used to measure the Q1 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na$^+$ ions followed by the 64 Cl$-$ ions.  Once again the average Q1 parameter is calculated and output to a
file called colvar

```plumed
q1: Q1 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2
q1_mean: MEAN ARG=q1 PERIODIC=NO
PRINT ARG=q1_mean FILE=colvar
```

If you simply want to examine the values of the Q1 parameters for each of the atoms in your system you can do so by exploiting the
command [DUMPATOMS](DUMPATOMS.md) as shown in the example below.  The following output file will output a file in an extended xyz format
called q1.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q1 parameter, columns
6-8 will contain the real parts of the director of the $q_{1m}$ vector while columns 9-11 will contain the imaginary parts of this director.

```plumed
q1: Q1 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q1 ARG=q1 FILE=q1.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells Q1 that it is safe not to calculate the Q1 parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the Q1 parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the fccubic parameter of the atoms
cc: Q1 ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the Q1 parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.
Below is an example where these deprecated keywords are used to calculate the histogram of Q1 parameters for the 64 atoms in a box of Lennard Jones print them
to a file called colvar:

```plumed
q1: Q1 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=q1.* FILE=colvar
```

The following example illustrates how you can use VSUM to calculate a global vector of $Q_{1m}$ values as follows:

$$
Q_{1m} = \sum_i \frac{q_{1m}(i)}{\sum_j \sigma(r_{ij})}
$$

where the sum runs over all the atoms.  You can then take these $Q_{1m}$ values and compute the following norm:

$$
s = \sqrt{ \sum_{m=-1}^1 Q_{1m}^{*} Q_{1m} }
$$

The VMEAN command that is also used in the input below performs a similar operations.  The only difference is that
we divide the sums in the first expression above by the number of atoms.

```plumed
q1: Q1 SPECIES=1-64 D_0=1.3 R_0=0.2 VMEAN VSUM
PRINT ARG=q1.* FILE=colvar
```

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q3
/*
Calculate 3rd order Steinhardt parameters.

The 3rd order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, $i$ is complex vector whose components are
calculated using the following formula:

$$
q_{3m}(i) = \sum_j \sigma( r_{ij} ) Y_{3m}(\mathbf{r}_{ij})
$$

where $Y_{3m}$ is one of the 3rd order spherical harmonics so $m$ is a number that runs from $-3$ to
$+3$.  The function $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  The parameters of this function should be set so that it the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html), the Steinhardt parameters can
be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

$$
Q_3(i) = \frac{1}{\sum_j \sigma(r_{ij}) } \sqrt{ \sum_{m=-3}^3 q_{3m}(i)^{*} q_{3m}(i) }
$$

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered.  The following input illustrates
how by averaging the value of this norm over all the atoms in the system you can measure the global degree of order in the system:

```plumed
q3: Q3 SPECIES=1-64 D_0=1.3 R_0=0.2
q3_mean: MEAN ARG=q3 PERIODIC=NO
PRINT ARG=q3_mean FILE=colvar
```

In the above input the rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
q3: Q3 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
q3_mean: MEAN ARG=q3 PERIODIC=NO
PRINT ARG=q3_mean FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atoms

The command below could be used to measure the Q3 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na$^+$ ions followed by the 64 Cl$-$ ions.  Once again the average Q3 parameter is calculated and output to a
file called colvar

```plumed
q3: Q3 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2
q3_mean: MEAN ARG=q3 PERIODIC=NO
PRINT ARG=q3_mean FILE=colvar
```

If you simply want to examine the values of the Q3 parameters for each of the atoms in your system you can do so by exploiting the
command [DUMPATOMS](DUMPATOMS.md) as shown in the example below.  The following output file will output a file in an extended xyz format
called q3.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q3 parameter, columns
6-12 will contain the real parts of the director of the $q_{3m}$ vector while columns 13-19 will contain the imaginary parts of this director.

```plumed
q3: Q3 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q3 ARG=q3 FILE=q3.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells Q3 that it is safe not to calculate the Q3 parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the Q3 parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the fccubic parameter of the atoms
cc: Q3 ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the Q3 parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.
Below is an example where these deprecated keywords are used to calculate the histogram of Q3 parameters for the 64 atoms in a box of Lennard Jones print them
to a file called colvar:

```plumed
q3: Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=q3.* FILE=colvar
```

The following example illustrates how you can use VSUM to calculate a global vector of $Q_{3m}$ values as follows:

$$
Q_{3m} = \sum_i \frac{q_{3m}(i)}{\sum_j \sigma(r_{ij})}
$$

where the sum runs over all the atoms.  You can then take these $Q_{3m}$ values and compute the following norm:

$$
s = \sqrt{ \sum_{m=-3}^3 Q_{3m}^{*} Q_{3m} }
$$

The VMEAN command that is also used in the input below performs a similar operations.  The only difference is that
we divide the sums in the first expression above by the number of atoms.

```plumed
q3: Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 VMEAN VSUM
PRINT ARG=q3.* FILE=colvar
```

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q4
/*
Calculate fourth order Steinhardt parameters.

The 4th order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, $i$ is complex vector whose components are
calculated using the following formula:

$$
q_{4m}(i) = \sum_j \sigma( r_{ij} ) Y_{4m}(\mathbf{r}_{ij})
$$

where $Y_{4m}$ is one of the 4th order spherical harmonics so $m$ is a number that runs from $-4$ to
$+4$.  The function $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  The parameters of this function should be set so that it the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html), the Steinhardt parameters can
be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

$$
Q_4(i) = \frac{1}{\sum_j \sigma(r_{ij}) } \sqrt{ \sum_{m=-4}^4 q_{4m}(i)^{*} q_{4m}(i) }
$$

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. The following input illustrates
how by averaging the value of this norm over all the atoms in the system you can measure the global degree of order in the system:

```plumed
q4: Q4 SPECIES=1-64 D_0=1.3 R_0=0.2
q4_mean: MEAN ARG=q4 PERIODIC=NO
PRINT ARG=q4_mean FILE=colvar
```

In the above input the rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
q4: Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
q4_mean: MEAN ARG=q4 PERIODIC=NO
PRINT ARG=q4_mean FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atoms

The command below could be used to measure the Q4 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na$^+$ ions followed by the 64 Cl$-$ ions.  Once again the average Q4 parameter is calculated and output to a
file called colvar

```plumed
q4: Q4 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2
q4_mean: MEAN ARG=q4 PERIODIC=NO
PRINT ARG=q4_mean FILE=colvar
```

If you simply want to examine the values of the Q4 parameters for each of the atoms in your system you can do so by exploiting the
command [DUMPATOMS](DUMPATOMS.md) as shown in the example below.  The following output file will output a file in an extended xyz format
called q4.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q4 parameter, columns
6-14 will contain the real parts of the director of the $q_{4m}$ vector while columns 15-23 will contain the imaginary parts of this director.

```plumed
q4: Q4 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q4 ARG=q4 FILE=q4.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells Q4 that it is safe not to calculate the Q4 parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the Q4 parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the fccubic parameter of the atoms
cc: Q4 ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the Q4 parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.
Below is an example where these deprecated keywords are used to calculate the histogram of Q4 parameters for the 64 atoms in a box of Lennard Jones print them
to a file called colvar:

```plumed
q4: Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=q4.* FILE=colvar
```

The following example illustrates how you can use VSUM to calculate a global vector of $Q_{4m}$ values as follows:

$$
Q_{4m} = \sum_i \frac{q_{4m}(i)}{\sum_j \sigma(r_{ij})}
$$

where the sum runs over all the atoms.  You can then take these $Q_{4m}$ values and compute the following norm:

$$
s = \sqrt{ \sum_{m=-4}^1 Q_{4m}^{*} Q_{4m} }
$$

The VMEAN command that is also used in the input below performs a similar operations.  The only difference is that
we divide the sums in the first expression above by the number of atoms.

```plumed
q4: Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 VMEAN VSUM
PRINT ARG=q4.* FILE=colvar
```

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q6
/*
Calculate sixth order Steinhardt parameters.

The 6th order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, $i$ is complex vector whose components are
calculated using the following formula:

$$
q_{6m}(i) = \sum_j \sigma( r_{ij} ) Y_{6m}(\mathbf{r}_{ij})
$$

where $Y_{6m}$ is one of the 6th order spherical harmonics so $m$ is a number that runs from $-6$ to
$+6$.  The function $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  The parameters of this function should be set so that it the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html), the Steinhardt parameters can
be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

$$
Q_6(i) = \frac{1}{\sum_j \sigma(r_{ij}) } \sqrt{ \sum_{m=-6}^6 q_{6m}(i)^{*} q_{6m}(i) }
$$

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. The following input illustrates
how by averaging the value of this norm over all the atoms in the system you can measure the global degree of order in the system:

```plumed
q6: Q6 SPECIES=1-64 D_0=1.3 R_0=0.2
q6_mean: MEAN ARG=q6 PERIODIC=NO
PRINT ARG=q6_mean FILE=colvar
```

In the above input the rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
q6: Q6 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
q6_mean: MEAN ARG=q6 PERIODIC=NO
PRINT ARG=q6_mean FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atoms

The command below could be used to measure the Q6 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na$^+$ ions followed by the 64 Cl$-$ ions.  Once again the average Q6 parameter is calculated and output to a
file called colvar

```plumed
q6: Q6 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2
q6_mean: MEAN ARG=q6 PERIODIC=NO
PRINT ARG=q6_mean FILE=colvar
```

If you simply want to examine the values of the Q6 parameters for each of the atoms in your system you can do so by exploiting the
command [DUMPATOMS](DUMPATOMS.md) as shown in the example below.  The following output file will output a file in an extended xyz format
called q6.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q6 parameter, columns
6-18 will contain the real parts of the director of the $q_{6m}$ vector while columns 19-31 will contain the imaginary parts of this director.

```plumed
q6: Q6 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q6 ARG=q6 FILE=q6.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells Q6 that it is safe not to calculate the Q6 parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the Q6 parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the fccubic parameter of the atoms
cc: Q6 ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the Q6 parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.
Below is an example where these deprecated keywords are used to calculate the histogram of Q6 parameters for the 64 atoms in a box of Lennard Jones print them
to a file called colvar:

```plumed
q6: Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=q6.* FILE=colvar
```

The following example illustrates how you can use VSUM to calculate a global vector of $Q_{6m}$ values as follows:

$$
Q_{6m} = \sum_i \frac{q_{6m}(i)}{\sum_j \sigma(r_{ij})}
$$

where the sum runs over all the atoms.  You can then take these $Q_{6m}$ values and compute the following norm:

$$
s = \sqrt{ \sum_{m=-6}^6 Q_{6m}^{*} Q_{6m} }
$$

The VMEAN command that is also used in the input below performs a similar operations.  The only difference is that
we divide the sums in the first expression above by the number of atoms.

```plumed
q6: Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 VMEAN VSUM
PRINT ARG=q6.* FILE=colvar
```

*/
//+ENDPLUMEDOC

class Steinhardt : public ActionShortcut {
private:
  static std::string getSymbol( int m );
  void createVectorNormInput( const std::string& ilab,
                              const std::string& olab,
                              int l,
                              const std::string& sep,
                              const std::string& vlab,
                              bool usegpu = false);
public:
  static void registerKeywords( Keywords& keys );
  explicit Steinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Steinhardt,"Q1")
PLUMED_REGISTER_ACTION(Steinhardt,"Q3")
PLUMED_REGISTER_ACTION(Steinhardt,"Q4")
PLUMED_REGISTER_ACTION(Steinhardt,"Q6")

void Steinhardt::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.addDeprecatedFlag("LOWMEM","");
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","scalar","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","scalar","the norm of the mean vector");
  keys.needsAction("GROUP");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("SPHERICAL_HARMONIC");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("MEAN");
  keys.needsAction("SUM");
  keys.setValueDescription("vector",
                           "the norms of the vectors of spherical harmonic coefficients");
  keys.addFlag("USEGPU",false,"run part of this calculation on the GPU");
}

Steinhardt::Steinhardt( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  {
    bool lowmem;
    parseFlag("LOWMEM",lowmem);
    if( lowmem ) {
      warning("LOWMEM flag is deprecated and is no longer required for this action");
    }
  }

  bool usegpu;
  parseFlag("USEGPU",usegpu);
  const std::string doUSEGPU = usegpu?" USEGPU":"";

  std::string sp_str;
  parse("SPECIES",sp_str);
  std::string specA;
  parse("SPECIESA",specA);
  std::string specB;
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true,
                                     getShortcutLabel(),
                                     sp_str,
                                     specA,
                                     specB,
                                     this );
  int l;
  std::string sph_input = getShortcutLabel() + "_sh: SPHERICAL_HARMONIC ARG="
                          + getShortcutLabel() + "_mat.x,"
                          + getShortcutLabel() + "_mat.y,"
                          + getShortcutLabel() + "_mat.z,"
                          + getShortcutLabel() + "_mat.w";

  if( getName()=="Q1" ) {
    sph_input +=" L=1";
    l=1;
  } else if( getName()=="Q3" ) {
    sph_input += " L=3";
    l=3;
  } else if( getName()=="Q4" ) {
    sph_input += " L=4";
    l=4;
  } else if( getName()=="Q6" ) {
    sph_input += " L=6";
    l=6;
  } else {
    plumed_merror("invalid input");
  }
  readInputLine( sph_input + doUSEGPU);

  // Input for denominator (coord)
  ActionWithValue* av = plumed.getActionSet()
                        .selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_denom_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT "
                 "ARG=" + getShortcutLabel() + "_mat.w,"
                 + getShortcutLabel() + "_denom_ones"
                 + doUSEGPU);
  readInputLine( getShortcutLabel() + "_sp: MATRIX_VECTOR_PRODUCT "
                 "ARG=" + getShortcutLabel() + "_sh.*,"
                 + getShortcutLabel() + "_denom_ones"
                 + doUSEGPU);

  // If we are doing VMEAN determine sum of vector components
  std::string snum;
  bool do_vmean;
  parseFlag("VMEAN",do_vmean);
  bool do_vsum;
  parseFlag("VSUM",do_vsum);
  if( do_vmean || do_vsum ) {
    auto makeString=[&](const std::string& numstr,
    const char realImg)->std::string{
      //realImg is "r" or "i"
      return getShortcutLabel() + "_" +realImg + "mn-" + numstr + ": CUSTOM "
      "ARG=" + getShortcutLabel() + "_sp." +realImg + "m-" + numstr + ","
      + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO";
    };
    // Divide all components by coordination numbers
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine(makeString(snum,'r'));
      // Imaginary part
      readInputLine(makeString(snum,'i'));
    }
  }

  if( do_vmean ) {
    auto makeString=[&](const std::string& numstr,
    const char realImg)->std::string{
      //realImg is "r" or "i"
      return getShortcutLabel() + "_" +realImg + "ms-" + numstr + ": MEAN "
      "ARG=" + getShortcutLabel() + "_" +realImg + "mn-" + numstr
      + " PERIODIC=NO";
    };
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine(makeString(snum,'r'));
      // Imaginary part
      readInputLine(makeString(snum,'i'));
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(),
                           getShortcutLabel() + "_vmean",
                           l,
                           "_",
                           "ms",
                           usegpu );
  }
  if( do_vsum ) {
    auto makeString=[&](const std::string& numstr,
    const std::string& realImg)->std::string{
      return getShortcutLabel() + "_" +realImg + "mz-" + numstr + ": SUM "
      "ARG=" + getShortcutLabel() + "_" + realImg + "mn-" + numstr
      + " PERIODIC=NO";
    };
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine(makeString(snum,"r"));
      // Imaginary part
      readInputLine(makeString(snum,"i"));
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(),
                           getShortcutLabel() + "_vsum",
                           l,
                           "_",
                           "mz",
                           usegpu );
  }

  // Now calculate the total length of the vector
  createVectorNormInput( getShortcutLabel() + "_sp",
                         getShortcutLabel() + "_norm",
                         l,
                         ".",
                         "m",
                         usegpu );
  // And take average
  readInputLine( getShortcutLabel() + ": CUSTOM "
                 "ARG=" + getShortcutLabel() + "_norm,"
                 + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(),
      getShortcutLabel(), "", this );
}

void Steinhardt::createVectorNormInput( const std::string& ilab,
                                        const std::string& olab,
                                        const int l,
                                        const std::string& sep,
                                        const std::string& vlab,
                                        const bool usegpu) {
  std::string norm_input = olab + "2: COMBINE PERIODIC=NO POWERS=";
  std::string snum = getSymbol( -l );
  std::string arg_inp = "";

  std::string arg_inp_real = ilab + sep + "r" + vlab + "-";
  std::string arg_inp_img = ilab + sep + "i" + vlab + "-";
  std::string comma="";
  for(int i=-l; i<=l; ++i) {
    snum = getSymbol( i );
    arg_inp += comma + arg_inp_real + snum + "," + arg_inp_img + snum;
    norm_input += comma + "2,2";
    comma=",";
  }
  readInputLine( norm_input + " ARG=" + arg_inp + (usegpu?" USEGPU":"") );
  readInputLine( olab + ": CUSTOM ARG=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO");
}

std::string Steinhardt::getSymbol( const int m ) {
  if( m<0 ) {
    std::string num;
    Tools::convert( -m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num;
    Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

} // namespace symfunc
} // namespace PLMD

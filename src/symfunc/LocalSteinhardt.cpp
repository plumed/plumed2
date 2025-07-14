/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVARF LOCAL_Q1
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_1 vector on the central atom and the q_3 vector on the atoms in the first coordination sphere.

The [Q1](Q1.md) command allows one to calculate one complex vector for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average [Q1](Q1.md) parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html) a variety of variations on the Steinhardt parameters have been
introduced to better describe nucleation. That page also shows how PLUMED provides you with flexibility that you can use to implement new combinations of the
Steinhardt parameters. However, the inputs that you would need to write to implement common symmetry functions are rather complex so we also provide shortcuts
like this one to help you compute CVs that have been widely used in the literature easily.

This particular shortcut allows you to compute the LOCAL_Q1 parameter for a particular atom, which is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

$$
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-1}^1 q_{1m}^{*}(i)q_{1m}(j) }{ \sum_j \sigma( r_{ij} ) }
$$

where $q_{1m}(i)$ and $q_{1m}(j)$ are the 1st order Steinhardt vectors calculated for atom $i$ and atom $j$ respectively and the asterisk denotes complex
conjugation.  The function
$\sigma( r_{ij} )$ is a switching function that acts on the distance between atoms $i$ and $j$.  The parameters of this function should be set
so that it the function is equal to one when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms $i$ and $j$ and thus measures the degree to which the orientations of these
adjacent atoms are correlated.

The following input shows how this works in practice.  This input calculates the average value of the LOCAL_Q1 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

```plumed
q1: Q1 SPECIES=1-64 D_0=1.3 R_0=0.2
lq1: LOCAL_Q1 SPECIES=q1 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
lq1_mean: MEAN ARG=lq1 PERIODIC=NO
PRINT ARG=lq1_mean FILE=colvar
```

The following input calculates the distribution of LOCAL_Q1 parameters at any given time and outputs this information to a file.

```plumed
q1: Q1 SPECIES=1-64 D_0=1.3 R_0=0.2
lq1: LOCAL_Q1 SPECIES=q1 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=lq1.* FILE=colvar
```

The following calculates the LOCAL_Q1 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

```plumed
q1a: Q1 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2
q1b: Q1 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2
w1: LOCAL_Q1 SPECIESA=q1a SPECIESB=q1b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=w1.* FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the $s_i$ parameter that is defined by the equation above for only those atoms that
lie in a certain part of the simulation box.

```plumed
# Calculate the Q1 parameters for all the atoms
q1: Q1 SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local_q1 parameters
lq1: LOCAL_Q1 SPECIES=q1 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=lq1,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average $s_i$ parameter for those atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$. By including the MASK
keyword in the LOCAL_Q1 line we reduce the number of $s_i$ values we have to compute using the expression above. However, we are still asking PLUMED to calculate the $q_{lm}$ vectors for many atoms that
will not contribute to the final averaged quantity. The documentation for [LOCAL_AVERAGE](LOCAL_AVERAGE.md) discusses how you can create a mask vector that can act upon
the [Q1](Q1.md) action here that will ensure that you are not calculating $q_{lm}$  parameters that are not required.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q3
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_3 vector on the central atom and the q_3 vector on the atoms in the first coordination sphere.

The [Q3](Q3.md) command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average [Q3](Q3.md) parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html) a variety of variations on the Steinhardt parameters have been
introduced to better describe nucleation. That page also shows how PLUMED provides you with flexibility that you can use to implement new combinations of the
Steinhardt parameters. However, the inputs that you would need to write to implement common symmetry functions are rather complex so we also provide shortcuts
like this one to help you compute CVs that have been widely used in the literature easily.

This particular shortcut allows you to compute the LOCAL_Q3 parameter for a particular atom, which is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

$$
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-3}^3 q_{3m}^{*}(i)q_{3m}(j) }{ \sum_j \sigma( r_{ij} ) }
$$

where $q_{3m}(i)$ and $q_{3m}(j)$ are the 3rd order Steinhardt vectors calculated for atom $i$ and atom $j$ respectively and the asterisk denotes complex
conjugation.  The function
$\sigma( r_{ij} )$ is a switching function that acts on the distance between atoms $i$ and $j$.  The parameters of this function should be set
so that it the function is equal to one when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms $i$ and $j$ and thus measures the degree to which the orientations of these
adjacent atoms are correlated.

The following input shows how this works in practice.  This input calculates the average value of the LOCAL_Q3 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

```plumed
q3: Q3 SPECIES=1-64 D_0=1.3 R_0=0.2
lq3: LOCAL_Q3 SPECIES=q3 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
lq3_mean: MEAN ARG=lq3 PERIODIC=NO
PRINT ARG=lq3.mean FILE=colvar
```

The following input calculates the distribution of LOCAL_Q3 parameters at any given time and outputs this information to a file.

```plumed
q3: Q3 SPECIES=1-64 D_0=1.3 R_0=0.2
lq3: LOCAL_Q3 SPECIES=q3 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=lq3.* FILE=colvar
```

The following calculates the LOCAL_Q3 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

```plumed
q3a: Q3 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2
q3b: Q3 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2
w3: LOCAL_Q3 SPECIESA=q3a SPECIESB=q3b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=w3.* FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the $s_i$ parameter that is defined by the equation above for only those atoms that
lie in a certain part of the simulation box.

```plumed
# Calculate the Q3 parameters for all the atoms
q3: Q3 SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local_q3 parameters
lq3: LOCAL_Q3 SPECIES=q3 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=lq3,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average $s_i$ parameter for those atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$. By including the MASK
keyword in the LOCAL_Q3 line we reduce the number of $s_i$ values we have to compute using the expression above. However, we are still asking PLUMED to calculate the $q_{lm}$ vectors for many atoms that
will not contribute to the final averaged quantity. The documentation for [LOCAL_AVERAGE](LOCAL_AVERAGE.md) discusses how you can create a mask vector that can act upon
the [Q3](Q3.md) action here that will ensure that you are not calculating $q_{lm}$ vectors that are not required.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q4
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_4 vector on the central atom and the q_4 vector on the atoms in the first coordination sphere.

The [Q4](Q4.md) command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average [Q4](Q4.md) parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html) a variety of variations on the Steinhardt parameters have been
introduced to better describe nucleation. That page also shows how PLUMED provides you with flexibility that you can use to implement new combinations of the
Steinhardt parameters. However, the inputs that you would need to write to implement common symmetry functions are rather complex so we also provide shortcuts
like this one to help you compute CVs that have been widely used in the literature easily.

This particular shortcut allows you to compute the LOCAL_Q4 parameter for a particular atom, which is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

$$
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-4}^4 q_{4m}^{*}(i)q_{4m}(j) }{ \sum_j \sigma( r_{ij} ) }
$$

where $q_{4m}(i)$ and $q_{4m}(j)$ are the 4th order Steinhardt vectors calculated for atom $i$ and atom $j$ respectively and the asterisk denotes complex
conjugation.  The function
$\sigma( r_{ij} )$ is a switching function that acts on the distance between atoms $i$ and $j$.  The parameters of this function should be set
so that it the function is equal to one when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms $i$ and $j$ and thus measures the degree to which the orientations of these
adjacent atoms are correlated.

The following input shows how this works in practice.  This input calculates the average value of the LOCAL_Q4 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

```plumed
q4: Q4 SPECIES=1-64 D_0=1.3 R_0=0.2
lq4: LOCAL_Q4 SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
lq4_mean: MEAN ARG=lq4 PERIODIC=NO
PRINT ARG=lq4_mean FILE=colvar
```

The following input calculates the distribution of LOCAL_Q4 parameters at any given time and outputs this information to a file.

```plumed
q4: Q4 SPECIES=1-64 D_0=1.3 R_0=0.2
lq4: LOCAL_Q4 SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=lq4.* FILE=colvar
```

The following calculates the LOCAL_Q4 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

```plumed
q4a: Q4 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2
q4b: Q4 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2
w4: LOCAL_Q4 SPECIESA=q4a SPECIESB=q4b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=w4.* FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the $s_i$ parameter that is defined by the equation above for only those atoms that
lie in a certain part of the simulation box.

```plumed
# Calculate the Q4 parameters for all the atoms
q4: Q4 SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local_q4 parameters
lq4: LOCAL_Q4 SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=lq4,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average $s_i$ parameter for those atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$. By including the MASK
keyword in the LOCAL_Q4 line we reduce the number of $s_i$ values we have to compute using the expression above. However, we are still asking PLUMED to calculate the $q_{lm}$ vectors for many atoms that
will not contribute to the final averaged quantity. The documentation for [LOCAL_AVERAGE](LOCAL_AVERAGE.md) discusses how you can create a mask vector that can act upon
the [Q4](Q4.md) action here that will ensure that you are not calculating $q_{lm}$ vectors that are not required.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q6
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_6 vector on the central atom and the q_6 vector on the atoms in the first coordination sphere.

The [Q6](Q6.md) command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average [Q6](Q6.md) parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

As discussed on [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html) a variety of variations on the Steinhardt parameters have been
introduced to better describe nucleation. That page also shows how PLUMED provides you with flexibility that you can use to implement new combinations of the
Steinhardt parameters. However, the inputs that you would need to write to implement common symmetry functions are rather complex so we also provide shortcuts
like this one to help you compute CVs that have been widely used in the literature easily.

This particular shortcut allows you to compute the LOCAL_Q6 parameter for a particular atom, which is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

$$
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-6}^6 q_{6m}^{*}(i)q_{6m}(j) }{ \sum_j \sigma( r_{ij} ) }
$$

where $q_{6m}(i)$ and $q_{6m}(j)$ are the 6th order Steinhardt vectors calculated for atom $i$ and atom $j$ respectively and the asterisk denotes complex
conjugation.  The function
$\sigma( r_{ij} )$ is a switching function that acts on the distance between atoms $i$ and $j$.  The parameters of this function should be set
so that it the function is equal to one when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms $i$ and $j$ and thus measures the degree to which the orientations of these
adjacent atoms are correlated.

The following input shows how this works in practice.  This input calculates the average value of the LOCAL_Q6 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

```plumed
q6: Q6 SPECIES=1-64 D_0=1.3 R_0=0.2
lq6: LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
lq6_mean: MEAN ARG=lq6 PERIODIC=NO
PRINT ARG=lq6_mean FILE=colvar
```

The following input calculates the distribution of LOCAL_Q6 parameters at any given time and outputs this information to a file.

```plumed
q6: Q6 SPECIES=1-64 D_0=1.3 R_0=0.2
lq6: LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1}
PRINT ARG=lq6.* FILE=colvar
```

The following calculates the LOCAL_Q6 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

```plumed
q6a: Q6 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2
q6b: Q6 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2
w6: LOCAL_Q6 SPECIESA=q6a SPECIESB=q6b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=w6.* FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the $s_i$ parameter that is defined by the equation above for only those atoms that
lie in a certain part of the simulation box.

```plumed
# Calculate the Q6 parameters for all the atoms
q6: Q6 SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local_q6 parameters
lq6: LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=lq6,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average $s_i$ parameter for those atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$. By including the MASK
keyword in the LOCAL_Q6 line we reduce the number of $s_i$ values we have to compute using the expression above. However, we are still asking PLUMED to calculate the $q_{lm}$ vectors for many atoms that
will not contribute to the final averaged quantity. The documentation for [LOCAL_AVERAGE](LOCAL_AVERAGE.md) discusses how you can create a mask vector that can act upon
the [Q6](Q6.md) action here that will ensure that you are not calculating $q_{lm}$ vectors that are not required.

*/
//+ENDPLUMEDOC

class LocalSteinhardt : public ActionShortcut {
private:
  std::string getSymbol( const int& m ) const ;
  std::string getArgsForStack( const int& l, const std::string& lab ) const;
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalSteinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q1")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q3")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q4")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q6")

void LocalSteinhardt::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","SPECIES","the label of the action that computes the Steinhardt parameters for which you would like to calculate local steinhardt parameters");
  keys.add("optional","SPECIESA","the label of the action that computes the Steinhardt parameters for which you would like to calculate local steinhardt parameters");
  keys.add("optional","SPECIESB","the label of the action that computes the Steinhardt parameters that you would like to use when calculating the loal steinhardt parameters");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("optional","MASK","the label/s for vectors that are used to determine which local steinhardt parameters to compute");
  keys.addDeprecatedFlag("LOWMEM","");
  keys.setValueDescription("vector","the values of the local steinhardt parameters for the input atoms");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("MATRIX_PRODUCT");
  keys.needsAction("GROUP");
  keys.needsAction("ONES");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("VSTACK");
  keys.needsAction("CONCATENATE");
  keys.needsAction("CUSTOM");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.addFlag("USEGPU",false,"run part of this calculation on the GPU");
}

std::string LocalSteinhardt::getSymbol( const int& m ) const {
  if( m<0 ) {
    std::string num;
    Tools::convert( -1*m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num;
    Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

std::string LocalSteinhardt::getArgsForStack( const int& l, const std::string& sp_lab ) const {
  std::string numstr;
  Tools::convert( l, numstr );
  std::string data_mat = " ARG=" + sp_lab + "_sp.rm-n" + numstr + ","
                         + sp_lab + "_sp.im-n" + numstr;
  for(int i=-l+1; i<=l; ++i) {
    numstr = getSymbol( i );
    data_mat += "," + sp_lab + "_sp.rm-" + numstr + ","
                + sp_lab + "_sp.im-" + numstr;
  }
  return data_mat;
}

LocalSteinhardt::LocalSteinhardt(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {

  bool usegpu;
  parseFlag("USEGPU",usegpu);
  const std::string doUSEGPU = usegpu?" USEGPU":"";

#define createLabel(name) const std::string name##Lab = getShortcutLabel()+"_"#name;
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  // Get the Q value
  int l;
  Tools::convert( getName().substr(7), l);
  // Create a vector filled with ones
  std::string twolplusone;
  Tools::convert( 2*(2*l+1), twolplusone );
  createLabel( uvec );
  readInputLine( uvecLab + ": ONES SIZE=" + twolplusone );
  // Read in species keyword
  std::string sp_str;
  parse("SPECIES",sp_str);
  std::string spa_str;
  parse("SPECIESA",spa_str);
  createLabel( dpmat );
  createLabel( cmap );
  createLabel( grp );
  if( sp_str.length()>0 ) {
    // Create a group with these atoms
    readInputLine( grpLab + ": GROUP ATOMS=" + sp_str );
    std::vector<std::string> sp_lab = Tools::getWords(sp_str, "\t\n ,");
    // This creates the stash to hold all the vectors
    createLabel( uvecs );
    createLabel( nmat );
    if( sp_lab.size()==1 ) {
      // The lengths of all the vectors in a vector
      readInputLine( nmatLab + ": OUTER_PRODUCT "
                     "ARG=" + sp_lab[0] + "_norm," + uvecLab);
      // The unormalised vectors
      readInputLine( uvecsLab + ": VSTACK" + getArgsForStack( l, sp_lab[0] ) );
    } else {
      createLabel( mags );
      std::string len_vec = magsLab + ": CONCATENATE ARG=" + sp_lab[0] + "_norm";
      for(unsigned i=1; i<sp_lab.size(); ++i) {
        len_vec += "," + sp_lab[i] + "_norm";
      }
      // This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = uvecsLab + ": CONCATENATE";
      for(unsigned i=0; i<sp_lab.size(); ++i) {
        std::string snum;
        Tools::convert( i+1, snum );
        concat_str += " MATRIX" + snum + "1=" + uvecsLab + snum;
        readInputLine( uvecsLab + snum + ": VSTACK" + getArgsForStack( l, sp_lab[i] ) );
      }
      // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
      readInputLine( nmatLab + ": OUTER_PRODUCT ARG=" + magsLab + "," + uvecLab);
      // The unormalised vectors
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    createLabel( vecs );
    readInputLine( vecsLab + ": CUSTOM ARG=" + uvecsLab + ","
                   + nmatLab + " FUNC=x/y PERIODIC=NO");
    // And transpose the matrix
    createLabel( vecsT );
    readInputLine( vecsTLab + ": TRANSPOSE ARG=" + vecsLab );
    std::string sw_str;
    parse("SWITCH",sw_str);
    std::string maskstr;
    parse("MASK",maskstr);
    if( maskstr.length()>0 ) {
      maskstr=" MASK=" + maskstr;
    }
    readInputLine( cmapLab + ": CONTACT_MATRIX GROUP=" + sp_str + " "
                   "SWITCH={" + sw_str + "}" + maskstr + doUSEGPU);
    // And the matrix of dot products
    readInputLine( dpmatLab + ": MATRIX_PRODUCT ARG=" + vecsLab + ","
                   + vecsTLab + " MASK=" + cmapLab + doUSEGPU);
  } else if( spa_str.length()>0 ) {
    // Create a group with these atoms
    readInputLine( grpLab + ": GROUP ATOMS=" + spa_str );
    std::string spb_str;
    parse("SPECIESB",spb_str);
    if( spb_str.length()==0 ) {
      plumed_merror("need both SPECIESA and SPECIESB in input");
    }
    const std::vector<std::string> sp_laba = Tools::getWords(spa_str, "\t\n ,");
    const std::vector<std::string> sp_labb = Tools::getWords(spb_str, "\t\n ,");
    createLabel(nmatA);
    createLabel(uvecsA);
    if( sp_laba.size()==1 ) {
      // The matrix that is used for normalising
      readInputLine( nmatALab + ": OUTER_PRODUCT "
                     "ARG=" +  sp_laba[0] + "_norm," + uvecLab);
      // The unormalised vectors
      readInputLine( uvecsALab + ": VSTACK" + getArgsForStack( l, sp_laba[0] ) );
    } else {
      createLabel(magsA);
      std::string len_vec = magsALab + ": CONCATENATE "
                            "ARG=" + sp_laba[0] + "_norm";
      for(unsigned i=1; i<sp_laba.size(); ++i) {
        len_vec += "," + sp_laba[i] + "_norm";
      }
      //  This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = uvecsALab + ": CONCATENATE";
      for(unsigned i=0; i<sp_laba.size(); ++i) {
        std::string snum;
        Tools::convert( i+1, snum );
        concat_str += " MATRIX" + snum + "1=" + uvecsALab + snum;
        readInputLine( uvecsALab + snum + ": VSTACK" + getArgsForStack( l, sp_laba[i] ) );
      }
      // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
      readInputLine( nmatALab + ": OUTER_PRODUCT ARG=" + magsALab + "," + uvecLab);
      // The unormalised vector
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    createLabel( vecsA );
    readInputLine( vecsALab + ": CUSTOM ARG=" + uvecsALab + ","
                   + nmatALab + " FUNC=x/y PERIODIC=NO");
    // Now do second matrix
    createLabel(nmatB);
    createLabel(uvecsBT);
    createLabel(uvecsB);
    if( sp_labb.size()==1 ) {
      readInputLine( nmatBLab + ": OUTER_PRODUCT ARG=" +  uvecLab + ","
                     + sp_labb[0] + "_norm");
      readInputLine( uvecsBTLab + ": VSTACK" + getArgsForStack( l, sp_labb[0] ) );
      readInputLine( uvecsBLab + ": TRANSPOSE ARG=" + uvecsBTLab);
    } else {
      createLabel(magsB);
      std::string len_vec = magsBLab + ": CONCATENATE ARG=" +  sp_labb[0] + "_norm";
      for(unsigned i=1; i<sp_labb.size(); ++i) {
        len_vec += "," + sp_labb[i] + "_norm";
      }
      //  This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = uvecsBLab + ": CONCATENATE";
      for(unsigned i=0; i<sp_labb.size(); ++i) {
        std::string snum;
        Tools::convert( i+1, snum );
        concat_str += " MATRIX1" + snum + "=" + uvecsBLab + snum;
        readInputLine( uvecsBTLab + snum + ": VSTACK" + getArgsForStack( l, sp_labb[i] ) );
        readInputLine( uvecsBLab + snum + ": TRANSPOSE ARG=" + uvecsBTLab + snum );
      }
      // And the normalising matrix
      readInputLine( nmatBLab + ": OUTER_PRODUCT ARG=" + uvecLab + "," + magsBLab);
      // The unormalised vectors
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    createLabel(vecsB);
    readInputLine( vecsBLab + ": CUSTOM ARG=" + uvecsBLab + ","
                   + nmatBLab + " FUNC=x/y PERIODIC=NO");
    std::string sw_str;
    parse("SWITCH",sw_str);
    std::string maskstr;
    parse("MASK",maskstr);
    if( maskstr.length()>0 ) {
      maskstr=" MASK=" + maskstr;
    }
    readInputLine( cmapLab + ": CONTACT_MATRIX GROUPA=" + spa_str
                   + " GROUPB=" + spb_str + " SWITCH={" + sw_str + "}" + maskstr + doUSEGPU);
    readInputLine( dpmatLab + ": MATRIX_PRODUCT ARG=" + vecsALab + "," + vecsBLab
                   + " MASK=" + cmapLab + doUSEGPU);
  }

  // Now create the product matrix
  createLabel( prod );
  readInputLine( prodLab + ": CUSTOM ARG=" + cmapLab + "," + dpmatLab + " FUNC=x*y PERIODIC=NO");
  // Now the sum of coordination numbers times the switching functions
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( cmapLab);
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  createLabel( ones );
  readInputLine( onesLab + ": ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT "
                 "ARG=" + prodLab +"," + getShortcutLabel() +"_ones" + doUSEGPU);
  // And just the sum of the coordination numbers
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT "
                 "ARG=" + cmapLab + "," + getShortcutLabel() +"_ones" + doUSEGPU);
  // And matheval to get the final quantity
  createLabel( av );
  readInputLine( avLab + ": CUSTOM ARG=" + getShortcutLabel() + ","
                 + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), avLab, "", this );
#undef createLabel
}

}
}

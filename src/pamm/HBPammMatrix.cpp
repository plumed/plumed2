/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "adjmat/AdjacencyMatrixBase.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/IFile.h"

//+PLUMEDOC MATRIX HBPAMM_MATRIX
/*
Adjacency matrix in which two electronegative atoms are adjacent if they are hydrogen bonded

This method allows you to calculate an adjacency matrix and use all the methods discussed in the documentation for
the [CONTACT_MATRIX](CONTACT_MATRIX.md) upon it to define CVs.  The $i,j$ element of the matrix that is calculated
by this action is one if there is a hydrogen bond connecting atom $i$ to atom $j$.  Furthermore, we determine whether
there is a hydrogen atom between these two atoms by using the PAMM technique that is discussed in the articles from the
bibliography below and in the documentation for the [PAMM](PAMM.md) action.

The example shown below illustrates how the method is used in practice

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
...
PRINT ARG=m FILE=colvar
```

The input above is outputting the full hbpamm matrix.  However, this action is perhaps more usefully used to investigate the connectivity between
a collection of water molecules in liquid water. Importantly, however,
the output matrix here is __not__ symmetric.  We thus calculate the number of hydrogen bonds that these atoms are donating using the
following input:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
...
ones: ONES SIZE=64
rsums: MATRIX_VECTOR_PRODUCT ARG=m,ones
DUMPATOMS ATOMS=1-192:3 ARG=rsums FILE=donors.xyz
```

To calculate the number of hydrogen bonds these atoms accept we would use the following input:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
...
mT: TRANSPOSE ARG=m
ones: ONES SIZE=64
rsums: MATRIX_VECTOR_PRODUCT ARG=mT,ones
DUMPATOMS ATOMS=1-192:3 ARG=rsums FILE=acceptors.xyz
```

To calculate the total number of hydorgen bonds these atoms participate in you would use an input like this:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
...
mT: TRANSPOSE ARG=m
hbmat: CUSTOM ARG=m,mT FUNC=x+y PERIODIC=NO
ones: ONES SIZE=64
rsums: MATRIX_VECTOR_PRODUCT ARG=hbmat,ones
DUMPATOMS ATOMS=1-192:3 ARG=rsums FILE=hbonds.xyz
```

If you want to investigate whether there are hydrogen bonds between two groups of molecules you can use an input like this:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUPA=1 GROUPB=2-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
...
PRINT ARG=m FILE=colvar
```

This input outputs a $1\times 63$ matrix in which the $1,i$th element tells you whether or not atom 1 donates a hydrogen bond
to the $i$th element in the group of 63 atoms that was specified using the ACCEPTORS keyword.  The $i,1$th element of the
transpose of this matrix tells you if the $i$th center donates a hydrogen bond to atom 1.

In general, it is better to use this action through the [HBPAMM_SA](HBPAMM_SA.md), [HBPAMM_SD](HBPAMM_SD.md) and [HBPAMM_SH](HBPAMM_SH.md)
keywords, which can be used to calculate the number of hydrogen bonds each donor, acceptor or hydrogen atom in your system participates in.

## Periodic boundary conditions

Notice that in all the inputs above the distances values that enter the pamm expressions are calculated in a way that takes the
periodic boundary conditions into account.  If you want to ignore the periodic boundary conditions you can use the NOPBC flag as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
m: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  NOPBC
...
```

## COMPONENTS flag

If you add the flag COMPONENTS to the input as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
c4: HBPAMM_MATRIX ...
  GROUP=1-192:3 GROUPC=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  COMPONENTS
...
```

then four matrices with the labels `c4.w`, `c4.x`, `c4.y` and `c4.z` are output by the action. The matrix with the label `c4.w` is the adjacency matrix
that would be output if you had not added the COMPONENTS flag. The $i,j$ component of the matrices `c4.x`, `c4.y` and `c4.z` contain the $x$, $y$ and $z$
components of the vector connecting atoms $j$ and $k$. Importantly, however, the components of these vectors are only stored in `c4.x`, `c4.y` and `c4.z`
if the elements of `c4.w` are non-zero. Using the COMPONENTS flag in this way ensures that you can use HBPAMM_MATRIX in tandem with many of the functionalities
that are part of the [symfunc module](module_symfunc.md).  Remember, however, that the $i,j$ element of the HBPAMM_MATRIX is only non-zero if atom $i$ donates
a hydrogen bond to atom $j$.  __You cannot use HBPAMM_MATRIX to identify the set of atoms that each atom is hydrogen bonded to.__

## The MASK keyword

You use the MASK keyword with HBPAMM_MATRIX in the same way that is used in [CONTACT_MATRIX](CONTACT_MATRIX.md).  This keyword thus expects a vector in input,
which tells HBOND_MATRIX that it is safe to not calculate certain rows of the output matrix.  An example where this keyword is used is shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
# The atoms that are of interest
ow: GROUP ATOMS=1-1650
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=ow CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculates cooordination numbers
cmap: HBPAMM_MATRIX ...
  GROUP=ow GROUPC=1650-3000
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  MASK=sphere
...
ones: ONES SIZE=1650
cc: MATRIX_VECTOR_PRODUCT ARG=cmap,ones
# Multiply coordination numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=prod,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average number of hydrogen bonds each of the atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$ donate.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SA
/*
Calculate the number of hydrogen bonds each acceptor participates in using the HBPamm method

This shortcut action allows you to calculate the number of hydrogen bonds each of the atoms specified using the SITES
keyword accepts from its neighbours. The number of hydrogen bonds that a particular site accepts is determined by using the
PAMM tehcnique that is discussed in the articles from the bibliography below and in the documentation for the [PAMM](PAMM.md) action

The following example shows how you can use this action to calculate how many hydrogen bonds each of the water moecules in
a box of water accepts from its neighbours.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
acceptors: HBPAMM_SA ...
  SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  REGULARISE=0.001 GAUSS_CUTOFF=6.25
...
DUMPATOMS ARG=acceptors ATOMS=1-192:3 FILE=acceptors.xyz
```

The output here is an xyz with five columns. As explained in the documentation for [DUMPATOMS](DUMPATOMS.md), the first four
columns are the usual columns that you would expect in an xyz file.  The fifth column then contains the number of hydrogen bonds
that have been accepted by each of the atoms.

In the example input above all the atoms specified using the SITE keyword can both accept and donate hydrogen bonds.  If the group
of atoms that can accept hydrogen bonds is different from the group of atoms that can donate hydrogen bonds then you use the ACCEPTORS and
DONORS keyword as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
acceptors: HBPAMM_SA ...
  ACCEPTORS=1-9:3 DONORS=1-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  HYDROGENS=2-192:3,3-192:3
...
DUMPATOMS ARG=acceptors ATOMS=1-9:3 FILE=acceptors.xyz
```

This input will still output an xyz file with five columns. The fifth column in this file will give the number of hydrogen bonds that each
of the three atoms that were specified using the ACCEPTORS keyword accept from the atoms that were specified by the DONORS keyword.  Notice also that the
as the REGULARISE and GAUSS_CUTOFF keywords were not used in this input these quantities are set to the default values that are given in the table below.

## The NOPBC flag

When you use the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) you need to calculate the distances between various pairs of atoms.  In all the inputs that have been provided
above we calculate these distances between pairs of atoms in a way that takes the periodic boundary conditions into account.  If you want to ignore the periodic
boundary conditions you can use the NOPBC flag as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
acceptors: HBPAMM_SA ...
   SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
   CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
   NOPBC
...
DUMPATOMS ARG=acceptors ATOMS=1-192:3 FILE=acceptors.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the number of hydrogen bonds.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average number of hydrogen bonds those atoms that lie in a certain part of the simulation box accept from their neighbors.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-192:3 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the number of hydrogen bonds that have been accepted
acceptors: HBPAMM_SA ...
   SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
   CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
   MASK=sphere
...
# Multiply number of hydrogen bonds by sphere vector
prod: CUSTOM ARG=acceptors,sphere FUNC=x*y PERIODIC=NO
# Total number of hydrogen bonds accepted for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average number of accepted hydrogen bonds for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the number of hydrogen bonds that are accepted by those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Optimisation details

If you expand the inputs above you will see that they all the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) action.  The matrix that this action computes is sparse
as each atom only forms hydroen bonds with a relatively small number of neighbors.  The vast majority of the elements in the [HBPAMM_MATRIX](HBPAMM_MATRIX.md)
are thsu zero.  To reduce the amount of memory that PLUMED requires PLUMED uses sparse matrix storage.  Consequently, whenever you calculate and store a contact
matrix only the elements of the matrix that are non-zero are stored.

We also use the sparsity of the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) to make the time required to compute the matrix scale linearly rather than quadratically with
the number of atoms. We know that element $i,j$ is only non-zero if the donor and acceptor atoms are within:

$$
r_c = \textrm{max}_i \left( \mu_i^{(ad)} + \sqrt(2 g_{cut}) | \sqrt{\lambda_{i,max}} v_\textrm{i,max}^{(ad)} | \right)
$$

where $i$ runs over all the PAMM kernels.  $g_{cut}$ is the gaussian cutoff that is specified using the REGULARIZE
keyword,  $\mu_i^{(ad)}$ is the mean acceptor-donor distance for the $i$th kernel, $\lambda_{i,max}$ is the largest
eigenvalue of the covariance for the $i$th kernel and $v_\textrm{i,max}^{(ad)}$ is the acceptor-donor distance component
of its corresponding eigenvector.  We can determine that many pairs of atoms are further appart than $r_c$ without computing the
distance between these atoms by using divide and conquer strategies such as linked lists and neighbour lists. Furthermore, we do not even
need to set a `D_MAX` parameter as we would have to if we were using [CONTACT_MATRIX](CONTACT_MATRIX.md) as the above cutoff can be determined
from the parameters of the Gaussian kernels directly.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SD
/*
Calculate the number of hydrogen bonds each donor participates in using the HBPamm method

This shortcut action allows you to calculate the number of hydrogen bonds each of the atoms specified using the SITES
keyword donates to its neighbours. The number of hydrogen bonds that a particular site donates is determined by using the
PAMM tehcnique that is discussed in the articles from the bibliography below and in the documentation for the [PAMM](PAMM.md) action

The following example shows how you can use this action to calculate how many hydrogen bonds each of the water moecules in
a box of water donates to its neighbours.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
donors: HBPAMM_SD ...
  SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  REGULARISE=0.001 GAUSS_CUTOFF=6.25
...
DUMPATOMS ARG=donors ATOMS=1-192:3 FILE=donors.xyz
```

The output here is an xyz with five columns. As explained in the documentation for [DUMPATOMS](DUMPATOMS.md), the first four
columns are the usual columns that you would expect in an xyz file.  The fifth column then contains the number of hydrogen bonds
that have been donated by each of the atoms.

In the example input above all the atoms specified using the SITE keyword can both accept and donate hydrogen bonds.  If the group
of atoms that can accept hydrogen bonds is different from the group of atoms that can donate hydrogen bonds then you use the ACCEPTORS and
DONORS keyword as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
donors: HBPAMM_SD ...
  DONORS=1-9:3 ACCEPTORS=1-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  HYDROGENS=2-192:3,3-192:3
...
DUMPATOMS ARG=donors ATOMS=1-9:3 FILE=donors.xyz
```

This input will still output an xyz file with five columns. The fifth column in this file will give the number of hydrogen bonds that each
of the three atoms that were specified using the DONORS keyword donate to atoms that were specified by the ACCEPTORS keyword.  Notice also that the
as the REGULARISE and GAUSS_CUTOFF keywords were not used in this input these quantities are set to the default values that are given in the table below.

## The NOPBC flag

When you use the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) you need to calculate the distances between various pairs of atoms.  In all the inputs that have been provided
above we calculate these distances between pairs of atoms in a way that takes the periodic boundary conditions into account.  If you want to ignore the periodic
boundary conditions you can use the NOPBC flag as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
donors: HBPAMM_SD ...
   SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
   CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
   NOPBC
...
DUMPATOMS ARG=donors ATOMS=1-192:3 FILE=donors.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the number of hydrogen bonds.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average number of hydrogen bonds those atoms that lie in a certain part of the simulation box donoate to their neighbors.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-192:3 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the number of hydrogen bonds that have been donated
donors: HBPAMM_SD ...
   SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
   CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
   MASK=sphere
...
# Multiply number of hydrogen bonds by sphere vector
prod: CUSTOM ARG=donors,sphere FUNC=x*y PERIODIC=NO
# Total number of hydrogen bonds donated for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average number of donated hydrogen bonds for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the number of hydrogen bonds that have been donated by those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Optimisation details

If you expand the inputs above you will see that they all the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) action.  The matrix that this action computes is sparse
as each atom only forms hydroen bonds with a relatively small number of neighbors.  The vast majority of the elements in the [HBPAMM_MATRIX](HBPAMM_MATRIX.md)
are thsu zero.  To reduce the amount of memory that PLUMED requires PLUMED uses sparse matrix storage.  Consequently, whenever you calculate and store a contact
matrix only the elements of the matrix that are non-zero are stored.

We also use the sparsity of the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) to make the time required to compute the matrix scale linearly rather than quadratically with
the number of atoms. We know that element $i,j$ is only non-zero if the donor and acceptor atoms are within:

$$
r_c = \textrm{max}_i \left( \mu_i^{(ad)} + \sqrt(2 g_{cut}) | \sqrt{\lambda_{i,max}} v_\textrm{i,max}^{(ad)} | \right)
$$

where $i$ runs over all the PAMM kernels.  $g_{cut}$ is the gaussian cutoff that is specified using the REGULARIZE
keyword,  $\mu_i^{(ad)}$ is the mean acceptor-donor distance for the $i$th kernel, $\lambda_{i,max}$ is the largest
eigenvalue of the covariance for the $i$th kernel and $v_\textrm{i,max}^{(ad)}$ is the acceptor-donor distance component
of its corresponding eigenvector.  We can determine that many pairs of atoms are further appart than $r_c$ without computing the
distance between these atoms by using divide and conquer strategies such as linked lists and neighbour lists. Furthermore, we do not even
need to set a `D_MAX` parameter as we would have to if we were using [CONTACT_MATRIX](CONTACT_MATRIX.md) as the above cutoff can be determined
from the parameters of the Gaussian kernels directly.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SH
/*
Calculate the number of hydrogen bonds each hydrogen participates in using the HBPamm method

This shortcut action allows you to calculate the number of hydrogen bonds each of the atoms specified using the SITES
keyword donates to its neighbours. The number of hydrogen bonds that a particular site donates is determined by using the
PAMM tehcnique that is discussed in the articles from the bibliography below and in the documentation for the [PAMM](PAMM.md) action

The following example shows how you can use this action to calculate how many hydrogen bonds each of the water moecules in
a box of water donates to its neighbours.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
hyd: HBPAMM_SH ...
  SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  REGULARISE=0.001 GAUSS_CUTOFF=6.25
...
DUMPATOMS ARG=hyd ATOMS=2-192:3,3-192:3 FILE=hydrogens.xyz
```

The output here is an xyz with five columns. As explained in the documentation for [DUMPATOMS](DUMPATOMS.md), the first four
columns are the usual columns that you would expect in an xyz file.  The fifth column then contains the number of hydrogen bonds
that each hydrogen atom participates in.

In the example input above all the atoms specified using the SITE keyword can both accept and donate hydrogen bonds.  If the group
of atoms that can accept hydrogen bonds is different from the group of atoms that can donate hydrogen bonds then you use the ACCEPTORS and
DONORS keyword as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
hyd: HBPAMM_SH ...
  ACCEPTORS=1-9:3 DONORS=1-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  HYDROGENS=2-192:3,3-192:3
...
DUMPATOMS ARG=hyd ATOMS=2-192:3,3-192:3 FILE=hydrogens.xyz
```

This input will still output an xyz file with five columns. The fifth column in this file will give the number of hydrogen bonds that each of the hydrogen
atoms that were specified in this input participate within.  For this atom for a hydrogen atom to particpiate in a hydrogen bond the second closest atom to
it must be either atom 1, 4 or 7 as these are the only atoms that can accept hydrogen bonds according to the above input.  Notice also that the
as the REGULARISE and GAUSS_CUTOFF keywords were not used in this input these quantities are set to the default values that are given in the table below.

## The NOPBC flag

When you use the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) you need to calculate the distances between various pairs of atoms.  In all the inputs that have been provided
above we calculate these distances between pairs of atoms in a way that takes the periodic boundary conditions into account.  If you want to ignore the periodic
boundary conditions you can use the NOPBC flag as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
hyd: HBPAMM_SH ...
  SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
  CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
  NOPBC
...
DUMPATOMS ARG=hyd ATOMS=2-192:3,3-192:3 FILE=hydrogens.xyz
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the number of hydrogen bonds.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average number of hydrogen bonds those hydrogen atoms that lie in a certain part of the simulation box partipipate in.

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-hbpamm/b3lyp.pamm
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if hydrogen i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=2-192:3,3-192:3 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the number of bonds each hydrogen participates in
hyd: HBPAMM_SH ...
   SITES=1-192:3 HYDROGENS=2-192:3,3-192:3
   CLUSTERS=regtest/pamm/rt-hbpamm/b3lyp.pamm
   MASK=sphere
...
# Multiply number of hydrogen bonds by sphere vector
prod: CUSTOM ARG=hyd,sphere FUNC=x*y PERIODIC=NO
# Total number of hydrogn bonds that the atoms in the sphere of interest participate in
numer: SUM ARG=prod PERIODIC=NO
# Number of hydrogen atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average number of hydrogen bonds each hydrogen atoms in sphere of interest participates in
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average number of hydrogen bonds each of the hydrogen atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$ participates in.

## Optimisation details

If you expand the inputs above you will see that they all the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) action.  The matrix that this action computes is sparse
as each atom only forms hydroen bonds with a relatively small number of neighbors.  The vast majority of the elements in the [HBPAMM_MATRIX](HBPAMM_MATRIX.md)
are thsu zero.  To reduce the amount of memory that PLUMED requires PLUMED uses sparse matrix storage.  Consequently, whenever you calculate and store a contact
matrix only the elements of the matrix that are non-zero are stored.

We also use the sparsity of the [HBPAMM_MATRIX](HBPAMM_MATRIX.md) to make the time required to compute the matrix scale linearly rather than quadratically with
the number of atoms. We know that element $i,j$ is only non-zero if the donor and acceptor atoms are within:

$$
r_c = \textrm{max}_i \left( \mu_i^{(ad)} + \sqrt(2 g_{cut}) | \sqrt{\lambda_{i,max}} v_\textrm{i,max}^{(ad)} | \right)
$$

where $i$ runs over all the PAMM kernels.  $g_{cut}$ is the gaussian cutoff that is specified using the REGULARIZE
keyword,  $\mu_i^{(ad)}$ is the mean acceptor-donor distance for the $i$th kernel, $\lambda_{i,max}$ is the largest
eigenvalue of the covariance for the $i$th kernel and $v_\textrm{i,max}^{(ad)}$ is the acceptor-donor distance component
of its corresponding eigenvector.  We can determine that many pairs of atoms are further appart than $r_c$ without computing the
distance between these atoms by using divide and conquer strategies such as linked lists and neighbour lists. Furthermore, we do not even
need to set a `D_MAX` parameter as we would have to if we were using [CONTACT_MATRIX](CONTACT_MATRIX.md) as the above cutoff can be determined
from the parameters of the Gaussian kernels directly.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace pamm {

class HBPammMatrix {
public:
  double regulariser;
  Tensor incoord_to_hbcoord;
  std::vector<double> weight;
  std::vector<Vector> centers;
  std::vector<Tensor> kmat;
  static void registerKeywords( Keywords& keys );
  void parseInput( adjmat::AdjacencyMatrixBase<HBPammMatrix>* action );
  HBPammMatrix& operator=( const HBPammMatrix& m ) {
    regulariser = m.regulariser;
    incoord_to_hbcoord = m.incoord_to_hbcoord;
    weight = m.weight;
    centers = m.centers;
    kmat = m.kmat;
    return *this;
  }
  static void calculateWeight( const HBPammMatrix& data, const adjmat::AdjacencyMatrixInput& input, adjmat::MatrixOutput output );
};

typedef adjmat::AdjacencyMatrixBase<HBPammMatrix> hbpmap;
PLUMED_REGISTER_ACTION(hbpmap,"HBPAMM_MATRIX")

void HBPammMatrix::registerKeywords( Keywords& keys ) {
  keys.use("GROUPC");
  keys.add("compulsory","ORDER","dah","the order in which the groups are specified in the input.  Can be dah (donor/acceptor/hydrogens), "
           "adh (acceptor/donor/hydrogens) or hda (hydrogens/donor/acceptors");
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the kernels for PAMM");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.add("compulsory","GAUSS_CUTOFF","6.25","the cutoff at which to stop evaluating the kernel function is set equal to sqrt(2*x)*(max(adc)+cov(adc))");
  keys.needsAction("PAMM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.addDOI("10.1063/1.4900655");
  keys.addDOI("10.1021/acs.jctc.7b00993");
}

void HBPammMatrix::parseInput( adjmat::AdjacencyMatrixBase<HBPammMatrix>* action ) {
  double DP2CUTOFF;
  action->parse("GAUSS_CUTOFF",DP2CUTOFF);
  std::string sorder;
  action->parse("ORDER",sorder);
  if( sorder=="dah" ) {
    incoord_to_hbcoord(0,0)=1;
    incoord_to_hbcoord(0,1)=-1;
    incoord_to_hbcoord(0,2)=0;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=1;
    incoord_to_hbcoord(1,2)=0;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=0;
    incoord_to_hbcoord(2,2)=1;
    action->log.printf("  GROUPA is list of donor atoms \n");
  } else if( sorder=="adh" ) {
    incoord_to_hbcoord(0,0)=-1;
    incoord_to_hbcoord(0,1)=1;
    incoord_to_hbcoord(0,2)=0;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=1;
    incoord_to_hbcoord(1,2)=0;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=0;
    incoord_to_hbcoord(2,2)=1;
    action->log.printf("  GROUPA is list of acceptor atoms \n");
  } else if( sorder=="hda" ) {
    incoord_to_hbcoord(0,0)=-1;
    incoord_to_hbcoord(0,1)=0;
    incoord_to_hbcoord(0,2)=1;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=0;
    incoord_to_hbcoord(1,2)=1;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=1;
    incoord_to_hbcoord(2,2)=0;
    action->log.printf("  GROUPA is list of hydrogen atoms \n");
  } else {
    plumed_error();
  }
  // Read in the regularisation parameter
  action->parse("REGULARISE",regulariser);

  // Read in the kernels
  double sqr2pi = sqrt(2*pi);
  double sqrt2pi3 = sqr2pi*sqr2pi*sqr2pi;
  std::string fname;
  action->parse("CLUSTERS", fname);
  double sfmax=0, ww;
  Vector cent;
  Tensor covar;
  IFile ifile;
  ifile.open(fname);
  ifile.allowIgnoredFields();
  while(true) {
    if( !ifile.scanField("height",ww) ) {
      break;
    }
    ifile.scanField("ptc",cent[0]);
    ifile.scanField("ssc",cent[1]);
    ifile.scanField("adc",cent[2]);
    ifile.scanField("sigma_ptc_ptc",covar[0][0]);
    ifile.scanField("sigma_ptc_ssc",covar[0][1]);
    ifile.scanField("sigma_ptc_adc",covar[0][2]);
    covar[1][0] = covar[0][1];
    ifile.scanField("sigma_ssc_ssc",covar[1][1]);
    ifile.scanField("sigma_ssc_adc",covar[1][2]);
    covar[2][0] = covar[0][2];
    covar[2][1] = covar[1][2];
    ifile.scanField("sigma_adc_adc",covar[2][2]);
    weight.push_back( ww / ( sqrt2pi3 * sqrt(covar.determinant()) ) );
    centers.push_back( cent );
    kmat.push_back( covar.inverse() );

    Vector eigval;
    Tensor eigvec;
    diagMatSym( covar, eigval, eigvec );
    unsigned ind_maxeval=0;
    double max_eval=eigval[0];
    for(unsigned i=1; i<3; ++i) {
      if( eigval[i]>max_eval ) {
        max_eval=eigval[i];
        ind_maxeval=i;
      }
    }
    double rcut = cent[2] + sqrt(2.0*DP2CUTOFF)*fabs(sqrt(max_eval)*eigvec(2,ind_maxeval));
    if( rcut > sfmax ) {
      sfmax = rcut;
    }
    ifile.scanField();
  }
  ifile.close();
  action->setLinkCellCutoff( false, sfmax );
}

void HBPammMatrix::calculateWeight( const HBPammMatrix& data, const adjmat::AdjacencyMatrixInput& input, adjmat::MatrixOutput output ) {
  Vector ddik, ddin, in_dists, hb_pamm_dists, hb_pamm_ders, real_ders;
  ddin = input.pos;
  in_dists[2] = ddin.modulo();
  if( in_dists[2]<epsilon ) {
    return;
  }

  output.val[0]=0;
  Vector disp, der, tmp_der;
  for(unsigned i=0; i<input.natoms; ++i) {
    Vector ddij( input.extra_positions[i][0], input.extra_positions[i][1], input.extra_positions[i][2] );
    in_dists[0] = ddij.modulo();
    ddik = input.pbc->distance( input.pos, ddij );
    in_dists[1] = ddik.modulo();
    if( in_dists[1]<epsilon ) {
      continue;
    }

    hb_pamm_dists = matmul( data.incoord_to_hbcoord, in_dists );
    disp = hb_pamm_dists - data.centers[0];
    der = matmul( data.kmat[0], disp );
    double vv = data.weight[0]*exp( -dotProduct( disp, der ) / 2. );
    der *= -vv;

    double denom = data.regulariser + vv;
    for(unsigned j=0; j<3; ++j) {
      hb_pamm_ders[j] = der[j];
    }
    for(unsigned k=1; k<data.weight.size(); ++k) {
      disp = hb_pamm_dists - data.centers[k];
      tmp_der = matmul( data.kmat[k], disp );
      double tval = data.weight[k]*exp( -dotProduct( disp, tmp_der ) / 2. );
      denom += tval;
      hb_pamm_ders += -tmp_der*tval;
    }
    double vf = vv / denom;
    output.val[0] += vf;
    if( fabs(vf)<epsilon ) {
      continue;
    }
    // Now get derivatives
    real_ders = matmul( der / denom - vf*hb_pamm_ders/denom, data.incoord_to_hbcoord );

    // And add the derivatives to the underlying atoms
    Vector d1 = -(real_ders[0]/in_dists[0])*ddij - (real_ders[2]/in_dists[2])*ddin;
    output.deriv[0] += d1[0];
    output.deriv[1] += d1[1];
    output.deriv[2] += d1[2];
    Vector d2 = -(real_ders[1]/in_dists[1])*ddik + (real_ders[2]/in_dists[2])*ddin;
    output.deriv[3] += d2[0];
    output.deriv[4] += d2[1];
    output.deriv[5] += d2[2];
    Vector d3 = (real_ders[0]/in_dists[0])*ddij + (real_ders[1]/in_dists[1])*ddik;
    output.deriv[6+i*3+0] = d3[0];
    output.deriv[6+i*3+1] = d3[1];
    output.deriv[6+i*3+2] = d3[2];
    Tensor vir = -(real_ders[0]/in_dists[0])*Tensor( ddij, ddij )
                 -(real_ders[1]/in_dists[1])*Tensor( ddik, ddik )
                 -(real_ders[2]/in_dists[2])*Tensor( ddin, ddin );
    output.deriv[6 + 3*input.natoms + 0] += vir[0][0];
    output.deriv[6 + 3*input.natoms + 1] += vir[0][1];
    output.deriv[6 + 3*input.natoms + 2] += vir[0][2];
    output.deriv[6 + 3*input.natoms + 3] += vir[1][0];
    output.deriv[6 + 3*input.natoms + 4] += vir[1][1];
    output.deriv[6 + 3*input.natoms + 5] += vir[1][2];
    output.deriv[6 + 3*input.natoms + 6] += vir[2][0];
    output.deriv[6 + 3*input.natoms + 7] += vir[2][1];
    output.deriv[6 + 3*input.natoms + 8] += vir[2][2];
  }
}

class HBPammShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  HBPammShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SD")
PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SA")
PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SH")

void HBPammShortcut::registerKeywords( Keywords& keys ) {
  adjmat::AdjacencyMatrixBase<HBPammMatrix>::registerKeywords( keys );
  keys.remove("GROUP");
  keys.remove("GROUPA");
  keys.remove("GROUPB");
  keys.remove("GROUPC");
  keys.remove("ORDER");
  keys.remove("COMPONENTS");
  keys.reset_style("NL_CUTOFF","hidden");
  keys.reset_style("NL_STRIDE","hidden");
  keys.add("atoms","SITES","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
           "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified"
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","DONORS","The list of atoms which can donate a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","ACCEPTORS","The list of atoms which can accept a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list, "
           "an index range or by using a \\ref GROUP");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("HBPAMM_MATRIX");
  keys.setValueDescription("vector","a vector specifiying the number of hydrogen bonds each of the specified atoms participates within");
}

HBPammShortcut::HBPammShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string mwords = getShortcutLabel() + "_mat: HBPAMM_MATRIX";
  if( getName()=="HBPAMM_SD" ) {
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUP=" + site_str;
    } else {
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPA=" + d_str;
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPB=" + a_str;
    }
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPC=" + h_str + " ORDER=dah";
  } else if( getName()=="HBPAMM_SA" ) {
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUP=" + site_str;
    } else {
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPA=" + a_str;
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPB=" + d_str;
    }
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPC=" + h_str + " ORDER=adh";
  } else if( getName()=="HBPAMM_SH" ) {
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPA=" + h_str + " ORDER=hda";
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUPB=" + site_str;
      mwords += " GROUPC=" + site_str;
    } else {
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPB=" + d_str;
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPC=" + a_str;
    }
  }
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  readInputLine( mwords + " " + convertInputLineToString() );
  ActionWithValue* mb=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( mb );
  std::string nsize;
  Tools::convert( (mb->copyOutput(getShortcutLabel() + "_mat"))->getShape()[1], nsize );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + nsize );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_ones");
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}

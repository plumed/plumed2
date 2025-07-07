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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/RMSD.h"
#include "tools/PDB.h"

namespace PLMD {
namespace colvar {

class RMSD : public Colvar {

  bool squared;
  bool nopbc;
  PLMD::RMSD myrmsd;
  std::vector<Vector> der;
public:
  explicit RMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

//+PLUMEDOC DCOLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.

One colvar that has been shown to be very successful in studying protein folding is the distance between the instantaneous configuration
and a reference configuration - often the structure of the folded state.  When the free energy of a protein is shown as a function
of this collective variable there is a minima for low values of the CV, which is due to the folded state of the protein.  There is
then a second minima at higher values of the CV, which is the minima corresponding to the unfolded state.  The aim of this colvar is
thus to calculate something like:

$$
d(X,X') = \vert X-X' \vert
$$

where $X$ is the instantaneous position of all the atoms in the system and
$X'$ is the positions of the atoms in some reference structure that is provided as input.
$d(X,X')$ thus measures the distance all the atoms have moved away from this reference configuration.
Oftentimes, it is only the internal motions of the structure - i.e. not the translations of the center of
mass or the rotations of the reference frame - that are interesting.  Hence, when calculating
the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some way. At present PLUMED provides two distinct ways
of performing this superposition.  The first method is applied when you use TYPE=SIMPLE in the input
line.  This instruction tells PLUMED that the root mean square deviation is to be calculated after the
positions of the geometric centers in the reference and instantaneous configurations are aligned.  In
other words $d(X,x')$ is to be calculated using:

$$
d(X,X') = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j} ( X _{i,\alpha}-com _\alpha(X)-{X'} _{i,\alpha}+com _\alpha(X') )^2 }
$$

with

$$
com_\alpha (X) = \sum_i \frac{w'_{i}}{\sum_j w'_j} X _{i,\alpha}
$$

and

$$
com_\alpha(X')= \sum_i  \frac{w'_{i}}{\sum_j w'_j}X' _{i,\alpha}
$$

Obviously, $com_\alpha(X)$ and  $com_\alpha(X')$  represent the positions of the center of mass in the reference
and instantaneous configurations if the weights $w'$ are set equal to the atomic masses.  If the weights are all set equal to
one, however, $com_\alpha(X)$ and  $com_\alpha(X')$ are the positions of the geometric centers.

An example input that can be used to calculate and print this RMSD distance is shown below:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=SIMPLE REFERENCE=regtest/basic/rt19/test0.pdb
PRINT ARG=rmsd0 FILE=colvar
```

Notice that there are sets of weights: $w'$ and $w$ in the formulas above. The first of these weights is used to calculate the position of the center of mass
(so it determines how the atoms are _aligned_).  The second set of weights, $w$ is used when calculating how far the atoms have been
_displaced_. These weights are assigned in the reference configuration that you provide as input (i.e. the appear in the input file
to this action that you set using REFERENCE=whatever.pdb). As you can see in the input above, this input consists of a simple pdb file
containing the set of atoms for which you want to calculate the RMSD displacement and their positions in the reference configuration.
It is important to note that the indices in this pdb need to be set correctly.  The indices in this file determine the indices of the
instantaneous atomic positions that are used by PLUMED when calculating this colvar.  As such if you want to calculate the RMSD distance
moved by the first, fourth, sixth and twenty eighth atoms in the MD codes input file then the indices of the corresponding reference positions in this pdb
file should be set equal to 1, 4, 6 and 28.

The pdb input file should also contain the weights $w$ and $w'$. These second of these sets of weights, $w'$, appears in the OCCUPANCY column (the first column after the coordinates).
These are the weights that are used to calculate the position of the center of mass.  The BETA column (the second column
after the Cartesian coordinates) is used to provide the $w$ values which are used in the the calculation of the displacement.
Please note that it is possible to use fractional values for beta and for the occupancy. However, we recommend you only do this when
you really know what you are doing however as the results can be rather strange.  A more common practise is to use different sets of atoms
for the alignment and the displacement.  You can do this by setting the $w$ and $w'$ values for all the atoms you wish to use for alignment only equal to zero and
one respectively and by setting the $w$ and $w'$ values for all the atoms you wish to use for displacement only to one and zero respectively.

In the PDB input files that you use for RMSD the atomic coordinates and box lengths should be in Angstroms unless
you are working with natural units.  If you are working with natural units then the coordinates
should be in your natural length unit.  You can find more details on the PDB file format [here](http://www.wwpdb.org/docs.html).
Please make sure your PDB file is correctly formatted.  More detail on the format for PDB files can be found in the documentation for the [PDB2CONSTANT](PDB2CONSTANT.md) action.

The following input uses a different method to calculate the RMSD distance as you can tell from the TYPE=OPTIMAL on the input line.  In addition, because we have added the
SQURED flag on the input line, we are calculating $d(X,X')^2$ rather than $d(X,X')$.  Calculating $d(X,X')^2$ is slightly less computationally expensive than computing $d(X,X')$ as
you avoid a square root operation.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=OPTIMAL SQUARED REFERENCE=regtest/basic/rt19/test0.pdb
PRINT ARG=rmsd0 FILE=colvar
```

In this case  the root mean square deviation is calculated after the positions of geometric centers in the reference and instantaneous configurations are aligned AND after
an optimal alignment of the two frames is performed so that motion due to rotation of the reference frame between the two structures is
removed.  The equation for $d(X,X')$ in this case is:

$$
d(X,X') = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j}[ X_{i,\alpha}-com_\alpha(X)- \sum_\beta M(X,X',w')_{\alpha,\beta}({X'}_{i,\beta}-com_\beta(X')) ]^2 }
$$

where $M(X,X',w')$ is the optimal alignment matrix which is calculated using the Kearsley algorithm that is described in the paper from the bibliography below.  Again different sets of
weights are used for the alignment ($w'$) and for the displacement calculations ($w$).
This gives a great deal of flexibility as it allows you to use a different sets of atoms (which may or may not overlap) for the alignment and displacement
parts of the calculation. This may be very useful when you want to calculate how a ligand moves about in a protein cavity as you can use the protein as a reference
system and do no alignment of the ligand.

(Note: when this form of RMSD is used to calculate the secondary structure variables ([ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md)
all the atoms in the segment are assumed to be part of both the alignment and displacement sets and all weights are set equal to one)

## Numerical derivatives

If you are calculating the RMSD distance from a single RMSD reference frame you calculate numerical derivatives of the RMSD distance as shown below for the simple metric:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=SIMPLE REFERENCE=regtest/basic/rt19/test0.pdb
rmsd0n: RMSD TYPE=SIMPLE NUMERICAL_DERIVATIVES REFERENCE=regtest/basic/rt19/test0.pdb
DUMPDERIVATIVES ARG=rmsd0,rmsd0n FILE=deriv
```

and as shown below for the OPTIMAL metric:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=OPTIMAL REFERENCE=regtest/basic/rt19/test0.pdb
rmsd0n: RMSD TYPE=OPTIMAL NUMERICAL_DERIVATIVES REFERENCE=regtest/basic/rt19/test0.pdb
DUMPDERIVATIVES ARG=rmsd0,rmsd0n FILE=deriv
```

In practice there is no reason to use the NUMERICAL_DERIVATIVES option unless you are testing the RMSD implementation in the way the above two inputs are doing by printing
a file that allows us to compare the analytic and numerical derivatives.  Evaluating the derivtives numerically is __much__ more expensive than calculating them analytically.

## Computing RMSD displacements

The $d(X,X')$ values that are calculated when you use the TYPE=SIMPLE and TYPE=OPTIMAL variants of RMSD are scalars. These scalars tell you the length of the vector of displacements,
$X - X'$, between the instantaneous and reference positions.  If you would like to access this vector of displacements instead of its length you can by using the DISPLACEMENT keyword
as shown below for TYPE=SIMPLE

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=SIMPLE DISPLACEMENT REFERENCE=regtest/basic/rt19/test0.pdb
PRINT ARG=rmsd0.disp,rmsd0.dist FILE=colvar
```

or as shown below for TYPE=OPTIMAL

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=OPTIMAL DISPLACEMENT REFERENCE=regtest/basic/rt19/test0.pdb
PRINT ARG=rmsd0.disp,rmsd0.dist FILE=colvar
```

The RMSD command for these inputs output two components

- `dist` - the length of displacement vector that is output when you don't use the DISPLACEMENT keyword
- `disp` - the 3N dimensional vector of atomic dispacements, where N is the number of atoms.

These vectors of displacements are used if you use the [PCAVARS](PCAVARS.md) action to compute the projection of the displacement on a particular reference vector.

## Computing multiple RMSD values

You can also define multiple reference configurations in the reference input as is done in the following example:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt39/all1.pdb
rmsd: RMSD TYPE=OPTIMAL REFERENCE=regtest/mapping/rt39/all1.pdb
PRINT ARG=rmsd FILE=colvar
```

The output from RMSD in this case is a vector that contains the RMSD distances from each of the reference configurations in your input file.  This feature is used in the
[PATH](PATH.md) shortcut.  Furthermore, you can use the DISPLACEMENT keyword when there are multiple reference configurations in the input file as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt39/all1.pdb
rmsd: RMSD TYPE=OPTIMAL DISPLACEMENT REFERENCE=regtest/mapping/rt39/all1.pdb
PRINT ARG=rmsd.disp,rmsd.dist FILE=colvar
```

The RMSD command here still outputs the `dist` and `disp` components but now `dist` is a vector and `disp` is a matrix.  This type of command is used to implement
the [GEOMETRIC_PATH](GEOMETRIC_PATH.md) shortcut.  For this command you need information on the distances from a set of reference configurations that can be found
in the `dist` component as well as the information on the displacement vectors between the instantaneous position and each of the reference configurations that is contained in the `dist` matrix.

Please note that there are a number of other methods for calculating the distance between the instantaneous configuration and a reference configuration that are available in plumed.

## A note on periodic boundary conditions

When periodic boundary conditions are used, the atoms should be
in the proper periodic image. This has been done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding molecules using a procedure
that is equivalent to that done in [WHOLEMOLECULES](WHOLEMOLECULES.md). Notice that
rebuilding is local to this action. This is different from [WHOLEMOLECULES](WHOLEMOLECULES.md)
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
rmsd0: RMSD TYPE=SIMPLE NOPBC REFERENCE=regtest/basic/rt19/test0.pdb
PRINT ARG=rmsd0 FILE=colvar
```

In that case you need to take care that atoms are in the correct
periodic image.

## A final thought

Notice that there are many other ways of calculating the distance from a particular reference structure. These normally work by computing the square root of the sum of the squares
of the differences in collective variable values. To compute these distances you use the functionality in the [refdist](module_refdist.md) module.

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(RMSD,"RMSD_SCALAR")

void RMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("RMSD");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
  keys.setValueDescription("scalar","the RMSD between the instantaneous structure and the reference structure that was input");
}

RMSD::RMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  squared(false),
  nopbc(false) {
  std::string reference;
  parse("REFERENCE",reference);
  std::string type;
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);
  parseFlag("NOPBC",nopbc);
  checkRead();

  addValueWithDerivatives();
  setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + reference );
  }
  myrmsd.set( pdb, type, true, true );

  std::vector<AtomNumber> atoms( pdb.getAtomNumbers() );
  requestAtoms( atoms );
  der.resize( atoms.size() );

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  with indices : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) {
      log<<"\n";
    }
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared) {
    log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  }
  if(nopbc) {
    log.printf("  without periodic boundary conditions\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }
}


// calculator
void RMSD::calculate() {
  if(!nopbc) {
    makeWhole();
  }
  double r=myrmsd.calculate( getPositions(), der, squared );

  setValue(r);
  for(unsigned i=0; i<getNumberOfAtoms(); i++) {
    setAtomsDerivatives( i, der[i] );
  }
  setBoxDerivativesNoPbc();
}

}
}




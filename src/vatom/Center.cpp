/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include <cmath>
#include <limits>

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM CENTER
/*
Calculate the center for a group of atoms, with arbitrary weights.

The position of the center ${r}_{\rm COM}$ given by:

$$
{r}_{\rm COM}=\frac{\sum_{i=1}^{n} {r}_i\ w_i }{\sum_i^{n} w_i}
$$

In these expressions $r_i$ indicates the position of atom $i$ and $w_i$ is the weight for atom $i$. The following input
shows how you can calculate the expressions for a set of atoms by using PLUMED:

```plumed
# This action calculates the position of the point on the line connecting atoms 1 and 10 that is
# twice as far atom 10 as it is from atom 1
c1: CENTER ATOMS=1,1,10
# this is another way of calculating this position
c1bis: CENTER ATOMS=1,10 WEIGHTS=2,1

DUMPATOMS ATOMS=c1,c1bis FILE=atoms.xyz
```

Notice that center's position is stored as [a virtual atom](specifying_atoms.md). The positions of the centers in the above input
used in the DUMPATOMS command by using the labels for the CENTER actions. Notice, furthermore, that the mass and charge of this new center
are equal to the sums of the masses and charges of the input atoms.

The input below shows how you can use the CENTER action in place of the [COM](COM.md) action to calculate the center of mass for a group of atoms.

```plumed
c: CENTER ATOMS=1-5 MASS
```

Center is more powerful than COM because you can use arbitrary vectors of weights as in the first example above or vector of weights that are calculated by
another action as has been done in the input below.

```plumed
fcc: FCCUBIC SPECIES=1-1000 SWITCH={RATIONAL D_0=3.0 R_0=1.5}
sfcc: MORE_THAN ARG=fcc SWITCH={RATIONAL R_0=0.5}
c: CENTER ATOMS=1-1000 WEIGHTS=sfcc
DUMPATOMS ATOMS=c FILE=atom.xyz
```

This input assumes you have a cluster of solid atoms in a liquid. The actions with labels `fcc` and `sfcc`
are used to differentiate between atoms in solid-like and liquid-like atoms. `sfcc` is thus a vector with
one element for each atom. These elements are equal to one if the environment around the corresponding atom
are solid like and zero if the environment around the atom is liquid like.

## A note on periodic boundary conditions

If you run with periodic boundary conditions
these are taken into account automatically when computing the center of mass. The way this is
handled is akin to the way molecules are rebuilt in the [WHOLEMOLECULES](WHOLEMOLECULES.md) command.
However, at variance to [WHOLEMOLECULES](WHOLEMOLECULES.md), the copies of the atomic positions in this
action are modified.  The global positions (i.e. those that are used in all other actions)
are not changed when the alignment is performed.

If you believe that PBC should not be applied when calculating the position fo the center of mass you can use the
NOPBC flag as shown below:

```plumed
c: CENTER ATOMS=1-100 NOPBC
```

An additional way of managing periodic boundary conditions is offered in CENTER by using the PHASES keyword as shown
in the example input below

```plumed
c: CENTER ATOMS=1-100 PHASES
```

The scaled value for the $x$ component of the position of the center is calculated from the scaled components of the input atoms, $x_i$,
using the following expression when the PHASES option is employed

$$
x_\textrm{com} = \frac{1}{2\pi} \arctan\left( \frac{\sum_{i=1}^n w_i \sin(2 \pi x_i ) }{ \sum_{i=1}^n w_i \cos(2 \pi x_i ) } \right)
$$

Similar, expressions are used to calculae the values of the scaled $y$ and $z$ components.  The final cartesian coordinates of the center are then computed
by multiplying these scaled components by the cell vectors.  Notice that by construction this center position is
not invariant with respect to rotations of the atoms at fixed cell lattice.
In addition, for symmetric Bravais lattices, it is not invariant with respect
to special symmetries. E.g., if you have an hexagonal cell, the center will
not be invariant with respect to rotations of 120 degrees.
On the other hand, it might make the treatment of PBC easier in difficult cases.

As an alternative to PHASES you can use the SAFE_PHASES flag as shown below:

```plumed
c: CENTER ATOMS=1-100 SAFE_PHASES
```

This option will use the method described above for the PHASES keyword if the cell coordinates are set.  If, however, the cell coordinates are not set the
position of the center will be calculated in the usual way.  We wouldn't recommend using this option in place of PHASES as additional computational overhead
is introduced.  Furthermore, you normally know if the cell parameters are not passed to PLUMED in advance of doing the calculation

*/
//+ENDPLUMEDOC

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The following example input shows how you can use this shortcut to calculate the center of
mass for atoms  1,2,3,4,5,6,7 and for atoms 15,20.  You can then compute the distance between
the two center of masses.

```plumed
c1: COM ATOMS=1-7
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
```

The center of mass is stored as a [virtual atom](specifying_atoms.md).  As you can see from the
above to refer to the position of the center of mass when specifying the atoms that should be used
when calculating some other action you use the lable of the COM action that computes the center of mass
of interest.

The COM command is a shortcut because it is a wrapper to [CENTER](CENTER.md). CENTER is more powerful than
comm as it allows you to use arbitrary weights in place of the masses.

## Periodic boundary conditions

If you run with periodic boundary conditions
these are taken into account automatically when computing the center of mass. By default the way this is
handled is akin to the way molecules are rebuilt in the [WHOLEMOLECULES](WHOLEMOLECULES.md) command.
However, at variance to [WHOLEMOLECULES](WHOLEMOLECULES.md), the copies of the atomic positions in this
action are modified.  The global positions (i.e. those that are used in all other actions)
are not changed when the alignment is performed.

If you believe that PBC should not be applied when calculating the position fo the center of mass you can use the
NOPBC flag as shown below:

```plumed
c1: COM ATOMS=1-7 NOPBC
```

Alternatively, if you would like to calculate the position of the center of mass using the PHASES method that is dicussed
in the documentation for [CENTER](CENTER.md) you can add the keyword PHASES as shwo below:

```plumed
c1: COM ATOMS=1-7 PHASES
```

## The mass and charge of the vatom

By default the mass and charge of the COM are calculated by adding the masses and charges of the input atoms together.  In other
words, if your COM is calculated from the positions of $N$ atoms, the mass of the vatom is calculated as:

$$
m_\textrm{com} = \sum_{i=1}^N m_i
$$

where $m_i$ is the mass of the $i$th input atom. The charge is then calculated as:

$$
q_\textrm{com} = \sum_{i=1}^N q_i
$$

where $q_i$ is the charge of the $i$th atom.  If for any reason you don't want the mass and charge of the VATOM to be set equal to these
values you can set them manually by using the `SET_CHARGE` and `SET_MASS` keywords as shown below:

```plumed
c1: COM ATOMS=1-7 SET_MASS=1 SET_CHARGE=-1
```

*/
//+ENDPLUMEDOC


class Center:
  public ActionWithVirtualAtom {
  std::vector<double> weights;
  std::vector<Tensor> dcenter_sin;
  std::vector<Tensor> dcenter_cos;
  std::vector<Tensor> deriv;
  bool isChargeSet_;
  bool isMassSet_;
  bool weight_mass;
  bool nopbc;
  bool first;
  bool phases;
public:
  explicit Center(const ActionOptions&ao);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Center,"CENTER_FAST")
PLUMED_REGISTER_ACTION(Center,"COM")

void Center::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  if( keys.getDisplayName()!="COM" ) {
    keys.setDisplayName("CENTER");
    keys.add("optional","WEIGHTS","Center is computed as a weighted average.");
    keys.addFlag("MASS",false,"If set center is mass weighted");
  }
  keys.add("optional","SET_CHARGE","Set the charge of the virtual atom to a given value.");
  keys.add("optional","SET_MASS","Set the mass of the virtual atom to a given value.");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("PHASES",false,"Compute center using trigonometric phases");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  isChargeSet_(false),
  isMassSet_(false),
  weight_mass(false),
  nopbc(false),
  first(true),
  phases(false) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) {
    error("at least one atom should be specified");
  }
  if( getName()!="COM") {
    parseVector("WEIGHTS",weights);
    parseFlag("MASS",weight_mass);
  }
  parseFlag("NOPBC",nopbc);
  parseFlag("PHASES",phases);
  double charge_=std::numeric_limits<double>::lowest();
  parse("SET_CHARGE",charge_);
  setCharge(charge_);
  if(charge_!=std::numeric_limits<double>::lowest()) {
    isChargeSet_=true;
  }
  double mass_=-1;
  parse("SET_MASS",mass_);
  setMass(mass_);
  if(mass_>0.) {
    isMassSet_=true;
  }
  if(mass_==0.) {
    error("SETMASS must be greater than 0");
  }
  if( getName()=="COM") {
    weight_mass=true;
  }
  checkRead();
  log.printf("  of atoms:");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) {
      log<<"\n";
    }
    log.printf(" %d",atoms[i].serial());
  }
  log<<"\n";
  if(weight_mass) {
    log<<"  mass weighted\n";
    if(weights.size()!=0) {
      error("WEIGHTS and MASS keywords cannot not be used simultaneously");
    }
  } else {
    if( weights.size()==0) {
      log<<" using the geometric center\n";
      weights.resize( atoms.size() );
      for(unsigned i=0; i<atoms.size(); i++) {
        weights[i] = 1.;
      }
    } else {
      log<<" with weights:";
      if( weights.size()!=atoms.size() ) {
        error("number of elements in weight vector does not match the number of atoms");
      }
      for(unsigned i=0; i<weights.size(); ++i) {
        if(i%25==0) {
          log<<"\n";
        }
        log.printf(" %f",weights[i]);
      }
      log.printf("\n");
    }
  }
  if(phases) {
    log<<"  Phases will be used to take into account PBC\n";
  } else if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms);
}

void Center::calculate() {
  Vector pos;
  const bool dophases=(getPbc().isSet() ? phases : false);

  if(!nopbc && !dophases) {
    makeWhole();
  }

  if( first ) {
    if( weight_mass ) {
      for(unsigned i=0; i<getNumberOfAtoms(); i++) {
        if(std::isnan(getMass(i))) {
          error(
            "You are trying to compute a CENTER or COM but masses are not known.\n"
            "        If you are using plumed driver, please use the --mc option"
          );
        }
      }
    }
    double mass(0.0);
    for(unsigned i=0; i<getNumberOfAtoms(); i++) {
      mass+=getMass(i);
    }
    if( chargesWereSet && !isChargeSet_) {
      double charge(0.0);
      for(unsigned i=0; i<getNumberOfAtoms(); i++) {
        charge+=getCharge(i);
      }
      setCharge(charge);
    } else if( !isChargeSet_ ) {
      setCharge(0.0);
    }
    if(!isMassSet_) {
      setMass(mass);
    }

    if( weight_mass ) {
      weights.resize( getNumberOfAtoms() );
      for(unsigned i=0; i<weights.size(); i++) {
        weights[i] = getMass(i) / mass;
      }
    } else {
      double wtot=0.0;
      for(unsigned i=0; i<weights.size(); i++) {
        wtot+=weights[i];
      }
      for(unsigned i=0; i<weights.size(); i++) {
        weights[i]=weights[i]/wtot;
      }
      first=false;
    }
  }

  deriv.resize(getNumberOfAtoms());

  if(dophases) {
    dcenter_sin.resize(getNumberOfAtoms());
    dcenter_cos.resize(getNumberOfAtoms());
    Vector center_sin;
    Vector center_cos;
    Tensor invbox2pi=2*pi*getPbc().getInvBox();
    Tensor box2pi=getPbc().getBox() / (2*pi);
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      double w=weights[i];

      // real to scaled
      const Vector scaled=matmul(getPosition(i),invbox2pi);
      const Vector ccos(
        w*std::cos(scaled[0]),
        w*std::cos(scaled[1]),
        w*std::cos(scaled[2])
      );
      const Vector csin(
        w*std::sin(scaled[0]),
        w*std::sin(scaled[1]),
        w*std::sin(scaled[2])
      );
      center_cos+=ccos;
      center_sin+=csin;
      for(unsigned l=0; l<3; l++)
        for(unsigned k=0; k<3; k++) {
          // k over real coordinates
          // l over scaled coordinates
          dcenter_sin[i][l][k]=ccos[l]*invbox2pi[k][l];
          dcenter_cos[i][l][k]=-csin[l]*invbox2pi[k][l];
        }
    }
    const Vector c(
      std::atan2(center_sin[0],center_cos[0]),
      std::atan2(center_sin[1],center_cos[1]),
      std::atan2(center_sin[2],center_cos[2])
    );

    // normalization is convenient for doing derivatives later
    for(unsigned l=0; l<3; l++) {
      double norm=1.0/(center_sin[l]*center_sin[l]+center_cos[l]*center_cos[l]);
      center_sin[l]*=norm;
      center_cos[l]*=norm;
    }

    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      Tensor dd;
      for(unsigned l=0; l<3; l++)
        for(unsigned k=0; k<3; k++) {
          // k over real coordinates
          // l over scaled coordinates
          dd[l][k]= (center_cos[l]*dcenter_sin[i][l][k] - center_sin[l]*dcenter_cos[i][l][k]);
        }
      // scaled to real
      deriv[i]=matmul(dd,box2pi);
    }
    setAtomsDerivatives(deriv);
    // scaled to real
    setPosition(matmul(c,box2pi));
  } else {
    for(unsigned i=0; i<getNumberOfAtoms(); i++) {
      double w=weights[i];
      pos+=w*getPosition(i);
      deriv[i]=w*Tensor::identity();
    }
    setPosition(pos);
    setAtomsDerivatives(deriv);
  }
}

}
}

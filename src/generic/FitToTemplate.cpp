/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "tools/Vector.h"
#include "tools/Matrix.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "tools/RMSD.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "core/PbcAction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"

#include <vector>
#include <string>
#include <memory>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC FIT_TO_TEMPLATE
/*
This action is used to align a molecule to a template.

This can be used to move the coordinates stored in plumed
so as to be aligned with a provided template in PDB format. Pdb should contain
also weights for alignment (see the format of PDB files used e.g. for [RMSD](RMSD.md)).
Make sure your PDB file is correctly formatted as explained in the documentation for [PDB2CONSTANT](PDB2CONSTANT.md).
Weights for displacement are ignored, since no displacement is computed here.
Notice that all atoms (not only those in the template) are aligned.
To see what effect try
the [DUMPATOMS](DUMPATOMS.md) directive to output the atomic positions.

Also notice that PLUMED propagate forces correctly so that you can add a bias on a CV computed
after alignment. For many CVs this has no effect, but in some case the alignment can
change the result. Examples are:
- [POSITION](POSITION.md) CV since it is affected by a rigid shift of the system.
- [DISTANCE](DISTANCE.md) CV with COMPONENTS. Since the alignment could involve a rotation (with TYPE=OPTIMAL) the actual components could be different
  from the original ones.
- [CELL](CELL.md) components for a similar reason.
- [DISTANCE](DISTANCE.md) from a [FIXEDATOM](FIXEDATOM.md), provided the fixed atom is introduced _after_ the FIT_TO_TEMPLATE action.

!!! caution ""

    The implementation of TYPE=OPTIMAL is available but should be considered in testing phase. Please report any
    strange behavior.

!!! caution "directive modifies stored positions"

    This directive modifies the stored position at the precise moment
    it is executed. This means that only collective variables
    which are below it in the input script will see the corrected positions.
    As a general rule, put it at the top of the input file. Also, unless you
    know exactly what you are doing, leave the default stride (1), so that
    this action is performed at every MD step.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecules using a procedure
that is equivalent to that done in [WHOLEMOLECULES](WHOLEMOLECULES.md) . Notice that
rebuilding is local to this action. This is different from [WHOLEMOLECULES](WHOLEMOLECULES.md)
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

## Examples

Align the atomic position to a template then print them.
The following example is only translating the system so as
to align the center of mass of a molecule to the one in the reference
structure `ref.pdb`:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
# dump coordinates before fitting, to see the difference:
DUMPATOMS FILE=dump-before.xyz ATOMS=1-20

# fit coordinates to ref.pdb template
# this is a "TYPE=SIMPLE" fit, so that only translations are used.
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=regtest/basic/rt63/align.pdb TYPE=SIMPLE

# dump coordinates after fitting, to see the difference:
DUMPATOMS FILE=dump-after.xyz ATOMS=1-20
```

The following example instead performs a rototranslational fit.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
# dump coordinates before fitting, to see the difference:
DUMPATOMS FILE=dump-before.xyz ATOMS=1-20

# fit coordinates to ref.pdb template
# this is a "TYPE=OPTIMAL" fit, so that rototranslations are used.
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=regtest/basic/rt63/align.pdb TYPE=OPTIMAL

# dump coordinates after fitting, to see the difference:
DUMPATOMS FILE=dump-after.xyz ATOMS=1-20
```

In the following example you see two completely equivalent way
to restrain an atom close to a position that is defined in the reference
frame of an aligned molecule. You could for instance use this command to calculate the
position of the center of mass of a ligand after having aligned the atoms to the reference
frame of the protein that is determined by aligning the atoms in the protein to the coordinates
provided in the file ref.pdb

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
# center of the ligand:
center: CENTER ATOMS=100-110

FIT_TO_TEMPLATE REFERENCE=regtest/basic/rt63/align.pdb TYPE=OPTIMAL

# place a fixed atom in the protein reference coordinates:
fix: FIXEDATOM AT=1.0,1.1,1.0

# take the distance between the fixed atom and the center of the ligand
d: DISTANCE ATOMS=center,fix

# apply a restraint
RESTRAINT ARG=d AT=0.0 KAPPA=100.0
```

Notice that you could have obtained an (almost) identical result by adding a fictitious
atom to `ref.pdb` with the serial number corresponding to the atom labelled `center` (there is no automatic way
to get it, but in this example it should be the number of atoms of the system plus one),
and properly setting the weights for alignment and displacement in [RMSD](RMSD.md).
There are two differences to be expected:
(ab) FIT_TO_TEMPLATE might be slower since it has to rototranslate all the available atoms and
(b) variables employing periodic boundary conditions (such as [DISTANCE](DISTANCE.md) without `NOPBC`, as in the example above)
  are allowed after FIT_TO_TEMPLATE, whereas [RMSD](RMSD.md) expects the issues related to the periodic boundary conditions to be already solved.
The latter means that before the [RMSD](RMSD.md) statement one should use [WRAPAROUND](WRAPAROUND.md) or [WHOLEMOLECULES](WHOLEMOLECULES.md) to properly place
the ligand.


*/
//+ENDPLUMEDOC


class FitToTemplate:
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithValue {
  std::string type;
  bool nopbc;
  std::vector<double> weights;
  std::vector<std::pair<std::size_t,std::size_t> > p_aligned;
  Vector center;
  Vector shift;
  // optimal alignment related stuff
  std::unique_ptr<PLMD::RMSD> rmsd;
  Tensor rotation;
  Matrix< std::vector<Vector> > drotdpos;
  // not used anymore (see notes below at doNotRetrieve())
  // std::vector<Vector> positions;
  std::vector<Vector> DDistDRef;
  std::vector<Vector> ddistdpos;
  std::vector<Vector> centeredpositions;
  Vector center_positions;
  // Copy of the box value
  Value* boxValue;
  PbcAction* pbc_action;
public:
  explicit FitToTemplate(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  bool actionHasForces() override {
    return true;
  }
  void calculate() override;
  void apply() override;
  unsigned getNumberOfDerivatives() override {
    plumed_merror("You should not call this function");
  };
};

PLUMED_REGISTER_ACTION(FitToTemplate,"FIT_TO_TEMPLATE")

void FitToTemplate::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.setValueDescription("scalar","the RMSD distance");
}

FitToTemplate::FitToTemplate(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  nopbc(false) {
  std::string reference;
  parse("REFERENCE",reference);
  type.assign("SIMPLE");
  parse("TYPE",type);

  parseFlag("NOPBC",nopbc);
// if(type!="SIMPLE") error("Only TYPE=SIMPLE is implemented in FIT_TO_TEMPLATE");

  checkRead();

  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + reference );
  }

  requestAtoms(pdb.getAtomNumbers());
  log.printf("  found %zu atoms in input \n",pdb.getAtomNumbers().size());
  log.printf("  with indices : ");
  for(unsigned i=0; i<pdb.getAtomNumbers().size(); ++i) {
    if(i%25==0) {
      log<<"\n";
    }
    log.printf("%d ",pdb.getAtomNumbers()[i].serial());
  }
  log.printf("\n");

  std::vector<Vector> myPositions=pdb.getPositions();
  weights=pdb.getOccupancy();
  std::vector<AtomNumber> aligned=pdb.getAtomNumbers();
  p_aligned.resize( aligned.size() );
  for(unsigned i=0; i<aligned.size(); ++i) {
    p_aligned[i] = getValueIndices( aligned[i] );
  }


  // normalize weights
  double n=0.0;
  for(unsigned i=0; i<weights.size(); ++i) {
    n+=weights[i];
  }
  if(n==0.0) {
    error("PDB file " + reference + " has zero weights. Please check the occupancy column.");
  }
  n=1.0/n;
  for(unsigned i=0; i<weights.size(); ++i) {
    weights[i]*=n;
  }

  // normalize weights for rmsd calculation
  std::vector<double> weights_measure=pdb.getBeta();
  n=0.0;
  for(unsigned i=0; i<weights_measure.size(); ++i) {
    n+=weights_measure[i];
  }
  n=1.0/n;
  for(unsigned i=0; i<weights_measure.size(); ++i) {
    weights_measure[i]*=n;
  }

  // subtract the center
  for(unsigned i=0; i<weights.size(); ++i) {
    center+=myPositions[i]*weights[i];
  }
  for(unsigned i=0; i<weights.size(); ++i) {
    myPositions[i]-=center;
  }

  if(type=="OPTIMAL" or type=="OPTIMAL-FAST" ) {
    rmsd=Tools::make_unique<RMSD>();
    rmsd->set(weights,weights_measure,myPositions,type,false,false);// note: the reference is shifted now with center in the origin
    log<<"  Method chosen for fitting: "<<rmsd->getMethod()<<" \n";
  }
  if(nopbc) {
    log<<"  Ignoring PBCs when doing alignment, make sure your molecule is whole!<n";
  }
  // register the value of rmsd (might be useful sometimes)
  addValue();
  setNotPeriodic();

  // I remove this optimization now in order to use makeWhole()
  // Notice that for FIT_TO_TEMPLATE TYPE=OPTIMAL a copy was made anyway
  // (due to the need to store position to propagate forces on rotational matrix later)
  // For FIT_TO_TEMPLATE TYPE=SIMPLE in principle we could use it and write an ad hoc
  // version of makeWhole that only computes the center. Too lazy to do it now.
  // In case we do it later, remember that uncommenting this line means that
  // getPositions will not work anymore! GB
  // doNotRetrieve();

  // this is required so as to allow modifyGlobalForce() to return correct
  // also for forces that are not owned (and thus not zeored) by all processors.
  pbc_action=plumed.getActionSet().selectWithLabel<PbcAction*>("Box");
  if( !pbc_action ) {
    error("cannot align box has not been set");
  }
  boxValue=pbc_action->copyOutput(0);
}


void FitToTemplate::calculate() {

  if(!nopbc) {
    makeWhole();
  }

  if (type=="SIMPLE") {
    Vector cc;

    for(unsigned i=0; i<p_aligned.size(); ++i) {
      cc+=weights[i]*getPosition(i);
    }

    shift=center-cc;
    setValue(shift.modulo());
    unsigned nat = getTotAtoms();
    for(unsigned i=0; i<nat; i++) {
      std::pair<std::size_t,std::size_t> a = getValueIndices( AtomNumber::index(i));
      Vector ato=getGlobalPosition(a);
      setGlobalPosition(a,ato+shift);
    }
  } else if( type=="OPTIMAL" or type=="OPTIMAL-FAST") {
    // specific stuff that provides all that is needed
    double r=rmsd->calc_FitElements( getPositions(), rotation,  drotdpos, centeredpositions, center_positions);
    setValue(r);
    unsigned nat = getTotAtoms();
    for(unsigned i=0; i<nat; i++) {
      std::pair<std::size_t,std::size_t> a = getValueIndices( AtomNumber::index(i));
      Vector ato=getGlobalPosition(a);
      setGlobalPosition(a,matmul(rotation,ato-center_positions)+center);
    }
// rotate box
    pbc_action->getPbc().setBox(matmul(pbc_action->getPbc().getBox(),transpose(rotation)));
  }
}

void FitToTemplate::apply() {
  auto nat=getTotAtoms();
  if (type=="SIMPLE") {
    Vector totForce;
    for(unsigned i=0; i<nat; i++) {
      std::pair<std::size_t,std::size_t> a = getValueIndices( AtomNumber::index(i));
      totForce+=getForce(a);
    }
    Tensor vv=Tensor(center,totForce);
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        boxValue->addForce( 3*i+j, vv(i,j) );
      }
    for(unsigned i=0; i<p_aligned.size(); ++i) {
      addForce( p_aligned[i], -totForce*weights[i]);
    }
  } else if ( type=="OPTIMAL" or type=="OPTIMAL-FAST") {
    Vector totForce;
    for(unsigned i=0; i<nat; i++) {
      std::pair<std::size_t,std::size_t> a = getValueIndices( AtomNumber::index(i));
      Vector f=getForce(a);
// rotate back forces
      Vector nf=matmul(transpose(rotation),f);
      addForce(a, nf-f);
// accumulate rotated c.o.m. forces - this is already in the non rotated reference frame
      totForce+=nf;
    }
    Tensor virial;
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        virial[i][j] = boxValue->getForce( 3*i+j );
      }
// notice that an extra Tensor(center,matmul(rotation,totForce)) is required to
// compute the derivatives of the rotation with respect to center
    Tensor ww=matmul(transpose(rotation),virial+Tensor(center,matmul(rotation,totForce)));
// rotate back virial
    virial=matmul(transpose(rotation),matmul(virial,rotation));

// now we compute the force due to alignment
    for(unsigned i=0; i<p_aligned.size(); i++) {
      Vector g;
      for(unsigned k=0; k<3; k++) {
// this could be made faster computing only the diagonal of d
        Tensor d=matmul(ww,RMSD::getMatrixFromDRot(drotdpos,i,k));
        g[k]=(d(0,0)+d(1,1)+d(2,2));
      }
// here is the extra contribution
      addForce( p_aligned[i], -g-weights[i]*totForce );
// here it the contribution to the virial
// notice that here we can use absolute positions since, for the alignment to be defined,
// positions should be in one well defined periodic image
      virial+=extProduct(getPosition(i),g);
    }
// finally, correction to the virial
    boxValue->clearInputForce();
    virial+=extProduct(matmul(transpose(rotation),center),totForce);
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        boxValue->addForce( 3*i+j, virial(i,j) );
      }
  }
}

}
}

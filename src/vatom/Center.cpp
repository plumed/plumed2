/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
\plumedfile
c1: COM ATOMS=1-7
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC VATOM CENTER
/*
Calculate the center for a group of atoms, with arbitrary weights.

The computed
center is stored as a virtual atom that can be accessed in
an atom list through the label for the CENTER action that creates it.
Notice that the generated virtual atom has charge equal to the sum of the
charges and mass equal to the sum of the masses. If used with the MASS flag,
then it provides a result identical to \ref COM.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding the molecule using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\note As an experimental feature, CENTER also supports a keyword PHASES.
This keyword finds the center of mass for sets of atoms that have been split by the period boundaries by computing scaled coordinates and average
trigonometric functions, similarly to \ref CENTER_OF_MULTICOLVAR.
Notice that by construction this center position is
not invariant with respect to rotations of the atoms at fixed cell lattice.
In addition, for symmetric Bravais lattices, it is not invariant with respect
to special symmetries. E.g., if you have an hexagonal cell, the center will
not be invariant with respect to rotations of 120 degrees.
On the other hand, it might make the treatment of PBC easier in difficult cases.

\par Examples

\plumedfile
# a point which is on the line connecting atoms 1 and 10, so that its distance
# from 10 is twice its distance from 1:
c1: CENTER ATOMS=1,1,10
# this is another way of stating the same:
c1bis: CENTER ATOMS=1,10 WEIGHTS=2,1

# center of mass among these atoms:
c2: CENTER ATOMS=2,3,4,5 MASS

d1: DISTANCE ATOMS=c1,c2

PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding the molecule using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
\plumedfile
c1: COM ATOMS=1-7
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC


class Center:
  public ActionWithVirtualAtom
{
  std::vector<double> weights;
  bool weight_mass, weight_charge;
  Value* val_weights;
  unsigned myx, myw, bufstart, nspace;
  std::vector<Tensor> deriv;
  std::vector<double> val_forces;
  std::vector<std::vector<double> > val_deriv;
  Tensor invbox2pi;
  bool nopbc;
  bool first;
  bool phases;
  bool dophases;
public:
  explicit Center(const ActionOptions&ao);
  void calculate();
  void apply() override;
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfDerivatives() const ;
  void setStashIndices( unsigned& nquants );
  void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize );
  void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  void gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const ;
  void transformFinalValueAndDerivatives( const std::vector<double>& buffer );
};

PLUMED_REGISTER_ACTION(Center,"CENTER")

void Center::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys); keys.remove("NUMERICAL_DERIVATIVES");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
  keys.addFlag("PHASES",false,"Compute center using trigonometric phases");
  keys.addFlag("MASS",false,"calculate the center of mass");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  weight_mass(false),
  weight_charge(false),
  val_weights(NULL),
  myx(0), myw(0), bufstart(0), nspace(1),
  nopbc(false),
  first(true),
  phases(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms); parseFlag("NOPBC",nopbc); parseFlag("PHASES",phases);
  bool usemass = false; parseFlag("MASS",usemass);
  if(atoms.size()==0) error("at least one atom should be specified");
  std::vector<std::string> str_weights; parseVector("WEIGHTS",str_weights);
  if( usemass ) {
      if( str_weights.size()>0 ) error("USEMASS is incompatible with WEIGHTS");
      str_weights.resize(1); str_weights[0]="@masses";
  }
  if( str_weights.size()==0 ) {
    log<<"  computing the geometric center of atoms:\n";
    weights.resize( atoms.size() );
    for(unsigned i=0; i<atoms.size(); i++) weights[i] = 1.;
  } else if( str_weights.size()==1 ) {
    if( str_weights[0]=="@masses" ) {
      weight_mass=true;
      log<<"  computing the center of mass of atoms:\n";
    } else if( str_weights[0]=="@charges" ) {
      weight_charge=true;
      log<<"  computing the center of charge of atoms:\n";
    } else {
      std::size_t dot=str_weights[0].find_first_of("."); unsigned nargs=0; std::vector<Value*> args;
      if( dot!=std::string::npos ) {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0].substr(0,dot) );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        action->interpretDataLabel( str_weights[0], this, nargs, args );
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0] );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        action->interpretDataLabel( str_weights[0], this, nargs, args );
      }
      if( args.size()!=1 ) error("should only have one value as input to WEIGHT");
      if( args[0]->getRank()!=1 || args[0]->getShape()[0]!=atoms.size() ) error("value input for WEIGHTS has wrong shape");
      val_weights = args[0]; std::vector<std::string> empty(1); empty[0] = (val_weights->getPntrToAction())->getLabel();
      if( (val_weights->getPntrToAction())->valuesComputedInChain() ) (val_weights->getPntrToAction())->addActionToChain( empty, this );
      log.printf("  atoms are weighted by values in vector labelled %s \n",val_weights->getName().c_str() );
    }
  } else {
    log<<" with weights:";
    if( str_weights.size()!=atoms.size() ) error("number of elements in weight vector does not match the number of atoms");
    weights.resize( atoms.size() );
    for(unsigned i=0; i<weights.size(); ++i) {
      if(i%25==0) log<<"\n";
      Tools::convert( str_weights[i], weights[i] ); log.printf(" %f",weights[i]);
    }
    log.printf("\n");
  }
  log.printf("  of atoms:");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i>0 && i%25==0) log<<"\n";
    log.printf("  %d",atoms[i].serial());
  }
  log<<"\n";
  checkRead();
  if( !nopbc && !phases ) {
    log<<"  PBC will be ignored\n";
  } else if( phases ) {
    nopbc=true; log<<"  phases will be used to take into account PBC\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms); deriv.resize( atoms.size() );
  if( val_weights ) {
    addDependency( val_weights->getPntrToAction() );
    val_deriv.resize(3); val_forces.resize( (val_weights->getPntrToAction())->getNumberOfDerivatives() );
    val_deriv[0].resize( (val_weights->getPntrToAction())->getNumberOfDerivatives() );
    val_deriv[1].resize( (val_weights->getPntrToAction())->getNumberOfDerivatives() );
    val_deriv[2].resize( (val_weights->getPntrToAction())->getNumberOfDerivatives() );
  }
  // And create task list
  for(unsigned i=0; i<atoms.size(); ++i) addTaskToList( i );
}

unsigned Center::getNumberOfDerivatives() const {
  if( val_weights ) return 3*getNumberOfAtoms() + (val_weights->getPntrToAction())->getNumberOfDerivatives();
  return 3*getNumberOfAtoms();
}

void Center::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  // Check that we have the mass information we need on the first step
  if( first ) {
    // Check if we have masses if we need them
    if( weight_mass ) {
        for(unsigned i=0; i<getNumberOfAtoms(); i++) {
          if(std::isnan(getMass(i))) {
            error(
              "You are trying to compute a CENTER or COM but masses are not known.\n"
              "        If you are using plumed driver, please use the --mc option"
            );
          }
        }
    // Check we have charges if we need them
    } else if ( weight_charge && !plumed.getAtoms().chargesWereSet() ) {
        error(
            "You are trying to compute a center of charnge but chargest are not known.\n"
            "        If you are using plumed driver, please use the --mc option"
          );
    }
    first=false;
  }
  dophases=(getPbc().isSet() ? phases : false);

  // Check if we need to make the whole thing
  if(!nopbc && !dophases) makeWhole();
  // Setup inverse box thing so we can do Berry phase
  invbox2pi=2*pi*getPbc().getInvBox();

  // Set mass for center
  double mass=0.;
  for(unsigned i=0; i<getNumberOfAtoms(); i++) mass+=getMass(i);
  setMass( mass );
  // Set charge for center
  if( plumed.getAtoms().chargesWereSet() ) {
    double charge=0.;
    for(unsigned i=0; i<getNumberOfAtoms(); i++) charge+=getCharge(i);
    setCharge(charge);
  } else setCharge(0.0);
}

void Center::setStashIndices( unsigned& nquants ) {
  if( phases ) {
    myx = nquants; myw = nquants + 6; nquants +=7;
  } else {
    myx = nquants; myw = nquants + 3; nquants += 4;
  }
}

void Center::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  bufstart = bufsize; if( !doNotCalculateDerivatives() ) nspace = 1 + getNumberOfDerivatives();
  if( phases ) bufsize += 7*nspace; else bufsize += 4*nspace;
  ActionWithValue::getSizeOfBuffer( nactive_tasks, bufsize );
}

void Center::calculate() {
  if( actionInChain() ) return;
  runAllTasks();
}

void Center::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  Vector pos = getPosition( task_index ); double w;
  if( weight_mass ) {
    w = getMass(task_index);
  } else if( weight_charge ) {
    if( !plumed.getAtoms().chargesWereSet() ) plumed_merror("cannot calculate center of charge if chrages are unset");
    w = getCharge(task_index);
  } else if( val_weights && actionInChain() ) {
    w = myvals.get( val_weights->getPositionInStream() );
  } else if( val_weights ) {
    w = val_weights->get(task_index);
  } else {
    plumed_dbg_assert( task_index<weights.size() );
    w = weights[task_index];
  }
  if( dophases ) {
    Vector stmp, ctmp, fpos = matmul(pos,invbox2pi); myvals.addValue( myw, w );
    for(unsigned j=0; j<3; ++j) {
      stmp[j] = sin( fpos[j] ); ctmp[j] = cos( fpos[j] );
      myvals.addValue( myx+j, w*stmp[j] ); myvals.addValue( myx+3+j, w*ctmp[j] );
    }
    if( !doNotCalculateDerivatives() ) {
      for(unsigned j=0; j<3; ++j) {
        for(unsigned k=0;k<3;++k) {
            myvals.addDerivative( myx+j, 3*task_index+k, w*ctmp[j]*invbox2pi[k][j] ); myvals.updateIndex( myx+j, 3*task_index+k );
            myvals.addDerivative( myx+3+j, 3*task_index+k, -w*stmp[j]*invbox2pi[k][j] ); myvals.updateIndex( myx+3+j, 3*task_index+k );
        }
      }
      if( val_weights ) {
        unsigned base = 3*getNumberOfAtoms();
        unsigned istrn = val_weights->getPositionInStream();
        for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
          unsigned kindex = myvals.getActiveIndex(istrn,k);
          double der = myvals.getDerivative( istrn, kindex );
          for(unsigned j=0; j<3; ++j) {
            myvals.addDerivative( myx+j, base+kindex, der*stmp[j] ); myvals.updateIndex( myx+j, base+kindex );
            myvals.addDerivative( myx+3+j, base+kindex, der*ctmp[j] ); myvals.updateIndex( myx+3+j, base+kindex );
          }
          myvals.addDerivative( myw, base+kindex, der ); myvals.updateIndex( myw, base+kindex );
        }
      }
    }
  } else {
    myvals.addValue( myw, w ); for(unsigned j=0; j<3; ++j) myvals.addValue( myx+j, w*pos[j] );
    if( !doNotCalculateDerivatives() ) {
      for(unsigned j=0; j<3; ++j) {
        myvals.addDerivative( myx+j, 3*task_index+j, w ); myvals.updateIndex( myx+j, 3*task_index+j );
      }
      if( val_weights ) {
        unsigned base = 3*getNumberOfAtoms();
        unsigned istrn = val_weights->getPositionInStream();
        for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
          unsigned kindex = myvals.getActiveIndex(istrn,k);
          double der = myvals.getDerivative( istrn, kindex );
          for(unsigned j=0; j<3; ++j) { myvals.addDerivative( myx+j, kindex, der*pos[j] ); myvals.updateIndex( myx+j, base+kindex ); }
          myvals.addDerivative( myw, base+kindex, der ); myvals.updateIndex( myw, base+kindex );
        }
      }
    }
  }
}

void Center::gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const {
  unsigned ntmp_vals = 4;
  if( dophases ) {
    unsigned sstart = bufstart, cstart = bufstart + 3*nspace;
    for(unsigned j=0; j<3; ++j) {
      buffer[sstart] += myvals.get(myx+j); sstart+=nspace;
      buffer[cstart] += myvals.get(myx+3+j); cstart+=nspace;
    }
    buffer[bufstart+6*nspace] += myvals.get( myw ); ntmp_vals = 7;
  } else {
    unsigned bstart = bufstart;
    for(unsigned j=0; j<3; ++j) { buffer[bstart] += myvals.get(myx+j); bstart+=nspace; }
    buffer[bufstart+3*nspace] += myvals.get( myw );
    ntmp_vals = 4;
  }
  if( !doNotCalculateDerivatives() ) {
    unsigned bstart = bufstart;
    for(unsigned j=0; j<ntmp_vals; ++j) {
      for(unsigned k=0; k<myvals.getNumberActive(myx+j); ++k) {
        unsigned kindex = myvals.getActiveIndex(myx+j,k);
        plumed_dbg_assert( bstart + 1 + kindex<buffer.size() );
        buffer[bstart + 1 + kindex] += myvals.getDerivative( myx+j, kindex );
      }
      bstart += nspace;
    }
  }
}

void Center::transformFinalValueAndDerivatives( const std::vector<double>& buffer ) {
  // Get final position
  if( dophases ) {
    double ww = buffer[bufstart + 6*nspace]; Vector fpos;
    Tensor box2pi=getPbc().getBox() / (2*pi); Vector stmp, ctmp; 
    for(unsigned i=0; i<3; ++i) {
      stmp[i]=buffer[bufstart+i*nspace]/ww; ctmp[i]=buffer[bufstart+(3+i)*nspace]/ww;
      fpos[i] = atan2( stmp[i], ctmp[i] );
    }
    setPosition( matmul(fpos,box2pi) );
    // And derivatives
    if( !doNotCalculateDerivatives() ) {
      double inv_weight = 1.0 / ww; Vector tander;
      for(unsigned j=0; j<3; ++j) {
        double tmp = stmp[j] / ctmp[j];
        tander[j] = 1.0 / ( 1 + tmp*tmp );
      }
      double sderv, cderv; Tensor dd;
      for(unsigned i=0; i<getNumberOfAtoms(); ++i ) {
        for(unsigned j=0; j<3; ++j) {
          for(unsigned k=0; k<3; ++k) {
              sderv = inv_weight*buffer[bufstart + 1 + j*nspace + 3*i + k]; cderv = inv_weight*buffer[bufstart + 1 + (3+j)*nspace + 3*i + k];
              dd(j,k) = tander[j]*( sderv/ctmp[j]  - stmp[j]*cderv/(ctmp[j]*ctmp[j]) );
          }
        }
        deriv[i]=matmul(dd,box2pi);
      }
      setAtomsDerivatives(deriv);
      if( val_weights ) {
        unsigned k=0; Vector val_dev;
        for(unsigned i=3*getNumberOfAtoms(); i<getNumberOfDerivatives(); ++i ) {
          double wder = buffer[bufstart + 1 + 6*nspace +i];
          for(unsigned j=0; j<3; ++j) {
            sderv = inv_weight*buffer[bufstart + 1 + j*nspace + i] - inv_weight*stmp[j]*wder;
            cderv = inv_weight*buffer[bufstart + 1 + (3+j)*nspace + i] - inv_weight*ctmp[j]*wder;
            val_dev[j] = tander[j]*( sderv/ctmp[j]  - stmp[j]*cderv/(ctmp[j]*ctmp[j]) );
          }
          Vector vvv = matmul( box2pi, val_dev );
          for(unsigned j=0; j<3; ++j) val_deriv[j][k] = vvv[j];
          k++;
        }
      }
    }
  } else {
    double ww = buffer[bufstart + 3*nspace]; Vector pos;
    for(unsigned i=0; i<3; ++i) pos[i]=buffer[bufstart + i*nspace]/ww;
    setPosition(pos);
    // And final derivatives
    if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0; i<getNumberOfAtoms(); ++i ) {
        for(unsigned j=0; j<3; ++j) {
          deriv[i](0,j) = buffer[bufstart + 1 + 3*i + j ] / ww;
          deriv[i](1,j) = buffer[bufstart + nspace + 1 + 3*i + j ] / ww;
          deriv[i](2,j) = buffer[bufstart + 2*nspace + 1 + 3*i +j ] / ww;
        }
      }
      setAtomsDerivatives(deriv);
      if( val_weights ) {
        unsigned k=0;
        for(unsigned i=3*getNumberOfAtoms(); i<getNumberOfDerivatives(); ++i ) {
          for(unsigned j=0; j<3; ++j) {
            val_deriv[j][k] = buffer[bufstart + 1 + j*nspace + i ] / ww - pos[j]*buffer[bufstart + 1 + 6*nspace + i] / ww;
          }
          k++;
        }
      }
    }
  }
}

void Center::apply() {
  Vector & f(atoms.getVatomForces(getIndex())); unsigned start;
  if( f.modulo2()>epsilon && val_weights ) {
    ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( val_weights->getPntrToAction() );
    ActionAtomistic* aat = dynamic_cast<ActionAtomistic*>( val_weights->getPntrToAction() );
    for(unsigned j=0; j<3; ++j) {
      for(unsigned k=0; k<val_deriv[j].size(); ++k) val_forces[k] = f[j]*val_deriv[j][k];
      start=0;
      if( aarg ) aarg->setForcesOnArguments( 0, val_forces, start );
      if( aat ) aat->setForcesOnAtoms( val_forces, start );
    }
  }
  // And apply the forces to the centers
  ActionWithVirtualAtom::apply();
}

}
}

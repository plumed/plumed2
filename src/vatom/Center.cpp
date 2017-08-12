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
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

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
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.


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


class Center:
  public ActionWithVirtualAtom
{
  std::vector<double> weights;
  bool weight_mass, weight_charge;
  bool nopbc;
  unsigned myx, myy, myz, myw, bufstart, nspace;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  explicit Center(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfDerivatives() const ;
  void setStashIndices( unsigned& nquants );
  void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize );
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  void gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const ;
  void transformFinalValueAndDerivatives( const std::vector<double>& buffer );
};

PLUMED_REGISTER_ACTION(Center,"CENTER")
PLUMED_REGISTER_SHORTCUT(Center,"CENTER")
PLUMED_REGISTER_SHORTCUT(Center,"COM")

void Center::shortcutKeywords( Keywords& keys ){
  keys.addFlag("MASS",false,"calculate the center of mass");
}

void Center::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ) {
  std::vector<std::string> input; 
  input.push_back(lab +":"); input.push_back("CENTER");
  if( words[0]=="COM" ) input.push_back("WEIGHTS=@masses"); 
  else if( keys.count("MASS") ) input.push_back("WEIGHTS=@masses");
  else plumed_error();
  for(unsigned i=1;i<words.size();++i) input.push_back(words[i]);
  actions.push_back( input ); 
}

void Center::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys); keys.remove("NUMERICAL_DERIVATIVES");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
                                "If WEIGHTS=@masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
                                "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
                                "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  weight_mass(false),
  weight_charge(false),
  nopbc(false),
  myx(0), myy(0), myz(0), myw(0), bufstart(0), nspace(1)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  std::vector<std::string> str_weights; parseVector("WEIGHTS",str_weights);
  if( str_weights.size()==0) {
      log<<"  computing the geometric center of atoms:\n";
      weights.resize( atoms.size() );
      for(unsigned i=0; i<atoms.size(); i++) weights[i] = 1.;
  } else if( str_weights.size()==1 ) {
      if( str_weights[0]=="@masses" ){
          weight_mass=true;
          log<<"  computing the center of mass of atoms:\n";
      } else if( str_weights[0]=="@charges" ){
          weight_charge=true;
          log<<"  computing the center of charge of atoms:\n";
      } else {
         error("not implemented yet");
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
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i>0 && i%25==0) log<<"\n";
    log.printf("  %d",atoms[i].serial());
  }
  log<<"\n";
  parseFlag("NOPBC",nopbc);
  checkRead();
  if(!nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms);
  // And create task list
  for(unsigned i=0; i<atoms.size(); ++i) addTaskToList( i ); 
}

unsigned Center::getNumberOfDerivatives() const {
  return 3*getNumberOfAtoms();
}

void Center::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  tflags.assign(tflags.size(),1);
}

void Center::setStashIndices( unsigned& nquants ) {
  myx = nquants; myy = nquants + 1; myz = nquants + 2; myw = nquants + 3; nquants += 4;
}

void Center::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ){
  bufstart = bufsize; if( !doNotCalculateDerivatives() ) nspace = 1 + 3*getNumberOfAtoms();
  bufsize += 4*nspace; ActionWithValue::getSizeOfBuffer( nactive_tasks, bufsize );
}

void Center::calculate() {
  if(!nopbc) makeWhole();
  runAllTasks();
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

void Center::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  Vector pos = getPosition( task_index ); double w;
  if( weight_mass ){
      w = getMass(task_index);
  } else if( weight_charge ){
      if( !plumed.getAtoms().chargesWereSet() ) plumed_merror("cannot calculate center of charge if chrages are unset");
      w = getCharge(task_index);
  } else {
      plumed_dbg_assert( task_index<weights.size() );
      w = weights[task_index];
  }
  myvals.addValue( myx, w*pos[0] ); myvals.addValue( myy, w*pos[1] ); 
  myvals.addValue( myz, w*pos[2] ); myvals.addValue( myw, w );
  if( !doNotCalculateDerivatives() ) {
      myvals.addDerivative( myx, 3*task_index+0, w ); myvals.updateIndex( myx, 3*task_index+0 );
      myvals.addDerivative( myy, 3*task_index+1, w ); myvals.updateIndex( myy, 3*task_index+1 );
      myvals.addDerivative( myz, 3*task_index+2, w ); myvals.updateIndex( myz, 3*task_index+2 );
  }
}

void Center::gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const {
  buffer[bufstart] += myvals.get( myx ); buffer[bufstart+nspace] += myvals.get( myy );
  buffer[bufstart+2*nspace] += myvals.get( myz ); buffer[bufstart+3*nspace] += myvals.get( myw );
  if( !doNotCalculateDerivatives() ) {
      unsigned bstart = bufstart;
      for(unsigned k=0;k<myvals.getNumberActive(myx);++k){
          unsigned kindex = myvals.getActiveIndex(myx,k); 
          plumed_dbg_assert( bstart + 1 + kindex<buffer.size() );
          buffer[bstart + 1 + kindex] += myvals.getDerivative( myx, kindex ); 
      }
      bstart += nspace;
      for(unsigned k=0;k<myvals.getNumberActive(myy);++k){
          unsigned kindex = myvals.getActiveIndex(myy,k); 
          plumed_dbg_assert( bstart + 1 + kindex<buffer.size() );
          buffer[bstart + 1 + kindex] += myvals.getDerivative( myy, kindex );
      }
      bstart += nspace;
      for(unsigned k=0;k<myvals.getNumberActive(myz);++k){
          unsigned kindex = myvals.getActiveIndex(myz,k); 
          plumed_dbg_assert( bstart + 1 + kindex<buffer.size() );
          buffer[bstart + 1 + kindex] += myvals.getDerivative( myz, kindex );
      }
  }
}

void Center::transformFinalValueAndDerivatives( const std::vector<double>& buffer ) {
  // Get final position
  double ww = buffer[bufstart + 3*nspace];
  Vector pos; pos[0]=buffer[bufstart]/ww; pos[1]=buffer[bufstart+nspace]/ww; pos[2]=buffer[bufstart+2*nspace]/ww;
  setPosition(pos);
  // And final derivatives
  if( !doNotCalculateDerivatives() ) {
      std::vector<Tensor> deriv(getNumberOfAtoms());
      for(unsigned i=0; i<getNumberOfAtoms(); ++i ){
          for(unsigned j=0;j<3;++j){
              deriv[i](0,j) = buffer[bufstart + 1 + 3*i + j ] / ww;
              deriv[i](1,j) = buffer[bufstart + nspace + 1 + 3*i + j ] / ww;
              deriv[i](2,j) = buffer[bufstart + 2*nspace + 1 + 3*i +j ] / ww;
          }
      }  
      setAtomsDerivatives(deriv);
  }
}

}
}

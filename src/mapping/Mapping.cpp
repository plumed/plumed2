/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceAtoms.h"
#include "reference/MetricRegister.h"
#include "tools/PDB.h"
#include "tools/Matrix.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/ActionRegister.h"

namespace PLMD {

namespace mapping {

class Mapping :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
private:
/// Are we calculating squared distances
  bool squared;
/// This holds all the reference information
  std::vector<ReferenceConfiguration*> myframes;
/// The forces on each of the derivatives (used in apply)
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  explicit Mapping(const ActionOptions&);
  ~Mapping();
/// Overload the virtual functions that appear in both ActionAtomistic and ActionWithArguments
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void lockRequests();
  void unlockRequests();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Turn on the tasks that are currently active
  void buildCurrentTaskList( std::vector<unsigned>& tflags ) const ;
/// Do the actual calculation
  void calculate();
/// This calculates the distance from the reference configuration of interest
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
/// Apply the forces
  void apply();
};

PLUMED_REGISTER_ACTION(Mapping,"EUCLIDEAN_DISSIMILARITIES_VECTOR")

// void Mapping::shortcutKeywords( Keywords& keys ) {
//   keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
// }

void Mapping::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("DISABLE_CHECKS",false,"disable checks on reference input structures.");
}

Mapping::Mapping(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  // Read the squared flag
  parseFlag("SQUARED",squared);
  // Read the input
  std::string mtype; parse("TYPE",mtype);
  bool skipchecks; parseFlag("DISABLE_CHECKS",skipchecks);

  // Open reference file
  std::string reference; parse("REFERENCE",reference);
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations
  bool do_read=true; 
  while (do_read) {
    // Read the pdb file
    PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    // Break if we are done
    if( !do_read ) break ;
    // Fix argument names
    expandArgKeywordInPDB( mypdb );
    // And add a task to the list
    addTaskToList( myframes.size() );
    // And read the frame
    myframes.push_back( metricRegister().create<ReferenceConfiguration>( mtype, mypdb ) );
  }
  fclose(fp);

  if(myframes.size()==0 ) error("no reference configurations were specified");
  log.printf("  found %u configurations in file %s\n",myframes.size(),reference.c_str() );
  std::vector<unsigned> shape(1); shape[0]=myframes.size(); addValue( shape ); setNotPeriodic();

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  for(unsigned i=0; i<myframes.size(); ++i) { myframes[i]->getAtomRequests( atoms, skipchecks ); myframes[i]->getArgumentRequests( args, skipchecks ); }
  requestAtoms( atoms ); std::vector<Value*> req_args;
  interpretArgumentList( args, req_args ); requestArguments( req_args );
  // Resize forces array
  if( getNumberOfAtoms()>0 ) {
    forcesToApply.resize( 3*getNumberOfAtoms() + 9 + getNumberOfArguments() );
  } else {
    forcesToApply.resize( getNumberOfArguments() );
  }
}

void Mapping::buildCurrentTaskList( std::vector<unsigned>& tflags ) const {
  tflags.assign(tflags.size(),1);
}

void Mapping::calculate(){
  plumed_dbg_assert( !done_over_stream && getFullNumberOfTasks()>0 ); 
  runAllTasks();
}

void Mapping::performTask( const unsigned& current, MultiValue& myvals ) const {
  ReferenceValuePack mypack( getNumberOfArguments(), getNumberOfAtoms(), myvals ); mypack.setValIndex( getPntrToOutput(0)->getPositionInStream() );
  unsigned nargs2=myframes[current]->getNumberOfReferenceArguments(); unsigned nat2=myframes[current]->getNumberOfReferencePositions();
  if( mypack.getNumberOfAtoms()!=nat2 || mypack.getNumberOfArguments()!=nargs2 ) mypack.resize( nargs2, nat2 );
  if( nat2>0 ) {
    ReferenceAtoms* myat2=dynamic_cast<ReferenceAtoms*>( myframes[current] ); plumed_dbg_assert( myat2 );
    for(unsigned i=0; i<nat2; ++i) mypack.setAtomIndex( i, myat2->getAtomIndex(i) );
  }
  double dd = myframes[current]->calculate( getPositions(), getPbc(), getArguments(), mypack, squared );
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), dd );
  if( !doNotCalculateDerivatives() && getNumberOfAtoms()>0 && !mypack.virialWasSet() ) {
    Tensor tvir; tvir.zero();
    for(unsigned i=0; i<mypack.getNumberOfAtoms(); ++i) tvir +=-1.0*Tensor( getPosition( mypack.getAtomIndex(i) ), mypack.getAtomDerivative(i) );
    mypack.addBoxDerivatives( tvir );
  }
}

unsigned Mapping::getNumberOfDerivatives() {
  unsigned nat=getNumberOfAtoms();
  if(nat>0) return 3*nat + 9 + getNumberOfArguments();
  return getNumberOfArguments();
}

void Mapping::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void Mapping::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void Mapping::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

void Mapping::apply() {
//  if( getForcesFromVessels( forcesToApply ) ) {
//    addForcesOnArguments( forcesToApply );
//    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
//  }
}


}
}

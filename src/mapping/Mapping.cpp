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
#include "reference/ReferenceAtoms.h"
#include "reference/MetricRegister.h"
#include "tools/PDB.h"
#include "tools/Matrix.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/ActionRegister.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

PLUMED_REGISTER_ACTION(Mapping,"EUCLIDEAN_DISSIMILARITIES_VECTOR")
PLUMED_REGISTER_SHORTCUT(Mapping,"DRMSD")
PLUMED_REGISTER_SHORTCUT(Mapping,"RMSD")
PLUMED_REGISTER_SHORTCUT(Mapping,"MULTI-RMSD")
PLUMED_REGISTER_SHORTCUT(Mapping,"TARGET")

void Mapping::shortcutKeywords( Keywords& keys ) {
  keys.add("compulsory","TYPE","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("optional","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.addFlag("NOPBC",false,"don't use PBC in calculation of DRMSD vectors");
}

void Mapping::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ){
  std::vector<std::string> thisact; thisact.push_back( lab + ":" );
  thisact.push_back( "EUCLIDEAN_DISSIMILARITIES_VECTOR" );
  for(unsigned i=1;i<words.size();++i) thisact.push_back( words[i] );
  if( words[0]=="DRMSD" ){
      if( keys.find("TYPE")!=keys.end() ) thisact.push_back( keys.find("TYPE")->first + "=" + keys.find("TYPE")->second );
      else thisact.push_back( "TYPE=DRMSD");

      if( keys.find("LOWER_CUTOFF")!=keys.end() ) thisact.push_back( keys.find("LOWER_CUTOFF")->first + "=" + keys.find("LOWER_CUTOFF")->second );
      else plumed_merror("LOWER_CUTOFF must be specified in DRMSD actions");

      if( keys.find("UPPER_CUTOFF")!=keys.end() ) thisact.push_back( keys.find("UPPER_CUTOFF")->first + "=" + keys.find("UPPER_CUTOFF")->second );
      else plumed_merror("UPPER_CUTOFF must be specified in DRMSD actions");

      if( keys.find("NOPBC")!=keys.end() ) thisact.push_back("NOPBC");
  } else if( words[0]=="MULTI-RMSD" ){
      if( keys.find("TYPE")!=keys.end() ) thisact.push_back( keys.find("TYPE")->first + "=" + keys.find("TYPE")->second ); 
      else thisact.push_back( "TYPE=MULTI-SIMPLE"); 
  } else if( words[0]=="TARGET" ){
      if( keys.find("TYPE")!=keys.end() ) thisact.push_back( keys.find("TYPE")->first + "=" + keys.find("TYPE")->second );
      else thisact.push_back( "TYPE=EUCLIDEAN" );
  } else if( words[0]=="RMSD" ){
      if( keys.find("TYPE")!=keys.end() ) thisact.push_back( keys.find("TYPE")->first + "=" + keys.find("TYPE")->second ); 
      else thisact.push_back( "TYPE=SIMPLE" );
  } else {
      plumed_assert( words[0]=="EUCLIDEAN_DISSIMILARITIES_VECTOR" );
      actions.push_back( words ); return;
  }
  actions.push_back( thisact );
}

void Mapping::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.add("optional","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("optional","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.addFlag("NOPBC",false,"don't use PBC in calculation of DRMSD vectors");
  keys.addFlag("DISABLE_CHECKS",false,"disable checks on reference input structures.");
}

Mapping::Mapping(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  // Read the squared flag
  if( keywords.exists("SQUARED") ) parseFlag("SQUARED",squared);
  // Read the input
  std::string mtype; parse("TYPE",mtype);
  bool skipchecks; parseFlag("DISABLE_CHECKS",skipchecks);

  // Read stuff for DRMSD 
  double lcutoff=0; parse("LOWER_CUTOFF",lcutoff);
  double ucutoff=std::numeric_limits<double>::max(); parse("UPPER_CUTOFF",ucutoff);
  bool nopbc; parseFlag("NOPBC",nopbc); std::vector<std::string> drmsd_remarks;
  if( nopbc || lcutoff>0 || ucutoff<std::numeric_limits<double>::max() ){
      std::string lstr; Tools::convert( lcutoff, lstr ); drmsd_remarks.push_back( "LOWER_CUTOFF=" + lstr );
      std::string ustr; Tools::convert( ucutoff, ustr ); drmsd_remarks.push_back( "UPPER_CUTOFF=" + ustr );
      if( nopbc ) drmsd_remarks.push_back( "NOPBC" );
  }

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
    if( myframes.size()>0 && !do_read ) break ;
    // Add remarks that have been read from the input line
    mypdb.addRemark( drmsd_remarks );
    // Fix argument names
    expandArgKeywordInPDB( mypdb );
    // And add a task to the list
    addTaskToList( myframes.size() );
    // And read the frame
    myframes.push_back( metricRegister().create<ReferenceConfiguration>( mtype, mypdb ) );
  }
  fclose(fp);

  if(myframes.size()==0 ) error("no reference configurations were specified");
  log.printf("  calculating distance from reference configurations using %s metric \n", mtype.c_str() );
  log.printf("  found %u configurations in file %s\n",myframes.size(),reference.c_str() );
  std::vector<unsigned> shape; 
  if( myframes.size()>1 ){ shape.resize(1); shape[0]=myframes.size(); addValue( shape ); }
  else addValueWithDerivatives( shape ); 
  setNotPeriodic();

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  for(unsigned i=0; i<myframes.size(); ++i) { myframes[i]->getAtomRequests( atoms, skipchecks ); myframes[i]->getArgumentRequests( args, skipchecks ); }
  requestAtoms( atoms ); std::vector<Value*> req_args;
  interpretArgumentList( args, req_args ); requestArguments( req_args, false );
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
  plumed_dbg_assert( !actionInChain() && getFullNumberOfTasks()>0 ); 
  runAllTasks();
}

double Mapping::calculateDistanceFromReference( const unsigned& current, ReferenceValuePack& mypack ) const {
  double dd = myframes[current]->calculate( getPositions(), getPbc(), getArguments(), mypack, squared );
  if( !doNotCalculateDerivatives() && getNumberOfAtoms()>0 && !mypack.virialWasSet() ) {
    Tensor tvir; tvir.zero();
    for(unsigned i=0; i<mypack.getNumberOfAtoms(); ++i) tvir +=-1.0*Tensor( getPosition( mypack.getAtomIndex(i) ), mypack.getAtomDerivative(i) );
    mypack.addBoxDerivatives( tvir ); 
  }
  return dd;
}

double Mapping::calculateDistanceBetweenReferenceAndThisPoint( const unsigned& current, const std::vector<Vector>& pos, 
                                                               const std::vector<double>& args, ReferenceValuePack& mypack ) const {
  double dd = myframes[current]->calc( pos, getPbc(), getArguments(), args, mypack, squared );
  if( !doNotCalculateDerivatives() && getNumberOfAtoms()>0 && !mypack.virialWasSet() ) {
    Tensor tvir; tvir.zero();
    for(unsigned i=0; i<mypack.getNumberOfAtoms(); ++i) tvir +=-1.0*Tensor( getPosition( mypack.getAtomIndex(i) ), mypack.getAtomDerivative(i) );
    mypack.addBoxDerivatives( tvir );
  } 
  return dd;  
}

void Mapping::performTask( const unsigned& current, MultiValue& myvals ) const {
  ReferenceValuePack mypack( getNumberOfArguments(), getNumberOfAtoms(), myvals ); mypack.setValIndex( getPntrToOutput(0)->getPositionInStream() );
  unsigned nargs2=myframes[current]->getNumberOfReferenceArguments(); unsigned nat2=myframes[current]->getNumberOfReferencePositions();
  if( mypack.getNumberOfAtoms()!=nat2 || mypack.getNumberOfArguments()!=nargs2 ) mypack.resize( nargs2, nat2 );
  if( nat2>0 ) {
    ReferenceAtoms* myat2=dynamic_cast<ReferenceAtoms*>( myframes[current] ); plumed_dbg_assert( myat2 );
    for(unsigned i=0; i<nat2; ++i) mypack.setAtomIndex( i, myat2->getAtomIndex(i) );
  }
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), calculateDistanceFromReference( current, mypack ) );
  // double dd = myframes[current]->calculate( getPositions(), getPbc(), getArguments(), mypack, squared );
  // myvals.setValue( getPntrToOutput(0)->getPositionInStream(), dd );
  // if( !doNotCalculateDerivatives() && getNumberOfAtoms()>0 && !mypack.virialWasSet() ) {
  //   Tensor tvir; tvir.zero();
  //   for(unsigned i=0; i<mypack.getNumberOfAtoms(); ++i) tvir +=-1.0*Tensor( getPosition( mypack.getAtomIndex(i) ), mypack.getAtomDerivative(i) );
  //   mypack.addBoxDerivatives( tvir );
  // }
}

unsigned Mapping::getNumberOfDerivatives() const {
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
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ){
      unsigned mm=0; setForcesOnArguments( forcesToApply, mm );
      setForcesOnAtoms( forcesToApply, getNumberOfArguments()  );
  }
}


}
}

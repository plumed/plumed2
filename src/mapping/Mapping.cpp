/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "Mapping.h"
#include "vesselbase/Vessel.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"
#include "tools/PDB.h"
#include "tools/Matrix.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace mapping {

void Mapping::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys ); 
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); 
  ActionAtomistic::registerKeywords( keys ); 
  vesselbase::ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","TYPE","OPTIMAL","the manner in which distances are calculated");
  keys.addFlag("DISABLE_CHECKS",false,"disable checks on reference input structures.");
}

Mapping::Mapping(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithArguments(ao),
ActionWithValue(ao),
ActionWithVessel(ao)
{
  // Read the input
  std::string mtype; parse("TYPE",mtype);
  std::string reference; parse("REFERENCE",reference);
  bool skipchecks; parseFlag("DISABLE_CHECKS",skipchecks);

  // Read the properties we require
  bool hasprop=false;
  if( keywords.exists("PROPERTY") ){
     hasprop=true; parseVector("PROPERTY",property);
     if(property.size()==0) error("no properties were specified");
  } else {
     property.resize(1); property[0]="sss";
  }

  // Open reference file
  FILE* fp=fopen(reference.c_str(),"r"); 
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations 
  bool do_read=true; 
  std::vector<AtomNumber> atoms_to_retrieve;
  std::vector<std::string> args_to_retrieve;
  while (do_read){
     PDB mypdb; 
     // Read the pdb file
     do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
     // Fix argument names
     expandArgKeywordInPDB( mypdb );
     if(do_read){
        // Create the reference configuration
        ReferenceConfiguration* mymsd;
        // If skipchecks are enabled metric types must be specified in the input file
        if(skipchecks) mymsd=metricRegister().create<ReferenceConfiguration>( "", mypdb );
        else mymsd=metricRegister().create<ReferenceConfiguration>( mtype, mypdb ); 
        // Get the atoms and arguments we require
        mymsd->getAtomRequests( atoms_to_retrieve, skipchecks );
        mymsd->getArgumentRequests( args_to_retrieve, skipchecks );

        // Now get the low dimensional projection
        std::vector<double> labelvals;
        if(hasprop){
           labelvals.resize( property.size() );
           for(unsigned i=0;i<property.size();++i) mymsd->parse( property[i], labelvals[i] );  
        } else {
           labelvals.resize(1);  
           labelvals[0]=static_cast<double>( frames.size() + 1 );
        }
        mymsd->checkRead();                       // Check that everything in the input has been read in
        frames.push_back( mymsd );                // Store the reference 
        low_dim.push_back( labelvals );           // Store the low-dimensional projection
     } else {
        break;
     } 
  }
  fclose(fp); 
  if(frames.size()==0 ) error("no reference configurations were specified");
  log.printf("  found %d configurations in file %s\n",frames.size(),reference.c_str() );

  // Create copies of all the frames (so we can store all derivatives stuff)
  unsigned nframes=frames.size();
  for(unsigned i=0;i<nframes;++i){ 
     frames.push_back( new FakeFrame( ReferenceConfigurationOptions("fake") ) ); 
  }
  fframes.resize( frames.size(), 0.0 ); dfframes.resize( frames.size(), 0.0 );
  if( atoms_to_retrieve.size()>0 ){ 
     requestAtoms( atoms_to_retrieve );
     for(unsigned i=0;i<frames.size();++i) frames[i]->setNumberOfAtoms( atoms_to_retrieve.size() );
  }
  if( args_to_retrieve.size()>0 ){
     std::vector<Value*> req_args;
     interpretArgumentList( args_to_retrieve, req_args ); requestArguments( req_args );
     for(unsigned i=0;i<frames.size();++i) frames[i]->setNumberOfArguments( args_to_retrieve.size() ); 
  }
  // Resize forces array
  forcesToApply.resize( 3*getNumberOfAtoms() + 9 + getNumberOfArguments() );
}

void Mapping::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
} 

Mapping::~Mapping(){
  for(unsigned i=0;i<frames.size();++i) delete frames[i];
}

unsigned Mapping::getPropertyIndex( const std::string& name ) const {
  for(unsigned i=0;i<property.size();++i){
     if( name==property[i] ) return i;
  }
  error("no property with name " + name + " found");
  return 0;
}

double Mapping::getLambda(){
  plumed_merror("lambda is not defined in this mapping type");
}

std::string Mapping::getArgumentName( unsigned& iarg ){
  if( iarg < getNumberOfArguments() ) return getPntrToArgument(iarg)->getName();
  unsigned iatom=iarg - getNumberOfArguments();
  std::string atnum; Tools::convert( getAbsoluteIndex(iatom).serial(),atnum);
  unsigned icomp=iatom%3;
  if(icomp==0) return "pos" + atnum + "x";
  if(icomp==1) return "pos" + atnum + "y";
  return "pos" + atnum + "z"; 
} 

double Mapping::calculateDistanceFunction( const unsigned& ifunc, const bool& squared ){
  // Calculate the distance
  double dd=frames[ifunc]->calculate( getPositions(), getPbc(), getArguments(), squared );

  // Transform distance by whatever
  fframes[ifunc]=transformHD( dd, dfframes[ifunc] );
  return fframes[ifunc];
}

void Mapping::calculateNumericalDerivatives( ActionWithValue* a ){
  if( getNumberOfAtoms()>0 ){
     ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ){
     Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
     for(unsigned j=0;j<getNumberOfComponents();++j){
        for(unsigned i=0;i<getNumberOfArguments();++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
     }
     calculateAtomicNumericalDerivatives( a, getNumberOfArguments() ); 
     for(unsigned j=0;j<getNumberOfComponents();++j){ 
        for(unsigned i=0;i<getNumberOfArguments();++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
     }
  }
}

void Mapping::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned frameno=ider*low_dim.size() + getCurrentTask();
  for(unsigned i=0;i<getNumberOfArguments();++i){
      accumulateDerivative( i, df*dfframes[frameno]*frames[getCurrentTask()]->getArgumentDerivative(i) );
  }
  if( getNumberOfAtoms()>0 ){
      Vector ader; Tensor tmpvir; tmpvir.zero();
      unsigned n=getNumberOfArguments(); 
      for(unsigned i=0;i<getNumberOfAtoms();++i){
          ader=frames[getCurrentTask()]->getAtomDerivative(i);
          accumulateDerivative( n, df*dfframes[frameno]*ader[0] ); n++;
          accumulateDerivative( n, df*dfframes[frameno]*ader[1] ); n++;
          accumulateDerivative( n, df*dfframes[frameno]*ader[2] ); n++;
          tmpvir += -1.0*Tensor( getPosition(i), ader );
      }
      Tensor vir; 
      if( !frames[getCurrentTask()]->getVirial( vir ) ) vir=tmpvir;
      accumulateDerivative( n, df*dfframes[frameno]*vir(0,0) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(0,1) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(0,2) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(1,0) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(1,1) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(1,2) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(2,0) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(2,1) ); n++;
      accumulateDerivative( n, df*dfframes[frameno]*vir(2,2) ); 
  }
}

void Mapping::apply(){
  if( getForcesFromVessels( forcesToApply ) ){
     addForcesOnArguments( forcesToApply );
     if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }
}

}
}



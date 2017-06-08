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
#include "Mapping.h"
#include "vesselbase/Vessel.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceAtoms.h"
#include "tools/PDB.h"
#include "tools/Matrix.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace mapping {

void Mapping::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  vesselbase::ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
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
  bool skipchecks; parseFlag("DISABLE_CHECKS",skipchecks);
  // Setup the object that does the mapping
  mymap = new PointWiseMapping( mtype, skipchecks );

  // Read the properties we require
  if( keywords.exists("PROPERTY") ) {
    std::vector<std::string> property;
    parseVector("PROPERTY",property);
    if(property.size()==0) error("no properties were specified");
    mymap->setPropertyNames( property, false );
  } else {
    std::vector<std::string> property(1);
    property[0]="spath";
    mymap->setPropertyNames( property, true );
  }

  // Open reference file
  std::string reference; parse("REFERENCE",reference);
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations
  bool do_read=true; std::vector<double> weights;
  unsigned nfram=0, wnorm=0., ww;
  while (do_read) {
    PDB mypdb;
    // Read the pdb file
    do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    // Fix argument names
    expandArgKeywordInPDB( mypdb );
    if(do_read) {
      mymap->readFrame( mypdb ); ww=mymap->getWeight( nfram );
      weights.push_back( ww );
      wnorm+=ww; nfram++;
    } else {
      break;
    }
  }
  fclose(fp);

  if(nfram==0 ) error("no reference configurations were specified");
  log.printf("  found %u configurations in file %s\n",nfram,reference.c_str() );
  for(unsigned i=0; i<weights.size(); ++i) weights[i] /= wnorm;
  mymap->setWeights( weights );

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  mymap->getAtomAndArgumentRequirements( atoms, args );
  requestAtoms( atoms ); std::vector<Value*> req_args;
  interpretArgumentList( args, req_args ); requestArguments( req_args );
  // Duplicate all frames (duplicates are used by sketch-map)
  // mymap->duplicateFrameList();
  // fframes.resize( 2*nfram, 0.0 ); dfframes.resize( 2*nfram, 0.0 );
  // plumed_assert( !mymap->mappingNeedsSetup() );
  // Resize all derivative arrays
  // mymap->setNumberOfAtomsAndArguments( atoms.size(), args.size() );
  // Resize forces array
  if( getNumberOfAtoms()>0 ) {
    forcesToApply.resize( 3*getNumberOfAtoms() + 9 + getNumberOfArguments() );
  } else {
    forcesToApply.resize( getNumberOfArguments() );
  }
}

void Mapping::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
}

Mapping::~Mapping() {
  delete mymap;
}

void Mapping::prepare() {
  if( mymap->mappingNeedsSetup() ) {
    // Get the arguments and atoms that are required
    std::vector<AtomNumber> atoms; std::vector<std::string> args;
    mymap->getAtomAndArgumentRequirements( atoms, args );
    requestAtoms( atoms ); std::vector<Value*> req_args;
    interpretArgumentList( args, req_args ); requestArguments( req_args );
    // Duplicate all frames (duplicates are used by sketch-map)
    //mymap->duplicateFrameList();
    // Get the number of frames in the path
    // unsigned nfram=getNumberOfReferencePoints();
    // fframes.resize( 2*nfram, 0.0 ); dfframes.resize( 2*nfram, 0.0 );
    // plumed_assert( !mymap->mappingNeedsSetup() );
    // Resize all derivative arrays
    // mymap->setNumberOfAtomsAndArguments( atoms.size(), args.size() );
    // Resize forces array
    if( getNumberOfAtoms()>0 ) {
      forcesToApply.resize( 3*getNumberOfAtoms() + 9 + getNumberOfArguments() );
    } else {
      forcesToApply.resize( getNumberOfArguments() );
    }
  }
}

unsigned Mapping::getPropertyIndex( const std::string& name ) const {
  return mymap->getPropertyIndex( name );
}

void Mapping::setPropertyValue( const unsigned& iframe, const unsigned& jprop, const double& property ) {
  mymap->setProjectionCoordinate( iframe, jprop, property );
}

double Mapping::getLambda() {
  plumed_merror("lambda is not defined in this mapping type");
}

std::string Mapping::getArgumentName( unsigned& iarg ) {
  if( iarg < getNumberOfArguments() ) return getPntrToArgument(iarg)->getName();
  unsigned iatom=iarg - getNumberOfArguments();
  std::string atnum; Tools::convert( getAbsoluteIndex(iatom).serial(),atnum);
  unsigned icomp=iatom%3;
  if(icomp==0) return "pos" + atnum + "x";
  if(icomp==1) return "pos" + atnum + "y";
  return "pos" + atnum + "z";
}

void Mapping::finishPackSetup( const unsigned& ifunc, ReferenceValuePack& mypack ) const {
  ReferenceConfiguration* myref=mymap->getFrame(ifunc); mypack.setValIndex(0);
  unsigned nargs2=myref->getNumberOfReferenceArguments(); unsigned nat2=myref->getNumberOfReferencePositions();
  if( mypack.getNumberOfAtoms()!=nat2 || mypack.getNumberOfArguments()!=nargs2 ) mypack.resize( nargs2, nat2 );
  if( nat2>0 ) {
    ReferenceAtoms* myat2=dynamic_cast<ReferenceAtoms*>( myref ); plumed_dbg_assert( myat2 );
    for(unsigned i=0; i<nat2; ++i) mypack.setAtomIndex( i, myat2->getAtomIndex(i) );
  }
}

double Mapping::calculateDistanceFunction( const unsigned& ifunc, ReferenceValuePack& myder, const bool& squared ) const {
  // Calculate the distance
  double dd = mymap->calcDistanceFromConfiguration( ifunc, getPositions(), getPbc(), getArguments(), myder, squared );
  // Transform distance by whatever
  double df, ff=transformHD( dd, df ); myder.scaleAllDerivatives( df );
  // And the virial
  if( getNumberOfAtoms()>0 && !myder.virialWasSet() ) {
    Tensor tvir; tvir.zero();
    for(unsigned i=0; i<myder.getNumberOfAtoms(); ++i) tvir +=-1.0*Tensor( getPosition( myder.getAtomIndex(i) ), myder.getAtomDerivative(i) );
    myder.addBoxDerivatives( tvir );
  }
  return ff;
}

ReferenceConfiguration* Mapping::getReferenceConfiguration( const unsigned& ifunc ) {
  return mymap->getFrame( ifunc );
}

void Mapping::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

void Mapping::apply() {
  if( getForcesFromVessels( forcesToApply ) ) {
    addForcesOnArguments( forcesToApply );
    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }
}

}
}



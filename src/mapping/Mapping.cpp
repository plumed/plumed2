/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

  // Read the properties we require
  bool ispath=false;
  if( keywords.exists("PROPERTY") ) {
    std::vector<std::string> propnames; parseVector("PROPERTY",propnames);
    if(propnames.size()==0) error("no properties were specified");
    for(unsigned i=0; i<propnames.size(); ++i) property.insert( std::pair<std::string,std::vector<double> >( propnames[i], std::vector<double>() ) );
  } else {
    property.insert( std::pair<std::string,std::vector<double> >( "spath", std::vector<double>() ) ); ispath=true;
  }

  // Open reference file
  std::string reference; parse("REFERENCE",reference);
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations
  bool do_read=true; unsigned nfram=0; double wnorm=0., ww;
  while (do_read) {
    // Read the pdb file
    PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    // Break if we are done
    if( !do_read ) break ;
    // Check for required properties
    if( !ispath ) {
      double prop;
      for(std::map<std::string,std::vector<double> >::iterator it=property.begin(); it!=property.end(); ++it) {
        if( !mypdb.getArgumentValue( it->first, prop ) ) error("pdb input does not have contain property named " + it->first );
        it->second.push_back(prop);
      }
    } else {
      property.find("spath")->second.push_back( myframes.size()+1 );
    }
    // Fix argument names
    expandArgKeywordInPDB( mypdb );
    // And read the frame
    myframes.emplace_back( metricRegister().create<ReferenceConfiguration>( mtype, mypdb ) );
    if( !mypdb.getArgumentValue( "WEIGHT", ww ) ) ww=1.0;
    weights.push_back( ww ); wnorm+=ww; nfram++;
  }
  fclose(fp);

  if(nfram==0 ) error("no reference configurations were specified");
  log.printf("  found %u configurations in file %s\n",nfram,reference.c_str() );
  for(unsigned i=0; i<weights.size(); ++i) weights[i] = weights[i]/wnorm;

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  for(unsigned i=0; i<myframes.size(); ++i) { myframes[i]->getAtomRequests( atoms, skipchecks ); myframes[i]->getArgumentRequests( args, skipchecks ); }
  std::vector<Value*> req_args; interpretArgumentList( args, req_args );
  if( req_args.size()>0 && atoms.size()>0 ) error("cannot mix atoms and arguments");
  if( req_args.size()>0 ) requestArguments( req_args );
  if( atoms.size()>0 ) {
    log.printf("  found %z atoms in input \n",atoms.size());
    log.printf("  with indices : ");
    for(unsigned i=0; i<atoms.size(); ++i) {
      if(i%25==0) log<<"\n";
      log.printf("%d ",atoms[i].serial());
    }
    log.printf("\n");
    requestAtoms( atoms );
  }
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
  mypack.setValIndex(0);
  unsigned nargs2=myframes[ifunc]->getNumberOfReferenceArguments();
  unsigned nat2=myframes[ifunc]->getNumberOfReferencePositions();
  if( mypack.getNumberOfAtoms()!=nat2 || mypack.getNumberOfArguments()!=nargs2 ) mypack.resize( nargs2, nat2 );
  if( nat2>0 ) {
    ReferenceAtoms* myat2=dynamic_cast<ReferenceAtoms*>( myframes[ifunc].get() ); plumed_dbg_assert( myat2 );
    for(unsigned i=0; i<nat2; ++i) mypack.setAtomIndex( i, myat2->getAtomIndex(i) );
  }
}

double Mapping::calculateDistanceFunction( const unsigned& ifunc, ReferenceValuePack& myder, const bool& squared ) const {
  // Calculate the distance
  double dd = myframes[ifunc]->calculate( getPositions(), getPbc(), getArguments(), myder, squared );
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
  return myframes[ifunc].get();
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
  if( getForcesFromVessels( forcesToApply ) ) {
    addForcesOnArguments( forcesToApply );
    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }
}

}
}



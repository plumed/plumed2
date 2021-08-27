/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "SetupReferenceBase.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace setup {

void SetupReferenceBase::registerKeywords( Keywords& keys ) {
  // ActionSetup is not registered as we don't want to remove LABEL
  // ActionWithValue is not registered as nothing useful is in that register
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys );
  keys.use("ARG"); keys.reset_style("ARG","optional");
}

SetupReferenceBase::SetupReferenceBase(const ActionOptions&ao):
Action(ao),
ActionSetup(ao),
ActionAtomistic(ao),
ActionWithArguments(ao),
ActionWithValue(ao),
hasatoms(false)
{
}

SetupReferenceBase::~SetupReferenceBase() {
  if( hasatoms ) { atoms.removeVirtualAtom( this ); } // atoms.removeGroup( getLabel() ); }
}

void SetupReferenceBase::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void SetupReferenceBase::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void SetupReferenceBase::calculateNumericalDerivatives( ActionWithValue* a ) {
  error("this should never be called");
}

AtomNumber SetupReferenceBase::getAtomNumber( const AtomNumber& anum ) const {
  for(unsigned i=0;i<mygroup.size();++i) {
      if( anum==mygroup[i] ) return myindices[i];
  }
  plumed_error(); return myindices[0];
}

void SetupReferenceBase::getNatomsAndNargs( unsigned& natoms, unsigned& nargs ) const {
  // Get the number of atoms
  natoms=0;
  for(unsigned i=0;i<myindices.size();++i) { 
      if( myindices[i].serial()>natoms ) natoms = myindices[i].serial();
  }
  // Get the number of arguments
  nargs=0; if( getNumberOfComponents()>0 ) nargs = getPntrToOutput(0)->getNumberOfValues();
}

void SetupReferenceBase::getAtomsFromReference( const unsigned& npos, std::vector<double>& masses, 
                                                std::vector<double>& charges, std::vector<Vector>& positions ) const {
  for(unsigned i=0;i<myindices.size();++i) {
      masses[npos + myindices[i].index()] = atoms.getVatomMass(mygroup[i]);
      charges[npos + myindices[i].index()] = atoms.getVatomCharge(mygroup[i]);
      positions[npos + myindices[i].index()] = atoms.getVatomPosition(mygroup[i]);
  }
}

void SetupReferenceBase::getReferenceConfiguration( std::vector<double>& ref ) const {
  unsigned n=0;
  if( getNumberOfComponents()>0 ) {
      Value* myval=getPntrToOutput(0);
      unsigned nvals = myval->getNumberOfValues();
      for(unsigned i=0;i<nvals;++i) { ref[n] = myval->get(i); n++; }
  }
  for(unsigned i=0;i<mygroup.size();++i) {
      Vector cpos( atoms.getVatomPosition(mygroup[i]) );
      for(unsigned k=0;k<3;++k) { ref[n] = cpos[k]; n++; } 
  }
}

void SetupReferenceBase::setReferenceConfiguration( std::vector<double>& ref ) {
  unsigned n=0;
  if( getNumberOfComponents()>0 ) {
      Value* myval=getPntrToOutput(0);
      unsigned nvals = myval->getNumberOfValues();
      for(unsigned i=0;i<nvals;++i) { myval->set( i, ref[n]) ; n++; }
  }
  for(unsigned i=0;i<mygroup.size();++i) {
      Vector cpos( atoms.getVatomPosition(mygroup[i]) );
      for(unsigned k=0;k<3;++k) { cpos[k] = ref[n]; n++; }
      atoms.setVatomPosition( mygroup[i], cpos ); 
  }
}

void SetupReferenceBase::displaceReferenceConfiguration( const double& val, const std::vector<double>& dir ) {
  unsigned n=0;
  if( getNumberOfComponents()>0 ) {
      Value* myval=getPntrToOutput(0);
      unsigned nvals = myval->getNumberOfValues();
      for(unsigned i=0;i<nvals;++i) { myval->set( i, myval->get(i) + val*dir[n]) ; n++; }
  }
  for(unsigned i=0;i<mygroup.size();++i) {
      Vector cpos( atoms.getVatomPosition(mygroup[i]) );
      for(unsigned k=0;k<3;++k) { cpos[k] = cpos[k] + val*dir[n+k]; }
      atoms.setVatomPosition( mygroup[i], cpos ); n += 3;
  }
}

std::string SetupReferenceBase::getArgName( const unsigned& k ) const {
  unsigned nn=0;
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      unsigned nvals = getPntrToArgument(i)->getNumberOfValues();
      if( k<nn+nvals ) {
          if( getPntrToArgument(i)->getRank()==0 ) return getPntrToArgument(i)->getName();
          if( getPntrToArgument(i)->getRank()==1 ) {
              std::string num; Tools::convert( k-nn+1, num );
              return getPntrToArgument(i)->getName() + "." + num;
          }
          if( getPntrToArgument(i)->getRank()==2 ) {
              unsigned rown = std::floor( (k-nn) / getPntrToArgument(i)->getShape()[1] );
              unsigned coln = k - nn - rown*getPntrToArgument(i)->getShape()[1]; 
              std::string rnum, cnum; Tools::convert( rown + 1, rnum ); Tools::convert( coln + 1,  cnum );
              return getPntrToArgument(i)->getName() + "." + rnum + "." + cnum;
          }
      }
      nn += nvals;
  }
  plumed_error(); return "";
}

}
}

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
#include "PointWiseMapping.h"
#include "FakeFrame.h"

namespace PLMD {

PointWiseMapping::PointWiseMapping( const std::string& type, const bool& checksoff ):
MultiReferenceBase(type,checksoff)
{
}

void PointWiseMapping::setPropertyNames( const std::vector<std::string>& prop, const bool isp ){
  property.resize( prop.size() ); ispath=isp;
  for(unsigned i=0;i<prop.size();++i) property[i]=prop[i];
}

void PointWiseMapping::readRestOfFrame(){
  plumed_assert( property.size()>0 );

  std::vector<double> labelvals;
  if( !ispath ){
      labelvals.resize( property.size() );
      for(unsigned i=0;i<property.size();++i) parse( property[i], labelvals[i] );
  } else {
      labelvals.resize(1);
      labelvals[0]=static_cast<double>( frames.size() ); 
  }
  low_dim.push_back( labelvals ); 
} 

void PointWiseMapping::getAtomAndArgumentRequirements( std::vector<AtomNumber>& atoms, std::vector<std::string>& args ){
  plumed_assert( atoms.size()==0 && args.size()==0 );
  for(unsigned i=0;i<frames.size();++i){
      frames[i]->getAtomRequests( atoms ); 
      frames[i]->getArgumentRequests( args );
  }
}

void PointWiseMapping::duplicateFrameList(){
  unsigned nframes=frames.size();
  for(unsigned i=0;i<nframes;++i){
     frames.push_back( new FakeFrame( ReferenceConfigurationOptions("fake") ) );
  }
}

void PointWiseMapping::setNumberOfAtomsAndArguments( const unsigned& natoms, const unsigned& nargs ){
  for(unsigned i=0;i<frames.size();++i){
      frames[i]->setNumberOfAtoms( natoms );
      frames[i]->setNumberOfArguments( nargs );
  }
}

unsigned PointWiseMapping::getPropertyIndex( const std::string& name ) const {
  for(unsigned i=0;i<property.size();++i){
     if( name==property[i] ) return i;
  }
  plumed_merror("no property with name " + name + " found");
  return 0;
}

}

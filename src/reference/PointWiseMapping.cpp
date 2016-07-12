/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "PointWiseMapping.h"
#include "Direction.h"
#include "FakeFrame.h"

namespace PLMD {

PointWiseMapping::PointWiseMapping( const std::string& type, const bool& checksoff ):
MultiReferenceBase(type,checksoff)
{
  ispath=false;
}

void PointWiseMapping::setPropertyNames( const std::vector<std::string>& prop, const bool isp ){
  property.resize( prop.size() ); ispath=isp;
  for(unsigned i=0;i<prop.size();++i) property[i]=prop[i];
}

void PointWiseMapping::readRestOfFrame(){
  Direction* tdir = dynamic_cast<Direction*>( getFrame( getNumberOfReferenceFrames() - 1 ) );
  if( tdir ) plumed_merror("cannot use directions in mapping");
  plumed_dbg_assert( property.size()>0 );

  std::vector<double> labelvals;
  if( !ispath ){
      labelvals.resize( property.size() );
      for(unsigned i=0;i<property.size();++i) parse( property[i], labelvals[i] );
  } else {
      labelvals.resize(1);
      labelvals[0]=static_cast<double>( frames.size() ); 
  }
  low_dim.push_back( labelvals ); 
  plumed_dbg_assert( low_dim.size()==getNumberOfReferenceFrames() );
}

void PointWiseMapping::clearRestOfData(){
  low_dim.resize(0);
} 

void PointWiseMapping::resizeRestOfFrame(){
  plumed_dbg_assert( property.size()>0 );
  std::vector<double> labelvals( property.size() );
  low_dim.push_back( labelvals );
  plumed_dbg_assert( low_dim.size()==getNumberOfReferenceFrames() );
}

void PointWiseMapping::duplicateFrameList(){
  unsigned nframes=frames.size();
  for(unsigned i=0;i<nframes;++i){
     frames.push_back( new FakeFrame( ReferenceConfigurationOptions("fake") ) );
  }
}

unsigned PointWiseMapping::getPropertyIndex( const std::string& name ) const {
  for(unsigned i=0;i<property.size();++i){
     if( name==property[i] ) return i;
  }
  plumed_merror("no property with name " + name + " found");
  return 0;
}

void PointWiseMapping::print( const std::string& method, const double & time, OFile& afile, 
                              const std::string& fmt, const double& lunits ){
  std::string descr2, descr="DESCRIPTION: results from %s analysis performed at time " + fmt +"\n";
  afile.printf(descr.c_str(), method.c_str(), time );
  if(fmt.find("-")!=std::string::npos){
     descr="REMARK WEIGHT=" + fmt + " %s=" + fmt + " "; descr2="%s=" + fmt;
  } else {
     // This ensures numbers are left justified (i.e. next to the equals sign
     std::size_t psign=fmt.find("%");
     plumed_assert( psign!=std::string::npos );
     descr="REMARK WEIGHT=%-" + fmt.substr(psign+1) + " %s=%-" + fmt.substr(psign+1) + " "; 
     descr2="%s=%-" + fmt.substr(psign+1);
  }
  for(unsigned i=0;i<frames.size();++i){
      afile.printf(descr.c_str(), frames[i]->getWeight(), property[0].c_str(), low_dim[i][0] );
      for(unsigned j=1;j<property.size();++j) afile.printf(descr2.c_str(), property[j].c_str(), low_dim[i][j]);
      afile.printf("\n"); 
      frames[i]->print( afile, fmt, lunits );
  }
} 

}

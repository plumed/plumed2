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
#include "MultiReferenceBase.h"
#include "MetricRegister.h"

namespace PLMD {

MultiReferenceBase::MultiReferenceBase( const std::string& type, const bool& checksoff ):
wasSet(false),
skipchecks(checksoff),
mtype(type)
{
  if(checksoff) plumed_assert( mtype.length()==0 );
}

MultiReferenceBase::~MultiReferenceBase(){
  for(unsigned i=0;i<frames.size();++i) delete frames[i];
}

void MultiReferenceBase::readFrame( PDB& mypdb ){
  wasSet=true; 
  // If skipchecks are enabled metric types must be specified in the input file
  ReferenceConfiguration* mymsd=metricRegister().create<ReferenceConfiguration>( mtype, mypdb );
  // Read in the weight we are allocating to this particular point in space
  double ww; mymsd->getWeight();
  // Save everything
  weights.push_back(ww); frames.push_back( mymsd );
  // Double check
  plumed_assert( weights.size()==frames.size() );
  // Do reading in derived class
  readRestOfFrame();
  // Check readin was succesfull
  mymsd->checkRead();
}

}

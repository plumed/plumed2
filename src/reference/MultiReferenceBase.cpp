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
#include "tools/Communicator.h"
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

void MultiReferenceBase::clearFrames(){
  for(unsigned i=0;i<frames.size();++i) delete frames[i];
  frames.resize(0); weights.resize(0);
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

void MultiReferenceBase::copyFrame( ReferenceConfiguration* frameToCopy, const double& weight ){
  // Create a reference configuration of the appropriate type
  ReferenceConfiguration* mymsd=metricRegister().create<ReferenceConfiguration>( frameToCopy->getName() );
  // Copy names of arguments and and indexes
  mymsd->setNamesAndAtomNumbers( frameToCopy->getAbsoluteIndexes(), frameToCopy->getArgumentNames() );
  // Copy reference positions, reference arguments and reference metric
  mymsd->setReference( frameToCopy->getReferencePositions(), frameToCopy->getReferenceArguments(), frameToCopy->getReferenceMetric() );
  // Easy bit - copy the weight and store the frame
  frames.push_back( mymsd ); weights.push_back( weight );
  plumed_dbg_assert( weight.size()==frames.size() );
  // This resizes the low dim array
  resizeRestOfFrame();
}

void MultiReferenceBase::calculateAllDistances( const Pbc& pbc, const std::vector<Value*> vals, Communicator& comm, Matrix<double>& distances ){
  distances=0.0;
  unsigned k=0, size=comm.Get_size(), rank=comm.Get_rank(); 
  for(unsigned i=1;i<frames.size();++i){
      for(unsigned j=0;j<i;++j){
          if( (k++)%size!=rank ) continue;         
          distances(i,j) = distances(j,i) = distance( pbc, vals, frames[i], frames[j], false );
      }
  }
//  comm.Sum( distances );
}

}

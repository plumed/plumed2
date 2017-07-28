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

MultiReferenceBase::~MultiReferenceBase() {
  for(unsigned i=0; i<frames.size(); ++i) delete frames[i];
}

void MultiReferenceBase::clearFrames() {
  for(unsigned i=0; i<frames.size(); ++i) delete frames[i];
  frames.resize(0);
  clearRestOfData();
}

void MultiReferenceBase::readFrame( PDB& mypdb ) {
  wasSet=true;
  // If skipchecks are enabled metric types must be specified in the input file
  ReferenceConfiguration* mymsd=metricRegister().create<ReferenceConfiguration>( mtype, mypdb );
  // Save everything
  frames.push_back( mymsd );
  // Do reading in derived class
  readRestOfFrame();
  // Check readin was succesfull
  mymsd->checkRead();
}

void MultiReferenceBase::getAtomAndArgumentRequirements( std::vector<AtomNumber>& atoms, std::vector<std::string>& args ) {
  plumed_assert( atoms.size()==0 && args.size()==0 );
  for(unsigned i=0; i<frames.size(); ++i) {
    frames[i]->getAtomRequests( atoms );
    frames[i]->getArgumentRequests( args );
  }
}

// void MultiReferenceBase::setNumberOfAtomsAndArguments( const unsigned& natoms, const unsigned& nargs ){
//   for(unsigned i=0;i<frames.size();++i){
//       frames[i]->setNumberOfAtoms( natoms );
//       frames[i]->setNumberOfArguments( nargs );
//   }
// }

void MultiReferenceBase::copyFrame( ReferenceConfiguration* frameToCopy ) {
  // Create a reference configuration of the appropriate type
  ReferenceConfiguration* mymsd=metricRegister().create<ReferenceConfiguration>( frameToCopy->getName() );
  // Copy names of arguments and and indexes
  mymsd->setNamesAndAtomNumbers( frameToCopy->getAbsoluteIndexes(), frameToCopy->getArgumentNames() );
  // Copy reference positions, reference arguments and reference metric
  mymsd->setReferenceConfig( frameToCopy->getReferencePositions(), frameToCopy->getReferenceArguments(), frameToCopy->getReferenceMetric() );
  // Copy weight
  mymsd->setWeight( frameToCopy->getWeight() );
  // Easy bit - copy the frame
  frames.push_back( mymsd );
  // This resizes the low dim array
  resizeRestOfFrame();
}

void MultiReferenceBase::setWeights( const std::vector<double>& weights ) {
  plumed_assert( weights.size()==frames.size() );
  for(unsigned i=0; i<weights.size(); ++i) frames[i]->setWeight( weights[i] );
}


void MultiReferenceBase::calculateAllDistances( const Pbc& pbc, const std::vector<Value*> & vals, Communicator& comm, Matrix<double>& distances, const bool& squared ) {
  distances=0.0;
  unsigned k=0, size=comm.Get_size(), rank=comm.Get_rank();
  for(unsigned i=1; i<frames.size(); ++i) {
    for(unsigned j=0; j<i; ++j) {
      if( (k++)%size!=rank ) continue;
      distances(i,j) = distances(j,i) = distance( pbc, vals, frames[i], frames[j], squared );
    }
  }
  comm.Sum( distances );
}

}

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
#ifndef __PLUMED_reference_MultiReferenceBase_h
#define __PLUMED_reference_MultiReferenceBase_h

#include "ReferenceConfiguration.h"
#include "tools/Matrix.h"

namespace PLMD {

class MultiReferenceBase {
private:
/// Everything has been set
  bool wasSet;
/// Skip all checking allows users to do really dodgy stuff :-)
  bool skipchecks;
/// The type of metric we are using
  std::string mtype;
protected:
/// These are the configurations that serve as references
  std::vector<ReferenceConfiguration*> frames;
/// Read something from the last frame
  template <class T>
  void parse(const std::string& key, T& val );
public:
  MultiReferenceBase( const std::string& type, const bool& checksoff );
/// Destructor deletes all polymorphic pointers
  virtual ~MultiReferenceBase();
/// Delete all the data in the reference object
  void clearFrames();
  virtual void clearRestOfData() {};
/// Read a frame from the input
  void readFrame( PDB& pdb );
/// Find what is required of us from the reference frames
  void getAtomAndArgumentRequirements( std::vector<AtomNumber>& atoms, std::vector<std::string>& args );
/// Finish setup of frames
//  void setNumberOfAtomsAndArguments( const unsigned& natoms, const unsigned& nargs );
/// Do additional reading required by derived class
  virtual void readRestOfFrame() {}
/// Do additional resizing required by derived class
  virtual void resizeRestOfFrame() {}
/// Return the size of the frames vector
  unsigned getNumberOfReferenceFrames() const ;
/// Calculate the distance from one of the reference points
  double calcDistanceFromConfiguration( const unsigned& ifunc, const std::vector<Vector>& pos, const Pbc& pbc,
                                        const std::vector<Value*>& arg, ReferenceValuePack& myder, const bool& squared ) const ;
/// Return the ith reference frame
  ReferenceConfiguration* getFrame( const unsigned& iframe );
/// Return a reference to all the reference frames
  std::vector<ReferenceConfiguration*>& getReferenceConfigurations();
/// Copy a reference configuration into the multi reference object
  void copyFrame( ReferenceConfiguration* frameToCopy );
/// Set the weight of the ith frame
  void setWeights( const std::vector<double>& ww );
/// Retrieve the weight of one of the frames
  double getWeight( const unsigned& ifram ) const ;
/// Calculate the distances between all the frames and store in a matrix
  void calculateAllDistances( const Pbc& pbc, const std::vector<Value*> & vals, Communicator& comm, Matrix<double>& distances, const bool& squared );
};

template <class T>
void MultiReferenceBase::parse(const std::string& key, T& val ) {
  frames[frames.size()-1]->parse(key,val);
}

inline
double MultiReferenceBase::calcDistanceFromConfiguration( const unsigned& ifunc, const std::vector<Vector>& pos, const Pbc& pbc,
    const std::vector<Value*>& arg, ReferenceValuePack& myder, const bool& squared ) const {
  return frames[ifunc]->calculate( pos, pbc, arg, myder, squared );
}

inline
unsigned MultiReferenceBase::getNumberOfReferenceFrames() const {
  return frames.size();
}

inline
double MultiReferenceBase::getWeight( const unsigned& ifram ) const {
  plumed_dbg_assert( ifram<frames.size() );
  return frames[ifram]->getWeight();
}

inline
ReferenceConfiguration* MultiReferenceBase::getFrame( const unsigned& iframe ) {
  plumed_dbg_assert( iframe<frames.size() );
  return frames[iframe];
}

inline
std::vector<ReferenceConfiguration*>& MultiReferenceBase::getReferenceConfigurations() {
  return frames;
}

}
#endif

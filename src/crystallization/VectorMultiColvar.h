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
#ifndef __PLUMED_crystallization_VectorMultiColvar_h
#define __PLUMED_crystallization_VectorMultiColvar_h

#include "tools/Matrix.h"
#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"

namespace PLMD {
namespace crystallization {

class VectorMultiColvar : public multicolvar::MultiColvarBase {
  friend class OrientationSphere;
  friend class VolumeGradientBase;
private:
/// Are we storing the director of the vector of the vector
  bool store_director;
/// How many components does the vector have
  unsigned ncomponents;
/// These are tempory vectors that are used to store values and directors
  std::vector<double> vv1, vv2;
protected:
/// Set the dimensionality of the vector
  void setVectorDimensionality( const unsigned& );
/// Used in vector average to add forces from vector the the forces from here
  void addForcesOnAtoms( const std::vector<double>& inforces );
public:
  static void registerKeywords( Keywords& keys );
  explicit VectorMultiColvar(const ActionOptions&);
  ~VectorMultiColvar() {}
/// The norm of a vector is not periodic
  virtual bool isPeriodic() { return false; }
/// Calculate the multicolvar
//  double doCalculation( const unsigned& taskIndex, multicolvar::AtomValuePack& myatoms ) const ;
/// This shouldn't do anything
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
/// Calculate the vector
  virtual void calculateVector( multicolvar::AtomValuePack& myatoms ) const=0;
/// Get the number of components in the vector
  unsigned getNumberOfComponentsInVector() const ;
/// Get the number of quantities we are calculating per step
  unsigned getNumberOfQuantities() const ;
/// Can we differentiate the orientation - yes we can the multicolvar is a vector
  bool hasDifferentiableOrientation() const { return true; }
///  This makes sure we are not calculating the director when we do LocalAverage
  virtual void doNotCalculateDirector();
/// This does normalizeing of vectors for storeDataVessel
  virtual void normalizeVector( std::vector<double>& vals ) const ;
  virtual void normalizeVectorDerivatives( MultiValue& myvals ) const ;
};

inline
unsigned VectorMultiColvar::getNumberOfComponentsInVector() const {
  return ncomponents;
}

inline
unsigned VectorMultiColvar::getNumberOfQuantities() const {
  return 2 + ncomponents;
}

}
}
#endif


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "multicolvar/MultiColvar.h"

namespace PLMD {
namespace crystallization {

class VectorMultiColvar : public multicolvar::MultiColvar {
friend class OrientationSphere;
friend class VolumeGradientBase;
private:
/// Are we storing the director of the vector of the vector
  bool store_director;
/// Used to make sure central atom position is only calculated
/// once when using orientation sphere
  bool firstcall;
/// How many components does the vector have
  unsigned ncomponents;
/// These are tempory vectors that are used to store values and directors
  std::vector<double> vv1, vv2;
protected:
/// Set the dimensionality of the vector
  void setVectorDimensionality( const unsigned&, const int& );
//  /// Add some value to the ith component of the vector
//    void addComponent( const unsigned&, const double& );
//  /// Get the ith component
//    double getComponent( const unsigned& ) const ;
//  /// Set the ith component
//    void setComponent( const unsigned&, const double& );
//  /// Add derivatives of ith component of vector with repect to jth atom
//    void addAtomsDerivative( const unsigned&, const unsigned&, const Vector& );
//  /// Add atomic derivatives to all components of matrix (note iatom is treated literally here - cf above)
//    void addAtomDerivativeToAllRealComponents( const unsigned& iatom, const std::vector<double>& vec, const Vector& avec );
//  /// Add derivatives of ith component of vector with respect to the box 
//    void addBoxDerivatives( const unsigned&, const Tensor& );
//  /// Add box derivatives to all components of matrix
//    void addBoxDerivativesToAllRealComponents( const std::vector<double>& vec, const Tensor& avec );
//  /// Add some value to the imaginary part of the ith component of the vector
//    void addImaginaryComponent( const unsigned&, const double& );
//  /// Get the ith component
//    double getImaginaryComponent( const unsigned& ) const ;
//  /// Set the ith component
//    void setImaginaryComponent( const unsigned&, const double& );
//  /// Add derivatives of the imaginary part of the ith component of vector with repect to jth atom
//    void addImaginaryAtomsDerivative( const unsigned&, const unsigned&, const Vector& );
//  /// Add atomic derivatives to all components of matrix (note iatom is treated literally here - cf above)
//    void addAtomDerivativeToAllImagComponents( const unsigned& iatom, const std::vector<double>& vec, const Vector& avec );
//  /// Add derivatives of the imaginary part of the ith component of vector with respect to the box 
//    void addImaginaryBoxDerivatives( const unsigned&, const Tensor& );
//  /// Add box derivatives to all components of matrix
//    void addBoxDerivativesToAllImagComponents( const std::vector<double>& vec, const Tensor& avec );
//  /// This can be used to accumulate derivative from a store of vectors
//    void accumulateDerivativesFromVector( const unsigned& ivec, const unsigned& base_cv_no, const double& weight, StoreVectorsVessel* vectors );
/// Used in vector average to add forces from vector the the forces from here
  void addForcesOnAtoms( const std::vector<double>& inforces );
public:
  static void registerKeywords( Keywords& keys );
  VectorMultiColvar(const ActionOptions&);
  ~VectorMultiColvar(){}
/// The norm of a vector is not periodic
  virtual bool isPeriodic(){ return false; }
/// Calculate the multicolvar
  double doCalculation( const unsigned& taskIndex, multicolvar::AtomValuePack& myatoms );
/// This shouldn't do anything
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ){ plumed_error(); }
/// Calculate the vector
  virtual void calculateVector( multicolvar::AtomValuePack& myatoms )=0;
/// Get the number of components in the vector
  unsigned getNumberOfComponentsInVector() const ;
/// Get the number of quantities we are calculating per step
  unsigned getNumberOfQuantities();
/// Create places to store the data
//  vesselbase::StoreDataVessel* buildDataStashes( const bool& allow_wcutoff, const double& wtol );
/// Get the vector
//  void getValueForTask( const unsigned& iatom, std::vector<double>& vals );
/// Used to accumulate values
//  void addWeightedValueDerivatives( const unsigned& iatom, const unsigned& base_cv_no, const double& weight, multicolvar::MultiColvarFunction* func );
/// Used for calculating weighted averages
//  void finishWeightedAverageCalculation( multicolvar::MultiColvarFunction* func );
/// Used in functions to add derivatives to the orientation vector
//  void addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
//                                        const std::vector<double>& der, multicolvar::MultiColvarFunction* func );
/// Can we differentiate the orientation - yes we can the multicolvar is a vector
  bool hasDifferentiableOrientation() const { return true; }
///  This makes sure we are not calculating the director when we do LocalAverage
  virtual void doNotCalculateDirector();
/// Used by ActionVolume
//  void copyElementsToBridgedColvar( multicolvar::BridgedMultiColvarFunction* func );
};

inline
unsigned VectorMultiColvar::getNumberOfComponentsInVector() const {
  return ncomponents; 
}

// inline
// void VectorMultiColvar::addComponent( const unsigned& icomp, const double& val ){
//   plumed_dbg_assert( icomp<ncomponents );
//   addElementValue( 5 + icomp, val );
// }
// 
// inline
// void VectorMultiColvar::setComponent( const unsigned& icomp, const double& val ){
//   plumed_dbg_assert( icomp<ncomponents );
//   setElementValue( 5 + icomp, val );
// } 
//   
// inline
// double VectorMultiColvar::getComponent( const unsigned& icomp ) const {
//   plumed_dbg_assert( icomp<ncomponents );
//   return getElementValue( 5 + icomp );
// } 
// 
// 
// inline
// void VectorMultiColvar::addAtomsDerivative( const unsigned& icomp, const unsigned& jatom, const Vector& der ){
//   plumed_dbg_assert( icomp<ncomponents && jatom<getNAtoms() );
//   MultiColvarBase::addAtomsDerivatives( 5 + icomp, current_atoms[jatom], der );
// }

// inline
// void VectorMultiColvar::addBoxDerivatives( const unsigned& icomp, const Tensor& vir ){
//   plumed_dbg_assert( icomp<ncomponents );
//   MultiColvarBase::addBoxDerivatives( 5 + icomp, vir );
// }
// 
// inline
// void VectorMultiColvar::addImaginaryComponent( const unsigned& icomp, const double& val ){
//   plumed_dbg_assert( icomp<ncomponents && complexvec );
//   addElementValue( 5 + ncomponents + icomp, val );
// }
// 
// inline
// void VectorMultiColvar::setImaginaryComponent( const unsigned& icomp, const double& val ){
//   plumed_dbg_assert( icomp<ncomponents && complexvec );
//   setElementValue( 5 + ncomponents + icomp, val );
// }
// 
// inline 
// double VectorMultiColvar::getImaginaryComponent( const unsigned& icomp ) const {
//   plumed_dbg_assert( icomp<ncomponents && complexvec );
//   return getElementValue( 5 + ncomponents + icomp );
// } 
// 
// inline
// void VectorMultiColvar::addImaginaryAtomsDerivative( const unsigned& icomp, const unsigned& jatom, const Vector& der){
//   plumed_dbg_assert( icomp<ncomponents && complexvec && jatom<getNAtoms() );
//   MultiColvarBase::addAtomsDerivatives( 5 + ncomponents + icomp, current_atoms[jatom], der );
// }
// 
// inline
// void VectorMultiColvar::addImaginaryBoxDerivatives( const unsigned& icomp, const Tensor& vir ){
//   plumed_dbg_assert( icomp<ncomponents && complexvec );
//   MultiColvarBase::addBoxDerivatives( 5 + ncomponents + icomp, vir ); 
// }

inline
unsigned VectorMultiColvar::getNumberOfQuantities(){
  return 2 + ncomponents;
}

}
}
#endif


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_reference_ReferenceAtoms_h
#define __PLUMED_reference_ReferenceAtoms_h

#include "ReferenceConfiguration.h"

namespace PLMD {

class Pbc;

/// \ingroup TOOLBOX
/// In many applications (e.g. paths, fields, property maps) it is necessary to calculate
/// the distance between two configurations.  These distances can be calculated in a variety of
/// different ways.  For instance, one can assert that the distance between the two configuration
/// is the distance one would have to move all the atoms to transform configuration 1 into configuration
/// 2. Alternatively, one could calculate the values of a large set of collective coordinates in the two
/// configurations and then calculate the Euclidean distances between these two points in the resulting
/// high-dimensional vector space.  Lastly, one can combine these two forms of distance calculation to calculate
/// a hybrid distance.  Plumed allows one to use all these forms of distance calculations and also to implement
/// new forms of distance.  You should inherit from this class if your distance involves reference atomic positions.
/// This class and \ref PLMD::ReferenceArguments mirror the functionalities in and \ref PLMD::ActionAtomistic
/// and \ref PLMD::ActionWithArguments respectively but for distances.

class ReferenceAtoms :
  virtual public ReferenceConfiguration
{
  friend class Direction;
  friend class SingleDomainRMSD;
  friend class MultiDomainRMSD;
  friend class ReferenceConfiguration;
private:
/// This flag tells us if the user has disabled checking of the input in order to
/// do fancy paths with weird inputs
  bool checks_were_disabled;
/// The atoms to be used to align the instantaneous atomic positions
/// to the reference configuration
  std::vector<double> align;
/// The atoms to be used to calculate the distance the atoms have moved
/// from the reference configuration
  std::vector<double> displace;
/// The positions of the atoms in the reference configuration
  std::vector<Vector> reference_atoms;
/// The indices of the atoms in the pdb file
  std::vector<AtomNumber> indices;
/// The indeces for setting derivatives
  std::vector<unsigned> atom_der_index;
protected:
/// Read in the atoms from the pdb file
  void readAtomsFromPDB( const PDB&, const bool allowblocks=false );
/// Add atom indices to list
  void setAtomIndices( const std::vector<AtomNumber>& atomnumbers );
/// Read a list of atoms from the pdb input file
  bool parseAtomList( const std::string&, std::vector<unsigned>& );
/// Get the position of the ith atom
  Vector getReferencePosition( const unsigned& iatom ) const ;
/// Add derivatives to iatom th atom in list
//  void addAtomicDerivatives( const unsigned& , const Vector& );
/// Get the atomic derivatives on the ith atom in the list
//  Vector retrieveAtomicDerivatives( const unsigned& ) const ;
/// Add derivatives to the viral
//  void addBoxDerivatives( const Tensor& );
/// This does the checks that are always required
  void singleDomainRequests( std::vector<AtomNumber>&, bool disable_checks );
public:
  explicit ReferenceAtoms( const ReferenceConfigurationOptions& ro );
/// This returns the number of reference atom positions
  unsigned getNumberOfReferencePositions() const override;
/// Get the reference positions
  const std::vector<Vector> & getReferencePositions() const override;
/// This allows us to use a single pos array with RMSD objects using different atom indexes
  unsigned getAtomIndex( const unsigned& ) const ;
/// Get the atoms required (additional checks are required when we have multiple domains)
  void getAtomRequests( std::vector<AtomNumber>&, bool disable_checks=false ) override;
/// Set the positions of the reference atoms
  virtual void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in )=0;
/// Return all atom indexes
  const std::vector<AtomNumber>& getAbsoluteIndexes() override;
/// This returns how many atoms there should be
  unsigned getNumberOfAtoms() const ;
/// Displace the positions of the reference atoms a bit
  void displaceReferenceAtoms( const double& weight, const std::vector<Vector>& dir );
/// Extract a displacement from a position in space
  virtual void extractAtomicDisplacement( const std::vector<Vector>& pos, std::vector<Vector>& direction ) const {
    plumed_error();
  }
/// Project the displacement on a vector
  virtual double projectAtomicDisplacementOnVector( const bool& normalized, const std::vector<Vector>& eigv, ReferenceValuePack& mypack ) const {
    plumed_error(); return 1;
  }
/// Get the vector of alignment weights
  const std::vector<double> & getAlign() const ;
/// Get the vector of displacement weights
  const std::vector<double> & getDisplace() const ;
};

inline
const std::vector<double> & ReferenceAtoms::getAlign() const {
  return align;
}

inline
const std::vector<double> & ReferenceAtoms::getDisplace() const {
  return displace;
}

inline
unsigned ReferenceAtoms::getNumberOfReferencePositions() const {
  plumed_dbg_assert( atom_der_index.size()==reference_atoms.size() );
  return reference_atoms.size();
}

inline
unsigned ReferenceAtoms::getNumberOfAtoms() const {
  return atom_der_index.size(); // reference_atoms.size();
}

inline
unsigned ReferenceAtoms::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_der_index.size() );
  plumed_dbg_assert( atom_der_index[iatom]<reference_atoms.size() );
  return atom_der_index[iatom];
}

inline
Vector ReferenceAtoms::getReferencePosition( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<reference_atoms.size() );
  return reference_atoms[iatom];
}

inline
const std::vector<Vector> & ReferenceAtoms::getReferencePositions() const {
  return reference_atoms;
}

inline
const std::vector<AtomNumber>& ReferenceAtoms::getAbsoluteIndexes() {
  return indices;
}


}
#endif


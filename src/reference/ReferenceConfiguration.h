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
#ifndef __PLUMED_reference_ReferenceConfiguration_h
#define __PLUMED_reference_ReferenceConfiguration_h

#include <vector>
#include <string>
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Tools.h"
#include "tools/Exception.h"
#include "ReferenceValuePack.h"
#include "tools/Matrix.h"

namespace PLMD {

class Value;
class Pbc;
class OFile;
class PDB;
class SetupMolInfo;

/// \ingroup TOOLBOX
/// Abstract base class for calculating the distance from a reference configuration.
/// A reference configuration can either have a particular set of atoms in a particular
/// given configuration or it can be that a particular set of colvars have a particular
/// set of values.  It could also be a combination of both.  To allow all the posible
/// permutations and in order make it easy to add new ways of calculating the distance
/// we have implemented this using polymorphism and multiple inheritance.

class Direction;

class ReferenceConfigurationOptions {
  friend class ReferenceConfiguration;
private:
  std::string tt;
public:
  explicit ReferenceConfigurationOptions( const std::string& type );
  bool usingFastOption() const ;
  std::string getMultiRMSDType() const ;
};

/// \ingroup INHERIT
/// Abstract base class for calculating the distance from a reference configuration.
/// A reference configuration can either have a particular set of atoms in a particular
/// given configuration or it can be that a particular set of colvars have a particular
/// set of values.  It could also be a combination of both.  To allow all the posible
/// permutations and in order make it easy to add new ways of calculating the distance
/// we have implemented this using polymorphism and multiple inheritance.  The following
/// provides \ref AddingAMetric "information" on how to implement a new method for
/// calculating the distance between a pair of configurations

class ReferenceConfiguration {
  friend class SingleDomainRMSD;
  friend double distance( const Pbc& pbc, const std::vector<Value*> & vals, ReferenceConfiguration*, ReferenceConfiguration*, const bool& squared );
private:
/// The name of this particular config
  std::string name;
/// A vector containing all the remarks from the pdb input
  std::vector<std::string> line;
/// These are used to do fake things when we copy frames
  std::vector<AtomNumber> fake_atom_numbers;
  std::vector<std::string> fake_arg_names;
/// These are use by the distance function above
  std::vector<Vector> fake_refatoms;
  std::vector<double> fake_refargs;
  std::vector<double> fake_metric;
protected:
/// Crash with an error
  void error(const std::string& msg);
public:
  explicit ReferenceConfiguration( const ReferenceConfigurationOptions& ro );
/// Destructor
  virtual ~ReferenceConfiguration();
/// Return the name of this metric
  std::string getName() const ;
///
  virtual unsigned getNumberOfReferencePositions() const ;
  virtual unsigned getNumberOfReferenceArguments() const ;
/// Retrieve the atoms that are required for this guy
  virtual void getAtomRequests( std::vector<AtomNumber>&, bool disable_checks=false ) {}
/// Retrieve the arguments that are required for this guy
  virtual void getArgumentRequests( std::vector<std::string>&, bool disable_checks=false ) {}
/// Do all local business for setting the configuration
  virtual void read( const PDB& )=0;
/// Calculate the distance from the reference configuration
  double calculate( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, ReferenceValuePack& myder, const bool& squared=false ) const ;
/// Calculate the distance from the reference configuration
  virtual double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& args,
                       ReferenceValuePack& myder, const bool& squared ) const=0;
/// Parse something from the pdb remarks
/// Copy derivatives from one frame to this frame
  void copyDerivatives( const ReferenceConfiguration* );
/// Get one of the referene arguments
  virtual double getReferenceArgument( const unsigned& i ) const { plumed_error(); return 0.0; }
/// These are overwritten in ReferenceArguments and ReferenceAtoms but are required here
/// to make PLMD::distance work
  virtual const std::vector<Vector>& getReferencePositions() const ;
  virtual const std::vector<double>& getReferenceArguments() const ;
  virtual const std::vector<double>& getReferenceMetric();
/// These are overwritten in ReferenceArguments and ReferenceAtoms to make frame copying work
  virtual const std::vector<AtomNumber>& getAbsoluteIndexes();
  virtual const std::vector<std::string>& getArgumentNames();
/// Extract a Direction giving you the displacement from some position
  void extractDisplacementVector( const std::vector<Vector>& pos, const std::vector<Value*>& vals,
                                  const std::vector<double>& arg, const bool& nflag,
                                  Direction& mydir ) const ;
/// Stuff for pca
  virtual bool pcaIsEnabledForThisReference() { return false; }
  double projectDisplacementOnVector( const Direction& mydir, const std::vector<Value*>& vals,
                                      const std::vector<double>& arg, ReferenceValuePack& mypack ) const ;
/// Stuff to setup pca
  virtual void setupPCAStorage( ReferenceValuePack& mypack ) { plumed_error(); }
/// Move the reference configuration by an ammount specified using a Direction
  void displaceReferenceConfiguration( const double& weight, Direction& dir );
};

inline
const std::vector<Vector>& ReferenceConfiguration::getReferencePositions() const {
  return fake_refatoms;
}

inline
const std::vector<double>& ReferenceConfiguration::getReferenceArguments() const {
  return fake_refargs;
}

inline
const std::vector<double>& ReferenceConfiguration::getReferenceMetric() {
  return fake_metric;
}

inline
const std::vector<AtomNumber>& ReferenceConfiguration::getAbsoluteIndexes() {
  return fake_atom_numbers;
}

inline
const std::vector<std::string>& ReferenceConfiguration::getArgumentNames() {
  return fake_arg_names;
}

inline
unsigned ReferenceConfiguration::getNumberOfReferencePositions() const {
  return 0;
}

inline
unsigned ReferenceConfiguration::getNumberOfReferenceArguments() const {
  return 0;
}

}
#endif

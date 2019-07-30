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
#ifndef __PLUMED_reference_ReferenceArguments_h
#define __PLUMED_reference_ReferenceArguments_h

#include "ReferenceConfiguration.h"
#include "tools/Matrix.h"

namespace PLMD {

/// \ingroup TOOLBOX
/// In many applications (e.g. paths, fields, property maps) it is necessary to calculate
/// the distance between two configurations.  These distances can be calculated in a variety of
/// different ways.  For instance, one can assert that the distance between the two configuration
/// is the distance one would have to move all the atoms to transform configuration 1 into configuration
/// 2. Alternatively, one could calculate the values of a large set of collective coordinates in the two
/// configurations and then calculate the Euclidean distances between these two points in the resulting
/// high-dimensional vector space.  Lastly, one can combine these two forms of distance calculation to calculate
/// a hybrid distance.  Plumed allows one to use all these forms of distance calculations and also to implement
/// new forms of distance.  You should inherit from this class if your distance involves reference colvar values.
/// This class and \ref PLMD::ReferenceAtoms mirror the functionalities in \ref PLMD::ActionWithArguments and
/// \ref PLMD::ActionAtomistic respectively but for distances.

class ReferenceArguments :
  virtual public ReferenceConfiguration
{
  friend class Direction;
  friend class ReferenceConfiguration;
private:
/// The weights for normed euclidean distance
  std::vector<double> weights, sqrtweight;
/// The N X N matrix we are using to calculate our Malanobius distance
  Matrix<double> metric;
  std::vector<double> trig_metric;
/// The values of the colvars in the reference configuration
  std::vector<double> reference_args;
/// The names of the arguments
  std::vector<std::string> arg_names;
/// The indices for setting derivatives
  std::vector<unsigned> arg_der_index;
protected:
/// Are we reading weights from input
  bool hasweights;
/// Are we calculating a Malanobius distance
  bool hasmetric;
/// Read in the atoms from the pdb file
  void readArgumentsFromPDB( const PDB& pdb );
/// Set the values of the colvars based on their current instantanous values (used in Analysis)
  void setReferenceArguments();
public:
  explicit ReferenceArguments( const ReferenceConfigurationOptions& ro );
/// Get the number of reference arguments
  unsigned getNumberOfReferenceArguments() const override;
/// Get the arguments required
  void getArgumentRequests( std::vector<std::string>&, bool disable_checks=false ) override;
/// Set the positions of the refernce arguments
  void setReferenceArguments( const std::vector<double>& arg_vals, const std::vector<double>& sigma );
/// Set the positions of the reference arguments
  void moveReferenceArguments( const std::vector<double>& arg_vals );
/// Get the value of the ith reference argument
  double getReferenceArgument( const unsigned& i ) const override;
/// Return all the reference arguments
  const std::vector<double>& getReferenceArguments() const override;
  const std::vector<double>& getReferenceMetric() override;
/// Return names
  const std::vector<std::string>& getArgumentNames() override;
/// Calculate the euclidean/malanobius distance the atoms have moved from the reference
/// configuration in CV space
  virtual double calculateArgumentDistance( const std::vector<Value*> & vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const ;
/// Displace the positions of the reference atoms
  void displaceReferenceArguments( const double& weight, const std::vector<double>& displace );
/// Extract the displacement from a position in a space
  virtual void extractArgumentDisplacement( const std::vector<Value*>& vals, const std::vector<double>& arg, std::vector<double>& dirout ) const ;
/// Project the displacement of the arguments on a vector
  double projectArgDisplacementOnVector( const std::vector<double>& eigv, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& mypack ) const ;
};

inline
double ReferenceArguments::getReferenceArgument( const unsigned& i ) const {
  plumed_dbg_assert( i<reference_args.size() );
  return reference_args[i];
}

inline
const std::vector<double>& ReferenceArguments::getReferenceArguments() const {
  return reference_args;
}

inline
const std::vector<std::string>& ReferenceArguments::getArgumentNames() {
  return arg_names;
}

inline
unsigned ReferenceArguments::getNumberOfReferenceArguments() const {
  return reference_args.size();
}

}
#endif


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
#ifndef __PLUMED_reference_ReferenceConfiguration_h
#define __PLUMED_reference_ReferenceConfiguration_h

#include <vector>
#include <string>
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Tools.h"
#include "tools/Exception.h"

namespace PLMD{

class Value;
class Pbc;
class OFile;
class PDB;

/// \ingroup TOOLBOX
/// Abstract base class for calculating the distance from a reference configuration.
/// A reference configuration can either have a particular set of atoms in a particular
/// given configuration or it can be that a particular set of colvars have a particular 
/// set of values.  It could also be a combination of both.  To allow all the posible
/// permutations and in order make it easy to add new ways of calculating the distance 
/// we have implemented this using polymorphism and multiple inheritance. 

class ReferenceConfigurationOptions {
friend class ReferenceConfiguration;
private:
  std::string tt;
public:
  ReferenceConfigurationOptions( const std::string& type );
  bool usingFastOption() const ;
  std::string getMultiRMSDType() const ;
};

class ReferenceConfiguration {
friend class SingleDomainRMSD;
private:
/// The name of this particular config
  std::string name;
/// A weight assigned to this particular frame
  double weight;
/// A vector containing all the remarks from the pdb input
  std::vector<std::string> line;
protected:
/// Derivatives wrt to the arguments
  std::vector<double> arg_ders;
/// The virial contribution has to be stored 
  bool virialWasSet;
  Tensor virial;
/// Derivatives wrt to the atoms
  std::vector<Vector> atom_ders;
/// Return the name of this metric
  std::string getName() const ;
/// Crash with an error
  void error(const std::string& msg);
/// Clear the derivatives 
  void clearDerivatives();
public:
  ReferenceConfiguration( const ReferenceConfigurationOptions& ro );
/// Destructor
  virtual ~ReferenceConfiguration();
/// Retrieve the atoms that are required for this guy
  virtual void getAtomRequests( std::vector<AtomNumber>&, bool disable_checks=false ){}
/// Retrieve the arguments that are required for this guy
  virtual void getArgumentRequests( std::vector<std::string>&, bool disable_checks=false ){}
/// Set the final number of arguments
  virtual void setNumberOfArguments( const unsigned& );
/// Set the final number of atoms
  virtual void setNumberOfAtoms( const unsigned& );
/// Set the reference configuration using a PDB 
  virtual void set( const PDB& );
/// Do all local business for setting the configuration 
  virtual void read( const PDB& )=0;
/// Return the weight for this frame
  double getWeight() const ;
/// Calculate the distance from the reference configuration
  double calculate( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const bool& squared=false );
/// Calculate the distance from the reference configuration
  virtual double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const bool& squared )=0;
/// Return the derivative wrt to the ith atom
  Vector getAtomDerivative( const unsigned& ) const ;
/// Return the derivative wrt to the ith argument
  double getArgumentDerivative( const unsigned& ) const ;
/// Return the derivatives of the distance wrt the cell vectors.  This returns false
/// for everything other than DRMSD as these sort of calculations have to be done 
/// separately when you use RMSD 
  bool getVirial( Tensor& virout ) const ;    
/// Parse something from the pdb remarks
  template<class T>
  bool parse( const std::string&key, T&t, bool ignore_missing=false );
/// Parse vector from the pdb remarks
  template<class T>
  bool parseVector( const std::string&key, std::vector<T>&t, bool ignore_missing=false );
/// Parse a flag
  void parseFlag(const std::string&key,bool&t);
/// Check that all the remarks in the pdb have been read in
  void checkRead();
/// Copy derivatives from one frame to this frame
  void copyDerivatives( const ReferenceConfiguration* );
/// Set the atom numbers and the argument names
  void setNamesAndAtomNumbers( const std::vector<AtomNumber>& numbers, const std::vector<Value*>& arg );
/// Set the reference structure (perhaps should also pass the pbc and align and displace )
  void setReference( const std::vector<Vector>& pos, const std::vector<Value*>& arg, const std::vector<double>& metric );
/// Print a pdb file containing the reference configuration
  void print( OFile& ofile, const double& time, const double& weight, const double& old_norm );
/// Get one of the referene arguments
  virtual double getReferenceArgument( const unsigned& i ){ plumed_error(); return 0.0; }
};

inline
Vector ReferenceConfiguration::getAtomDerivative( const unsigned& ider ) const {
  plumed_dbg_assert( ider<atom_ders.size() );
  return atom_ders[ider];
}

inline
double ReferenceConfiguration::getArgumentDerivative( const unsigned& ider ) const {
  plumed_dbg_assert( ider<arg_ders.size() );
  return arg_ders[ider];
}

inline
double ReferenceConfiguration::getWeight() const {
  return weight;
}

template<class T>
bool ReferenceConfiguration::parse(const std::string&key, T&t, bool ignore_missing ){
  bool found=Tools::parse(line,key,t);
  if(!ignore_missing && !found) error(key + " is missing"); 
  return found;
}

template<class T>
bool ReferenceConfiguration::parseVector(const std::string&key,std::vector<T>&t, bool ignore_missing){
  bool found=Tools::parseVector(line,key,t);
  if(!ignore_missing && !found)  error(key + " is missing");
  return found;
}

}
#endif

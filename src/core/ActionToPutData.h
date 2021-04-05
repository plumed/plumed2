/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#ifndef __PLUMED_core_ActionToPutData_h
#define __PLUMED_core_ActionToPutData_h

#include "ActionWithValue.h"
#include "DataPassingObject.h"

namespace PLMD {

class ActionToPutData : 
public ActionWithValue
{
friend class Atoms;
private:
/// Action has been set
  bool wasset;
/// Do the domains need to be summed
  bool sum_domains;
/// Are we not applying forces on this values
  bool noforce;
/// Is this quantity scattered over the domains
  bool scattered;
/// Is this quantity fixed
  bool fixed;
/// Is the the first step
  bool firststep;
/// Can we set data at the current time
  bool dataCanBeSet;
/// Have the original forces been scaled
  bool wasscaled;
/// This is the list of forces that must be scaled
  std::vector<ActionToPutData*> forces_to_scale;
/// This holds the pointer that we getting data from
  std::unique_ptr<DataPassingObject> mydata;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionToPutData(const ActionOptions&ao);
/// This resets the stride in the collection object
  void setStride( const unsigned& sss );
/// Check if the value has been set
  bool hasBeenSet() const ;
/// Override clear the input data 
  void clearDerivatives( const bool& force ){}
/// Do not add chains to setup actions
  bool canChainFromThisAction() const { return false; }
/// Have to override this to deal with wildcards
  void interpretDataLabel( const std::string& mystr, Action* myuser, unsigned& nargs, std::vector<Value*>& args );
/// The number of derivatives
  unsigned getNumberOfDerivatives() const { return 0; }
/// Do we need to sum this over all the domains
  bool sumOverDomains() const ;
/// Is this quantity scattered over the domains
  bool collectFromDomains() const;
/// Do we always need to collect the atoms from all domains
  bool collectAllFromDomains() const;
/// Set the unit for this quantity
  void setUnit( const double& u );
/// Set the unit of the force on this quantity
  void setForceUnit( const double& u );
/// Set the memory that holds the value
  void set_value(void* val );
/// Set the memory that holds the force
  void set_force(void* val );
/// Share the data from the holder when the data is distributed over domains
  void share( const unsigned& j, const unsigned& k );
  void share( const std::set<AtomNumber>&index, const std::vector<unsigned>& i );
/// Get the data to share
  void wait();
/// Actually set the values for the output
  void calculate(){ firststep=false; wasscaled=false; }
  void apply();
  void rescaleForces( const double& alpha );
/// For replica exchange
  void writeBinary(std::ostream&o);
  void readBinary(std::istream&i);
};

inline
bool ActionToPutData::sumOverDomains() const {
  return sum_domains;
}

inline
bool ActionToPutData::hasBeenSet() const {
  return wasset;
}

inline
bool ActionToPutData::collectFromDomains() const {
  if( scattered && fixed ) return firststep;
  return scattered;
}

inline
bool ActionToPutData::collectAllFromDomains() const {
  if( scattered && fixed ) return firststep;
  return false;
}

}
#endif

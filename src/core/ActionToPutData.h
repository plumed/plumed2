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
private:
/// Action has been set
  bool wasset;
/// Do the domains need to be summed
  bool sum_domains;
/// This holds the pointer that we getting data from
  std::unique_ptr<DataPassingObject> mydata;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionToPutData(const ActionOptions&ao);
/// Set the periodicity of the variable
  void set_domain( const bool& periodic, const std::string& min, const std::string& max );
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
/// Set the unit of the energy
  void setUnit( const double& u );
/// Set the memory that holds the value
  void set_value(void* val );
/// Set the memory that holds the force
  void set_force(void* val );
/// Get the data to share
  void wait();
/// Actually set the values for the output
  void calculate(){}
  void apply();
/// For replica exchange
  void writeBinary(std::ostream&o);
  void readBinary(std::istream&i);
};

inline
bool ActionToPutData::sumOverDomains() const {
  return sum_domains;
}

}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_ValueFromMDCode_h
#define __PLUMED_core_ValueFromMDCode_h

#include <memory>
#include "Value.h"
#include "tools/Communicator.h"

namespace PLMD {

class ValueFromMDCode { 
friend class Atoms;
private:
/// Creates an ValueFromMDCodeTyped<T> object such that sizeof(T)==n
  static std::unique_ptr<ValueFromMDCode> create(unsigned n, const std::string& name, const std::vector<unsigned>& shape );
/// Returns the data array from the value so that Atoms can do an allgather on the data
  std::vector<double>& getDataToGather();
protected:
/// Is the value not time dependent
  bool fixed;
/// This tell us if we need to collect this value
  bool collect;  
/// This is set true once the variable has been obtained
  bool set;
/// Is the value collected from the domain decomposition
  bool domain_decomposed;
/// The value that we are creating based on the data from the MD code
  std::unique_ptr<Value> value;
public:
  explicit ValueFromMDCode(const std::string& name, const std::vector<unsigned>& shape);
/// Setup the periodicity of the value
  void setupPeriodicity( const bool& isperiodic, const std::string& min, const std::string& max );
/// Set up a pointer to the value array in the MD code
  virtual void setv(void*m)=0;
/// Set up the pointer to the output forces for the value
  virtual void setf(void*f)=0;
/// Update the forces on this value
  virtual void updateForces()=0;
/// Gather the value that we have collected
  virtual void gather()=0;
/// Value does not change with time
  bool isFixed() const ;
/// Value must be collected from domain decomposition
  bool collectFromDomains() const ;
/// Get the name of the input value
  const std::string& getName() const ;
/// Return a pointer to the value
  Value* getPntrToValue();
};

inline
bool ValueFromMDCode::isFixed() const {
  return fixed;
}

inline
bool ValueFromMDCode::collectFromDomains() const {
  return domain_decomposed;
}

inline
const std::string&  ValueFromMDCode::getName() const {
  return value->getName();
}

inline
std::vector<double>& ValueFromMDCode::getDataToGather() {
  return value->data;
}

inline
Value* ValueFromMDCode::getPntrToValue() {
  return value.get();
}

}

#endif

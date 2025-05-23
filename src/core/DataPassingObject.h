/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#ifndef __PLUMED_core_DataPassingObject_h
#define __PLUMED_core_DataPassingObject_h

#include <vector>
#include <memory>
#include <set>
#include "tools/AtomNumber.h"
#include "tools/TypesafePtr.h"
#include "Value.h"

namespace PLMD {

class DataPassingObject {
protected:
/// The start of the data in the input pointer
  unsigned start;
/// The spacing between values in the input arrays
  unsigned stride;
/// The units of the quantity
  double unit;
/// The units of the force on this quantity
  double funit;
/// The backup value of the quantity (used if the value is passed directly)
  bool hasbackup;
  double bvalue;
public:
  static std::unique_ptr<DataPassingObject> create(unsigned n);
  explicit DataPassingObject() : start(0), stride(1), unit(1), funit(1), hasbackup(false), bvalue(0) {}
  /// Virtual destructor, just to allow inheritance.
  virtual ~DataPassingObject() {}
/// Convert what comes from the MD code to a double
  virtual double MD2double(const TypesafePtr & m) const=0;
///
  void setStart( const unsigned& s ) {
    start=s;
  }
/// Set the stride to use when getting data from the input array
  void setStride( const unsigned& s ) {
    stride=s;
  }
/// Set the unit for the value
  void setUnit( const double& u ) {
    unit=u;
  }
/// Set the unit for the force
  void setForceUnit( const double& u ) {
    funit=u;
  }
/// This is used when you want to save the passed object to a double variable in PLUMED rather than the pointer
/// this can be used even when you don't pass a pointer from the MD code
  virtual void saveValueAsDouble( const TypesafePtr & val )=0;
/// Set the pointer to the value
  virtual void setValuePointer( const TypesafePtr & val, const std::vector<std::size_t>& shape, const bool& isconst )=0;
/// Set the pointer to the force
  virtual void setForcePointer( const TypesafePtr & val, const std::vector<std::size_t>& shape )=0;
/// This gets the data in the pointer and passes it to the output value
  virtual void share_data( std::vector<double>& values ) const = 0;
/// Share the data and put it in the value from sequential data
  virtual void share_data( const unsigned& j, const unsigned& k, Value* value )=0;
/// Share the data and put it in the value from a scattered data
  virtual void share_data( const std::vector<AtomNumber>&index, const std::vector<unsigned>& i, Value* value )=0;
/// Pass the force from the value to the output value
  virtual void add_force( Value* vv )=0;
  virtual void add_force( const std::vector<int>& index, Value* value )=0;
  virtual void add_force( const std::vector<AtomNumber>& index, const std::vector<unsigned>& i, Value* value )=0;
/// Rescale the forces that were passed
  virtual void rescale_force( const unsigned& n, const double& factor, Value* value )=0;
/// This transfers everything to the output
  virtual void setData( Value* value )=0;
};

}
#endif

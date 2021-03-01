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
#ifndef __PLUMED_core_ActionToFetchData_h
#define __PLUMED_core_ActionToFetchData_h

#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "ActionPilot.h"

namespace PLMD {

class OutputDataObject {
public:
  static std::unique_ptr<OutputDataObject> create(unsigned n);
/// Set the pointer to the output
  virtual void setPointer( void* outval )=0;
/// This transfers everything to the output
  virtual void setData( const std::vector<double>& data )=0;  
};

class ActionToFetchData : 
public ActionPilot,
public ActionWithArguments 
{
private:
/// What do you want to collect to pass back to python
  enum {val,deriv,force} gtype;
/// This holds the pointer that we are setting 
  std::unique_ptr<OutputDataObject> mydata;
/// This temporarily holds the data so it can be passed out
  std::vector<double> data;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionToFetchData(const ActionOptions&ao);
/// Get the rank of the output
  void get_rank( long* rank );
/// Get the shape of the output
  void get_shape( long* dims );
/// Set the memory that holds the output
  void set_memory(void* val );
/// Actually set the values for the output
  void calculate();
  void apply(){}
};

}
#endif

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
#ifndef __PLUMED_core_DataFetchingObject_h
#define __PLUMED_core_DataFetchingObject_h

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>

namespace PLMD {

class ActionSet;
class PlumedMain;
class ActionWithValue;
class Value;

class DataFetchingObject {
protected:
/// Pointers to the various actions required by the grabber
  std::vector<ActionWithValue*> myactions;
/// The values required by the user
  std::vector<Value*> myvalues;
/// A copy of the plumed main object
  PlumedMain & plumed;
public:
  static std::unique_ptr<DataFetchingObject> create(unsigned n, PlumedMain& p);
/// A constructor so that we can create the plumed main object
  explicit DataFetchingObject(PlumedMain&p);
  virtual ~DataFetchingObject() {}
///
  bool activate() const ;
/// Return the rank required for a particular key
  static void get_rank( const ActionSet& a, const std::string& key, const std::string& type, long* rank );
/// Return the shape required for a particular key
  static void get_shape( const ActionSet& a, const std::string& key, const std::string& type, long* dims );
/// Find the action that calculates a particular value
  static ActionWithValue* findAction( const ActionSet& a, const std::string& key );
/// Set the pointer to the data
  virtual void setData( const std::string& key, const std::string& type, void* outval )=0;
/// After calc has been performed grab all the data and put it in the relevant arrays
  virtual void finishDataGrab()=0;
};

}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "DataFetchingObject.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "Action.h"
#include "ActionWithValue.h"
#include "Value.h"

namespace PLMD {

template <class T>
class DataFetchingObjectTyped : public DataFetchingObject {
private:
/// A map containing the data we are grabbing
  std::map<std::string,T*> data;
public:
  explicit DataFetchingObjectTyped(PlumedMain&plumed);
  ~DataFetchingObjectTyped() {}
  void setData( const std::string& key, const std::string& type, void* outval ) override;
  void finishDataGrab() override;
};

std::unique_ptr<DataFetchingObject> DataFetchingObject::create(unsigned n, PlumedMain& p) {
  if(n==sizeof(double)) {
    return std::unique_ptr<DataFetchingObjectTyped<double>>(new DataFetchingObjectTyped<double>(p));
  } else  if(n==sizeof(float)) {
    return std::unique_ptr<DataFetchingObjectTyped<float>>(new DataFetchingObjectTyped<float>(p));
  }
  std::string pp; Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

DataFetchingObject::DataFetchingObject(PlumedMain&p):
  plumed(p)
{
}

bool DataFetchingObject::activate() const {
  for(unsigned j=0; j<myactions.size(); ++j) myactions[j]->activate();
  if( myactions.size()>0 ) return true;
  return false;
}

ActionWithValue* DataFetchingObject::findAction( const ActionSet& a, const std::string& key ) {
  std::string aname = key; std::size_t dot = key.find(".");
  if( dot!=std::string::npos ) aname = key.substr(0,dot);
  return a.selectWithLabel<ActionWithValue*>( aname );
}

void DataFetchingObject::get_rank( const ActionSet& a, const std::string& key, const std::string& type, long* dims ) {
  plumed_assert( Tools::getWords(key,"\t\n ,").size()==1 );
  plumed_massert( key.find("*")==std::string::npos, "cannot use wildcards in python interface");

  // Find the appropriate action and store value containing quantity of interest
  ActionWithValue* myv = findAction( a, key );
  Value* val = myv->copyOutput( key );

  // Now work out what we are returning for this action
  if( type=="" ) {
    // Return a single value in this case
    dims[0]=1;
  } else if( type=="derivatives" ) {
    plumed_merror("not yet implemented");
  } else if( type=="forces" ) {
    plumed_merror("not yet implemented");
  } else {
    plumed_merror("invalid type specifier");
  }
}

void DataFetchingObject::get_shape( const ActionSet& a, const std::string& key, const std::string& type, long* dims ) {
  plumed_assert( Tools::getWords(key,"\t\n ,").size()==1 );
  plumed_massert( key.find("*")==std::string::npos, "cannot use wildcards in python interface");

  // Find the appropriate action and store value containing quantity of interest
  ActionWithValue* myv = findAction( a, key );
  Value* val = myv->copyOutput( key );

  // Now work out what we are returning for this action
  if( type=="" ) {
    // Return a single value in this case
    dims[0]=1;
  } else if( type=="derivatives" ) {
    plumed_merror("not yet implemented");
  } else if( type=="forces" ) {
    plumed_merror("not yet implemented");
  } else {
    plumed_merror("invalid type specifier");
  }
}

template <class T>
DataFetchingObjectTyped<T>::DataFetchingObjectTyped(PlumedMain&p):
  DataFetchingObject(p)
{
}

template <class T>
void DataFetchingObjectTyped<T>::setData( const std::string& key, const std::string& type, void* outval ) {
  plumed_assert( Tools::getWords(key,"\t\n ,").size()==1 );
  plumed_massert( key.find("*")==std::string::npos, "cannot use wildcards in python interface");
  plumed_massert( !data.count(key + " " + type), "already collecting this data elsewhere");
  // Add the space to store the data to the data map
  T* f=static_cast<T*>(outval);
  data.insert(std::pair<std::string,T*>(key + " " + type,f));

  // Find the appropriate action and store value containing quantity of interest
  ActionWithValue* myv = DataFetchingObject::findAction( plumed.getActionSet(), key );
  // Store the action if not already stored
  bool found=false;
  for(const auto & p : myactions) {
    if( p->getLabel()==myv->getLabel() ) { found=true; break; }
  }
  if( !found ) myactions.push_back( myv );
  // Store the value
  myvalues.push_back( myv->copyOutput( key ) );
}

template <class T>
void DataFetchingObjectTyped<T>::finishDataGrab() {
  // Run over all values and collect data
  for(const auto & p : myvalues ) {
    T* val = static_cast<T*>( data.find(p->getName() + " ")->second );
    if( data.find(p->getName() + " ")!=data.end() ) {
      val[0] = static_cast<T>( p->get() );
    }
    if( data.find(p->getName() + " derivatives")!=data.end() ) {
      plumed_merror("not implemented yet");
    }
    if( data.find(p->getName() + " forces")!=data.end() ) {
      plumed_merror("not implemented yet");
    }
  }
}

}


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_reference_MetricRegister_h
#define __PLUMED_reference_MetricRegister_h

#include <string>
#include <cstring>
#include <vector>
#include <map>
#include "tools/Exception.h"
#include "tools/Tools.h"
#include "tools/PDB.h"
#include "ReferenceConfiguration.h"

namespace PLMD {

class PDB;

class MetricRegister {
private:
/// Pointer to a function which, given the type for a ReferenceConfiguration, creates it
  typedef std::unique_ptr<ReferenceConfiguration> (*creator_pointer)(const ReferenceConfigurationOptions&);
/// The set of possible distribution functions we can work with
  std::map<std::string,creator_pointer> m;
public:
/// The destructor
  ~MetricRegister();
/// Add a new metric to the register of metrics
  void add( std::string type, creator_pointer );
/// Remove a metric from the register of metrics
  void remove(creator_pointer f);
/// Verify if a particular metric type is present in the register
  bool check(std::string type);
/// Create a reference configuration and don't set a point of reference
  template <class T>
  std::unique_ptr<T> create( const std::string& type );
/// Create a reference configuration and set the point of reference from the pdb
  template <class T>
  std::unique_ptr<T> create( const std::string& type, const PDB& pdb );
};

MetricRegister& metricRegister();

#define PLUMED_REGISTER_METRIC(classname,type) \
  static class classname##RegisterMe{ \
    static std::unique_ptr<ReferenceConfiguration> create(const PLMD::ReferenceConfigurationOptions&ro){return std::unique_ptr<ReferenceConfiguration>( new classname(ro) );} \
  public: \
    classname##RegisterMe(){PLMD::metricRegister().add(type,create);}; \
    ~classname##RegisterMe(){PLMD::metricRegister().remove(create);}; \
  } classname##RegisterMeObject;

template <class T>
std::unique_ptr<T> MetricRegister::create( const std::string& type ) {
  std::string ftype;
  if( type.find("MULTI-")!=std::string::npos ) {
    ftype="MULTI";
  } else {
    std::size_t dash=type.find("-FAST"); // We must remove the fast label
    ftype=type.substr(0,dash);
  }
  plumed_massert( check(ftype), "metric " + ftype + " does not exist" );
  ReferenceConfigurationOptions ro( type );
// put immediately the result in a safe pointer
  std::unique_ptr<ReferenceConfiguration> conf( m[ftype]( ro ) );
// try conversion
  T*ptr=dynamic_cast<T*>( conf.get() );
// if this throws, the unique_ptr conf is deleted.
// Notice that with the original version of the code (2.4) an error here
// would have lead to a memory leak.
  if(!ptr ) plumed_merror( type + " metric is not valid in this context");
// release ownership in order to transfer it to returned pointer
  conf.release();
// notice that I should pass ptr here rather than conf.release(),
// since the type is different
  return std::unique_ptr<T>(ptr);
}

template <class T>
std::unique_ptr<T> MetricRegister::create( const std::string& type, const PDB& pdb ) {
  std::string rtype;
  if( type.length()==0 ) {
    rtype=pdb.getMtype();
    plumed_massert(rtype.length()>0, "TYPE not specified in pdb input file");
  } else {
    rtype=type;
  }
  std::unique_ptr<T> confout( create<T>( rtype ) );
  confout->read( pdb );
  return confout;
}

}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_VesselRegister_h
#define __PLUMED_vesselbase_VesselRegister_h

#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include "tools/Exception.h"
#include "tools/Keywords.h"

namespace PLMD {
namespace vesselbase {

class Vessel;
class VesselOptions;

class VesselRegister {
private:
/// Pointer to a function which, given the keyword for a distribution function, creates it
  typedef std::unique_ptr<Vessel>(*creator_pointer)(const VesselOptions&);
/// Pointer to the function that reserves the keyword for the distribution
  typedef void(*keyword_pointer)(Keywords&);
/// The set of possible distribution functions we can work with
  std::map<std::string,creator_pointer> m;
/// Map action to a function which documents the related object
  std::map<std::string,keyword_pointer> mk;
/// A vector of function pointers - this is used to create the documentation
  Keywords keywords;
public:
/// The destructor
  ~VesselRegister();
/// Add a new distribution function option to the register of distribution functions
  void add(std::string keyword,creator_pointer,keyword_pointer k,keyword_pointer ik);
/// Remove a distribution function from the register of distribution functions
  void remove(creator_pointer f);
/// Verify if a distribution keyword is present in the register
  bool check(std::string keyname);
/// Create a distribution function of the specified type
  std::unique_ptr<Vessel> create(std::string keyword, const VesselOptions&da);
/// Return the keywords
  Keywords getKeywords();
};

VesselRegister& vesselRegister();

#define PLUMED_REGISTER_VESSEL(classname,keyword) \
  static class classname##RegisterMe{ \
    static std::unique_ptr<PLMD::vesselbase::Vessel> create(const PLMD::vesselbase::VesselOptions&da){return std::unique_ptr<classname>( new classname(da) );} \
  public: \
    classname##RegisterMe(){PLMD::vesselbase::vesselRegister().add(keyword,create,classname::reserveKeyword,classname::registerKeywords);} \
    ~classname##RegisterMe(){PLMD::vesselbase::vesselRegister().remove(create);} \
  } classname##RegisterMeObject;

}
}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#ifndef __PLUMED_analysis_LandmarkRegister_h
#define __PLUMED_analysis_LandmarkRegister_h

#include <string>
#include <cstring>
#include <vector>
#include <map>
#include "LandmarkSelectionBase.h"

namespace PLMD{

class PDB;

namespace analysis{

class LandmarkRegister{
private:
/// Pointer to a function which, given the type for a ReferenceConfiguration, creates it
  typedef LandmarkSelectionBase*(*creator_pointer)(const LandmarkSelectionOptions&);
/// The set of possible landmark selection algorithms we can work with
  std::map<std::string,creator_pointer> m;
public:
/// The destructor
  ~LandmarkRegister();
/// Add a new landmark selection style to the register of landmark selectors
  void add( std::string type, creator_pointer );
/// Remove a landmark selection style from the register of metrics
  void remove(creator_pointer f);
/// Verify if a landmark selection style is present in the register
  bool check(std::string type);
/// Create a landmark selection object
  LandmarkSelectionBase* create( const LandmarkSelectionOptions& lo );
};

LandmarkRegister& landmarkRegister();

#define PLUMED_REGISTER_LANDMARKS(classname,type) \
  static class classname##RegisterMe{ \
    static LandmarkSelectionBase * create(const LandmarkSelectionOptions&lo){return new classname(lo);} \
  public: \
    classname##RegisterMe(){landmarkRegister().add(type,create);}; \
    ~classname##RegisterMe(){landmarkRegister().remove(create);}; \
  } classname##RegisterMeObject;

}
}
#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#ifndef __PLUMED_multicolvar_MultiColvar_h
#define __PLUMED_multicolvar_MultiColvar_h

#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/SwitchingFunction.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for creating distributions of colvars and functions
thereof, whtin it there is \ref AddingAMultiColvar "information" as to how to go implementing these types of actions.
*/

class MultiColvar : public MultiColvarBase {
private:
/// Do we want lots of details in the output
  bool verbose_output;
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms, std::vector<AtomNumber>& all_atoms );
protected:
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( const std::string& str1, const std::string& str2, int& natoms, std::vector<AtomNumber>& all_atoms );
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string & key, int& natoms, std::vector<AtomNumber>& all_atoms );
/// Read two group keywords
  void readTwoGroups( const std::string& key1, const std::string& key2, std::vector<AtomNumber>& all_atoms );
/// Read three groups
  void readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3, 
                        const bool& allow2, const bool& no_third_dim_accum, std::vector<AtomNumber>& all_atoms );
/// Add a collective variable
  void addColvar( const std::vector<unsigned>& newatoms );
public:
  explicit MultiColvar(const ActionOptions&);
  ~MultiColvar(){}
  static void registerKeywords( Keywords& keys );
/// Get the position of atom iatom
  const Vector & getPosition(unsigned) const;
};

}
}

#endif

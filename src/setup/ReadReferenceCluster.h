/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#ifndef __PLUMED_setup_ReadReferenceCluster_h
#define __PLUMED_setup_ReadReferenceCluster_h

#include "SetupReferenceBase.h"

namespace PLMD {
namespace setup {

class ReadReferenceCluster: public SetupReferenceBase {
private:
  std::vector<std::string> read_args;
public:
  static void registerKeywords( Keywords& keys );
  static std::string convertFileToLine( const std::string& reference, const unsigned& number, const std::vector<std::string>& names );
  explicit ReadReferenceCluster(const ActionOptions&ao);
  std::string getArgName( const unsigned& k ) const ;
};

}
}
#endif

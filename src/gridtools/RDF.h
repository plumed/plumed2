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
#ifndef __PLUMED_gridtools_RDF_h
#define __PLUMED_gridtools_RDF_h

#include "core/ActionShortcut.h"

namespace PLMD {
namespace gridtools {

class RDF : public ActionShortcut {
public:
  static void createX2ReferenceObject( const std::string& lab,  const std::string& grid_setup, const bool& calc_dens, const bool& no_average, ActionShortcut* action );
  static void registerKeywords( Keywords& keys );
  static void getDistanceMatrixShape( const std::string& lab, ActionShortcut* action, std::vector<std::string>& shape_str );
  explicit RDF(const ActionOptions&ao);
};

}
}
#endif

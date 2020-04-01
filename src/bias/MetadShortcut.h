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
#ifndef __PLUMED_bias_MetadShortcut_h
#define __PLUMED_bias_MetadShortcut_h

#include "core/ActionShortcut.h"

namespace PLMD {
namespace bias {

class MetadShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MetadShortcut(const ActionOptions&);
  static void createMetadBias( const std::string& lab, const std::string& pacestr, const std::vector<std::string>& args,
                               const std::vector<std::string>& gmin, const std::vector<std::string>& gmax, const std::vector<std::string>& grid_nbins,
                               const std::string& weight_str, const std::string& truncflag1, const std::string& truncflag2, ActionShortcut* act );
};

}
}
#endif

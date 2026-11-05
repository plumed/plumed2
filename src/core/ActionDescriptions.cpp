/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionDescriptions.h"

#include <utility>

namespace PLMD {
const std::map<std::string,std::string> & getActionDescriptions() {
  static const std::map<std::string,std::string> descriptions=[] {
    // Insert rows in source order because duplicate registration names keep
    // the first summary, regardless of map initializer-list semantics.
    static const std::pair<const char*,const char*> rows[] = {
#include "ActionDescriptions.inc"
    };
    std::map<std::string,std::string> result;
    for(const auto& row : rows) {
      result.emplace(row.first,row.second);
    }
    return result;
  }();
  return descriptions;
}
} // namespace PLMD

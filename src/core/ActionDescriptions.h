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
#ifndef __PLUMED_core_ActionDescriptions_h
#define __PLUMED_core_ActionDescriptions_h
#include <map>
#include <string>

namespace PLMD {
/// Return the source-derived one-line descriptions embedded in this kernel.
/// The table follows the PLUMEDOC summary extraction used by gen_json.
/// Duplicate action names retain the first description in source-file order.
/// \return A kernel-owned map with static lifetime. An absent action has no
/// embedded description; dynamically loaded actions need not appear here.
const std::map<std::string,std::string> & getActionDescriptions();
} // namespace PLMD

#endif // __PLUMED_core_ActionDescriptions_h

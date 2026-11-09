/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_core_ModuleLoader_h
#define __PLUMED_core_ModuleLoader_h

#include <string>
#include <vector>

namespace PLMD {
class DLLoader;

/// Kernel-local discovery and activation of separately installed modules.
/// Native dependencies are loaded by the system loader. Each caller retains
/// its own handles, so explicit activation does not leak to other objects.
class ModuleLoader {
public:
  class Context {
    const DLLoader* previous;
  public:
    explicit Context(const DLLoader& loader);
    ~Context();
    Context(const Context&)=delete;
    Context& operator=(const Context&)=delete;
  };
  static void inheritContext(DLLoader& loader);
  static void autoload(DLLoader& loader);
  static void load(const std::string& name, DLLoader& loader);
  static bool available(const std::string& name, const DLLoader& loader);
  static std::vector<std::string> modules(const DLLoader& loader);
};
}
#endif

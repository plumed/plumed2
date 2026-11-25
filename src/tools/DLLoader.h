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
#ifndef __PLUMED_tools_DLLoader_h
#define __PLUMED_tools_DLLoader_h

#include <vector>
#include <string>

namespace PLMD {

/// \ingroup TOOLBOX
/// Class taking care of dynamic loading.
/// It contains wrappers to the dlopen() routine.
/// It is designed so that when an object of this class goes
/// out of scope all the libraries loaded by it are unloaded. In this
/// manner, loaded libraries are automatically unloaded at the end of
/// execution. Libraries are loaded with RTDL_LOCAL option, which
/// means that they are not accessible from outside. Still, if they
/// contain self-registering classes, they will register themselves
/// to the ActionRegister object.
class DLLoader {
  std::vector<void*> handles;
  /// Deleted copy constructor
  DLLoader(const DLLoader&) = delete;
  /// Deleted assignment
  DLLoader&operator=(const DLLoader&) = delete;
public:
  /// Default constructor
  DLLoader();
  /// Cleanup
  ~DLLoader();
  /// Load a library, returning its handle
  void* load(const std::string&);
  /// Returns true if the dynamic loader is available (on some systems it may not).
  static bool installed();
  /// RAII helper for promoting RTLD_LOCAL loaded objects to RTLD_GLOBAL
  class EnsureGlobalDLOpen {
    void* handle_=nullptr;
  public:
    /// makes sure that object defining ptr is globally available
    explicit EnsureGlobalDLOpen(const void* symbol) noexcept;
    /// dlclose the dlopened object
    ~EnsureGlobalDLOpen();
    ///Confevert a const reference to a
    template<typename T> EnsureGlobalDLOpen(const T&p) noexcept
      : EnsureGlobalDLOpen(reinterpret_cast<const void*>(p)) {}
  };

  /// Returns true if a PLUMED library is available in the global namespace.
  /// It does so by looking for the presence of the C interface.
  /// It will detect any kernel that is available in the global namespece,
  /// not just the one from which this call is issued. This is useful to
  /// detect possible conflicts in advance.
  static bool isPlumedGlobal();
  const std::vector<void*> & getHandles() const noexcept;
};

} // namespace PLMD

#endif

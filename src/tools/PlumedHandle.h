/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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
#ifndef __PLUMED_tools_PlumedHandle_h
#define __PLUMED_tools_PlumedHandle_h
#include "core/PlumedMainInitializer.h"
#include <memory>

namespace PLMD
{

class PlumedMain;

/**
Tiny local class to load a PLUMED kernel.

Maps command to either a loaded PLUMED kernel or to the present one.
It is a simplified version of the interface located at wrapper/Plumed.h.
Differences are:
- It does not use the `PLUMED_KERNEL` env var. Indeed, it would not make sense to use it,
  since this class is meant to load different kernels.
- It does not implement interoperability with C/FORTRAN interfaces.
- It does not implement global versions (e.g. PLMD::Plumed::gcmd).
- It does not implement PLMD::Plumed::ginitialized. If it cannot be loaded, it crashes in its constructor.
  This will make sure that once constructed the object is usable.

The mechanism for loading the kernel is anyway very similar to the one in wrapper/Plumed.c.
In particular, it can load both kernels from PLUMED <=2.4 and >=2.5, and it
tries to load the `libplumed.so` object if the `libplumedKernel.so` object does not load correctly.
It can also be created without passing any kernel path. In that case it refers to the current one
(the one to which this class belongs).

The following syntax creates a handle referring to the current kernel
\verbatim
PlumedHandle p;
// Alternatively:
// auto p=PlumedHandle();
p.cmd("init");
\endverbatim

The following syntax instead creates a handle referring to a loaded kernel
\verbatim
PlumedHandle p(PlumedHandle::dlopen("/path/to/libplumedkernel.so");
// Alternatively:
// auto p=PlumedHandle::dlopen("/path/to/libplumedkernel.so");
p.cmd("init");
\endverbatim

Notice that if there are problems loading the kernel an exception is thrown.
Thus, once constructed the object is guaranteed to be functional.

*/
class PlumedHandle {
/// Automatically dlclosing auxiliary class.
/// Just used to make sure that handle is dlclosed correctly.
  class DlHandle {
    void *handle=nullptr;
  public:
/// Default construct as nullptr
    DlHandle() {}
/// Construct from a void*
    DlHandle(void*h): handle(h) {}
/// Destructor will call dlclose if necessary
    ~DlHandle();
/// Covertible to void* so that it can be used directly
    operator void*() const {
      return handle;
    }
  };
/// Pointer to PlumedMain.
/// Used when using the current kernel in order to avoid unneeded indirections.
  std::unique_ptr<PlumedMain> local;
/// Pointer to dlsym handle used to open the kernel.
/// Null when using current kernel.
  DlHandle handle;
/// Pointer to symbol table.
/// Used for kernels>=2.5. We store it here since it is needed
/// in constructor to initialize create_/cmd_/finalize_.
/// Later on we might use the additional version information that it carries.
  plumed_symbol_table_type* const symbol_=nullptr;
/// Pointer to create function.
/// Used when kernel is dlopened.
  const plumed_create_pointer create_=nullptr;
/// Pointer to cmd function.
/// Used when kernel is dlopened.
  const plumed_cmd_pointer cmd_=nullptr;
/// Pointer to finalize function.
/// Used when kernel is dlopened.
  const plumed_finalize_pointer finalize_=nullptr;
/// Pointer to the plumed object.
/// Used when kernel is dlopened.
  void* const p=nullptr;
/// Constructor using the path to a kernel.
/// I keep it private to avoid confusion wrt the
/// similar constructor of PLMD::Plumed that accepts a string (conversion from FORTRAN).
  explicit PlumedHandle(const char* path);
public:
/// Default constructor.
/// Maps following commands to the current kernel.
  PlumedHandle();
/// Construct a PlumedHandle given the path to a kernel.
/// It just uses the private constructor PlumedHandle(const char* path).
  static PlumedHandle dlopen(const char* path);
/// Destructor.
/// In case a kernel was dlopened, it dlcloses it.
/// I make it virtual for future extensibility, though this is not necessary now.
  virtual ~PlumedHandle();
/// Move constructor.
  PlumedHandle(PlumedHandle &&) = default;
/// Execute cmd.
  void cmd(const char*key,const void*ptr=nullptr);
};

}
#endif

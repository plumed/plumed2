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
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

// The NVIDIA OpenACC runtime, statically linked into the plugin with
// -static-nvidia, holds a weak reference to acc_register_library and calls it
// at startup if it resolves (OpenACC profiling interface).
// When PLUMED is built with GCC and -fopenmp, libgomp is already loaded in the
// process before the plugin. libgomp exports an acc_register_library stub that
// aborts with "libgomp: TODO", and the dynamic linker binds the weak reference
// to it, so every run crashes at the first OpenACC call.
// A hidden definition here is bound when the plugin is linked, so libgomp's
// stub is never reached. As a consequence, OpenACC profiling tools that rely on
// acc_register_library being found are not registered by the plugin.
// The function has C linkage, so the enclosing namespace does not affect the
// symbol name; it is there to follow the PLUMED convention that every source
// file lives in namespace PLMD.
namespace PLMD {
extern "C" __attribute__((visibility("hidden")))
void acc_register_library(void*, void*, void*) {}
}

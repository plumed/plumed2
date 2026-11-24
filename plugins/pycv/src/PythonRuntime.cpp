/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 of Toni Giorgino

The pycv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The pycv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "plumed/tools/DLLoader.h"
#include <Python.h>

namespace PLMD {}

// Promote the embedding library without initializing another interpreter.
// A Python host already exports its interpreter symbols; dynamic lookup in
// the action library selects those ahead of the embedding library.
namespace {
auto runtime = PLMD::DLLoader::EnsureGlobalDLOpen(&Py_Initialize);
}

extern "C" void plumed_pycv_runtime_anchor() {}

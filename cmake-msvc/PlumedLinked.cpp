/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2026 The plumed team
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
/*
  Wrapper of the plumed executable linked to the shared kernel
  (CMake build with PLUMED_SHARED_KERNEL).

  The executable gets the kernel functions from plumed_symbol_table_reexport(),
  which the kernel exports (cmake-msvc/plumedKernel.def). Importing the symbol
  table itself (as src/wrapper/PlumedStatic.cpp does) would not work on Windows,
  since a data symbol can only be imported with __declspec(dllimport).
*/
#define __PLUMED_WRAPPER_REEXPORT_SYMBOL_TABLE 0
#define __PLUMED_WRAPPER_CXX 0
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#include "wrapper/Plumed.h"

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
#ifndef __PLUMED_gemmi_cif_h
#define __PLUMED_gemmi_cif_h

// The CIF parser of gemmi (grammar gemmi::cif::rules, built on PEGTL).
// This is the only file of PLUMED that includes it: see README.md in this directory.
// The upstream headers are compiled with PLUMED's warning flags, so they
// are treated as system headers (this covers the headers they include).
#if defined(__GNUC__)
#pragma GCC system_header
#endif

#ifdef __PLUMED_HAS_SYSTEM_GEMMI
#include <gemmi/cif.hpp>
// A gemmi installed in the system lives in the global namespaces gemmi and tao:
// PLUMED refers to them as PLMD::gemmi and PLMD::tao in both cases
namespace PLMD {
namespace gemmi = ::gemmi;
namespace tao = ::tao;
}
#else
// The shipped copy is in the namespaces PLMD::gemmi and PLMD::tao (see import.sh),
// so that it cannot clash with another copy of gemmi or PEGTL in the same program
#include "include/gemmi/cif.hpp"
#endif

#endif

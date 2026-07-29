/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023-2026 of Amr Alhossary and the KEnRef Authors.

The kenref module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The kenref module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------

This module is only compiled when PLUMED is configured with --enable-kenref,
in which case it links against the external KEnRef core library
(https://github.com/Smith-Group/KEnRef), which is distributed separately under
the BSD 3-Clause License. That licence is permissive and imposes no additional
restriction on this module or on PLUMED; see LICENSE_CORE.txt in the KEnRef
distribution for its text.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*
 * KEnRefBias_setup.cpp — frozen frame forwarder.
 *
 * The KEnRefBias constructor (one-time model + sub-indexing + driver setup) is hosted in the KEnRef
 * repository (src/plumedinterface/KEnRefBias_setup.cpp) so it can evolve with the KEnRef model
 * abstraction WITHOUT re-pushing this fork. It is compiled here, within PLUMED's build. The repo's
 * include path is supplied by `pkg-config --cflags kenref_plumed` (see this module's Makefile).
 */
#include "plumedinterface/KEnRefBias_setup.cpp"

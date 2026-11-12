# /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   Copyright (c) 2011-2016 The plumed team
#   (see the PEOPLE file at the root of the distribution for a list of names)
#
#   See http://www.plumed.org for more information.
#
#   This file is part of plumed, version 2.
#
#   plumed is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   plumed is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
"""Locate the optional kernel bundled alongside the Python interface."""

import os
import sys


def _bundled_kernel_path():
    """Return the installed private kernel path, if this is a bundled wheel."""
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runtime")
    if sys.platform == "win32":
        # as installed by CMake: the DLL next to plumed.exe
        path = os.path.join(root, "bin", "plumedKernel.dll")
    else:
        suffix = ".dylib" if sys.platform == "darwin" else ".so"
        path = os.path.join(root, "lib", "libplumedKernel" + suffix)
    return path if os.path.isfile(path) else None

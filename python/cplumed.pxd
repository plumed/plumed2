#/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
# This create cython wrappers to the main bits of the PLUMED libraray
#

from libcpp cimport bool

cdef extern from "exceptions_handler.h":
     cdef void exceptions_handler()

# Some of these functions are noexcept.
# We anyway use except + in case this changes later.
cdef extern from "Plumed.h" namespace "PLMD":
     cdef cppclass Plumed:
         Plumed() except +exceptions_handler
         void cmd_shaped "cmd" (const char*key, const void*val, const size_t* shape) except +exceptions_handler
         void cmd(const char*key, const void*val, size_t nelem) except +exceptions_handler
         void cmd(const char*key, const void*val) except +exceptions_handler
# see https://stackoverflow.com/questions/42610108/is-overloading-broken-in-cppclass-cython-c-definitions
         void cmd_int "cmd" (const char*key, int val) except +exceptions_handler
         void cmd_float "cmd" (const char*key, double val) except +exceptions_handler
         void cmd_mpi "cmd" (const char*key, const void*val) except +exceptions_handler
         void cmd(const char*key) except +exceptions_handler
         bool valid() except +exceptions_handler
         @staticmethod
         Plumed dlopen(const char*path) except +exceptions_handler
         @staticmethod
         Plumed makeValid() except +exceptions_handler

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

cdef extern from "Plumed.h":
    ctypedef struct plumed:
        pass # ignore content
    ctypedef struct plumed_safeptr:
        const void* ptr
        size_t nelem
        size_t* shape
        size_t flags
        # ignore other 
    ctypedef struct plumed_nothrow_handler:
        void* ptr
        void (*handler)(void*,int,const char*,const void*)
    ctypedef struct plumed_error:
        int code
        const char* what
        void* nested
        # ignore other members
    void plumed_cmd_safe_nothrow(plumed p,const char*key,plumed_safeptr safe,plumed_nothrow_handler nothrow)
    void plumed_error_set(void*ptr,int code,const char*what,const void* opt)
    void plumed_error_init(plumed_error* error)
    void plumed_error_finalize(plumed_error error)
    plumed plumed_create()
    plumed plumed_create_dlopen(const char*path)
    plumed plumed_create_invalid()
    void plumed_finalize(plumed p)
    int plumed_valid(plumed p)

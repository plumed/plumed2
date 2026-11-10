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
# This is a cython wrapper for the main parts of the PLUMED interface - the constructor and cmd
# The main purpose of this is to convert the python types to C types that PLUMED understands
#

# cython: binding=True, language_level=3

cimport cplumed  # This imports information from pxd file - including contents of this file here causes name clashes

import array
import re
import gzip
import math
import os
import sys
import warnings
import types

import cython
from cython.operator import dereference

if sys.version_info < (3,):
    raise ImportError("PLUMED >=2.6 only supports Python 3")

try:
     import numpy as np
     HAS_NUMPY=True
except ImportError:
     HAS_NUMPY=False

# flags used to pass type information
DEF type_void          = 0x10000*1
DEF type_nullptr       = 0x10000*2
DEF type_integral      = 0x10000*3
DEF type_real          = 0x10000*4

DEF type_value         = 0x2000000*1
DEF type_pointer       = 0x2000000*2
DEF type_const_pointer = 0x2000000*3

DEF type_nocopy        = 0x10000000
# Ordinary Python values cannot hold native MPI communicator storage. Keep
# recognized mpi4py communicators and the legacy null path unannotated.
DEF type_not_mpi_comm  = 0x20000000

from importlib import import_module as _import_module

# Resolve against the actual installed package name. Cython can otherwise reuse
# the default package's modules when differently named installations coexist.
_exceptions = _import_module("._exceptions", __package__)
PlumedError = _exceptions.PlumedError
LeptonError = _exceptions.LeptonError
_bundled_kernel_path = _import_module("._runtime", __package__)._bundled_kernel_path

cdef _decode_path(cplumed.plumed_error_filesystem_path path):
     # the kernel copies the native path (std::filesystem::path::native()), as
     # numbytes bytes without a terminator: char on Unix, wchar_t on Windows
     if not path.ptr:
         return None
     raw = (<char*>path.ptr)[:path.numbytes]
     if sys.platform == "win32":
         return raw.decode("utf-16-le", errors="replace")
     return os.fsdecode(raw)

cdef class Plumed:
     cdef cplumed.plumed c_plumed
     #cdef int initialized #not used
     def __cinit__(self):
         # this is guaranteed to be called once
         # we use it to make sure c_plumed is initialized correctly
         # in all other places we always finalize before resetting
         self.c_plumed=cplumed.plumed_create_invalid()
     def __dealloc__(self):
         # this is guaranteed to be called once
         # we use it to make sure c_plumed is finalized correctly
        cplumed.plumed_finalize(self.c_plumed)
     def __init__(self,kernel=None):
         cdef bytes py_kernel
         cdef char* ckernel
         if kernel is None and os.environ.get("PLUMED_KERNEL"):
            cplumed.plumed_finalize(self.c_plumed)
            self.c_plumed=cplumed.plumed_create()
            if not cplumed.plumed_valid(self.c_plumed):
                 raise RuntimeError("Error loading PLUMED_KERNEL at path " + os.environ["PLUMED_KERNEL"])
         else:
            if kernel is None:
                kernel=_bundled_kernel_path()
            if kernel is None:
                cplumed.plumed_finalize(self.c_plumed)
                self.c_plumed=cplumed.plumed_create()
                if not cplumed.plumed_valid(self.c_plumed):
                    raise RuntimeError("PLUMED not available, check your PLUMED_KERNEL environment variable")
            else:
                kernel=os.fspath(kernel)
                if not isinstance(kernel,str):
                    raise TypeError("kernel must be a path string or path-like object")
                py_kernel= kernel.encode() + b'\x00'  # Explicitly add null terminator
                ckernel = py_kernel
                cplumed.plumed_finalize(self.c_plumed)
                self.c_plumed=cplumed.plumed_create_dlopen(ckernel)
                if not cplumed.plumed_valid(self.c_plumed):
                    raise RuntimeError("Error loading PLUMED kernel at path " + kernel)
         api=array.array('i',[0])
         self.cmd( "getApiVersion", api)
         if(api[0]>=10):
           self.cmd("setNestedExceptions",1)
         self.cmd( "setRealPrecision", 8)
     def finalize(self):
         """ Explicitly finalize a Plumed object.

             Can be used in cases where one wants to make sure the Plumed object is finalized
             (so that all output files are flushed and all calculations are finalized) but there is
             a dangling reference to that Plumed object. Notice that after this call the self object
             will be invalid so that using cmd will raise an exception.

             It is also called by __exit__ in order to allow the following usage:
             ````
             with plumed.Plumed() as p:
                 p.cmd("init")
                 # ETC

             # p object will be finalized when exiting from this context
             ````
         """
         cplumed.plumed_finalize(self.c_plumed)
         self.c_plumed=cplumed.plumed_create_invalid()
     def __enter__(self):
         return self
     def __exit__(self, type, value, traceback):
        self.finalize()
     cdef raise_exception(self,cplumed.plumed_error error):
        if error.nested:
          try:
            self.raise_exception(dereference(<cplumed.plumed_error*> error.nested))
          except Exception as e:
            raise_from = e
        else:
          raise_from = None
        msg = error.what
        # messages from the system (e.g. FormatMessage on Windows) might not be UTF-8
        what = msg.decode("utf-8", errors="replace")

        path1 = _decode_path(error.path1)
        path2 = _decode_path(error.path2)

        # this map is from cython doc
        if error.code>=20300 and error.code<20400: # PLMD::Plumed::ExceptionTypeError
           raise TypeError(what) from raise_from
        elif error.code>=11400 and error.code<11500: # std::bad_alloc
           raise MemoryError(what) from raise_from
        elif error.code>=11000 and error.code<11100: # std::bad_typeid
           raise TypeError(what) from raise_from
        elif error.code>=11100 and error.code<11200: # std::bad_cast
           raise TypeError(what) from raise_from
        elif error.code>=10110 and error.code<10115: # std::domain_error
           raise ValueError(what) from raise_from
        elif error.code>=10105 and error.code<10110: # std::invalid_argument
           raise ValueError(what) from raise_from
        elif error.code==10229: # filesystem::filesystem_error
           raise OSError(error.error_code,what,path1,None,path2) from raise_from
        elif error.code>=10230 and error.code<10240: # std::ios_base::failure
           # Unfortunately, in standard C++ we have no way of distinguishing EOF
           # from other errors here; be careful with the exception mask
           raise IOError(what) from raise_from
        elif error.code>=10120 and error.code<10125: # std::out_of_range
           raise IndexError(what) from raise_from
        elif error.code>=10210 and error.code<10215: # std::overflow_error
           raise OverflowError(what) from raise_from
        elif error.code>=10205 and error.code<10210: # std::range_error
           raise ArithmeticError(what) from raise_from
        elif error.code>=10215 and error.code<10220: # std::underflow_error
           raise ArithmeticError(what) from raise_from
        elif error.code>=19900 and error.code<20000: # Lepton
           raise LeptonError(what) from raise_from
        elif error.code>=20000 and error.code<30000: # Plumed
           raise PlumedError(what) from raise_from
        else:
           raise RuntimeError(what) from raise_from
     cdef cmd_low_level(self,const char* ckey, const void* val,size_t nelem, size_t* shape,size_t flags):
         cdef cplumed.plumed_safeptr safe
         cdef cplumed.plumed_nothrow_handler nothrow
         cdef cplumed.plumed_error error
         safe.ptr=val
         safe.nelem=nelem
         safe.shape=shape
         safe.flags=flags
         cplumed.plumed_error_init(&error)
         nothrow.ptr=&error
         nothrow.handler=cplumed.plumed_error_set
         # see https://github.com/plumed/plumed2/pull/1129#issuecomment-2410867829
         with cython.nogil:
            cplumed.plumed_cmd_safe_nothrow(self.c_plumed,ckey,safe,nothrow)
         if(error.code):
           try:
             self.raise_exception(error)
           finally:
             cplumed.plumed_error_finalize(error)
     cdef cmd_ndarray_double(self, ckey, val):
         cdef double [:] abuffer = val.ravel()
         cdef size_t ashape[5]
         shape=val.shape
         assert len(shape)<5
         for i in range(len(shape)):
            ashape[i]=shape[i]
         ashape[len(shape)]=0
         self.cmd_low_level(ckey,&abuffer[0], 0, & ashape[0], sizeof(abuffer[0]) + type_real +  type_pointer + type_not_mpi_comm)
     cdef cmd_ndarray_int(self, ckey, val):
         cdef int [:] abuffer = val.ravel()
         cdef size_t ashape[5]
         shape=val.shape
         assert len(shape)<5
         for i in range(len(shape)):
            ashape[i]=shape[i]
         ashape[len(shape)]=0
         self.cmd_low_level(ckey,&abuffer[0], 0, & ashape[0], sizeof(abuffer[0]) + type_integral +  type_pointer + type_not_mpi_comm)
     cdef cmd_ndarray_longlong(self, ckey, val):
         cdef long long [:] abuffer = val.ravel()
         cdef size_t ashape[5]
         shape=val.shape
         assert len(shape)<5
         for i in range(len(shape)):
            ashape[i]=shape[i]
         ashape[len(shape)]=0
         self.cmd_low_level(ckey,&abuffer[0], 0, & ashape[0], sizeof(abuffer[0]) + type_integral +  type_pointer + type_not_mpi_comm)
     cdef cmd_array_double(self, ckey, val):
         cdef double [:] abuffer = val
         self.cmd_low_level(ckey,&abuffer[0], len(abuffer), NULL, sizeof(abuffer[0]) + type_real +  type_pointer + type_not_mpi_comm)
     cdef cmd_array_int(self, ckey, val):
         cdef int [:] abuffer = val
         self.cmd_low_level(ckey,&abuffer[0], len(abuffer), NULL, sizeof(abuffer[0]) + type_integral +  type_pointer + type_not_mpi_comm)
     cdef cmd_array_long(self, ckey, val):
         cdef long [:] abuffer = val
         self.cmd_low_level(ckey,&abuffer[0], len(abuffer), NULL, sizeof(abuffer[0]) + type_integral +  type_pointer + type_not_mpi_comm)
     cdef cmd_float(self, ckey, double val ):
         self.cmd_low_level(ckey,&val, 1, NULL, sizeof(val) + type_real +  type_value + type_not_mpi_comm)
     cdef cmd_int(self, ckey, int val):
         self.cmd_low_level(ckey,&val, 1, NULL, sizeof(val) + type_integral +  type_value + type_not_mpi_comm)
     cdef cmd_mpi(self, ckey, val):
         import mpi4py.MPI as MPI
         cdef size_t comm_addr = MPI._addressof(val)
         self.cmd_low_level(ckey,<void*>comm_addr, 0, NULL, type_void +  type_const_pointer)
     def cmd( self, key, val=None ):
         cdef bytes py_bytes = key.encode() + b'\x00'  # Explicitly add null terminator
         cdef char* ckey = py_bytes
         cdef char* cval
         if val is None :
            self.cmd_low_level(ckey,NULL,0,NULL,type_nullptr)
            return
         if isinstance(val, int ):
            self.cmd_int(ckey, val)
            return
         if isinstance(val, float ) :
            self.cmd_float(ckey, val)
            return
         if HAS_NUMPY and isinstance(val, np.ndarray) :
            # See https://numpy.org/doc/stable/user/basics.types.html
            if( val.dtype==np.double):
               self.cmd_ndarray_double(ckey, val)
            elif( val.dtype==np.intc ) :
               self.cmd_ndarray_int(ckey, val)
            elif( val.dtype==np.int_ ) :
               # np.int_ is 64 bit (if not np.intc), also on Windows with numpy>=2,
               # where long is 32 bit: pass it as long long
               self.cmd_ndarray_longlong(ckey, val)
            else :
               raise ValueError("ndarrays should be type np.double, np.intc, or np.int_")
            return
         if isinstance(val, array.array) :
            if( (val.typecode=="d" or val.typecode=="f") and val.itemsize==8):
               self.cmd_array_double(ckey, val)
            elif( (val.typecode=="i" or val.typecode=="I") ) :
               self.cmd_array_int(ckey, val)
            elif( (val.typecode=="l" or val.typecode=="L") ) :
               self.cmd_array_long(ckey, val)
            else :
               raise ValueError("arrays should be type double (size=8), int, or long")
            return
         if isinstance(val, str ) :
            py_bytes = val.encode() + b'\x00'  # Explicitly add null terminator
            cval = py_bytes
            # assume sizeof(char)=1
            self.cmd_low_level(ckey,cval,len(py_bytes), NULL,1 + type_integral + type_const_pointer + type_nocopy + type_not_mpi_comm)
            return
         if 'mpi4py' in sys.modules:
            import mpi4py.MPI as MPI
            if isinstance(val, MPI.Comm):
              self.cmd_mpi(ckey, val)
              return
         raise ValueError("Unknown value type ({})".format(str(type(val))))

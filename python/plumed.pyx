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

cimport cplumed  # This imports information from pxd file - including contents of this file here causes name clashes

from cpython cimport array
import array

try:
     import numpy as np
     HAS_NUMPY=True
except ImportError:
     HAS_NUMPY=False

cdef class Plumed:
     cdef cplumed.Plumed c_plumed
     def __cinit__(self,kernel=None):
         cdef bytes py_kernel
         cdef char* ckernel
         if kernel is None:
            self.c_plumed=cplumed.Plumed.makeValid()
            if not self.c_plumed.valid():
                 raise RuntimeError("PLUMED not available, check your PLUMED_KERNEL environment variable")
         else:
            py_kernel= kernel.encode()
            ckernel = py_kernel
            self.c_plumed=cplumed.Plumed.dlopen(ckernel)
            if not self.c_plumed.valid():
                 raise RuntimeError("Error loading PLUMED kernel at path " + kernel)
         cdef int pres = 8
         self.c_plumed.cmd( "setRealPrecision", <void*>&pres )
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
         self.c_plumed=cplumed.Plumed()
     def __enter__(self):
         return self
     def __exit__(self, type, value, traceback):
        self.finalize()
     def cmd_ndarray_real(self, ckey, val):
         cdef double [:] abuffer = val.ravel()
         self.c_plumed.cmd( ckey, <void*>&abuffer[0])
     def cmd_ndarray_int(self, ckey, val):
         cdef long [:] abuffer = val.ravel()
         self.c_plumed.cmd( ckey, <void*>&abuffer[0])
     cdef cmd_float(self, ckey, double val ):
         self.c_plumed.cmd( ckey, <void*>&val )
     cdef cmd_int(self, ckey, int val):
         self.c_plumed.cmd( ckey, <void*>&val)
     def cmd( self, key, val=None ):
         cdef bytes py_bytes = key.encode()
         cdef char* ckey = py_bytes
         cdef char* cval 
         cdef array.array ar
         if val is None :
            self.c_plumed.cmd( ckey, NULL )
         elif isinstance(val, (int,long) ):
            if key=="getDataRank" :
               raise ValueError("when using cmd with getDataRank option value must a size one ndarray")
            self.cmd_int(ckey, val)
         elif isinstance(val, float ) :
            if key=="getBias" :
               raise ValueError("when using cmd with getBias option value must be a size one ndarray")
            self.cmd_float(ckey, val) 
         elif HAS_NUMPY and isinstance(val, np.ndarray) : 
            if( val.dtype=="float64" ):
               self.cmd_ndarray_real(ckey, val)
            elif( val.dtype=="int64" ) : 
               self.cmd_ndarray_int(ckey, val)
            else :
               raise ValueError("ndarrys should be float64 or int64")
         elif isinstance(val, array.array) : 
            if( (val.typecode=="d" or val.typecode=="f") and val.itemsize==8): 
               ar = val
               self.c_plumed.cmd( ckey, <void*> ar.data.as_voidptr)
            elif( (val.typecode=="i" or val.typecode=="I") ) :
               ar = val
               self.c_plumed.cmd( ckey, <void*> ar.data.as_voidptr)
            else :
               raise ValueError("ndarrays should be double (size=8) or int")
         elif isinstance(val, basestring ) :
              py_bytes = val.encode()
              cval = py_bytes 
              self.c_plumed.cmd( ckey, <void*>cval )
         else :
            raise ValueError("Unknown value type ({})".format(str(type(val))))

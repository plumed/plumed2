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
import numpy as np
cimport numpy as np

cdef class Plumed:
     cdef cplumed.plumed c_plumed
     def __cinit__(self):
         self.c_plumed = cplumed.plumed_create()   #new cplumed.Plumed()
         cdef int pres = 8
         cplumed.plumed_cmd(self.c_plumed, "setRealPrecision", <void*>&pres )  
     def __dealloc__(self): 
         cplumed.plumed_finalize(self.c_plumed)

     def cmd_ndarray_real(self, ckey, val):
         cdef double [:] abuffer = val.ravel()
         cplumed.plumed_cmd(self.c_plumed, ckey, <void*>&abuffer[0])
     def cmd_ndarray_int(self, ckey, val):
         cdef long [:] abuffer = val.ravel()
         cplumed.plumed_cmd(self.c_plumed, ckey, <void*>&abuffer[0])
     cdef cmd_float(self, ckey, double val ):
         cplumed.plumed_cmd(self.c_plumed, ckey, <void*>&val )
     cdef cmd_int(self, ckey, int val):
         cplumed.plumed_cmd(self.c_plumed, ckey, <void*>&val)

     def cmd( self, key, val=None ):
         cdef bytes py_bytes = key.encode()
         cdef char* ckey = py_bytes
         cdef char* cval 
         cdef np.int_t[:] ibuffer
         cdef np.float64_t[:] dbuffer
         if val is None :
            cplumed.plumed_cmd( self.c_plumed, ckey, NULL )
         elif isinstance(val, (int,long) ):
            if key=="getDataRank" :
               raise ValueError("when using cmd with getDataRank option value must a size one ndarray")
            self.cmd_int(ckey, val)
         elif isinstance(val, float ) :
            if key=="getBias" :
               raise ValueError("when using cmd with getBias option value must be a size one ndarray")
            self.cmd_float(ckey, val) 
         elif isinstance(val, np.ndarray) : 
            if( val.dtype=="float64" ):
               self.cmd_ndarray_real(ckey, val)
            elif( val.dtype=="int64" ) : 
               self.cmd_ndarray_int(ckey, val)
            else :
               raise ValueError("ndarrys should be float64 or int64")
         elif isinstance(val, basestring ) :
              py_bytes = val.encode()
              cval = py_bytes 
              cplumed.plumed_cmd( self.c_plumed, ckey, <void*>cval )
         else :
            raise ValueError("Unknown value type ({})".format(str(type(val))))

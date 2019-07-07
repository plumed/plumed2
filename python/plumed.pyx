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

# cython: binding=True

cimport cplumed  # This imports information from pxd file - including contents of this file here causes name clashes

from cpython cimport array
import array
import re
import gzip

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


class FormatError(Exception):
    """Custom error reported by read_as_pandas.
    """
    pass

def read_as_pandas(file,chunksize=None,usecols=None,skiprows=None,nrows=None):
    """Import a plumed data file as a pandas dataset.

       file : either path to a file or open filed object

       chunksize : int, optional
           Return an iterable object.
       usecols : list-like or callable, optional
           Directly passed to pandas.
       skiprows : list-like, int or callable, optional
           Directly passed to pandas.
       nrows : int, optional
           Directly passed to pandas.

       Returns
       -------
       DataFrame (when chunksize is not provided) or iterable TextFileReader (when chunksize is provided).

       Comments
       --------

       Gzipped files are supported and automatically detected when a file name ends with '.gz'.

       `pandas` module is imported the first time this function is used. Since importing `pandas` is quite slow,
       the first call to this function will be significantly slower than the following ones.
       Following calls should be fast. The overall speed is comparable or better to loading with `numpy.loadtxt`.

       Examples
       --------

       colvar=plumed.read_as_pandas("COLVAR")
       print(colvar) # show the datasheet

       colvar=plumed.read_as_pandas("COLVAR",usecols=[0,4])
       print(colvar) # show the datasheet

       colvar=plumed.read_as_pandas("COLVAR",usecols=["time","distance"])
       print(colvar) # show the datasheet

       for chunk in plumed.read_as_pandas("COLVAR",chunksize=10):
           print(chunk) # show the datasheet. actually here you should process the chunk

       Limitations
       -----------

       1. Constants are presently ignored. As a consequence it is not possible to retrieve
       information such as ranges of variables.

       2. Text variables are not converted using standard PLUMED conversions (e.g. `pi` converted to
       the value of pi and arithmetic operations resolved).

       3. Only the initial header is read, which implies that files resulting from concatenating
       datasets with a different number of columns or different column names will not
       be read correctly.

       Issue 1 could be solved parsing the `#! SET` lines in python.
       Issue 2 could be solved calling `PLMD::Tools::convert`, which
       could be easily done through the `plumed_cmd` interface (TODO).
 
       Alternatively, all issues might be solved using `PLMD::IFile` for reading,
       which could be useful but possibly more complicated to implement.
    """
# importing pandas is pretty slow, so we only do it when needed
    import pandas as pd
# allow passing a string
    if isinstance(file,str):
        file=open(file)
# take care of gzipped files
    if re.match(".*\.gz",file.name):
        file = gzip.open(file.name,'rt')
# read first line
    line = file.readline()
    columns = line.split()
# check header
    if len(columns)<2:
        raise FormatError("Error reading PLUMED file "+file.name + ". Not enough columns")
    if columns[0] != "#!" or columns[1] != "FIELDS":
        raise FormatError("Error reading PLUMED file" +file.name + ". Columns: "+columns[0]+" "+columns[1])
# set column names
    columns = columns[2:]
# read the rest of the file
# notice that if chunksize was provided the result will be an iterable object
    return pd.read_csv(file, delim_whitespace=True, comment="#", header=None,names=columns,
                       usecols=usecols,skiprows=skiprows,nrows=nrows,chunksize=chunksize)

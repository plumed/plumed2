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

import array
import re
import gzip
import math
import sys
import warnings
import re
import types

if sys.version_info < (3,):
    raise ImportError("PLUMED 2.6 only supports Python 3")

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
     def cmd_array_real(self, ckey, val):
         cdef double [:] abuffer = val
         self.c_plumed.cmd( ckey, <void*>&abuffer[0])
     def cmd_array_int(self, ckey, val):
         cdef long [:] abuffer = val
         self.c_plumed.cmd( ckey, <void*>&abuffer[0])
     cdef cmd_float(self, ckey, double val ):
         self.c_plumed.cmd( ckey, <void*>&val )
     cdef cmd_int(self, ckey, int val):
         self.c_plumed.cmd( ckey, <void*>&val)
     def cmd( self, key, val=None ):
         cdef bytes py_bytes = key.encode()
         cdef char* ckey = py_bytes
         cdef char* cval
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
               self.cmd_array_real(ckey, val)
            elif( (val.typecode=="i" or val.typecode=="I") ) :
               self.cmd_array_int(ckey, val)
            else :
               raise ValueError("ndarrays should be double (size=8) or int")
         elif isinstance(val, str ) :
              py_bytes = val.encode()
              cval = py_bytes
              self.c_plumed.cmd( ckey, <void*>cval )
         else :
            raise ValueError("Unknown value type ({})".format(str(type(val))))


class FormatError(Exception):
    """Custom error reported by read_as_pandas.
    """
    pass

def _fix_file(file,mode):
    """Internal utility: returns a file open with mode.

       Takes care of opening file (if it receives a string)
       and or unzipping (if the file has ".gz" suffix).
    """
# allow passing a string
    if isinstance(file,str):
        file=open(file,mode)
# takes care of gzipped files
    if re.match(".*\.gz",file.name):
        file = gzip.open(file.name,mode)
    return file

def _build_convert_function(kernel=None):
    """Internal utility: returns a function that can be used for conversions.

       kernel : Plumed instance or str
           The object used to perform conversion.
           Pass a string to load a Plumed() instance giving the
           path to the libplumedKernel library, or pass None
           to load the default Plumed() instance.

       In case of failure, it writes a warning and returns None.

       Notice that this function will store a reference to the passed Plumed() object,
       thus potentially increasing its lifetime.
    """
    try:
# if necessary, load a kernel
        if not isinstance(kernel,Plumed):
            kernel=Plumed(kernel=kernel)
    except Exception:
        warnings.warn("cannot load PLUMED instance, conversions will not be available")
        return None
    try:
# define a function to convert data
        def convert_func(a):
            r=array.array('d',[float('nan')])
            convert_func.kernel.cmd("convert "+str(a),r)
            if math.isnan(r[0]):
               return a
            return r[0];
        convert_func.kernel=kernel
# check if convert_func is working correctly
        if (convert_func("pi")=="pi"):
            warnings.warn("PLUMED instance seems to have a non-working convert cmd, conversions do not work and will be disabled")
            return None
# set convert
        return convert_func
    except Exception:
        warnings.warn("PLUMED instance is too old, conversions do not work and will be disabled")
        return None

class Constants(list):
   """Custom class used to store plumed constants.
   """
   def __init__(self,l,kernel=None,convert=None):
       if(isinstance(l,dict)):
           for k in l:
              self.append((k,l[k]))
       else:
           self.extend(l)
       for i in range(len(self)):
           if(len(self[i])==2):
               if not convert:
                   convert=_build_convert_function(kernel)
               if convert:
                   self[i]=(self[i][0],convert(self[i][1]),str(self[i][1]))
               else:
                   self[i]=(self[i][0],self[i][1],str(self[i][1]))
           elif(len(self[i])==3):
               pass
           else:
               raise ValueError("plumed.Constants should be initialized with a list of 2- or 3-plets")

def read_as_pandas(file_or_path,enable_constants=True,enable_conversion=True,kernel=None,chunksize=None,usecols=None,skiprows=None,nrows=None,index_col=None):
    """Import a plumed data file as a pandas dataset.

       file_or_path : str or file
           Either string containing the path of the file or an already opened file object.

       enable_constants : str or boolean, optional (default is True)
           If 'columns', constants are read and added as constant columns.
           If 'metadata' or True, constants are read and stored as metadata.
           If 'no' or False, constants are not read at all.
       enable_conversion : str or boolean, optional (default is True)
           If 'constant' or True, only constants are converted.
           If 'all', all data are converted. Might be slow and probably useless.
           If 'no' or False, no data are converted.
       kernel : str or Plumed, optional
           The Plumed kernel used for conversions. If a string, it is interpreted
           as the path to a kernel. If None, the default Plumed loading procedure is used
           (with PLUMED_KERNEL env val). If an existing Plumed object, a pointer is stored
           and this object is used for conversion.

       chunksize : int, optional
           Return an iterable object. Useful to process large files in chunks.
       usecols : list-like or callable, optional
           Directly passed to pandas.
       skiprows : list-like, int or callable, optional
           Directly passed to pandas.
       nrows : int, optional
           Directly passed to pandas.
       index_col : int, str, sequence of int / str, or False, default None
           Directly passed to pandas.

       Returns
       -------
       By default, it returns a special subclass of pandas.DataFrame that includes
       metadata with constant values in an attribute named `plumed_constants`.
       If using `enable_constants='no'` or `enable_constants='columns'`,
       it returns a plain pandas.DataFrame.

       If `chunksize` is provided, it returns a special subclass of pandas.io.parsers.TextFileReader
       that can be iterated in order to read a file in chunks. Every iteration returns an object
       equivalent to the one that would have been returned with a call to
       read_pandas with chunksize=None (that is: either a pandas.DataFrame
       or a subclass of it).

       Comments
       --------

       Gzipped files are supported and automatically detected when a file name ends with '.gz'.

       `pandas` module is imported the first time this function is used. Since importing `pandas` is quite slow,
       the first call to this function will be significantly slower than the following ones.
       Following calls should be faster. The overall speed is comparable or better to loading with `numpy.loadtxt`.

       Examples
       --------

       colvar=plumed.read_as_pandas("COLVAR")
       print(colvar) # show the datasheet
       print(colvar.plumed_constants) # show the constant columns

       colvar=plumed.read_as_pandas("COLVAR",usecols=[0,4])
       print(colvar) # show the datasheet
       print(colvar.plumed_constants) # show the constant columns

       colvar=plumed.read_as_pandas("COLVAR",usecols=["time","distance"])
       print(colvar) # show the datasheet
       print(colvar.plumed_constants) # show the constant columns

       colvar=plumed.read_as_pandas("COLVAR",enable_constants='columns')
       print(colvar) # this dataframe will contain extra columns with the constants

       for chunk in plumed.read_as_pandas("COLVAR",chunksize=10):
           print(chunk) # show the datasheet. actually here you should process the chunk
           print(chunk.plumed_constants) # show the constant columns

       Limitations
       -----------

       Only the initial header is read, which implies that files resulting from concatenating
       datasets with a different number of columns or different column names will not
       be read correctly and that only constants set at the beginning of the file will be considered.

       This issues might be solved using `PLMD::IFile` for reading,
       which could be useful but possibly a bit complicated to implement.
    """

# importing pandas is pretty slow, so we only do it when needed
    import pandas as pd

# special classes used to attach metadata
# they are defined inside this function since they need pandas to be imported
# see https://pandas.pydata.org/pandas-docs/stable/development/extending.html
    class PlumedSeries(pd.Series):
        @property
        def _constructor(self):
            return PlumedSeries
        @property
        def _constructor_expanddim(self):
            return PlumedDataFrame

    class PlumedDataFrame(pd.DataFrame):
        _metadata=["plumed_constants"]
        @property
        def _constructor(self):
            return PlumedDataFrame
        @property
        def _constructor_sliced(self):
            return PlumedSeries

# auxiliary function to process a dataframe
# it is defined here since it requires PlumedDataFrame to be defined
    def process_dataframe(df,enable_constants,constants,convert_all):
        if convert_all: df=df.applymap(convert_all)
        if enable_constants=='columns':
            for c in constants: df[c[0]]=c[1]
        if enable_constants=='metadata':
            df=PlumedDataFrame(df)
            df.plumed_constants=Constants(constants)
        return df

# process arguments:
    if enable_conversion is True:
       enable_conversion='constants'
    if enable_conversion is False:
       enable_conversion='no'
    if enable_constants is True:
       enable_constants='metadata'
    if enable_constants is False:
       enable_constants='no'

# check arguments:
    if not (enable_conversion=='no' or enable_conversion=='constants' or enable_conversion=='all'):
        raise ValueError("enable_conversion not valid")
    if not (enable_constants=='no' or enable_constants=='metadata' or enable_constants=='columns'):
        raise ValueError("enable_conversion not valid")

# conversions functions:
    convert=None
    convert_all=None
# only create them if needed
    if (enable_conversion=='constants' and enable_constants) or enable_conversion=='all':
        convert=_build_convert_function(kernel)
# if necessary, set convert_all
        if enable_conversion=='all': convert_all=convert
         
# handle file
    file_or_path=_fix_file(file_or_path,'rt')

# read first line
    line = file_or_path.readline()
    columns = line.split()

# check header
    if len(columns)<2:
        raise FormatError("Error reading PLUMED file "+file_or_path.name + ". Not enough columns")
    if columns[0] != "#!" or columns[1] != "FIELDS":
        raise FormatError("Error reading PLUMED file" +file_or_path.name + ". Columns: "+columns[0]+" "+columns[1])

# read column names
    columns = columns[2:]

# read constants
    constants=[]
    if enable_constants!='no':
        while True:
            pos=file_or_path.tell()
            line = file_or_path.readline()
            file_or_path.seek(pos)
            if not line:
                break
            sets = line.split()
            if len(sets) < 4:
                break
            if sets[0]!="#!" or sets[1]!="SET":
                break
            if(convert):
                v=convert(sets[3])
            else:
                v=sets[3]
# name / value / string
            constants.append((sets[2],v,sets[3]))
            file_or_path.readline() # read again to go to next line

# read the rest of the file
# notice that if chunksize was provided the result will be an iterable object
    df=pd.read_csv(file_or_path, delim_whitespace=True, comment="#", header=None,names=columns,
                    usecols=usecols,skiprows=skiprows,nrows=nrows,chunksize=chunksize,index_col=index_col)

    if chunksize is None:
# just perform conversions and attach constants to the dataframe
        return process_dataframe(df,enable_constants,constants,convert_all)
    else:
# declare an alternate class that is iterable to read the file in chunks
        class TextFileReader(type(df)):
            """Subclass of pandas.io.TestFileReader, needed for storing constants"""
# some information (constant values and conversion function)
# should be stored in the class to be used while iterating on it
            def __init__(self,reader,enable_constants,constants,convert_all):
                self.TextFileReader=reader
                self.enable_constants=enable_constants
                self.constants=constants
                self.convert_all=convert_all
            def __next__(self):
# override __next__
                df=self.TextFileReader.__next__()
                return process_dataframe(df,self.enable_constants,self.constants,self.convert_all)
        return TextFileReader(df,enable_constants,constants,convert_all)

def write_pandas(df,file_or_path=None):
    """Save a pandas dataframe as a PLUMED file.

       df: pandas dataframe or derived class
           the dataframe. If it contains a list attribute `plumed_constants`, this is
           interpreted as a list of constants and written with `SET` lines.

       file_or_path: str, file, or None (default is None)
           path to the file to be written, or already opened file object.
           If None, stdout is used.

       Examples
       --------

       colvar=plumed.read_as_colvar("COLVAR")
       colvar["distance"]=colvar["distance"]*2
       plumed.write_pandas(colvar)
      
    """
# importing pandas is pretty slow, so we only do it when needed
    import pandas as pd
# check if there is an index. if so, write it as an additional field
    has_index=hasattr(df.index,'name') and df.index.name is not None
# check if there is a multi-index
    has_mindex=(not has_index) and hasattr(df.index,'names') and df.index.names[0] is not None
# writing multiindex is currently not supported
    if has_mindex:
        raise TypeError("Writing dataframes with MultiIndexes is not supported at this time")
# handle file
    if file_or_path is None:
        file_or_path=sys.stdout
    file_or_path=_fix_file(file_or_path,'wt')
# write header
    file_or_path.write("#! FIELDS")
    if has_index:
        file_or_path.write(" "+str(df.index.name))
    for n in df.columns:
        file_or_path.write(" "+str(n))
    file_or_path.write("\n")
# write constants
    if hasattr(df,"plumed_constants") and isinstance(df.plumed_constants,Constants):
        for c in df.plumed_constants:
# notice that string constants are written (e.g. pi) rather than the numeric ones (e.g. 3.14...)
            file_or_path.write("#! SET "+c[0]+" "+c[2]+"\n")
# write data
    for i in range(df.shape[0]):
        if has_index:
            file_or_path.write(" "+str(df.index[i]))
        for j in df.columns:
            file_or_path.write(" "+str(df[j][i]))
        file_or_path.write("\n")

def _guessplumedroot(kernel=None):
    """Guess plumed root.

       kernel: path to the plumed kernel

       In case the Plumed object cannot be created, try to launch a `plumed` executable and obtain the root
       dir from there.
    """
    try:
        import tempfile
        log=tempfile.mkstemp()[1]
        with Plumed(kernel) as p:
            p.cmd("setLogFile",log)
            p.cmd("init")
        i=0
        root=""
        with open(log) as fin:
            for line in fin:
                i=i+1
                if re.match("PLUMED: Root: ",line):
                    root=re.sub("PLUMED: Root: ","",line).rstrip("\n")
                    break
        if len(root)>0:
            return root
    except:
        pass
    # alternative solution, search for a plumed executable in the path
    import subprocess
    return subprocess.check_output(["plumed","--no-mpi","info","--root"]).decode("utf-8").rstrip()

def _readvimdict(plumedroot=None,kernel=None):
    """Read VIM dictionary given the path to PLUMED root and return (dictionary,doc).

       If plumedroot is not given, it is guessed by launching `plumed` executable.
       The dictionary is structured as follows:
       - The keys are the names of the actions (e.g. "RESTRAINT").
       - The values are dictionaries structured as follows:
         - The keys are the available options.
         - The value is a string describing the option type.
       For instance `dictionary["RESTRAINT"]["NUMERICAL_DERIVATIVES"]=="(flag)"`.

       The doc is a dictionary structured as follows:
       - The keys are the names of the actions (e.g. "RESTRAINT").
       - The values are docstrings for the corresponding actions.
    """
    if plumedroot is None:
        plumedroot=_guessplumedroot(kernel)
    syntax_file=plumedroot + "/vim/syntax/plumed.vim"
    help_dir=plumedroot + "/vim/help"
# local dictionary, read from VIM
    plumedDictionary={}
    pattern = re.compile("^let b:plumedDictionary\[.*$")
    with open(syntax_file) as fin:
        for line in fin:
            if pattern.match(line):
                line=re.sub("^let b:","",line)
                exec(line,{'__builtins__': None},{'plumedDictionary':plumedDictionary})
    ret={}
    doc={}
    for action in plumedDictionary:
        ret[action]={}
        doc[action]={}
# read documentation
        with open(help_dir + "/" + action + ".txt") as fin:
            thisdoc=""
            for line in fin:
                thisdoc+=line
                if(line.rstrip("\n")=="****************************************"):
                    thisdoc="Create action " + action + "\n"
            doc[action]=thisdoc
# remove LaTex stuff
            doc[action]=doc[action].replace("\\f$","").replace("\\","")
# read dictionary
        for opt in plumedDictionary[action]:
# skip label (it is added automatically)
            if opt["menu"] != "(label)":                
                ret[action][re.sub("=$","",opt["word"])]=opt["menu"]
    return ret,doc

def _create_functions(dictionary,*,doc=None,append_underscores=False):
    """Create functions given dictionary and, optionally, documentation.

       Given a dictionary a doc produced with _readvimdict, it returns a dictionary
       containing the functions corresponding to each action. The functions are stored
       in string that should be then evaluated with exec. For each action (say, "RESTRAINT")
       we define both a function `def RESTRAINT and a docstring `RESTRAINT.__doc__`.
       These functions are only using builtins and functions from the _format_tools dictionary.
    """
    functions={}
    for action in dictionary:
        # skip actions with incorrect name
        if re.match(".*-.*",action):
            continue
        fname=action
        if append_underscores:
            fname+="__"
        string=""
        string+="def " + fname
        string+="(self,LABEL=\"\",verbatim=None"
        if len(dictionary[action])>0:
            string+=",*"
        for w in dictionary[action]:
            # skip arguments with incorrect name
            if re.match(".*-.*",w):
                continue
            if dictionary[action][w]=="(flag)" :
                string+="," + w + "=False"
            else:
                string+="," + w + "=None"
        string+=",**kwargs):\n"
        string+="  ret=\"\"\n"
        string+="  ret+=_format_label(self,LABEL)\n"
        string+="  ret+=\"" + action + "\"\n"
        string+="  retlist=[]\n"
        for w in dictionary[action]:
            # skip arguments with incorrect name
            if re.match(".*-.*",w):
                continue
            t=dictionary[action][w]
            if t=="(flag)" :
                string+="  retlist.append(_format_flag(self,\"" + w + "\"," + w + "))\n"
            elif t=="(option)" :
                string+="  retlist.append(_format_opt(self,\"" + w + "\"," + w + "))\n"
            elif t=="(numbered)":
                string+="  retlist.append(_format_numbered(self,\"" + w + "\"," + w + "))\n"
            else:
                raise TypeError("error parsing dictionary, unknown type "+t)
# now process kwargs to allow numbered arguments
        string+="  for arg in kwargs:\n"
        string+="      import re\n"
        string+="      allowed=[]\n"
        for x in dictionary[action]:
          string+="      allowed.append(\"" + x + "\")\n"
        string+="      if not re.sub(\"[0-9]*$\",\"\",arg) in allowed:\n"
        string+="         raise TypeError(\"unknown arg \" + arg)\n"
        string+="      retlist.append(_format_anything(self,arg,kwargs[arg]))\n"
# sorting is necessary to make the line reproducible in regtests
        string+="  retlist.sort()\n"
        string+="  for x in retlist:\n"
        string+="      if(len(x)>0):\n"
        string+="          ret+=x\n"
        string+="  ret+=_format_verbatim(self,verbatim)\n"
        string+="  return _format_return(self,ret)\n"
# set docstring
        if doc is not None:
           string+=fname + ".__doc__ = \"\"\"\n"
           string+=doc[action]
           string+="\"\"\"\n"
        functions[fname]=string
    return functions

# formatting tools

class _numbered():
    def __init__(self,*args):
# this is experimental:
# it allows calling _numbered(arg1,arg2,arg3)
# it is however on purpose not implemented in the numbered() method below
        if(len(args)==1):
            if isinstance(args[0],dict):
                self.arg=args[0]
            elif hasattr(args[0],'__iter__'):
                self.arg={}
                i=0
                for x in args[0]:
                     self.arg[i]=x
                     i+=1
            else:
                raise TypeError("when calling numbered with 1 argument, it should be a list/tuple/dictionary")
        else:
            self.arg={}
            i=0
            for x in args:
                 self.arg[i]=x
                 i+=1

class _replicas():
    def __init__(self,*args):
# this is experimental:
# it allows calling _replicas(arg1,arg2,arg3)
# it is however on purpose not implemented in the replicas() method below
        if(len(args)==1):
          if hasattr(args[0],'__iter__'):
            self.arg=args[0]
          else:
            raise TypeError("when calling replicas with 1 argument, it should be a list/tuple")
        else:
          self.arg=args

## tool to format at strings (optional)

def _format_at_one_residue(builder,name,residue,chain):
      if isinstance(chain,int):
        chain=str(chain)
      digit=False
      for i in chain:
        if i.isdigit():
           digit=True
      if digit:
        chain=chain+"_"
      if isinstance(residue,int) or isinstance(residue,str):
        return "@" + name + "-" + chain + str(residue)
      else:
        assert False
    
def _format_at_one_chain(builder,name,residue,chain):
      res=""
      if hasattr(residue,'__iter__') and not isinstance(residue,str):
        for x in residue:
              res+=builder._separator + _format_at_one_residue(builder,name,x,chain)
      else:
        res+=builder._separator + _format_at_one_residue(builder,name,residue,chain)
              
      return res

def _format_at(builder,name,residue,chain=""):
      res=""
      if hasattr(chain,'__iter__') and not isinstance(chain,str):
         for x in chain:
             res+=_format_at_one_chain(builder,name,residue,x)
      else:
         res+=_format_at_one_chain(builder,name,residue,chain)
      return res[len(builder._separator):]

class _at():
    def __init__(self,builder):
      import weakref
      self._builder=weakref.ref(builder)
      _at_global=["mdatoms","allatoms","water","nucleic","protein","water","ions","hydrogens","nonhydrogens"]
      for x in _at_global:
        exec("self." + x + "=\"@" + x + "\"",None,{"self":self})
      _at_residue=["phi","psi","omega","chi1","alpha","beta","gamma","delta","epsilon","zeta","v0","v1","v2","v3","v4","chi","back","sugar","base","lcs"]
      for x in _at_residue:
            ldict={}
            exec("def " + x + "(self,residue,chain=\"\"):\n          return _format_at(self._builder(),\"" + x + "\",residue,chain)\n",{"_format_at":_format_at},ldict)
            exec("self." + x + " = types.MethodType( ldict['" + x + "'],self)",None,{"self":self, "types":types, "ldict":ldict})
    def __call__(self,name,residue,chain=""):
       return _format_at(self._builder(),name,residue,chain)

## end of tool to format at strings

def _fix_braces(builder,arg,comma):
    """ Fix braces.

        If comma is true, consider comma as a separator. Otherwise, only consider space as a separator.
    """
    always=not builder._minimize_braces
    # recursively remove braces to find non-matching ones
    tmp=arg
    go=1
    while go>0:
      (tmp,go)=re.subn("{[^{}]*}","",tmp)
    if re.match(".*[{}]",tmp):
        raise ValueError("option \"" + arg + "\" contains nonmatching braces")
    # @replicas: beginning strings are not treated.
    # the reason is that these strings are expected to be already fixed by possibly adding braces after the :
    if re.match("@replicas:",arg):
        return arg
    if always or arg=="" or (not comma and re.match(".*[ {}]",arg)) or ( comma and re.match(".*[ {},]",arg) ):
        return "{" + arg + "}"
    return arg

def _format_single(builder,arg,level=0):
    """Format a single argument.

       Level keeps track of recursions.
    """
    import numbers
    if builder._pre_format is not None:
       arg=builder._pre_format(arg)
# only import if already loaded
# this is to avoid slow startup times.
    if builder._enable_mda_groups and 'MDAnalysis' in sys.modules:
       import MDAnalysis
       if isinstance(arg,MDAnalysis.core.groups.AtomGroup):
         arg=arg.indices+1
    if isinstance(arg,str):
        return re.sub("[\n\t]"," ",arg)
    if isinstance(arg,numbers.Number):
        return str(arg)
    if isinstance(arg,dict):
        raise TypeError("options cannot be a dictionary")
    if isinstance(arg,_replicas):
        if level>1:
            raise TypeError("@replica syntax only allowed for scalar or rank 1 vectors")
        return "@replicas:" + _fix_braces(builder,_format_single(builder,arg.arg,level+1),comma=True)
    if hasattr(arg,'__iter__'):
        string=""
        for x in arg:
            string+=builder._separator + _fix_braces(builder,_format_single(builder,x,level+1),comma=True)
        return string[len(builder._separator):]
    raise TypeError("options should be string, number or iterable")

def _format_numbered(builder,name,arg):
    """Format a numbered argument."""
    if isinstance(arg,_numbered):
        arg=arg.arg
        ret=""
        for x in arg:
            if not isinstance(x,int):
                raise TypeError("numbered types should have integer keys")
            ret+=_format_opt(builder,name+str(x),arg[x])
        return ret
    return _format_opt(builder,name,arg)

def _format_flag(builder,name,arg):
    """Format a flag."""
    if arg is None:
        return ""
    if isinstance(arg,bool):
        if arg:
            return " " + name
        else:
            return ""
    raise TypeError(name + " should be of bool type")

def _format_label(builder,label):
    """Format a label."""
    if label is None:
        return ""
    if isinstance(label,str):
        if len(label)==0:
            return ""
        return label+": "
    raise TypeError("label should be of str type")

def _format_opt(builder,name,arg):
    """Format an option."""
    if arg is None:
        return ""
    if isinstance(arg,_replicas):
        string="@replicas:" + _fix_braces(builder,_format_single(builder,arg.arg),comma=False)
    else:
        string=_fix_braces(builder,_format_single(builder,arg),comma=False)
    return " " + name + "=" + string

def _format_return(builder,ret):
    """Format the return statement"""
    if not builder._toplumed is None :
        try:
          builder._toplumed.cmd("readInputLine",ret)
        except:
          builder.last_failure=ret
          raise
    if not builder._tofile is None:
        builder._tofile.write(ret + "\n")
    builder.history.append(ret)
    if builder._post_format is not None:
        return builder._post_format(ret)
    return ret + "\n"

def _format_verbatim(builder,verbatim):
  if not verbatim is None:
    if(len(verbatim)>0):
      return " "+re.sub("[\n\t]"," ",verbatim)
  return ""

def _format_anything(builder,name,arg):
   """Choose format based on arg type"""
   ret=""
   if name == "verbatim":
       ret+=_format_verbatim(builder,arg)
   elif isinstance(arg,bool) : 
       ret+=_format_flag(builder,name,arg)
   elif isinstance(arg,_numbered):
       ret+=_format_numbered(builder,name,arg)
   else:
       ret+=_format_opt(builder,name,arg)
   return ret

# this functions are made available for the exec commands
_format_tools={
    "_format_numbered":_format_numbered,
    "_format_flag":_format_flag,
    "_format_label":_format_label,
    "_format_opt":_format_opt,
    "_format_return":_format_return,
    "_format_verbatim":_format_verbatim,
    "_format_anything":_format_anything,
}

class InputBuilder:
    """Object used to construct plumed input files.

       An instance of this object can be used to construct plumed input files.
       Check the constructor to see all the available options.
    """
    def __init__(self,
                 *,
                 tofile=None,
                 toplumed=None,
                 kernel=None,
                 append_underscores=False,
                 comma_separator=False,
                 minimize_braces=True,
                 enable_at=True,
                 load_dict=True,
                 enable_mda_groups=True,
                 post_format=None,
                 pre_format=None):
        """Constructor.

           Parameters
           ----------
           tofile:
             PLUMED input is also forwarded to this file. Useful for preparing input files.
           toplumed:
             PLUMED input is also forwarded to this PLUMED object. Useful for testing interactively the input.
           kernel:
             Path to PLUMED kernel. By default, a new PLUMED object is created with default kernel (PLUMED_KERNEL).
           append_underscores:
             Append the two underscores to method names to avoid autocompletion problems in older ipython versions.
           comma_separator:
             Use comma as a separator rather than a space.
           minimize_braces:
             Minimize the number of braces. Setting this to False triggers a bug in the PLUMED parser
             (fixed in 2.4.5 and 2.5.1).
           enable_at:
             Enable `ib.at.chi(1,"A")` syntax.
           load_dict:
             Load full dictionary. Set to false to make initialization faster, at the price of loosing autocompletion.
           enable_mda_groups:
             Enable MDAnalysis groups. Notice that MDAnalysis is not explicitly imported. Only set to false if you
             have a non-working MDAnalysis module installed.
        """
        self._toplumed=toplumed
        self._tofile=tofile
        if comma_separator:
          self._separator=","
        else:
          self._separator=" "
        self._minimize_braces=minimize_braces
        self._append_underscores=append_underscores
        # history of plumed lines
        self.history=[]
        # last line leading to a failure in a Plumed object
        self.last_failure=""
        self._has_dict=False
        self._enable_mda_groups=enable_mda_groups
        self._pre_format=pre_format
        self._post_format=post_format
        if load_dict:
            # stored for debugging:
            self._vimdict,self._doc=_readvimdict(kernel=kernel)
            # stored for debugging:
            self._functions=_create_functions(self._vimdict,doc=self._doc,append_underscores=append_underscores)
            for action in self._functions:
                ldict={}
                try:
                    # create free-standing functions
                    exec(self._functions[action],_format_tools,ldict)
                except:
                    print("ERROR interpreting " + action)
                    print(self._functions[action])
                    raise
                # create the object method
                exec("self." + action + " = types.MethodType( ldict['" + action + "'],self)",None,{"self":self, "types":types, "ldict":ldict})
            self._has_dict=True

        if enable_at:
            self.at=_at(self)

    def __call__(self,action,LABEL="",verbatim=None,**kwargs):
        ret=""
        ret+=_format_label(self,LABEL)
        ret+=action
        retlist=[]
        for arg in sorted(kwargs):
            retlist.append(_format_anything(self,arg,kwargs[arg]))
        retlist.sort()
        for x in retlist:
            if(len(x)>0):
                ret+=x
        ret+=_format_verbatim(self,verbatim)
        return _format_return(self,ret)

    def __getattr__(self, name):
       if self._has_dict:
           class _callme:
               def __init__(self,builder,name):
                   self._builder=builder
                   self._name=name
                   self.__doc__=getattr(self._builder,name).__doc__
               def __call__(self,LABEL="",verbatim=None,**kwargs):
                   func=getattr(self._builder,self._name)
                   return func(LABEL,verbatim,**kwargs)
           if self._append_underscores and not re.match(".*__$",name):
                name__=name+"__"
                if name__ in self._functions.keys():
                    return _callme(self,name__)
           if not self._append_underscores and re.match(".*__$",name):
                name__=re.sub("__$","",name)
                if name__ in self._functions.keys():
                    return _callme(self,name__)
           raise AttributeError("unknown method " + name)
       else:
           class _callme:
               def __init__(self,builder,name):
                   self._builder=builder
                   self._name=name
                   self.__doc__="dynamic method printing " + name
               def __call__(self,LABEL="",verbatim=None,**kwargs):
                   return self._builder(self._name,LABEL,verbatim,**kwargs)
           if re.match(".*[^_]__$",name):
               name=name[:-2]
           return _callme(self,name)

    def verbatim(self,line):
        """Create an arbitrary line."""
        return _format_return(self,re.sub("[\n\t]"," ",line))

    def numbered(self,arg):
        """Shortcut for numbered syntax.

        Accepts either a list/tuple or a dictionary with integer keys.
        """
        return _numbered(arg)

    def replicas(self,arg):
        """Shortcut for replica syntax.

        Accepts a list/tuple.
        """
        return _replicas(arg)




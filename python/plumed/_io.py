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
"""Readers, writers and HILLS utilities for PLUMED data files."""

import array
import gzip
import math
import re
import sys
import warnings

from ._plumed import Plumed

try:
    import numpy as np
except ImportError:
    pass


class FormatError(Exception):
    """Custom error reported by read_as_pandas."""

    pass


def _fix_file(file, mode):
    """Internal utility: returns a file open with mode.

    Takes care of opening file (if it receives a string)
    and or unzipping (if the file has ".gz" suffix).
    """
    # allow passing a string
    if isinstance(file, str):
        file = open(file, mode)
    # takes care of gzipped files
    if re.match(r".*\.gz", file.name):
        file = gzip.open(file.name, mode)
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
        if not isinstance(kernel, Plumed):
            kernel = Plumed(kernel=kernel)
    except Exception:
        warnings.warn("cannot load PLUMED instance, conversions will not be available")
        return None
    try:
        # define a function to convert data
        def convert_func(a):
            r = array.array("d", [float("nan")])
            convert_func.kernel.cmd("convert " + str(a), r)
            if math.isnan(r[0]):
                return a
            return r[0]

        convert_func.kernel = kernel
        # check if convert_func is working correctly
        if convert_func("pi") == "pi":
            warnings.warn(
                "PLUMED instance seems to have a non-working convert cmd, conversions do not work and will be disabled"
            )
            return None
        # set convert
        return convert_func
    except Exception:
        warnings.warn(
            "PLUMED instance is too old, conversions do not work and will be disabled"
        )
        return None


class Constants(list):
    """Custom class used to store plumed constants."""

    def __init__(self, l, kernel=None, convert=None):
        if isinstance(l, dict):
            for k in l:
                self.append((k, l[k]))
        else:
            self.extend(l)
        for i in range(len(self)):
            if len(self[i]) == 2:
                if not convert:
                    convert = _build_convert_function(kernel)
                if convert:
                    self[i] = (self[i][0], convert(self[i][1]), str(self[i][1]))
                else:
                    self[i] = (self[i][0], self[i][1], str(self[i][1]))
            elif len(self[i]) == 3:
                pass
            else:
                raise ValueError(
                    "plumed.Constants should be initialized with a list of 2- or 3-plets"
                )


def read_as_pandas(
    file_or_path,
    enable_constants=True,
    enable_conversion=True,
    kernel=None,
    chunksize=None,
    usecols=None,
    skiprows=None,
    nrows=None,
    index_col=None,
):
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
    read_as_pandas with chunksize=None (that is: either a pandas.DataFrame
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
        _metadata = ["plumed_constants"]

        @property
        def _constructor(self):
            return PlumedDataFrame

        @property
        def _constructor_sliced(self):
            return PlumedSeries

    # auxiliary function to process a dataframe
    # it is defined here since it requires PlumedDataFrame to be defined
    def process_dataframe(df, enable_constants, constants, convert_all):
        if convert_all:
            df = df.applymap(convert_all)
        if enable_constants == "columns":
            for c in constants:
                df[c[0]] = c[1]
        if enable_constants == "metadata":
            df = PlumedDataFrame(df)
            df.plumed_constants = Constants(constants)
        return df

    # process arguments:
    if enable_conversion is True:
        enable_conversion = "constants"
    if enable_conversion is False:
        enable_conversion = "no"
    if enable_constants is True:
        enable_constants = "metadata"
    if enable_constants is False:
        enable_constants = "no"

    # check arguments:
    if not (
        enable_conversion == "no"
        or enable_conversion == "constants"
        or enable_conversion == "all"
    ):
        raise ValueError("enable_conversion not valid")
    if not (
        enable_constants == "no"
        or enable_constants == "metadata"
        or enable_constants == "columns"
    ):
        raise ValueError("enable_conversion not valid")

    # conversions functions:
    convert = None
    convert_all = None
    # only create them if needed
    if (
        enable_conversion == "constants" and enable_constants
    ) or enable_conversion == "all":
        convert = _build_convert_function(kernel)
        # if necessary, set convert_all
        if enable_conversion == "all":
            convert_all = convert

    # handle file
    file_or_path = _fix_file(file_or_path, "rt")

    # read first line
    line = file_or_path.readline()
    columns = line.split()

    # check header
    if len(columns) < 2:
        raise FormatError(
            "Error reading PLUMED file " + file_or_path.name + ". Not enough columns"
        )
    if columns[0] != "#!" or columns[1] != "FIELDS":
        raise FormatError(
            "Error reading PLUMED file"
            + file_or_path.name
            + ". Columns: "
            + columns[0]
            + " "
            + columns[1]
        )

    # read column names
    columns = columns[2:]

    # read constants
    constants = []
    if enable_constants != "no":
        while True:
            pos = file_or_path.tell()
            line = file_or_path.readline()
            file_or_path.seek(pos)
            if not line:
                break
            sets = line.split()
            if len(sets) < 4:
                break
            if sets[0] != "#!" or sets[1] != "SET":
                break
            if convert:
                v = convert(sets[3])
            else:
                v = sets[3]
            # name / value / string
            constants.append((sets[2], v, sets[3]))
            file_or_path.readline()  # read again to go to next line

    # read the rest of the file
    # notice that if chunksize was provided the result will be an iterable object
    df = pd.read_csv(
        file_or_path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=columns,
        usecols=usecols,
        skiprows=skiprows,
        nrows=nrows,
        chunksize=chunksize,
        index_col=index_col,
    )

    if chunksize is None:
        # just perform conversions and attach constants to the dataframe
        return process_dataframe(df, enable_constants, constants, convert_all)
    else:
        # declare an alternate class that is iterable to read the file in chunks
        class TextFileReader(type(df)):
            """Subclass of pandas.io.TestFileReader, needed for storing constants"""

            # some information (constant values and conversion function)
            # should be stored in the class to be used while iterating on it
            def __init__(self, reader, enable_constants, constants, convert_all):
                self.TextFileReader = reader
                self.enable_constants = enable_constants
                self.constants = constants
                self.convert_all = convert_all

            def __next__(self):
                # override __next__
                df = self.TextFileReader.__next__()
                return process_dataframe(
                    df, self.enable_constants, self.constants, self.convert_all
                )

        return TextFileReader(df, enable_constants, constants, convert_all)


def write_pandas(df, file_or_path=None):
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
    import pandas as pd  # noqa: F401 - Preserve lazy dependency validation.

    # check if there is an index. if so, write it as an additional field
    has_index = hasattr(df.index, "name") and df.index.name is not None
    # check if there is a multi-index
    has_mindex = (
        (not has_index) and hasattr(df.index, "names") and df.index.names[0] is not None
    )
    # writing multiindex is currently not supported
    if has_mindex:
        raise TypeError(
            "Writing dataframes with MultiIndexes is not supported at this time"
        )
    # handle file
    if file_or_path is None:
        file_or_path = sys.stdout
    file_or_path = _fix_file(file_or_path, "wt")
    # write header
    file_or_path.write("#! FIELDS")
    if has_index:
        file_or_path.write(" " + str(df.index.name))
    for n in df.columns:
        file_or_path.write(" " + str(n))
    file_or_path.write("\n")
    # write constants
    if hasattr(df, "plumed_constants") and isinstance(df.plumed_constants, Constants):
        for c in df.plumed_constants:
            # notice that string constants are written (e.g. pi) rather than the numeric ones (e.g. 3.14...)
            file_or_path.write("#! SET " + c[0] + " " + c[2] + "\n")
    # write data
    for i in range(df.shape[0]):
        if has_index:
            file_or_path.write(" " + str(df.index[i]))
        for j in df.columns:
            file_or_path.write(" " + str(df[j][i]))
        file_or_path.write("\n")


def hills_time_average(
    hills, tofile=None, *, t0=None, t1=None, frac0=None, frac1=None, inplace=False
):
    """Compute a time-averaged hills file.

    hills: dataframe containing a hills file, or file, or filename

    tofile: write the resulting dataframe to a file or filename

    t0: initial time for averaging

    t1: final time for averaging

    frac0: initial fraction for averaging

    frac1: final fraction for averaging

    inplace: pass true to avoid copying the dataframe (for very large dataframes)

    Returns
    -------

    if inplace=False and tofile=None (default): a modified copy of the dataframe
    if inplace=True or tofile!=None: None

    This tool takes as an input a pandas dataframe read from a HILLS file and
    returns an equivalent dataframe where Gaussian heights have been scaled
    with a windowing function so as to effectively result in a time-averaged
    potential. By default, the average is taken along to entire list based on the
    value of the time column. The averaging windows can be specified in two ways:
    - using the time field, from t0 to t1. If t0 (t1) is omitted, its default is
      to be set to np.min(hills.time) (np.max(hills.time))
    - using the index of the hill, from frac0 to frac1 specified as a fraction
      of the entire length; that is, frac0=0.5 implies avering from the mid of the
      simulation. If frac0 (frac1) is omitted, its default is to be set to 0.0 (1.0)

    If tofile is specified, nothing is returned and the resulting dataframe is written
    on file.

    Examples
    --------

    plumed.hills_time_average("HILLS","HILLSOUT") # time average along entire file

    plumed.hills_time_average("HILLS","HILLSOUT1",t0=100,t1=200) # time average between 100 and 200 ps

    plumed.hills_time_average("HILLS","HILLSOUT2",frac0=0.5) # time average in second half of the file

    df=plumed.hills_time_average("HILLS") # time average and return dataframe

    df=plumed.read_as_pandas("HILLS") # first read hills file
    df=plumed.hills_time_average(df) # then perform time average


    """
    import pandas as pd

    if not isinstance(hills, pd.DataFrame):
        hills = read_as_pandas(hills)
    elif not inplace:
        hills = hills.copy()

    use_time = t0 is not None or t1 is not None
    use_index = frac0 is not None or frac1 is not None

    if use_time and use_index:
        raise ValueError("cannot use simultaneously time and fraction")

    if not use_index and not use_time:
        use_time = True

    if use_time:
        if t0 is None:
            t0 = np.min(hills.time)
        if t1 is None:
            t1 = np.max(hills.time)
        w = np.clip((t1 - hills.time) / (t1 - t0), 0.0, 1.0)
    else:
        if frac0 is None:
            frac0 = 0
        if frac1 is None:
            frac1 = 1.0
        i0 = int(len(hills) * frac0)
        i1 = int(len(hills) * frac1)
        w = np.hstack(
            (np.ones(i0), np.linspace(1.0, 0.0, i1 - i0), np.zeros(len(hills) - i1))
        )
    hills.height *= w

    if tofile is not None:
        write_pandas(hills, tofile)
    elif not inplace:
        return hills

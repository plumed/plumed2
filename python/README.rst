Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

The package supports CPython 3.11 through 3.15. The Linux and macOS wheels
include a basic PLUMED kernel and the ``plumed`` command. Linux wheels target
the manylinux2014 baseline (glibc 2.17). The macOS deployment target is 10.15
for Intel and at least 11.0 for Apple Silicon; the wheel tag records the actual
minimum, which also depends on the Python version. The bundled kernel is used
by default. You
can select another kernel for a process with ``PLUMED_KERNEL`` or for one
object with ``Plumed(kernel=...)``.

To override the bundled kernel for all new PLUMED objects, set `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple PLUMED versions on your system with environment modules, setting ``PLUMED_KERNEL`` selects that external version instead of the bundled one.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This allows you to mix different PLUMED versions in the same Python script.

Source builds bundle a kernel by default. Set ``PLUMED_BUILD_KERNEL=0`` when building the package to produce a wrapper-only distribution.

Package layout
--------------

The public API remains available through ``import plumed``. Its compiled
native interface is the private ``plumed._plumed`` extension. InputBuilder,
metadata parsing, data-file IO and HILLS averaging are implemented in ordinary
Python modules, with their existing names reexported from ``plumed``. The native
command path, including argument marshalling and GIL release, remains compiled.
Bundled executables, libraries and resources live under ``plumed/runtime``. No separate top-level
``plumed_runtime`` package is installed. Wrapper-only installations use the same
Python package layout without bundled runtime resources.

``runtime`` is a private installation prefix: executables are in ``bin``,
shared libraries in ``lib``, headers in ``include/plumed``, and PLUMED resources
in ``lib/plumed``. The wheel builder uses ``make install`` to populate this prefix, including
Fortran sources, wrapper archives, build integration files, shell tools, JSON
and Vim metadata, and Bash completion. Python installation is disabled in this
native build because the wheel separately supplies the Python package. The
wheel build enables ``__PLUMED_PIP_BUNDLE`` for resource discovery and applies
platform loader repair afterward. This is a normal PLUMED installation, entirely inside
the Python package. The bundled shared kernel discovers its resources and
headers relative to its library location, including when loaded directly from
Python or a native host. ``PLUMED_ROOT``, ``PLUMED_INCLUDEDIR`` and
``PLUMED_HTMLDIR`` retain their normal installed meanings. HTML manuals are not
bundled; ``plumed info --user-doc`` falls back to the online manual.

Pandas text and HDF5 tables
--------------------------

``read_as_pandas`` and ``write_pandas`` support both PLUMED text tables and
PLUMED HDF5 tables. Install the optional pandas and HDF5 dependencies with
``pip install 'plumed[hdf5]'`` (or install ``h5py`` separately)::

    import plumed

    data = plumed.read_as_pandas("COLVAR.h5")
    plumed.write_pandas(data, "COLVAR-copy.h5")
    plumed.write_pandas(data, "COLVAR")  # convert back to text

The filename extensions ``.h5`` and ``.hdf5`` select HDF5 automatically.
Use ``format='hdf5'`` for anonymous binary streams or other filenames;
``format='text'`` explicitly selects text. Open ``h5py.File`` objects are also
accepted. Supplied files/streams remain open. Internally opened files are
closed after use; when stopping chunked reading early, close the reader or
use a context manager::

    with plumed.read_as_pandas("HILLS.h5", chunksize=10000) as reader:
        for chunk in reader:
            print(chunk.plumed_constants)

Column selection, row skipping/limits, and named indexes work for both formats.
HDF5 row skipping counts data rows from zero across blocks. HDF5 chunks keep
block boundaries, including changes in constants. Unchunked reading combines
columns across blocks with missing values where necessary. Integer columns
with missing values use pandas nullable integers; fields that change between
integer and floating-point types use objects to preserve integer precision.
When constants change, use ``enable_constants='columns'`` to retain their values per row,
``enable_constants=False`` to discard them, or chunked reading for per-block
metadata. Default metadata mode rejects changing constants in one DataFrame.

The HDF5 layout is the native PLUMED table format, compatible with
``plumed h5convert``; pandas' ``read_hdf``/``to_hdf`` use a different layout.
Floating-point columns use double precision, integer columns use signed
64-bit integers, and string columns contain UTF-8 bytes. Named indexes become
ordinary fields. MultiIndexes, missing integer/string values, and unsupported
column types raise errors. Constants retain their text expressions, such as
``pi``, when writing a DataFrame that carries ``plumed_constants`` metadata.
Numeric HDF5 constants also retain their types when conversion is disabled.

The default ``enable_conversion=True`` evaluates constants with the selected
PLUMED kernel, just as for text tables: for example, a constant written as
``pi`` is returned with its numeric value. Its original expression is retained
in ``plumed_constants`` for writing it back. Ordinary HDF5 numeric columns are
read directly with their stored types; they do not need expression conversion.

Use ``enable_conversion=False`` to leave constants as stored and read without
loading a PLUMED kernel::

    data = plumed.read_as_pandas("COLVAR.h5", enable_conversion=False)

This option is useful for a wrapper-only installation without an available
kernel. It is optional for HDF5 reading; the default examples above use the
bundled or otherwise selected kernel. Writing HDF5 files needs no kernel.

PLUMED itself must be built with ``--enable-libhdf5`` to read/write these files
in native tools and simulations. Python HDF5 I/O uses ``h5py`` independently
of whether the selected PLUMED kernel was built with HDF5 support.

Data-file ownership
-------------------

``read_as_pandas`` and ``write_pandas`` close files they open, including gzip
streams and error paths. Streams supplied by the caller, including stdout,
remain open. A chunked reader keeps its internally opened file alive until
exhaustion, an error, or ``close()``. Use ``with plumed.read_as_pandas(...,
chunksize=...) as reader:`` when stopping iteration early.

InputBuilder metadata
---------------------

With kernel API version 11 or newer, ``InputBuilder`` reads metadata directly
from the selected kernel. It requires no executable, source checkout, or Vim
files. The kernel generates JSON from its live action registry, including
hidden keywords and the exact numbered-keyword classification. The same
serializer serves ``plumed gen_json``; its optional ``--actions`` argument
still accepts a legacy action-summary file.

Older kernels continue to use their installation's Vim dictionary and help
files. Keep those files installed when selecting an older kernel.

With ``InputBuilder(toplumed=p)``, metadata is read from that existing PLUMED
object, including actions it has already loaded. A successful ``ib.LOAD(...)``
refreshes the builder's methods and help automatically. If an action is loaded
directly with ``p.cmd(...)``, call ``ib.refresh()`` afterward. The supplied
object remains owned by the caller. Refresh is unavailable with older kernels
using static Vim files, or when ``load_dict=False``. For an attached older kernel,
select its Vim installation through ``kernel=...`` or ``PLUMED_ROOT`` as before.

Vim files
---------

``plumed gen_vim --outdir DIR`` generates both Vim syntax plugins and all
per-action help files in one invocation. It reads two static templates from
``<PLUMED_ROOT>/vim/templates``; use ``--template-dir DIR`` to select another
location. The templates are installed resources, while the action help comes
from the kernel's keyword registrations. An invalid template leaves existing
syntax and help files in place; successful regeneration removes stale help.

Maintainer release notes
------------------------

The ``python-package.yml`` workflow builds and validates the source distribution on pull requests and publishes it when a ``v*`` release tag is pushed. Before the first automated release, configure PyPI Trusted Publishing for the ``plumed`` project with GitHub owner ``plumed``, repository ``plumed2``, workflow ``python-package.yml``, and the GitHub environment ``pypi``. For a TestPyPI rehearsal, register the same workflow and repository with the ``testpypi`` environment on TestPyPI, then run the workflow manually with ``publish_to_testpypi`` enabled. No upload token is stored in GitHub Actions.

CHANGES: See the PLUMED documentation.

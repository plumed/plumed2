Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and make sure the file `libplumedKernel.so` (or `libplumedKernel.dylib`) is available on your system.

You should then make sure the library is found setting the environment variable `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple plumed versions on your system using tcl environment modules, this should be taken care automatically
by the plumed module.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This will allow you to mix different plumed versions in the same python script.

CHANGES: See the PLUMED documentation.

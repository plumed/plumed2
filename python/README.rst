Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and set the environment variable
`PLUMED_KERNEL` to point to the file `libplumedKernel.so` (or `libplumedKernel.dylib`). Something like this should work::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

CHANGES: See the PLUMED documentation.

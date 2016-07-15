\defgroup internal-blas Internal BLAS

Internal implementation of BLAS, imported from GROMACS.

The module in src/blas contains an internal implementation
of BLAS routines which is automatically imported from GROMACS
using the src/blas/import.sh script. This set of routines
is compiled when __PLUMED_HAS_EXTERNAL_BLAS is not defined
and allow PLUMED to be used when installed BLAS libraries
are not available. Notice that the import script
creates a blas.h file with function declarations which
are used also when installed blas are employed. This is
done because there are blas installation written in FORTRAN
that do not provide header files.

Since files are automatically generated, do not edit them directly.
In case you need PLUMED specific modifications
please do it by modifying the import script.

Within the PLUMED doxygen (this page) the available
macros are listed but not documented. Have a look
at the corresponding documentation at http://www.netlib.org/blas


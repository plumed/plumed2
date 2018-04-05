/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib BLAS library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference BLAS implementation.

The reference BLAS implementation is available from http://www.netlib.org/blas 

BLAS does not come with a formal named "license", but a general statement that 

"The reference BLAS is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified BLAS files as the original netlib versions, so do what you want with them.
However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_blas_real_h
#define __PLUMED_blas_real_h
    /*! \brief Double precision accuracy */
#define PLUMED_GMX_DOUBLE_EPS   2.2204460492503131e-16

    /*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_DOUBLE_MAX   1.7976931348623157e+308

    /*! \brief Minimum double precision value */
#define PLUMED_GMX_DOUBLE_MIN   2.2250738585072014e-308

    /*! \brief Single precision accuracy */
#define PLUMED_GMX_FLOAT_EPS    1.19209290e-07F

    /*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_FLOAT_MAX    3.40282346E+38F

    /*! \brief Minimum single precision value */
#define PLUMED_GMX_FLOAT_MIN    1.175494351E-38F

#if defined(F77_NO_UNDERSCORE)
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower
#else
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower ## _
#endif

#endif

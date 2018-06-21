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
#ifndef __PLUMED_blas_def_external_h
#define __PLUMED_blas_def_external_h
#define plumed_blas_dasum PLUMED_BLAS_F77_FUNC(dasum,DASUM)
#define plumed_blas_daxpy PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)
#define plumed_blas_dcopy PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)
#define plumed_blas_ddot PLUMED_BLAS_F77_FUNC(ddot,DDOT)
#define plumed_blas_dgemm PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)
#define plumed_blas_dgemv PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)
#define plumed_blas_dger PLUMED_BLAS_F77_FUNC(dger,DGER)
#define plumed_blas_dnrm2 PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)
#define plumed_blas_drot PLUMED_BLAS_F77_FUNC(drot,DROT)
#define plumed_blas_dscal PLUMED_BLAS_F77_FUNC(dscal,DSCAL)
#define plumed_blas_dswap PLUMED_BLAS_F77_FUNC(dswap,DSWAP)
#define plumed_blas_dsymv PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)
#define plumed_blas_dsyr2 PLUMED_BLAS_F77_FUNC(dsyr2,DSYR2)
#define plumed_blas_dsyr2k PLUMED_BLAS_F77_FUNC(dsyr2k,DSYR2K)
#define plumed_blas_dtrmm PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)
#define plumed_blas_dtrmv PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)
#define plumed_blas_dtrsm PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)
#define plumed_blas_idamax PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)
#define plumed_blas_sasum PLUMED_BLAS_F77_FUNC(sasum,SASUM)
#define plumed_blas_saxpy PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)
#define plumed_blas_scopy PLUMED_BLAS_F77_FUNC(scopy,SCOPY)
#define plumed_blas_sdot PLUMED_BLAS_F77_FUNC(sdot,SDOT)
#define plumed_blas_sgemm PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)
#define plumed_blas_sgemv PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)
#define plumed_blas_sger PLUMED_BLAS_F77_FUNC(sger,SGER)
#define plumed_blas_snrm2 PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)
#define plumed_blas_srot PLUMED_BLAS_F77_FUNC(srot,SROT)
#define plumed_blas_sscal PLUMED_BLAS_F77_FUNC(sscal,SSCAL)
#define plumed_blas_sswap PLUMED_BLAS_F77_FUNC(sswap,SSWAP)
#define plumed_blas_ssymv PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)
#define plumed_blas_ssyr2 PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)
#define plumed_blas_ssyr2k PLUMED_BLAS_F77_FUNC(ssyr2k,SSYR2K)
#define plumed_blas_strmm PLUMED_BLAS_F77_FUNC(strmm,STRMM)
#define plumed_blas_strmv PLUMED_BLAS_F77_FUNC(strmv,STRMV)
#define plumed_blas_strsm PLUMED_BLAS_F77_FUNC(strsm,STRSM)
#define plumed_blas_isamax PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)
#endif

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
#ifndef __PLUMED_blas_def_internal_h
#define __PLUMED_blas_def_internal_h
/** \ingroup internal-blas */
#define plumed_blas_dasum PLMD::blas::PLUMED_BLAS_F77_FUNC(dasum,DASUM)
/** \ingroup internal-blas */
#define plumed_blas_daxpy PLMD::blas::PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)
/** \ingroup internal-blas */
#define plumed_blas_dcopy PLMD::blas::PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)
/** \ingroup internal-blas */
#define plumed_blas_ddot PLMD::blas::PLUMED_BLAS_F77_FUNC(ddot,DDOT)
/** \ingroup internal-blas */
#define plumed_blas_dgemm PLMD::blas::PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)
/** \ingroup internal-blas */
#define plumed_blas_dgemv PLMD::blas::PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)
/** \ingroup internal-blas */
#define plumed_blas_dger PLMD::blas::PLUMED_BLAS_F77_FUNC(dger,DGER)
/** \ingroup internal-blas */
#define plumed_blas_dnrm2 PLMD::blas::PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)
/** \ingroup internal-blas */
#define plumed_blas_drot PLMD::blas::PLUMED_BLAS_F77_FUNC(drot,DROT)
/** \ingroup internal-blas */
#define plumed_blas_dscal PLMD::blas::PLUMED_BLAS_F77_FUNC(dscal,DSCAL)
/** \ingroup internal-blas */
#define plumed_blas_dswap PLMD::blas::PLUMED_BLAS_F77_FUNC(dswap,DSWAP)
/** \ingroup internal-blas */
#define plumed_blas_dsymv PLMD::blas::PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)
/** \ingroup internal-blas */
#define plumed_blas_dsyr2 PLMD::blas::PLUMED_BLAS_F77_FUNC(dsyr2,DSYR2)
/** \ingroup internal-blas */
#define plumed_blas_dsyr2k PLMD::blas::PLUMED_BLAS_F77_FUNC(dsyr2k,DSYR2K)
/** \ingroup internal-blas */
#define plumed_blas_dtrmm PLMD::blas::PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)
/** \ingroup internal-blas */
#define plumed_blas_dtrmv PLMD::blas::PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)
/** \ingroup internal-blas */
#define plumed_blas_dtrsm PLMD::blas::PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)
/** \ingroup internal-blas */
#define plumed_blas_idamax PLMD::blas::PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)
/** \ingroup internal-blas */
#define plumed_blas_sasum PLMD::blas::PLUMED_BLAS_F77_FUNC(sasum,SASUM)
/** \ingroup internal-blas */
#define plumed_blas_saxpy PLMD::blas::PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)
/** \ingroup internal-blas */
#define plumed_blas_scopy PLMD::blas::PLUMED_BLAS_F77_FUNC(scopy,SCOPY)
/** \ingroup internal-blas */
#define plumed_blas_sdot PLMD::blas::PLUMED_BLAS_F77_FUNC(sdot,SDOT)
/** \ingroup internal-blas */
#define plumed_blas_sgemm PLMD::blas::PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)
/** \ingroup internal-blas */
#define plumed_blas_sgemv PLMD::blas::PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)
/** \ingroup internal-blas */
#define plumed_blas_sger PLMD::blas::PLUMED_BLAS_F77_FUNC(sger,SGER)
/** \ingroup internal-blas */
#define plumed_blas_snrm2 PLMD::blas::PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)
/** \ingroup internal-blas */
#define plumed_blas_srot PLMD::blas::PLUMED_BLAS_F77_FUNC(srot,SROT)
/** \ingroup internal-blas */
#define plumed_blas_sscal PLMD::blas::PLUMED_BLAS_F77_FUNC(sscal,SSCAL)
/** \ingroup internal-blas */
#define plumed_blas_sswap PLMD::blas::PLUMED_BLAS_F77_FUNC(sswap,SSWAP)
/** \ingroup internal-blas */
#define plumed_blas_ssymv PLMD::blas::PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)
/** \ingroup internal-blas */
#define plumed_blas_ssyr2 PLMD::blas::PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)
/** \ingroup internal-blas */
#define plumed_blas_ssyr2k PLMD::blas::PLUMED_BLAS_F77_FUNC(ssyr2k,SSYR2K)
/** \ingroup internal-blas */
#define plumed_blas_strmm PLMD::blas::PLUMED_BLAS_F77_FUNC(strmm,STRMM)
/** \ingroup internal-blas */
#define plumed_blas_strmv PLMD::blas::PLUMED_BLAS_F77_FUNC(strmv,STRMV)
/** \ingroup internal-blas */
#define plumed_blas_strsm PLMD::blas::PLUMED_BLAS_F77_FUNC(strsm,STRSM)
/** \ingroup internal-blas */
#define plumed_blas_isamax PLMD::blas::PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)
#endif

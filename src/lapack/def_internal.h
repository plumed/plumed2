/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib LAPACK library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference LAPACK implementation.

The reference LAPACK implementation is available from http://www.netlib.org/lapack 

LAPACK does not come with a formal named "license", but a general statement saying:

"The reference LAPACK is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified LAPACK files as the original netlib versions, so do what you want with them.

However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_lapack_def_internal_h
#define __PLUMED_lapack_def_internal_h
/** \ingroup internal-lapack */
#define plumed_lapack_dbdsdc PLMD::lapack::PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)
/** \ingroup internal-lapack */
#define plumed_lapack_dgetf2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgetf2,DGETF2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlamrg PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarnv PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarnv,DLARNV)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd0 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd0,DLASD0)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasda PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasda,DLASDA)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq6,DLASQ6)
/** \ingroup internal-lapack */
#define plumed_lapack_dorgl2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorgl2,DORGL2)
/** \ingroup internal-lapack */
#define plumed_lapack_dbdsqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dbdsqr,DBDSQR)
/** \ingroup internal-lapack */
#define plumed_lapack_dgetrf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgetrf,DGETRF)
/** \ingroup internal-lapack */
#define plumed_lapack_dgetri PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgetri,DGETRI)
/** \ingroup internal-lapack */
#define plumed_lapack_dgetrs PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgetrs,DGETRS)
/** \ingroup internal-lapack */
#define plumed_lapack_dtrtri PLMD::lapack::PLUMED_BLAS_F77_FUNC(dtrtri,DTRTRI)
/** \ingroup internal-lapack */
#define plumed_lapack_dtrti2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dtrti2,DTRTI2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlange PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlange,DLANGE)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarrbx PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarrbx,DLARRBX)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd1 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd1,DLASD1)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasdq PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasr,DLASR)
/** \ingroup internal-lapack */
#define plumed_lapack_dorglq PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorglq,DORGLQ)
/** \ingroup internal-lapack */
#define plumed_lapack_dormtr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dormtr,DORMTR)
/** \ingroup internal-lapack */
#define plumed_lapack_dgebd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgebd2,DGEBD2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlabrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlabrd,DLABRD)
/** \ingroup internal-lapack */
#define plumed_lapack_dlanst PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlanst,DLANST)
/** \ingroup internal-lapack */
#define plumed_lapack_dlansy PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlansy,DLANSY)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarrex PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarrex,DLARREX)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd2,DLASD2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasdt PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasdt,DLASDT)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasrt PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasrt2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasrt2,DLASRT2)
/** \ingroup internal-lapack */
#define plumed_lapack_ilasrt2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(ilasrt2,ILASRT2)
/** \ingroup internal-lapack */
#define plumed_lapack_dorgqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorgqr,DORGQR)
/** \ingroup internal-lapack */
#define plumed_lapack_dstebz PLMD::lapack::PLUMED_BLAS_F77_FUNC(dstebz,DSTEBZ)
/** \ingroup internal-lapack */
#define plumed_lapack_dsteqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dsteqr,DSTEQR)
/** \ingroup internal-lapack */
#define plumed_lapack_dgebrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)
/** \ingroup internal-lapack */
#define plumed_lapack_dlacpy PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)
/** \ingroup internal-lapack */
#define plumed_lapack_dlapy2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarrfx PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarrfx,DLARRFX)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd3 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd3,DLASD3)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaset PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaset,DLASET)
/** \ingroup internal-lapack */
#define plumed_lapack_dlassq PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)
/** \ingroup internal-lapack */
#define plumed_lapack_dorm2l PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorm2l,DORM2L)
/** \ingroup internal-lapack */
#define plumed_lapack_dstegr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dstegr,DSTEGR)
/** \ingroup internal-lapack */
#define plumed_lapack_ssteqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(ssteqr,SSTEQR)
/** \ingroup internal-lapack */
#define plumed_lapack_dgelq2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgelq2,DGELQ2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlae2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaev2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaev2,DLAEV2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlar1vx PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlar1vx,DLAR1VX)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarrvx PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarrvx,DLARRVX)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd4 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd4,DLASD4)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq1 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq1,DLASQ1)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasv2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasv2,DLASV2)
/** \ingroup internal-lapack */
#define plumed_lapack_dorm2r PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorm2r,DORM2R)
/** \ingroup internal-lapack */
#define plumed_lapack_dstein PLMD::lapack::PLUMED_BLAS_F77_FUNC(dstein,DSTEIN)
/** \ingroup internal-lapack */
#define plumed_lapack_dgelqf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgelqf,DGELQF)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaebz PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaebz,DLAEBZ)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarf,DLARF)
/** \ingroup internal-lapack */
#define plumed_lapack_dlartg PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd5 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd5,DLASD5)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq2,DLASQ2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq3 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq3,DLASQ3)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaswp PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)
/** \ingroup internal-lapack */
#define plumed_lapack_dormbr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)
/** \ingroup internal-lapack */
#define plumed_lapack_dsterf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dsterf,DSTERF)
/** \ingroup internal-lapack */
#define plumed_lapack_dgeqr2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgeqr2,DGEQR2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaed6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaed6,DLAED6)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarfb PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)
/** \ingroup internal-lapack */
#define plumed_lapack_dlaruv PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlaruv,DLARUV)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd6,DLASD6)
/** \ingroup internal-lapack */
#define plumed_lapack_dlatrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlatrd,DLATRD)
/** \ingroup internal-lapack */
#define plumed_lapack_dorml2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorml2,DORML2)
/** \ingroup internal-lapack */
#define plumed_lapack_dstevr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dstevr,DSTEVR)
/** \ingroup internal-lapack */
#define plumed_lapack_dsytrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(dsytrd,DSYTRD)
/** \ingroup internal-lapack */
#define plumed_lapack_dsyevr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dsyevr,DSYEVR)
/** \ingroup internal-lapack */
#define plumed_lapack_dormql PLMD::lapack::PLUMED_BLAS_F77_FUNC(dormql,DORMQL)
/** \ingroup internal-lapack */
#define plumed_lapack_dormqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dormqr,DORMQR)
/** \ingroup internal-lapack */
#define plumed_lapack_dorgbr PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorgbr,DORGBR)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq5 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq5,DLASQ5)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd8 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd8,DLASD8)
/** \ingroup internal-lapack */
#define plumed_lapack_dlascl PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarft PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)
/** \ingroup internal-lapack */
#define plumed_lapack_dlagts PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlagts,DLAGTS)
/** \ingroup internal-lapack */
#define plumed_lapack_dgesdd PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgesdd,DGESDD)
/** \ingroup internal-lapack */
#define plumed_lapack_dsytd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dsytd2,DSYTD2)
/** \ingroup internal-lapack */
#define plumed_lapack_dormlq PLMD::lapack::PLUMED_BLAS_F77_FUNC(dormlq,DORMLQ)
/** \ingroup internal-lapack */
#define plumed_lapack_dorg2r PLMD::lapack::PLUMED_BLAS_F77_FUNC(dorg2r,DORG2R)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasq4 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasq4,DLASQ4)
/** \ingroup internal-lapack */
#define plumed_lapack_dlasd7 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlasd7,DLASD7)
/** \ingroup internal-lapack */
#define plumed_lapack_dlas2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlas2,DLAS2)
/** \ingroup internal-lapack */
#define plumed_lapack_dlarfg PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)
/** \ingroup internal-lapack */
#define plumed_lapack_dlagtf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dlagtf,DLAGTF)
/** \ingroup internal-lapack */
#define plumed_lapack_dgeqrf PLMD::lapack::PLUMED_BLAS_F77_FUNC(dgeqrf,DGEQRF)
/** \ingroup internal-lapack */
#define plumed_lapack_sbdsdc PLMD::lapack::PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)
/** \ingroup internal-lapack */
#define plumed_lapack_sgetf2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)
/** \ingroup internal-lapack */
#define plumed_lapack_slamrg PLMD::lapack::PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)
/** \ingroup internal-lapack */
#define plumed_lapack_slarnv PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarnv,SLARNV)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd0 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd0,SLASD0)
/** \ingroup internal-lapack */
#define plumed_lapack_slasda PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasda,SLASDA)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq6,SLASQ6)
/** \ingroup internal-lapack */
#define plumed_lapack_sorgl2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorgl2,SORGL2)
/** \ingroup internal-lapack */
#define plumed_lapack_sbdsqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sbdsqr,SBDSQR)
/** \ingroup internal-lapack */
#define plumed_lapack_sgetrf PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgetrf,SGETRF)
/** \ingroup internal-lapack */
#define plumed_lapack_sgetri PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgetri,SGETRI)
/** \ingroup internal-lapack */
#define plumed_lapack_sgetrs PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgetrs,SGETRS)
/** \ingroup internal-lapack */
#define plumed_lapack_strtri PLMD::lapack::PLUMED_BLAS_F77_FUNC(strtri,STRTRI)
/** \ingroup internal-lapack */
#define plumed_lapack_strti2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(strti2,STRTI2)
/** \ingroup internal-lapack */
#define plumed_lapack_slange PLMD::lapack::PLUMED_BLAS_F77_FUNC(slange,SLANGE)
/** \ingroup internal-lapack */
#define plumed_lapack_slarrbx PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarrbx,SLARRBX)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd1 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd1,SLASD1)
/** \ingroup internal-lapack */
#define plumed_lapack_slasdq PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)
/** \ingroup internal-lapack */
#define plumed_lapack_slasr PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasr,SLASR)
/** \ingroup internal-lapack */
#define plumed_lapack_sorglq PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorglq,SORGLQ)
/** \ingroup internal-lapack */
#define plumed_lapack_sormtr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sormtr,SORMTR)
/** \ingroup internal-lapack */
#define plumed_lapack_sgebd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgebd2,SGEBD2)
/** \ingroup internal-lapack */
#define plumed_lapack_slabrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(slabrd,SLABRD)
/** \ingroup internal-lapack */
#define plumed_lapack_slanst PLMD::lapack::PLUMED_BLAS_F77_FUNC(slanst,SLANST)
/** \ingroup internal-lapack */
#define plumed_lapack_slansy PLMD::lapack::PLUMED_BLAS_F77_FUNC(slansy,SLANSY)
/** \ingroup internal-lapack */
#define plumed_lapack_slarrex PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarrex,SLARREX)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd2,SLASD2)
/** \ingroup internal-lapack */
#define plumed_lapack_slasdt PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasdt,SLASDT)
/** \ingroup internal-lapack */
#define plumed_lapack_slasrt PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)
/** \ingroup internal-lapack */
#define plumed_lapack_slasrt2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasrt2,SLASRT2)
/** \ingroup internal-lapack */
#define plumed_lapack_sorgqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorgqr,SORGQR)
/** \ingroup internal-lapack */
#define plumed_lapack_sstebz PLMD::lapack::PLUMED_BLAS_F77_FUNC(sstebz,SSTEBZ)
/** \ingroup internal-lapack */
#define plumed_lapack_sgebrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)
/** \ingroup internal-lapack */
#define plumed_lapack_slacpy PLMD::lapack::PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)
/** \ingroup internal-lapack */
#define plumed_lapack_slapy2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)
/** \ingroup internal-lapack */
#define plumed_lapack_slarrfx PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarrfx,SLARRFX)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd3 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd3,SLASD3)
/** \ingroup internal-lapack */
#define plumed_lapack_slaset PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaset,SLASET)
/** \ingroup internal-lapack */
#define plumed_lapack_slassq PLMD::lapack::PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)
/** \ingroup internal-lapack */
#define plumed_lapack_sorm2l PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorm2l,SORM2L)
/** \ingroup internal-lapack */
#define plumed_lapack_sstegr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sstegr,SSTEGR)
/** \ingroup internal-lapack */
#define plumed_lapack_sgelq2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgelq2,SGELQ2)
/** \ingroup internal-lapack */
#define plumed_lapack_slae2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slae2,SLAE2)
/** \ingroup internal-lapack */
#define plumed_lapack_slaev2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaev2,SLAEV2)
/** \ingroup internal-lapack */
#define plumed_lapack_slar1vx PLMD::lapack::PLUMED_BLAS_F77_FUNC(slar1vx,SLAR1VX)
/** \ingroup internal-lapack */
#define plumed_lapack_slarrvx PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarrvx,SLARRVX)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd4 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd4,SLASD4)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq1 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq1,SLASQ1)
/** \ingroup internal-lapack */
#define plumed_lapack_slasv2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasv2,SLASV2)
/** \ingroup internal-lapack */
#define plumed_lapack_sorm2r PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorm2r,SORM2R)
/** \ingroup internal-lapack */
#define plumed_lapack_sstein PLMD::lapack::PLUMED_BLAS_F77_FUNC(sstein,SSTEIN)
/** \ingroup internal-lapack */
#define plumed_lapack_sgelqf PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgelqf,SGELQF)
/** \ingroup internal-lapack */
#define plumed_lapack_slaebz PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaebz,SLAEBZ)
/** \ingroup internal-lapack */
#define plumed_lapack_slarf PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarf,SLARF)
/** \ingroup internal-lapack */
#define plumed_lapack_slartg PLMD::lapack::PLUMED_BLAS_F77_FUNC(slartg,SLARTG)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd5 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd5,SLASD5)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq2,SLASQ2)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq3 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq3,SLASQ3)
/** \ingroup internal-lapack */
#define plumed_lapack_slaswp PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)
/** \ingroup internal-lapack */
#define plumed_lapack_sormbr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)
/** \ingroup internal-lapack */
#define plumed_lapack_ssterf PLMD::lapack::PLUMED_BLAS_F77_FUNC(ssterf,SSTERF)
/** \ingroup internal-lapack */
#define plumed_lapack_sgeqr2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgeqr2,SGEQR2)
/** \ingroup internal-lapack */
#define plumed_lapack_slaed6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaed6,SLAED6)
/** \ingroup internal-lapack */
#define plumed_lapack_slarfb PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)
/** \ingroup internal-lapack */
#define plumed_lapack_slaruv PLMD::lapack::PLUMED_BLAS_F77_FUNC(slaruv,SLARUV)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd6 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd6,SLASD6)
/** \ingroup internal-lapack */
#define plumed_lapack_slatrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(slatrd,SLATRD)
/** \ingroup internal-lapack */
#define plumed_lapack_sorml2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorml2,SORML2)
/** \ingroup internal-lapack */
#define plumed_lapack_sstevr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sstevr,SSTEVR)
/** \ingroup internal-lapack */
#define plumed_lapack_ssytrd PLMD::lapack::PLUMED_BLAS_F77_FUNC(ssytrd,SSYTRD)
/** \ingroup internal-lapack */
#define plumed_lapack_ssyevr PLMD::lapack::PLUMED_BLAS_F77_FUNC(ssyevr,SSYEVR)
/** \ingroup internal-lapack */
#define plumed_lapack_sormql PLMD::lapack::PLUMED_BLAS_F77_FUNC(sormql,SORMQL)
/** \ingroup internal-lapack */
#define plumed_lapack_sormqr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sormqr,SORMQR)
/** \ingroup internal-lapack */
#define plumed_lapack_sorgbr PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorgbr,SORGBR)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq5 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq5,SLASQ5)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd8 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd8,SLASD8)
/** \ingroup internal-lapack */
#define plumed_lapack_slascl PLMD::lapack::PLUMED_BLAS_F77_FUNC(slascl,SLASCL)
/** \ingroup internal-lapack */
#define plumed_lapack_slarft PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarft,SLARFT)
/** \ingroup internal-lapack */
#define plumed_lapack_slagts PLMD::lapack::PLUMED_BLAS_F77_FUNC(slagts,SLAGTS)
/** \ingroup internal-lapack */
#define plumed_lapack_sgesdd PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgesdd,SGESDD)
/** \ingroup internal-lapack */
#define plumed_lapack_ssytd2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(ssytd2,SSYTD2)
/** \ingroup internal-lapack */
#define plumed_lapack_sormlq PLMD::lapack::PLUMED_BLAS_F77_FUNC(sormlq,SORMLQ)
/** \ingroup internal-lapack */
#define plumed_lapack_sorg2r PLMD::lapack::PLUMED_BLAS_F77_FUNC(sorg2r,SORG2R)
/** \ingroup internal-lapack */
#define plumed_lapack_slasq4 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasq4,SLASQ4)
/** \ingroup internal-lapack */
#define plumed_lapack_slasd7 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slasd7,SLASD7)
/** \ingroup internal-lapack */
#define plumed_lapack_slas2 PLMD::lapack::PLUMED_BLAS_F77_FUNC(slas2,SLAS2)
/** \ingroup internal-lapack */
#define plumed_lapack_slarfg PLMD::lapack::PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)
/** \ingroup internal-lapack */
#define plumed_lapack_slagtf PLMD::lapack::PLUMED_BLAS_F77_FUNC(slagtf,SLAGTF)
/** \ingroup internal-lapack */
#define plumed_lapack_sgeqrf PLMD::lapack::PLUMED_BLAS_F77_FUNC(sgeqrf,SGEQRF)
#endif

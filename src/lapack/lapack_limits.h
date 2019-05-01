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
#ifndef __PLUMED_lapack_lapack_limits_h
#define __PLUMED_lapack_lapack_limits_h

#define DSTEBZ_BLOCKSIZE  1

#define DORGBR_BLOCKSIZE    32
#define DORGBR_MINBLOCKSIZE 2
#define DORGBR_CROSSOVER    128

#define DGEQRF_BLOCKSIZE    32
#define DGEQRF_MINBLOCKSIZE 2  
#define DGEQRF_CROSSOVER    128

#define DORGQR_BLOCKSIZE    32
#define DORGQR_MINBLOCKSIZE 2
#define DORGQR_CROSSOVER    128

#define DORMLQ_BLOCKSIZE    32
#define DORMLQ_MINBLOCKSIZE 2
#define DORMLQ_CROSSOVER    128

#define DORMQL_BLOCKSIZE    32
#define DORMQL_MINBLOCKSIZE 2
#define DORMQL_CROSSOVER    128

#define DSYTRD_BLOCKSIZE    32
#define DSYTRD_MINBLOCKSIZE 2
#define DSYTRD_CROSSOVER    128

#define DGEBRD_BLOCKSIZE    32
#define DGEBRD_MINBLOCKSIZE 2
#define DGEBRD_CROSSOVER    128

#define DORMQR_BLOCKSIZE    32
#define DORMQR_MINBLOCKSIZE 2
#define DORMQR_CROSSOVER    128

#define DGELQF_BLOCKSIZE    32
#define DGELQF_MINBLOCKSIZE 2  
#define DGELQF_CROSSOVER    128

#define DGETRF_BLOCKSIZE    64
#define DGETRF_MINBLOCKSIZE 2

#define DGETRI_BLOCKSIZE    64
#define DGETRI_MINBLOCKSIZE 2

#define DTRTRI_BLOCKSIZE    64

#define DBDSDC_SMALLSIZE 25

#define DBDSQR_MAXITR 6



#endif /* _LAPACK_LIMITS_H_ */

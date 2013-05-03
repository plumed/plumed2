/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib LAPACK library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the Gromacs developers, and
not the reference LAPACK implementation.

The reference LAPACK implementation is available from http://www.netlib.org/lapack 

LAPACK does not come with a formal named "license", but a general statement saying:

"The reference LAPACK is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of Gromacs is GPL, we think it's only fair to give you the same rights to
our modified LAPACK files as the original netlib versions, so do what you want with them.

However, be warned that we have only tested that they to the right thing in the cases used
in Gromacs (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgetrf,SGETRF)(int *m,
	int *n,
	float *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int mindim,jb;
  int i,j,k,l;
  int iinfo;
  float minusone = -1.0;
  float one = 1.0;

  if(*m<=0 || *n<=0)
    return;

  *info = 0;

  mindim = (*m < *n) ? *m : *n;

  if(DGETRF_BLOCKSIZE>=mindim) {

    /* unblocked code */
    PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)(m,n,a,lda,ipiv,info);

  } else {

    /* blocked case */

    for(j=1;j<=mindim;j+=DGETRF_BLOCKSIZE) {
      jb = ( DGETRF_BLOCKSIZE < (mindim-j+1)) ? DGETRF_BLOCKSIZE : (mindim-j+1);
      /* factor diag. and subdiag blocks and test for singularity */
      k = *m-j+1;
      PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)(&k,&jb,&(a[(j-1)*(*lda)+(j-1)]),lda,&(ipiv[j-1]),&iinfo);
      
      if(*info==0 && iinfo>0)
	*info = iinfo + j - 1;

      /* adjust pivot indices */
      k = (*m < (j+jb-1)) ? *m : (j+jb-1);
      for(i=j;i<=k;i++)
	ipiv[i-1] += j - 1;

      /* Apply to columns 1 throughj j-1 */
      k = j - 1;
      i = j + jb - 1;
      l = 1;
      PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(&k,a,lda,&j,&i,ipiv,&l);
      if((j+jb)<=*n) {
	/* Apply to cols. j+jb through n */
	k = *n-j-jb+1;
	i = j+jb-1;
	l = 1;
	PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(&k,&(a[(j+jb-1)*(*lda)+0]),lda,&j,&i,ipiv,&l);
	/* Compute block row of U */
	k = *n-j-jb+1;
	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left","Lower","No transpose","Unit",&jb,&k,&one,
	       &(a[(j-1)*(*lda)+(j-1)]),lda,&(a[(j+jb-1)*(*lda)+(j-1)]),lda);

	if((j+jb)<=*m) {
	  /* Update trailing submatrix */
	  k = *m-j-jb+1;
	  i = *n-j-jb+1;
	  PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose","No transpose",&k,&i,&jb,&minusone,
		 &(a[(j-1)*(*lda)+(j+jb-1)]),lda,
		 &(a[(j+jb-1)*(*lda)+(j-1)]),lda,&one,
		 &(a[(j+jb-1)*(*lda)+(j+jb-1)]),lda);
	}

      }
    }
  }
}
}
}

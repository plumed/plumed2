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
#include <ctype.h>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slatrd,SLATRD)(const char *  uplo,
       int  *   n,
       int  *   nb,
       float * a,
       int *    lda,
       float * e,
       float * tau,
       float * w,
       int *    ldw)
{
  int i,iw;
  int ti1,ti2,ti3;
  float one,zero,minusone,alpha;
  const char ch=toupper(*uplo);

  one=1.0;
  minusone=-1.0;
  zero=0.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n;i>=(*n-*nb+1);i--) {
      iw = i -*n + *nb;
      
      if(i<*n) {
	ti1 = *n-i;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(a[ i*(*lda) + 0]),lda,&(w[iw*(*ldw)+(i-1)]),
	       ldw,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(w[ iw*(*ldw) + 0]),ldw,&(a[i*(*lda)+(i-1)]),
	       lda,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
      }

      if(i>1) {
	/*  Generate elementary reflector H(i) to annihilate
	 *              A(1:i-2,i) 
	 */
	ti1 = i-1;
	ti2 = 1;

	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i-2)]),&(a[(i-1)*(*lda)+0]),&ti2,&(tau[i-2]));
      
	e[i-2] = a[(i-1)*(*lda)+(i-2)];
	a[(i-1)*(*lda)+(i-2)] = 1.0;

	/* Compute W(1:i-1,i) */
	ti1 = i-1;
	ti2 = 1;

	/* BLAS */
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("U",&ti1,&one,a,lda,&(a[(i-1)*(*lda)+0]),&ti2,&zero,
	       &(w[(iw-1)*(*ldw)+0]),&ti2);
	if(i<*n) {
	  ti1 = i-1;
	  ti2 = *n-i;
	  ti3 = 1;
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[iw*(*ldw)+0]),ldw,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[i*(*lda)+0]),lda,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[i*(*lda)+0]),lda,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[iw*(*ldw)+0]),ldw,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	}
      
	ti1 = i-1;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-2]),&(w[(iw-1)*(*ldw)+0]),&ti2);
      
	alpha = -0.5*tau[i-2]*PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(w[(iw-1)*(*ldw)+0]),&ti2,
				    &(a[(i-1)*(*lda)+0]),&ti2);
      
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+0]),&ti2,&(w[(iw-1)*(*ldw)+0]),&ti2);

      }
    }
  } else {
    /* lower */
    for(i=1;i<=*nb;i++) {

      ti1 = *n-i+1;
      ti2 = i-1;
      ti3 = 1;
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(a[ i-1 ]),lda,&(w[ i-1 ]),
	       ldw,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(w[ i-1 ]),ldw,&(a[ i-1 ]),
	       lda,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);

      if(i<*n) {
	ti1 = *n - i;
	ti2 = (*n < i+2 ) ? *n : (i+2);
	ti3 = 1;
	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+(ti2-1)]),&ti3,&(tau[i-1]));
	e[i-1] = a[(i-1)*(*lda)+(i)];
	a[(i-1)*(*lda)+(i)] = 1.0;
	
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("L",&ti1,&one,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),&ti2,
	       &zero,&(w[(i-1)*(*ldw)+i]),&ti2);
	ti1 = *n - i;
	ti2 = i-1;
	ti3 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[ i ]),ldw,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[ i ]),lda,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[ i ]),lda,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[ i ]),ldw,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);

	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-1]),&(w[(i-1)*(*ldw)+i]),&ti3);
	alpha = -0.5*tau[i-1]*PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(w[(i-1)*(*ldw)+i]),&ti3,
				   &(a[(i-1)*(*lda)+i]),&ti3);
	
	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti3,&(w[(i-1)*(*ldw)+i]),&ti3);
      }
    }
  }
  return;
}
	


  
}
}

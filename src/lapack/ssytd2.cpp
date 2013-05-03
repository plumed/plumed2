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
#include <math.h>

#include "simple.h"

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssytd2,SSYTD2)(const char *    uplo,
	int *     n,
	float *  a,
	int *     lda,
	float *  d,
	float *  e,
	float *  tau,
	int *     info)
{
  float minusone,zero;
  float taui,alpha,tmp;
  int ti1,ti2,ti3,i;
  const char ch=toupper(*uplo);

  zero = 0.0;
  minusone = -1.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n-1;i>=1;i--) {

      ti1 = 1;
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i,&(a[i*(*lda)+(i-1)]),&(a[i*(*lda)+0]),&ti1,&taui);
      e[i-1] = a[i*(*lda) + (i-1)];
      if(fabs(taui)>PLUMED_GMX_FLOAT_MIN) {
	a[i*(*lda)+(i-1)] = 1.0;
      
	ti1 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("U",&i,&taui,a,lda,&(a[i*(*lda)+0]),&ti1,&zero,tau,&ti1);

	tmp = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&i,tau,&ti1,&(a[i*(*lda)+0]),&ti1);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&i,&alpha,&(a[i*(*lda)+0]),&ti1,tau,&ti1);

	PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)("U",&i,&minusone,&(a[i*(*lda)+0]),&ti1,tau,&ti1,a,lda);

	a[i*(*lda)+(i-1)] = e[i-1]; 

      }
      d[i] = a[i*(*lda)+i];
      tau[i-1] = taui;
    }
    d[0] = a[0];
    
  } else {
    /* lower */

    for(i=1;i<*n;i++) {

      ti1 = *n - i;
      ti2 = ( *n < i+2) ? *n : i+2;
      ti3 = 1;
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+ti2-1]),&ti3,&taui);

      e[i-1] = a[(i-1)*(*lda) + (i)];

      if(fabs(taui)>PLUMED_GMX_FLOAT_MIN) {
	a[(i-1)*(*lda)+(i)] = 1.0;
      
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)(uplo,&ti1,&taui,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),
	       &ti2,&zero,&(tau[i-1]),&ti2);
	
	tmp = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(tau[i-1]),&ti2,&(a[(i-1)*(*lda)+i]),&ti2);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2);

	PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)(uplo,&ti1,&minusone,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2,
	       &(a[(i)*(*lda)+i]),lda);

	a[(i-1)*(*lda)+(i)] = e[i-1]; 

      }
      d[i-1] = a[(i-1)*(*lda)+i-1];
      tau[i-1] = taui;
    }
    d[*n-1] = a[(*n-1)*(*lda)+(*n-1)];
 
  }
  return;
}
}
}

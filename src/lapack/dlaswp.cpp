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

While the rest of Gromacs is LGPL, we think it's only fair to give you the same rights to
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
#include "lapack.h"

/* LAPACK */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(int *n,
	double *a,
	int *lda,
	int *k1,
	int *k2,
	int *ipiv,
	int *incx)
{
  int ix0,i1,i2,inc,n32;
  int ix,i,j,ip,k;
  double temp;

  if(*incx>0) {
    ix0 = *k1 - 1;
    i1 = *k1 - 1;
    i2 = *k2;
    inc = 1;
  } else if(*incx<0) {
    ix0 = *incx * (1- *k2);
    i1 = *k2 - 1;
    i2 = *k1;
    inc = -1;
  } else
    return;

  n32 = *n / 32;
  
  n32 *= 32;


  if(n32!=0) {
    for(j=0;j<n32;j+=32) {
      ix = ix0;
      for(i=i1;i<i2;i+=inc,ix+=*incx) {
	ip = ipiv[ix] - 1;
	if(ip != i) {
	  for(k=j;k<j+32;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	  }
	}
      }
    }
  }
  if(n32!=*n) {
    ix = ix0;
    for(i=i1;i<i2;i+=inc,ix+=*incx) {
      ip = ipiv[ix] - 1;
      if(ip != i) {
	for(k=n32;k<*n;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	}
      }
    }
  }
  return;
}
}
}

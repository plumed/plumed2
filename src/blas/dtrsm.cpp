/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib BLAS library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the Gromacs developers, and
not the reference BLAS implementation.

The reference BLAS implementation is available from http://www.netlib.org/blas 

BLAS does not come with a formal named "license", but a general statement that 

"The reference BLAS is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of Gromacs is GPL, we think it's only fair to give you the same rights to
our modified BLAS files as the original netlib versions, so do what you want with them.
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
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)(const char * side,
       const char * uplo,
       const char * transa,
       const char * diag,
       int *  m__,
       int *  n__,
       double *alpha__,
       double *a,
       int *  lda__,
       double *b,
       int *  ldb__)
{
  const char xside  = toupper(*side);
  const char xuplo  = toupper(*uplo);
  const char xtrans = toupper(*transa);
  const char xdiag  = toupper(*diag);
  int i,j,k;
  double temp;

  
  int m = *m__;
  int n = *n__;
  int lda = *lda__;
  int ldb = *ldb__;
  double alpha = *alpha__;

  if(n<=0)
    return;
  
  if(fabs(alpha)<PLUMED_GMX_DOUBLE_MIN) { 
    for(j=0;j<n;j++)
      for(i=0;i<m;i++)
	b[j*(ldb)+i] = 0.0;
    return;
  }

  if(xside=='L') {
    /* left side */
    if(xtrans=='N') {
      /* No transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  }
	  for(k=m-1;k>=0;k--) {
	    if(fabs(b[j*(ldb)+k])>PLUMED_GMX_DOUBLE_MIN) {
	      if(xdiag=='N')
		b[j*(ldb)+k] /= a[k*(lda)+k];
	      for(i=0;i<k;i++)
		b[j*(ldb)+i] -= b[j*(ldb)+k]*a[k*(lda)+i];
	    }
	  }
	}
      } else {
	/* lower */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<m;k++) {
	    if(fabs(b[j*(ldb)+k])>PLUMED_GMX_DOUBLE_MIN) {
	      if(xdiag=='N')
		b[j*(ldb)+k] /= a[k*(lda)+k];
	      for(i=k+1;i<m;i++)
		b[j*(ldb)+i] -= b[j*(ldb)+k]*a[k*(lda)+i];
	    }
	  }
	}
      }
    } else {
      /* Transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  for(i=0;i<m;i++) {
	    temp = alpha * b[j*(ldb)+i];
	    for(k=0;k<i;k++)
	      temp -= a[i*(lda)+k] * b[j*(ldb)+k];
	    if(xdiag=='N')
		temp /= a[i*(lda)+i];
	    b[j*(ldb)+i] = temp;
	  }
	}
      } else {
	/* lower */
	for(j=0;j<n;j++) {
	  for(i=m-1;i>=0;i--) {
	    temp = alpha * b[j*(ldb)+i];
	    for(k=i+1;k<m;k++)
	      temp -= a[i*(lda)+k] * b[j*(ldb)+k];
	    if(xdiag=='N')
		temp /= a[i*(lda)+i];
	    b[j*(ldb)+i] = temp;
	  }
	}
      }
    }
  } else {
    /* right side */
    if(xtrans=='N') {
      /* No transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<j;k++) {
	    if(fabs(a[j*(lda)+k])>PLUMED_GMX_DOUBLE_MIN) {
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= a[j*(lda)+k]*b[k*(ldb)+i];
	    }
	  }
	  if(xdiag=='N') {
	    temp = 1.0/a[j*(lda)+j];
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= temp;
	  }
	}
      } else {
	/* lower */
	for(j=n-1;j>=0;j--) {
	  if(fabs(alpha)>PLUMED_GMX_DOUBLE_MIN)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=j+1;k<n;k++) {
	    if(fabs(a[j*(lda)+k])>PLUMED_GMX_DOUBLE_MIN) {
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= a[j*(lda)+k]*b[k*(ldb)+i];
	    }
	  }
	  if(xdiag=='N') {
	    temp = 1.0/a[j*(lda)+j];
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= temp;
	  }
	}
      }
    } else {
      /* Transpose */
      if(xuplo=='U') {
	/* upper */
	for(k=n-1;k>=0;k--) {
	  if(xdiag=='N') {
	    temp = 1.0/a[k*(lda)+k];
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= temp;
	  }
	  for(j=0;j<k;j++) {
	    if(fabs(a[k*(lda)+j])>PLUMED_GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(fabs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      } else {
	/* lower */
	for(k=0;k<n;k++) {
	  if(xdiag=='N') {
	    temp = 1.0/a[k*(lda)+k];
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= temp;
	  }
	  for(j=k+1;j<n;j++) {
	    if(fabs(a[k*(lda)+j])>PLUMED_GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(fabs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      }      
    }
  }    
}
}
}

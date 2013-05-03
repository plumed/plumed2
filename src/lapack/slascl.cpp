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

#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slascl,SLASCL)(const char *type,
                        int *kl,
                        int *ku,
                        float *cfrom,
                        float *cto,
                        int *m,
                        int *n,
                        float *a,
                        int *lda,
                        int *info)
{
  const char ch=toupper(*type);
  int i,j,k,l,k1,k2,k3,k4;
  int done=0;
  float minval,smlnum,bignum;
  float cfromc, ctoc, cfrom1, cto1, mul;

  if(*n<=0 || *m<=0)
    return;

  minval = PLUMED_GMX_FLOAT_MIN;
  smlnum = minval / PLUMED_GMX_FLOAT_EPS;
  bignum = 1.0 / smlnum;

  cfromc = *cfrom;
  ctoc   = *cto;

  while(!done) {
    
    cfrom1 = cfromc * smlnum;
    cto1   = ctoc / bignum;

    if(fabs(cfrom1)>fabs(ctoc) && fabs(ctoc)>PLUMED_GMX_FLOAT_MIN) {
      mul = smlnum;
      done = 0;
      cfromc = cfrom1;
    } else if(fabs(cto1)>fabs(cfromc)) {
      mul = bignum;
      done = 0;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      done = 1;
    }

    switch(ch) {
    case 'G': 
      /* Full matrix */
      for(j=0;j<*n;j++)
	for(i=0;i<*m;i++)
	  a[j*(*lda)+i] *= mul;
      break;

    case 'L': 
      /* Lower triangular matrix */
      for(j=0;j<*n;j++)
	for(i=j;i<*m;i++)
	  a[j*(*lda)+i] *= mul;
      break;

    case 'U': 
      /* Upper triangular matrix */
      for(j=0;j<*n;j++) {
	k = (j < (*m-1)) ? j : (*m-1);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'H': 
      /* Upper Hessenberg matrix */
      for(j=0;j<*n;j++) {
	k = ((j+1) < (*m-1)) ? (j+1) : (*m-1);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'B': 
      /* Symmetric band matrix, lower bandwidth KL, upper KU,
       * only the lower half stored.
       */
      k3 = *kl;
      k4 = *n - 1;
      for(j=0;j<*n;j++) {
	k = (k3 < (k4-j)) ? k3 : (k4-j);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'Q': 
      /* Symmetric band matrix, lower bandwidth KL, upper KU,
       * only the upper half stored.
       */
      k1 = *ku;
      k3 = *ku;
      for(j=0;j<*n;j++) {
	k = ((k1-j) > 0) ? (k1-j) : 0;
	for(i=k;i<=k3;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'Z': 
      /* Band matrix, lower bandwidth KL, upper KU. */

      k1 = *kl + *ku;
      k2 = *kl;
      k3 = 2*(*kl) + *ku;
      k4 = *kl + *ku - 1 + *m;
      for(j=0;j<*n;j++) {
	k = ((k1-j) > k2) ? (k1-j) : k2;
	l = (k3 < (k4-j)) ? k3 : (k4-j);
	for(i=k;i<=l;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    default:
      *info = -1;
      return;
    }
  } /* finished */

  *info = 0;
  return;
}
}
}

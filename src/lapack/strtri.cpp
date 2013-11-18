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
#include <math.h>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "simple.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(strtri,STRTRI)(const char *uplo,
	const char *diag, 
	int *n,
	float *a, 
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__3, i__4, i__5;
    int j, jb, nb, nn;
    float c_b18 = 1.;
    float c_b22 = -1.;

    int upper;
    int nounit;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    *info = 0;
    upper = (*uplo=='U' || *uplo=='u');
    nounit = (*diag=='N' || *diag=='n');

    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    if (*n == 0) {
	return;
    }

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (fabs(a[*info + *info * a_dim1])<PLUMED_GMX_FLOAT_MIN) {
		return;
	    }
	}
	*info = 0;
    }

    nb = DTRTRI_BLOCKSIZE;
    if (nb <= 1 || nb >= *n) {

	PLUMED_BLAS_F77_FUNC(strti2,STRTI2)(uplo, diag, n, &a[a_offset], lda, info);
    } else {

	if (upper) {

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
		i__4 = nb, i__5 = *n - j + 1;
		jb = (i__4<i__5) ? i__4 : i__5;

		i__4 = j - 1;
		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda);
		i__4 = j - 1;
		PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], 
			lda);

		PLUMED_BLAS_F77_FUNC(strti2,STRTI2)("Upper", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	} else {

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
		i__1 = nb, i__4 = *n - j + 1;
		jb = (i__1<i__4) ? i__1 : i__4;
		if (j + jb <= *n) {

		    i__1 = *n - j - jb + 1;
		    PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda);
		    i__1 = *n - j - jb + 1;
		    PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda);
		}

		PLUMED_BLAS_F77_FUNC(strti2,STRTI2)("Lower", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	}
    }
    return;
}


}
}

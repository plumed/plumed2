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
PLUMED_BLAS_F77_FUNC(sstein,SSTEIN)(int *n, 
	float *d__, 
	float *e, 
	int *m, 
	float *w, 
	int *iblock,
	int *isplit, 
	float *z__,
	int *ldz, 
	float *work,
	int *iwork, 
	int *ifail,
	int *info)
{
    int z_dim1, z_offset, i__1, i__2, i__3;
    float d__2, d__3, d__4, d__5;

    int i__, j, b1, j1, bn;
    float xj, scl, eps, sep, nrm, tol;
    int its;
    float xjm, ztr, eps1;
    int jblk, nblk;
    int jmax;

    int iseed[4], gpind, iinfo;
    float ortol;
    int indrv1, indrv2, indrv3, indrv4, indrv5;
    int nrmchk;
    int blksiz;
    float onenrm, dtpcrt, pertol;
    int c__2 = 2;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    --e;
    --w;
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;

    *info = 0;

    xjm = 0.0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ifail[i__] = 0;
    }

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -4;
    } else if (*ldz < (*n)) {
	*info = -9;
    } else {
	i__1 = *m;
	for (j = 2; j <= i__1; ++j) {
	    if (iblock[j] < iblock[j - 1]) {
		*info = -6;
		break;
	    }
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
		*info = -5;
		break;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    if (*n == 0 || *m == 0) {
	return;
    } else if (*n == 1) {
	z__[z_dim1 + 1] = 1.;
	return;
    }

    eps = PLUMED_GMX_FLOAT_EPS;

    for (i__ = 1; i__ <= 4; ++i__) {
	iseed[i__ - 1] = 1;
    }

    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;

    j1 = 1;
    i__1 = iblock[*m];
    for (nblk = 1; nblk <= i__1; ++nblk) {

	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = isplit[nblk - 1] + 1;
	}
	bn = isplit[nblk];
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    continue;
	}
	gpind = b1;

	onenrm = fabs(d__[b1]) + fabs(e[b1]);
	d__3 = onenrm;
	d__4 = fabs(d__[bn]) + fabs(e[bn - 1]);
	onenrm = (d__3>d__4) ? d__3 : d__4;
	i__2 = bn - 1;
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
	  d__4 = onenrm;
	  d__5 = fabs(d__[i__]) + fabs(e[i__ - 1]) + fabs(e[i__]);
	    onenrm = (d__4>d__5) ? d__4 : d__5;
	}
	ortol = onenrm * .001;

	dtpcrt = sqrt(.1 / blksiz);

	jblk = 0;
	i__2 = *m;
	for (j = j1; j <= i__2; ++j) {
	    if (iblock[j] != nblk) {
		j1 = j;
		break;
	    }
	    ++jblk;
	    xj = w[j];

	    if (blksiz == 1) {
		work[indrv1 + 1] = 1.;
		goto L120;
	    }

	    if (jblk > 1) {
		eps1 = fabs(eps * xj);
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

	    PLUMED_BLAS_F77_FUNC(slarnv,SLARNV)(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
	    i__3 = blksiz - 1;
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
	    i__3 = blksiz - 1;
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

	    tol = 0.;
	    PLUMED_BLAS_F77_FUNC(slagtf,SLAGTF)(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

	    d__2 = eps;
	    d__3 = fabs(work[indrv4 + blksiz]);
	    scl = blksiz * onenrm * ((d__2>d__3) ? d__2 : d__3) / PLUMED_BLAS_F77_FUNC(sasum,SASUM)(&blksiz, &work[
		    indrv1 + 1], &c__1);
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&blksiz, &scl, &work[indrv1 + 1], &c__1);

	    PLUMED_BLAS_F77_FUNC(slagts,SLAGTS)(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

	    if (jblk == 1) {
		goto L90;
	    }
	    if (fabs(xj - xjm) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i__ = gpind; i__ <= i__3; ++i__) {
		    ztr = -PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
			    i__ * z_dim1], &c__1);
		    PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &
			    work[indrv1 + 1], &c__1);
		}
	    }

L90:
	    jmax = PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)(&blksiz, &work[indrv1 + 1], &c__1);
	    nrm = fabs(work[indrv1 + jmax]);

	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L110;

L100:
	    ++(*info);
	    ifail[*info] = j;

L110:
	    scl = 1. / PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(&blksiz, &work[indrv1 + 1], &c__1);
	    jmax = PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)(&blksiz, &work[indrv1 + 1], &c__1);
	    if (work[indrv1 + jmax] < 0.) {
		scl = -scl;
	    }
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L120:
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z__[i__ + j * z_dim1] = 0.;
	    }
	    i__3 = blksiz;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
	    }

	    xjm = xj;
	}
    }

    return;

}


}
}

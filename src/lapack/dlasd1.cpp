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
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd1,DLASD1)(int *nl, 
	int *nr, 
	int *sqre, 
	double *d__, 
	double *alpha, 
	double *beta, 
	double *u, 
	int *ldu, 
	double *vt, 
	int *ldvt, 
	int *idxq, 
	int *iwork, 
	double *work, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    double d__1, d__2;

    int i__, k, m, n, n1, n2, iq, iz, iu2, ldq, idx, ldu2, ivt2, 
	    idxc, idxp, ldvt2;
    int isigma;
    double orgnrm;
    int coltyp;
    int c__0 = 0;
    double one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --idxq;
    --iwork;
    --work;

    *info = 0;

    if (*nl < 1) {
	*info = -1;
    } else if (*nr < 1) {
	*info = -2;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    n = *nl + *nr + 1;
    m = n + *sqre;


    ldu2 = n;
    ldvt2 = m;

    iz = 1;
    isigma = iz + m;
    iu2 = isigma + n;
    ivt2 = iu2 + ldu2 * n;
    iq = ivt2 + ldvt2 * m;

    idx = 1;
    idxc = idx + n;
    coltyp = idxc + n;
    idxp = coltyp + n;

    d__1 = fabs(*alpha);
    d__2 = fabs(*beta);
    orgnrm = (d__1>d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(d__[i__]) > orgnrm) {
	    orgnrm = fabs(d__[i__]);
	}
    }
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    PLUMED_BLAS_F77_FUNC(dlasd2,DLASD2)(nl, nr, sqre, &k, &d__[1], &work[iz], alpha, beta, &u[u_offset], 
	    ldu, &vt[vt_offset], ldvt, &work[isigma], &work[iu2], &ldu2, &
	    work[ivt2], &ldvt2, &iwork[idxp], &iwork[idx], &iwork[idxc], &
	    idxq[1], &iwork[coltyp], info);

    ldq = k;
    PLUMED_BLAS_F77_FUNC(dlasd3,DLASD3)(nl, nr, sqre, &k, &d__[1], &work[iq], &ldq, &work[isigma], &u[
	    u_offset], ldu, &work[iu2], &ldu2, &vt[vt_offset], ldvt, &work[
	    ivt2], &ldvt2, &iwork[idxc], &iwork[coltyp], &work[iz], info);
    if (*info != 0) {
	return;
    }
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = k;
    n2 = n - k;
    PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}
}
}

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

#include "simple.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasr,DLASR)(const char *side, 
       const char *pivot, 
       const char *direct, 
       int *m,
       int *n, 
       double *c__, 
       double *s, 
       double *a, 
       int *lda)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, info;
    double temp;
    double ctemp, stemp;

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;

    if (*m == 0 || *n == 0) {
	return;
    }
    if (*side=='L' || *side=='l') {

	if (*pivot=='V' || *pivot=='v') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='T' || *pivot=='t') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
			}
		    }
		}
	    }
	} else if (*pivot=='B' || *pivot=='b') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    }
	}
    } else if (*side=='R' || *side=='r') {

	if (*pivot=='V' || *pivot=='v') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='T' || *pivot=='t') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='B' || *pivot=='b') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || fabs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    }
	}
    }

    return;

}


}
}

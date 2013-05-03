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
PLUMED_BLAS_F77_FUNC(dsyr2,DSYR2)(const char *    uplo,
                      int *     n__,
                      double *  alpha__,
                      double *  x,
                      int *     incx__,
                      double *  y,
                      int *     incy__,
                      double *  a,
                      int *     lda__)
{
    int kx,ky,ix,iy,jx,jy,j,i;
    double temp1,temp2;
    const char ch=toupper(*uplo);
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    float alpha = *alpha__;
    
    
    if(n<=0 || fabs(alpha)<PLUMED_GMX_DOUBLE_MIN || incx==0 || incy==0 ||
       (ch != 'U' && ch != 'L'))
        return;
    
    jx = jy = kx = ky = 0;
    
    /* init start points for non-unit increments */
    if(incx!=1 || incy!=1) {
        if(incx>0)
            kx = 1;
        else
            kx = 1 - (n - 1)*(incx);
        if(incy>0)
            ky = 1;
        else
            ky = 1 - (n - 1)*(incy);
        
        jx = kx;
        jy = ky;
    }
    
    if(ch == 'U') {
        /* Data in upper part of A */
        if(incx==1 && incy==1) {
            /* Unit increments for both x and y */
            for(j=1;j<=n;j++) {
                if( fabs(x[j-1])>PLUMED_GMX_DOUBLE_MIN  || fabs(y[j-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=1;i<=j;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( fabs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN || fabs(y[jy-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[jy-1];
                    temp2 = alpha * x[jx-1];
                    ix = kx;
                    iy = ky;
                    for(i=1;i<=j;i++) {
                        a[(j-1)*(lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
            }
        }
    } else {
        /* Data in lower part of A */
        if(incx==1 && incy==1) {
            /* Unit increments for both x and y */
            for(j=1;j<=n;j++) {
                if( fabs(x[j-1])>PLUMED_GMX_DOUBLE_MIN  || fabs(y[j-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=j;i<=n;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( fabs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN || fabs(y[jy-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[jy-1];
                    temp2 = alpha * x[jx-1];
                    ix = jx;
                    iy = jy;
                    for(i=j;i<=n;i++) {
                        a[(j-1)*(lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
            }
        }
    }
    
    return;
}
}
}

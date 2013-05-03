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
#include "blas.h"


namespace PLMD{
namespace blas{
float
PLUMED_BLAS_F77_FUNC(sdot,SDOT)(int *n_arg,
                    float *dx,
                    int *incx_arg,
                    float *dy,
                    int *incy_arg)
{
    int i,ix,iy,m;
    int n=*n_arg;
    int incx = *incx_arg;
    int incy = *incy_arg;
    float t1;
    
    if(n<=0)
        return 0.0;
    
    t1 = 0.0;
    
    if(incx!=1 || incy!=1) {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*incx;
        if(incy<0)
            iy = (1-n)*incy;
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            t1 += dx[ix] * dy[iy];
        
        return t1;
        
    } else {
        
        m = n%5;
        
        for(i=0;i<m;i++)
            t1 += dx[i] * dy[i];
        
        /* unroll */
        for(i=m;i<n;i+=5) 
            t1  =  t1 + dx[i] * dy[i]   
                +    dx[i+1] * dy[i+1] 
                +    dx[i+2] * dy[i+2] 
                +    dx[i+3] * dy[i+3]   
                +    dx[i+4] * dy[i+4];   
        
        return t1;
    }
}


}
}

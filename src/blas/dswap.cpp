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
void
PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(int *n__,
                      double *dx,
                      int *incx__,
                      double *dy,
                      int *incy__)
{
  int i,ix,iy;
  double d1,d2,d3;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  
  if(n<=0)
    return;

  if(incx==1 && incy==1) {
    for(i=0;i<(n-3);i+=3) {
      d1      = dx[i];
      d2      = dx[i+1];
      d3      = dx[i+2];
      dx[i]   = dy[i];
      dx[i+1] = dy[i+1];
      dx[i+2] = dy[i+2];
      dy[i]   = d1;
      dy[i+1] = d2;
      dy[i+2] = d3;
    }
    /* continue with last i value */
    for(;i<n;i++) {
      d1      = dx[i];
      dx[i]   = dy[i];
      dy[i]   = d1;
    }

  } else {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = incx * (1 - n);
    if(incy<0)
      iy = incy * (1 - n);

    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      d1     = dx[ix];
      dx[ix] = dy[iy];
      dy[iy] = d1;
    }
  }
  return;
}
 
}
}

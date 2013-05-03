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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(int   *n,
                        double *alpha,
                        double *x,
                        int    *incx,
                        double *tau)
{
  double xnorm,t;
  int    ti1,knt,j;
  double minval,safmin,rsafmn,beta;

  if(*n<=1) {
    *tau = 0;
    return;
  }

  ti1 = *n-1;

  xnorm = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(&ti1,x,incx);

  if(fabs(xnorm)<PLUMED_GMX_DOUBLE_MIN) {
    *tau = 0.0;
  } else {

    t = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(alpha,&xnorm);

    if(*alpha<0)
      beta = t;
    else
      beta = -t;

    minval = PLUMED_GMX_DOUBLE_MIN;
    
    safmin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS) / PLUMED_GMX_DOUBLE_EPS;

        
    if(fabs(beta)<safmin) {

      knt = 0;
      rsafmn = 1.0 / safmin;
      
      while(fabs(beta)<safmin) {
	knt++;
	ti1 = *n-1;
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&ti1,&rsafmn,x,incx);
	beta *= rsafmn;
	*alpha *= rsafmn;
      }
      
      /* safmin <= beta <= 1 now */
      ti1 = *n-1;
      xnorm = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(&ti1,x,incx);
      t = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(alpha,&xnorm);
      
      if(*alpha<0)
	beta = t;
      else
	beta = -t;
      
      *tau = (beta-*alpha)/beta;

      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&ti1,&t,x,incx);
   
      *alpha = beta;
      for(j=0;j<knt;j++)
	*alpha *= safmin;
    } else {
      *tau = (beta-*alpha)/beta;
      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&ti1,&t,x,incx);
      *alpha = beta;
    }
  }
   
  return;
}
}
}

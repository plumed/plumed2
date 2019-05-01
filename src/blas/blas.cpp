/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib BLAS library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference BLAS implementation.

The reference BLAS implementation is available from http://www.netlib.org/blas 

BLAS does not come with a formal named "license", but a general statement that 

"The reference BLAS is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified BLAS files as the original netlib versions, so do what you want with them.
However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#if ! defined (__PLUMED_HAS_EXTERNAL_BLAS)
#include <cmath>
#include "blas.h"

namespace PLMD{
namespace blas{
double
PLUMED_BLAS_F77_FUNC(dasum,DASUM)(int *n__, 
                      double *dx, 
                      int *incx__)
{
    int i__1, i__2;
    
    int i__, m, mp1;
    double dtemp;
    int nincx;
    
    int n = *n__;
    int incx = *incx__;
    
    --dx;
    
    dtemp = 0.;
    if (n <= 0 || incx <= 0) {
        return 0.0;
    }
    if (incx != 1) {
        nincx = n * incx;
        i__1 = nincx;
        i__2 = incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dtemp += std::abs(dx[i__]);
        }
        return dtemp;
    }
    
    m = n % 6;
    if (m != 0) {
        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dtemp += std::abs(dx[i__]);
        }
        if (n < 6) {
            return dtemp;
        }
    }
    mp1 = m + 1;
    i__2 = n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
        dtemp = dtemp + std::abs(dx[i__]) + std::abs(dx[i__ + 1]) + 
        std::abs(dx[i__ + 2]) + std::abs(dx[i__+ 3]) + std::abs(dx[i__ + 4]) +
        std::abs(dx[i__ + 5]);
    }
    return dtemp;
}


}
}
#include "blas.h"


namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(int   *   n_arg,
                      double *   da_arg,
                      double *   dx,
                      int *      incx_arg,
                      double *   dy,
                      int *      incy_arg)
{
  int i,ix,iy;
  int n=*n_arg;
  double da=*da_arg;
  int incx = *incx_arg;
  int incy = *incy_arg;

  if (n<=0)
    return;

  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*incx;
    if(incy<0)
      iy = (1-n)*incy;
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) 
      dy[iy] += da*dx[ix];

    return;

  } else {

    /* unroll */
    
    for(i=0;i<(n-4);i+=4) {
      dy[i]   += da*dx[i];
      dy[i+1] += da*dx[i+1];
      dy[i+2] += da*dx[i+2];
      dy[i+3] += da*dx[i+3];
    }
    /* continue with current value of i */
    for(;i<n;i++)
      dy[i]   += da*dx[i];
  }
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(int *n__,
                      double *dx,
                      int *incx__,
                      double *dy,
                      int *incy__)
{
    int i,ix,iy;

    int n= *n__;
    int incx = *incx__;
    int incy = *incy__;
    

    if(incx!=1 || incy!=1) {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*(incx);
        if(incy<0)
            iy = (1-n)*(incy);
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            dy[iy] = dx[ix];
        
        return;
        
    } else {
        
        /* unroll */
        
        for(i=0;i<(n-8);i+=8) {
            dy[i]   = dx[i];
            dy[i+1] = dx[i+1];
            dy[i+2] = dx[i+2];
            dy[i+3] = dx[i+3];
            dy[i+4] = dx[i+4];
            dy[i+5] = dx[i+5];
            dy[i+6] = dx[i+6];
            dy[i+7] = dx[i+7];
        }
        /* continue with current value of i */
        for(;i<n;i++)
            dy[i] = dx[i];
    }
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
double
PLUMED_BLAS_F77_FUNC(ddot,DDOT)(int *n_arg,
                    double *dx,
                    int *incx_arg,
                    double *dy,
                    int *incy_arg)
{
    int i,ix,iy,m;
    int n=*n_arg;
    int incx = *incx_arg;
    int incy = *incy_arg;
    double t1;
    
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
#include <cctype>
#include <cmath>

#include "real.h"

#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)(const char *transa,
                      const char *transb,
                      int *m__,
                      int *n__,
                      int *k__,
                      double *alpha__,
                      double *a,
                      int *lda__,
                      double *b,
                      int *ldb__,
                      double *beta__,
                      double *c,
                      int *ldc__)
{
  const char tra=std::toupper(*transa);
  const char trb=std::toupper(*transb);
  double temp;
  int i,j,l;

  int m = *m__;
  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  double alpha = *alpha__;
  double beta  = *beta__;
  
  if(m==0 || n==0 || (( std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN || k==0) && std::abs(beta-1.0)<PLUMED_GMX_DOUBLE_EPS))
    return;

  if(std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN) {
    if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) {
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] = 0.0;
    } else {
      /* nonzero beta */
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] *= beta;
    }
    return;
  }

  if(trb=='N') {
    if(tra=='N') {
      
      for(j=0;j<n;j++) {
	if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( std::abs(b[ j*(ldb) + l ])>PLUMED_GMX_DOUBLE_MIN) {
	    temp = alpha * b[ j*(ldb) + l ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
    } else {
      /* transpose A, but not B */
      for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[j*(ldb)+l];
	  if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
      }
    }
  } else {
    /* transpose B */
    if(tra=='N') {

      /* transpose B, but not A */

      for(j=0;j<n;j++) {
	if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( std::abs(b[ l*(ldb) + j ])>PLUMED_GMX_DOUBLE_MIN) {
	    temp = alpha * b[ l*(ldb) + j ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
 
    } else {
      /* Transpose both A and B */
       for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[l*(ldb)+j];
	  if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
       }
    }
  }
}
}
}
#include <cctype>
#include <cmath>

#include "real.h"

#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)(const char *trans, 
       int *m__,
       int *n__,
       double *alpha__,
       double *a,
       int *lda__,
       double *x,
       int *incx__,
       double *beta__,
       double *y,
       int *incy__)
{
  const char ch=std::toupper(*trans);
  int lenx,leny,kx,ky;
  int i,j,jx,jy,ix,iy;
  double temp;

  int m = *m__;
  int n = *n__;
  double alpha = *alpha__;
  double beta = *beta__;
  int incx = *incx__;
  int incy = *incy__;
  int lda = *lda__;
  
  if(n<=0 || m<=0 || (std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN && std::abs(beta-1.0)<PLUMED_GMX_DOUBLE_EPS))
    return;

  if(ch=='N') {
    lenx = n;
    leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  
   if(incx>0)
    kx = 1;
  else
    kx = 1 - (lenx -1)*(incx);

  if(incy>0)
    ky = 1;
  else
    ky = 1 - (leny -1)*(incy);
 
  if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS) {
    if(incy==1) {
      if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	for(i=0;i<leny;i++)
	  y[i] = 0.0;
      else
	for(i=0;i<leny;i++)
	  y[i] *= beta;
    } else {
      /* non-unit incr. */
      iy = ky;
      if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) 
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] = 0.0;
      else
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] *= beta;
    }
  }
  
  if(std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN)
    return;
  
  if(ch=='N') {
    jx = kx;
    if(incy==1) {
      for(j=1;j<=n;j++,jx+=incx) 
	if(std::abs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  for(i=1;i<=m;i++)
	    y[i-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jx+=incx) 
	if(std::abs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  iy = ky;
	  for(i=1;i<=m;i++,iy+=incy)
	    y[iy-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    }
  } else {
    /* transpose */
    jy = ky;
    if(incx==1) {
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	for(i=1;i<=m;i++)
	  temp += a[(j-1)*(lda)+(i-1)] * x[i-1];
	y[jy-1] += alpha * temp;
      }
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	ix = kx;
	for(i=1;i<=m;i++,ix+=incx)
	  temp += a[(j-1)*(lda)+(i-1)] * x[ix-1];
	y[jy-1] += alpha * temp;
      }
    }
  }
} 
   
}
}
#include <cmath>

#include "real.h"

#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dger,DGER)(int *m__,
                    int *n__,
                    double *alpha__,
                    double *x,
                    int *incx__,
                    double *y,
                    int *incy__,
                    double *a,
                    int *lda__)
{
    int ix,kx,jy;
    int i,j;
    double temp;
    
    
    int m = *m__;
    int n = *n__;
    int incx = *incx__;
    int incy = *incy__;
    int lda = *lda__;
    double alpha = *alpha__;
    
    if(m<=0 || n<=0 || std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN)
        return;
    
    if(incy>0)
        jy = 0;
    else
        jy = incy * (1 - n);
    
    if(incx==1) {
        for(j=0;j<n;j++,jy+=incy)
            if(std::abs(y[jy])>PLUMED_GMX_DOUBLE_MIN) {
                temp = alpha * y[jy];
                for(i=0;i<m;i++)
                    a[j*(lda)+i] += temp*x[i];
            }
    } else {
        /* non-unit incx */
        if(incx>0) 
            kx = 0;
        else
            kx = incx * (1 - m);
        
        for(j=0;j<n;j++,jy+=incy) {
            if(std::abs(y[jy])>PLUMED_GMX_DOUBLE_MIN) {
                temp = alpha * y[jy];
                ix = kx;
                for(i=0;i<m;i++,ix+=incx)
                    a[j*(lda)+i] += temp*x[ix];
            }
        }
    }
        return;
}
}
}
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
double
PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(int  *     n__,
                      double *    x,
                      int    *    incx__)
{
    int ix,max_ix;
    double ssq,scale,absxi,t;
    
    int n = *n__;
    int incx = *incx__;
    
    if(n<1 || incx<1)
        return 0;
    else if (n==1) {
        t = x[0];
        if(t>=0)
            return t;
        else 
            return -t;
    }
    
    scale = 0.0;
    ssq   = 1.0;
    
    max_ix = 1+(n-1)*(incx);
    for(ix=1;ix<=max_ix;ix+=incx) {
        t = x[ix-1];
        if(std::abs(t)>PLUMED_GMX_DOUBLE_MIN) {
            absxi = (t>=0) ? t : (-t);
            if(scale<absxi) {
                t = scale/absxi;
                t = t*t;
                ssq = ssq*t + 1.0;
                scale = absxi;
            } else {
                t = absxi/scale;
                ssq += t*t;
            }
        }
    }
    return scale*std::sqrt(ssq);
    
}


 
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(drot,DROT)(int *n__,
      double *dx,
      int *incx__,
      double *dy,
      int *incy__,
      double *c__,
      double *s__)
{
  int i,ix,iy;
  double dtemp;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  double c = *c__;
  double s = *s__;
  
  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*(incx);
    if(incy<0)
      iy = (1-n)*(incy);
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      dtemp  = (c) * dx[ix] + (s) * dy[iy];
      dy[iy] = (c) * dy[iy] - (s) * dx[ix];
      dx[ix] = dtemp;
    }

    return;

  } else {

    /* unit increments */   
    for(i=0;i<n;i++) {
      dtemp = (c) * dx[i] + (s) * dy[i];
      dy[i] = (c) * dy[i] - (s) * dx[i];
      dx[i] = dtemp;      
    }

  }
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(int  *    n__,
                      double *   fact__,
                      double *   dx,
                      int    *   incx__)
{
    int nincx,i;

    int n = *n__;
    double fact = *fact__;
    int incx = *incx__;
    
    if(n<=0 || incx<=0)
        return;
    
    if(incx==1) {
        /* Unrool factor 5 */
        for(i=0;i<(n-5);i+=5) {
            dx[i]   *= fact;
            dx[i+1] *= fact;
            dx[i+2] *= fact;
            dx[i+3] *= fact;
            dx[i+4] *= fact;
        }    
        /* continue with current value of i */
        for(;i<n;i++)
            dx[i]   *= fact;
        
        return;
    } else {
        /* inc != 1 */
        nincx = n * (incx);
        for (i=0;i<nincx;i+=incx)
            dx[i] *= fact;
        
        return;
    } 
    
}
}
}
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
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)(const char *uplo,
       int *n__,
       double *alpha__,
       double *a,
       int *lda__,
       double *x,
       int *incx__,
       double *beta__,
       double *y,
       int *incy__)
{
    const char ch=std::toupper(*uplo);
    int kx,ky,i,j,ix,iy,jx,jy;
    double temp1,temp2;
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    double alpha = *alpha__;
    double beta  = *beta__;
    
    if(n<=0 || incx==0 || incy==0)
        return;
    
    if(incx>0)
        kx = 1;
    else
        kx = 1 - (n -1)*(incx);
    
    if(incy>0)
        ky = 1;
    else
        ky = 1 - (n -1)*(incy);
    
    if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS) {
        if(incy==1) {
            if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) 
                for(i=1;i<=n;i++)
                    y[i-1] = 0.0;
            else
                for(i=1;i<=n;i++)
                    y[i-1] *= beta;
        } else {
            /* non-unit incr. */
            iy = ky;
            if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) 
                for(i=1;i<=n;i++) {
                    y[iy-1] = 0.0;
                    iy += incy;
                }
                    else
                        for(i=1;i<=n;i++) {
                            y[iy-1] *= beta;
                            iy += incy;
                        }
        }
    }
        
        if(std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN) 
            return;
        
        if(ch=='U') {
            if(incx==1 && incy==1) {
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[j-1];
                    temp2 = 0.0;
                    for(i=1;i<j;i++) {
                        y[i-1] += temp1*a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[i-1];
                    }
                    y[j-1] += temp1*a[(j-1)*(lda)+(j-1)] + alpha *temp2;
                }
            } else {
                /* non-unit incr. */
                jx = kx;
                jy = ky;
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[jx-1];
                    temp2 = 0.0;
                    ix = kx;
                    iy = ky;
                    for(i=1;i<j;i++) {
                        y[iy-1] += temp1 * a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[ix-1];
                        ix += incx;
                        iy += incy;
                    }
                    y[jy-1] += temp1*a[(j-1)*(lda)+(j-1)] + alpha*temp2;
                    jx += incx;
                    jy += incy;
                }
            }
        } else {
            /* lower */
            if(incx==1 && incy==1) {
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[j-1];
                    temp2 = 0.0;
                    y[j-1] += temp1 * a[(j-1)*(lda)+(j-1)];
                    for(i=j+1;i<=n;i++) {
                        y[i-1] += temp1*a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[i-1];
                    }
                    y[j-1] += alpha *temp2;
                }
            } else {
                /* non-unit incr. */
                jx = kx;
                jy = ky;
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[jx-1];
                    temp2 = 0.0;
                    y[jy-1] += temp1 * a[(j-1)*(lda)+(j-1)];
                    ix = jx;
                    iy = jy;
                    for(i=j+1;i<=n;i++) {
                        ix += incx;
                        iy += incy;
                        y[iy-1] += temp1 * a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[ix-1];
                    }
                    y[jy-1] += alpha*temp2;
                    jx += incx;
                    jy += incy;
                }
            }
        }
        return;
}    
}
}
#include <cctype>
#include <cmath>

#include "real.h"

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
    const char ch=std::toupper(*uplo);
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    float alpha = *alpha__;
    
    
    if(n<=0 || std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN || incx==0 || incy==0 ||
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
                if( std::abs(x[j-1])>PLUMED_GMX_DOUBLE_MIN  || std::abs(y[j-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=1;i<=j;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN || std::abs(y[jy-1])>PLUMED_GMX_DOUBLE_MIN ) {
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
                if( std::abs(x[j-1])>PLUMED_GMX_DOUBLE_MIN  || std::abs(y[j-1])>PLUMED_GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=j;i<=n;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>PLUMED_GMX_DOUBLE_MIN || std::abs(y[jy-1])>PLUMED_GMX_DOUBLE_MIN ) {
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
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(dsyr2k,DSYR2K)(const char *uplo, 
	const char *trans,
	int *n__,
	int *k__,
	double *alpha__,
	double *a,
	int *lda__,
	double *b,
	int *ldb__,
	double *beta__,
	double *c,
	int *ldc__)
{
  char ch1,ch2;
  int i,j,l;
  double temp1,temp2;

  
  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  double alpha = *alpha__;
  double beta  = *beta__;
  
  ch1 = std::toupper(*uplo);
  ch2 = std::toupper(*trans);

  if(n==0 || ( ( std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN || k==0 ) && std::abs(beta-1.0)<PLUMED_GMX_DOUBLE_EPS))
    return;

  if(std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN ) {
    if(ch1=='U') {
      if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    } else {
      /* lower */
      if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    }
    return;
  }

  if(ch2=='N') {
    if(ch1=='U') {
      for(j=1;j<=n;j++) {
	if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	  for(i=1;i<=j;i++)
	     c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS)
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( std::abs(a[(l-1)*(lda)+(j-1)])>PLUMED_GMX_DOUBLE_MIN ||
	      std::abs(b[(l-1)*(ldb)+(j-1)])>PLUMED_GMX_DOUBLE_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=1;i<=j;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    } else {
      /* lower */
      for(j=1;j<=n;j++) {
	if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(std::abs(beta-1.0)>PLUMED_GMX_DOUBLE_EPS)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( std::abs(a[(l-1)*(lda)+(j-1)])>PLUMED_GMX_DOUBLE_MIN ||
	      std::abs(b[(l-1)*(ldb)+(j-1)])>PLUMED_GMX_DOUBLE_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=j;i<=n;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    }
  } else {
    /* transpose */
    if(ch1=='U') {
      for(j=1;j<=n;j++) 
	for(i=1;i<=j;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    } else {
      /* lower */
      for(j=1;j<=n;j++) 
	for(i=j;i<=n;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(std::abs(beta)<PLUMED_GMX_DOUBLE_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    }
  }
  return;
}
}
}
#include <cmath>

#include "real.h"

#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)(const char *side, 
       const char *uplo, 
       const char *transa, 
       const char *diag, 
       int *m__, 
       int *n__, 
       double *alpha__, 
       double *a, 
       int *lda__, 
       double *b, 
       int *ldb__)
{
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    int m = *m__;
    int n = *n__;
    int lda = *lda__;
    int ldb = *ldb__;
    double alpha = *alpha__;
    
    /* Local variables */
    int i__, j, k;
    double temp;
    int lside;
    int upper;
    int nounit;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = (*side=='L' || *side=='l');

    nounit = (*diag=='N' || *diag=='n');
    upper = (*uplo=='U' || *uplo=='u');

    if (n == 0) {
	return;
    }
    if (std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN) {
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
	    }
	}
	return;
    }
    if (lside) {
	if (*transa=='N' || *transa=='n') {
	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (k = 1; k <= i__2; ++k) {
			if (std::abs(b[k + j * b_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = m; k >= 1; --k) {
			if (std::abs(b[k + j * b_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (*transa=='N' || *transa=='n') {

	    if (upper) {
		for (j = n; j >= 1; --j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (std::abs(a[k + j * a_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__2 = n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (std::abs(a[k + j * a_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (std::abs(a[j + k * a_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (std::abs(temp-1.0)>PLUMED_GMX_DOUBLE_EPS) {
			i__2 = m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    } else {
		for (k = n; k >= 1; --k) {
		    i__1 = n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (std::abs(a[j + k * a_dim1])>PLUMED_GMX_DOUBLE_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (std::abs(temp-1.0)>PLUMED_GMX_DOUBLE_EPS) {
			i__1 = m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
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
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)(const char *uplo, 
       const char *trans,
       const char *diag, 
       int *n__, 
       double *a, 
       int *lda__, 
       double *x, 
       int *incx__)
{
    int a_dim1, a_offset, i__1, i__2;

    int i__, j, ix, jx, kx;
    double temp;
    int nounit;
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    if (n == 0) {
	return;
    }

    nounit = (*diag=='n' || *diag=='N');

    if (incx <= 0) {
	kx = 1 - (n - 1) * incx;
    } else {
	kx = 1;
    }

    if (*trans=='N' || *trans=='n') {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (std::abs(x[j])>PLUMED_GMX_DOUBLE_MIN) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (std::abs(x[jx])>PLUMED_GMX_DOUBLE_MIN) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix += incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += incx;
		}
	    }
	} else {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    if (std::abs(x[j])>PLUMED_GMX_DOUBLE_MIN) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		kx += (n - 1) * incx;
		jx = kx;
		for (j = n; j >= 1; --j) {
		    if (std::abs(x[jx])>PLUMED_GMX_DOUBLE_MIN) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix -= incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= incx;
		}
	    }
	}
    } else {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx + (n - 1) * incx;
		for (j = n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx -= incx;
		}
	    }
	} else {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx += incx;
		}
	    }
	}
    }

    return;

}


}
}
#include <cctype>
#include <cmath>

#include "real.h"
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
  const char xside  = std::toupper(*side);
  const char xuplo  = std::toupper(*uplo);
  const char xtrans = std::toupper(*transa);
  const char xdiag  = std::toupper(*diag);
  int i,j,k;
  double temp;

  
  int m = *m__;
  int n = *n__;
  int lda = *lda__;
  int ldb = *ldb__;
  double alpha = *alpha__;

  if(n<=0)
    return;
  
  if(std::abs(alpha)<PLUMED_GMX_DOUBLE_MIN) { 
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  }
	  for(k=m-1;k>=0;k--) {
	    if(std::abs(b[j*(ldb)+k])>PLUMED_GMX_DOUBLE_MIN) {
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<m;k++) {
	    if(std::abs(b[j*(ldb)+k])>PLUMED_GMX_DOUBLE_MIN) {
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<j;k++) {
	    if(std::abs(a[j*(lda)+k])>PLUMED_GMX_DOUBLE_MIN) {
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
	  if(std::abs(alpha)>PLUMED_GMX_DOUBLE_MIN)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=j+1;k<n;k++) {
	    if(std::abs(a[j*(lda)+k])>PLUMED_GMX_DOUBLE_MIN) {
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
	    if(std::abs(a[k*(lda)+j])>PLUMED_GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(std::abs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
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
	    if(std::abs(a[k*(lda)+j])>PLUMED_GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(std::abs(alpha-1.0)>PLUMED_GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      }      
    }
  }    
}
}
}
#include <cmath>
#include "blas.h"

namespace PLMD{
namespace blas{
int
PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)(int *n__,
                        double *dx,
                        int *incx__)
{
  int i,ix,idxmax;
  double dmax,tmp;

  int n    = *n__;
  int incx = *incx__;
  
  if(n<1 || incx<=0)
    return -1;

  if(n==1)
    return 1;

  dmax = std::abs(dx[0]);
  idxmax = 1;

  if(incx==1) {
    for(i=1;i<n;i++) {
      tmp = std::abs(dx[i]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = i+1;
      }
    }
  } else {
    /* Non-unit increments */
    ix = incx; /* this is really 0 + an increment */
    for(i=1;i<n;i++,ix+=incx) {
      tmp = std::abs(dx[ix]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = ix+1;
      }
    }    
  }
  return idxmax;
}
}
}
#include <cmath>
#include "blas.h"

namespace PLMD{
namespace blas{
int
PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)(int *n__,
       float *dx,
       int *incx__)
{
  int i,ix,idxmax;
  float dmax,tmp;

  int n    = *n__;
  int incx = *incx__;
  
  if(n<1 || incx<=0)
    return -1;

  if(n==1)
    return 1;

  dmax = std::abs(dx[0]);
  idxmax = 1;

  if(incx==1) {
    for(i=1;i<n;i++) {
      tmp = std::abs(dx[i]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = i+1;
      }
    }
  } else {
    /* Non-unit increments */
    ix = incx; /* this is really 0 + an increment */
    for(i=1;i<n;i++,ix+=incx) {
      tmp = std::abs(dx[ix]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = ix+1;
      }
    }    
  }
  return idxmax;
}
}
}
#include <cmath>
#include "blas.h"

namespace PLMD{
namespace blas{
float
PLUMED_BLAS_F77_FUNC(sasum,SASUM)(int *n__, 
       float *dx, 
       int *incx__)
{
    int i__1, i__2;
    
    int i__, m, mp1;
    float dtemp;
    int nincx;
    
    int n = *n__;
    int incx = *incx__;
    
    
    --dx;
    
    dtemp = 0.;
    if (n <= 0 || incx <= 0) {
        return 0.0;
    }
    if (incx != 1) {
        nincx = n * incx;
        i__1 = nincx;
        i__2 = incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dtemp += std::abs(dx[i__]);
        }
        return dtemp;
    }
    
    m = n % 6;
    if (m != 0) {
        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dtemp += std::abs(dx[i__]);
        }
        if (n < 6) {
            return dtemp;
        }
    }
    mp1 = m + 1;
    i__2 = n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
        dtemp = dtemp + std::abs(dx[i__]) + std::abs(dx[i__ + 1]) + 
        std::abs(dx[i__ + 2]) + std::abs(dx[i__+ 3]) + std::abs(dx[i__ + 4]) +
        std::abs(dx[i__ + 5]);
    }
    return dtemp;
}


}
}
#include "blas.h"


namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(int   *   n_arg,
                      float *   da_arg,
                      float *   dx,
                      int *      incx_arg,
                      float *   dy,
                      int *      incy_arg)
{
  int i,ix,iy;
  int n=*n_arg;
  float da=*da_arg;
  int incx = *incx_arg;
  int incy = *incy_arg;

  if (n<=0)
    return;

  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*incx;
    if(incy<0)
      iy = (1-n)*incy;
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) 
      dy[iy] += da*dx[ix];

    return;

  } else {

    /* unroll */
    
    for(i=0;i<(n-4);i+=4) {
      dy[i]   += da*dx[i];
      dy[i+1] += da*dx[i+1];
      dy[i+2] += da*dx[i+2];
      dy[i+3] += da*dx[i+3];
    }
    /* continue with current value of i */
    for(;i<n;i++)
      dy[i]   += da*dx[i];
  }
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(int *n__,
                      float *dx,
                      int *incx__,
                      float *dy,
                      int *incy__)
{
    int i,ix,iy;

    int n= *n__;
    int incx = *incx__;
    int incy = *incy__;
    
    if(incx!=1 || incy!=1) 
    {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*(incx);
        if(incy<0)
            iy = (1-n)*(incy);
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            dy[iy] = dx[ix];
        
        return;
        
    } else {
        
        /* unroll */
        
        for(i=0;i<(n-8);i+=8) {
            dy[i]   = dx[i];
            dy[i+1] = dx[i+1];
            dy[i+2] = dx[i+2];
            dy[i+3] = dx[i+3];
            dy[i+4] = dx[i+4];
            dy[i+5] = dx[i+5];
            dy[i+6] = dx[i+6];
            dy[i+7] = dx[i+7];
        }
        /* continue with current value of i */
        for(;i<n;i++)
            dy[i] = dx[i];
    }
}
}
}
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
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)(const char *transa,
       const char *transb,
       int *m__,
       int *n__,
       int *k__,
       float *alpha__,
       float *a,
       int *lda__,
       float *b,
       int *ldb__,
       float *beta__,
       float *c,
       int *ldc__)
{
  const char tra=std::toupper(*transa);
  const char trb=std::toupper(*transb);
  float temp;
  int i,j,l;

  int m   = *m__;
  int n   = *n__;
  int k   = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  float alpha = *alpha__;
  float beta  = *beta__;
  
  if(m==0 || n==0 || (( std::abs(alpha)<PLUMED_GMX_FLOAT_MIN || k==0) && std::abs(beta-1.0)<PLUMED_GMX_FLOAT_EPS))
    return;

  if(std::abs(alpha)<PLUMED_GMX_FLOAT_MIN) {
    if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) {
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] = 0.0;
    } else {
      /* nonzero beta */
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] *= beta;
    }
    return;
  }

  if(trb=='N') {
    if(tra=='N') {
      
      for(j=0;j<n;j++) {
	if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( std::abs(b[ j*(ldb) + l ])>PLUMED_GMX_FLOAT_MIN) {
	    temp = alpha * b[ j*(ldb) + l ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
    } else {
      /* transpose A, but not B */
      for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[j*(ldb)+l];
	  if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
      }
    }
  } else {
    /* transpose B */
    if(tra=='N') {

      /* transpose B, but not A */

      for(j=0;j<n;j++) {
	if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( std::abs(b[ l*(ldb) + j ])>PLUMED_GMX_FLOAT_MIN) {
	    temp = alpha * b[ l*(ldb) + j ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
 
    } else {
      /* Transpose both A and B */
       for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[l*(ldb)+j];
	  if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
       }
    }
  }
}
}
}
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)(const char *trans, 
                      int *m__,
                      int *n__,
                      float *alpha__,
                      float *a,
                      int *lda__,
                      float *x,
                      int *incx__,
                      float *beta__,
                      float *y,
                      int *incy__)
{
  const char ch=std::toupper(*trans);
  int lenx,leny,kx,ky;
  int i,j,jx,jy,ix,iy;
  float temp;

  int m = *m__;
  int n = *n__;
  float alpha = *alpha__;
  float beta = *beta__;
  int incx = *incx__;
  int incy = *incy__;
  int lda = *lda__;
  
  if(n<=0 || m<=0 || (std::abs(alpha)<PLUMED_GMX_FLOAT_MIN && std::abs(beta-1.0)<PLUMED_GMX_FLOAT_EPS))
    return;

  if(ch=='N') {
    lenx = n;
    leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  
   if(incx>0)
    kx = 1;
  else
    kx = 1 - (lenx -1)*(incx);

  if(incy>0)
    ky = 1;
  else
    ky = 1 - (leny -1)*(incy);
 
  if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS) {
    if(incy==1) {
      if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	for(i=0;i<leny;i++)
	  y[i] = 0.0;
      else
	for(i=0;i<leny;i++)
	  y[i] *= beta;
    } else {
      /* non-unit incr. */
      iy = ky;
      if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) 
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] = 0.0;
      else
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] *= beta;
    }
  }
  
  if(std::abs(alpha)<PLUMED_GMX_FLOAT_MIN)
    return;
  
  if(ch=='N') {
    jx = kx;
    if(incy==1) {
      for(j=1;j<=n;j++,jx+=incx) 
	if( std::abs(x[jx-1])>PLUMED_GMX_FLOAT_MIN) {
	  temp = alpha * x[jx-1];
	  for(i=1;i<=m;i++)
	    y[i-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jx+=incx) 
	if( std::abs(x[jx-1])>PLUMED_GMX_FLOAT_MIN) {
	  temp = alpha * x[jx-1];
	  iy = ky;
	  for(i=1;i<=m;i++,iy+=incy)
	    y[iy-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    }
  } else {
    /* transpose */
    jy = ky;
    if(incx==1) {
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	for(i=1;i<=m;i++)
	  temp += a[(j-1)*(lda)+(i-1)] * x[i-1];
	y[jy-1] += alpha * temp;
      }
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	ix = kx;
	for(i=1;i<=m;i++,ix+=incx)
	  temp += a[(j-1)*(lda)+(i-1)] * x[ix-1];
	y[jy-1] += alpha * temp;
      }
    }
  }
} 
   
}
}
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(sger,SGER)(int *m__,
                    int *n__,
                    float *alpha__,
                    float *x,
                    int *incx__,
                    float *y,
                    int *incy__,
                    float *a,
                    int *lda__)
{
    int ix,kx,jy;
    int i,j;
    float temp;
    
    int m = *m__;
    int n = *n__;
    int incx = *incx__;
    int incy = *incy__;
    int lda = *lda__;
    float alpha = *alpha__;
    
    if(m<=0 || n<=0 || std::abs(alpha)<PLUMED_GMX_FLOAT_MIN)
        return;
    
    if(incy>0)
        jy = 0;
    else
        jy = incy * (1 - n);
    
    if(incx==1) {
        for(j=0;j<n;j++,jy+=incy)
            if(std::abs(y[jy])>PLUMED_GMX_FLOAT_MIN) {
                temp = alpha * y[jy];
                for(i=0;i<m;i++)
                    a[j*(lda)+i] += temp*x[i];
            }
    } else {
        /* non-unit incx */
        if(incx>0) 
            kx = 0;
        else
            kx = incx * (1 - m);
        
        for(j=0;j<n;j++,jy+=incy) {
            if(std::abs(y[jy])>PLUMED_GMX_FLOAT_MIN) {
                temp = alpha * y[jy];
                ix = kx;
                for(i=0;i<m;i++,ix+=incx)
                    a[j*(lda)+i] += temp*x[ix];
            }
        }
    }
        return;
}
}
}
#include <cmath>


#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
float
PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(int  *     n__,
                      float *    x,
                      int    *    incx__)
{
    int ix,max_ix;
    float ssq,scale,absxi,t;
    
    int n = *n__;
    int incx = *incx__;
    
    if(n<1 || incx<1)
        return 0;
    else if (n==1) {
        t = x[0];
        if(t>=0)
            return t;
        else 
            return -t;
    }
    
    scale = 0.0;
    ssq   = 1.0;
    
    max_ix = 1+(n-1)*(incx);
    for(ix=1;ix<=max_ix;ix+=incx) {
        t = x[ix-1];
        if(std::abs(t)>PLUMED_GMX_FLOAT_MIN) {
            absxi = (t>=0) ? t : (-t);
            if(scale<absxi) {
                t = scale/absxi;
                t = t*t;
                ssq = ssq*t + 1.0;
                scale = absxi;
            } else {
                t = absxi/scale;
                ssq += t*t;
            }
        }
    }
    return scale*std::sqrt(ssq);
    
}


 
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(srot,SROT)(int *n__,
                    float *dx,
                    int *incx__,
                    float *dy,
                    int *incy__,
                    float *c__,
                    float *s__)
{
  int i,ix,iy;
  float dtemp;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  float c = *c__;
  float s = *s__;
  
  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*(incx);
    if(incy<0)
      iy = (1-n)*(incy);
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      dtemp  = (c) * dx[ix] + (s) * dy[iy];
      dy[iy] = (c) * dy[iy] - (s) * dx[ix];
      dx[ix] = dtemp;
    }

    return;

  } else {

    /* unit increments */   
    for(i=0;i<n;i++) {
      dtemp = (c) * dx[i] + (s) * dy[i];
      dy[i] = (c) * dy[i] - (s) * dx[i];
      dx[i] = dtemp;      
    }

  }
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(int  *n__,
                      float *fact__,
                      float *dx,
                      int   *incx__)
{
    int nincx,i;

    int n = *n__;
    float fact = *fact__;
    int incx = *incx__;
    
    if(n<=0 || incx<=0)
        return;
    
    if(incx==1) {
        /* Unrool factor 5 */
        for(i=0;i<(n-5);i+=5) {
            dx[i]   *= fact;
            dx[i+1] *= fact;
            dx[i+2] *= fact;
            dx[i+3] *= fact;
            dx[i+4] *= fact;
        }    
        /* continue with current value of i */
        for(;i<n;i++)
            dx[i]   *= fact;
        
        return;
    } else {
        /* inc != 1 */
        nincx = n * (incx);
        for (i=0;i<nincx;i+=incx)
            dx[i] *= fact;
        
        return;
    } 
    
}
}
}
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(int *n__,
                      float *dx,
                      int *incx__,
                      float *dy,
                      int *incy__)
{
  int i,ix,iy;
  float d1,d2,d3;

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
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)(const char *uplo,
                      int *n__,
                      float *alpha__,
                      float *a,
                      int *lda__,
                      float *x,
                      int *incx__,
                      float *beta__,
                      float *y,
                      int *incy__)
{
    const char ch=std::toupper(*uplo);
    int kx,ky,i,j,ix,iy,jx,jy;
    float temp1,temp2;
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    float alpha = *alpha__;
    float beta  = *beta__;
    
    if(n<=0 || incx==0 || incy==0)
        return;
    
    if(incx>0)
        kx = 1;
    else
        kx = 1 - (n -1)*(incx);
    
    if(incy>0)
        ky = 1;
    else
        ky = 1 - (n -1)*(incy);
    
    if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS) {
        if(incy==1) {
            if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) 
                for(i=1;i<=n;i++)
                    y[i-1] = 0.0;
            else
                for(i=1;i<=n;i++)
                    y[i-1] *= beta;
        } else {
            /* non-unit incr. */
            iy = ky;
            if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) 
                for(i=1;i<=n;i++) {
                    y[iy-1] = 0.0;
                    iy += incy;
                }
                    else
                        for(i=1;i<=n;i++) {
                            y[iy-1] *= beta;
                            iy += incy;
                        }
        }
    }
        
        if(std::abs(alpha)<PLUMED_GMX_FLOAT_MIN) 
            return;
        
        if(ch=='U') {
            if(incx==1 && incy==1) {
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[j-1];
                    temp2 = 0.0;
                    for(i=1;i<j;i++) {
                        y[i-1] += temp1*a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[i-1];
                    }
                    y[j-1] += temp1*a[(j-1)*(lda)+(j-1)] + alpha *temp2;
                }
            } else {
                /* non-unit incr. */
                jx = kx;
                jy = ky;
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[jx-1];
                    temp2 = 0.0;
                    ix = kx;
                    iy = ky;
                    for(i=1;i<j;i++) {
                        y[iy-1] += temp1 * a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[ix-1];
                        ix += incx;
                        iy += incy;
                    }
                    y[jy-1] += temp1*a[(j-1)*(lda)+(j-1)] + alpha*temp2;
                    jx += incx;
                    jy += incy;
                }
            }
        } else {
            /* lower */
            if(incx==1 && incy==1) {
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[j-1];
                    temp2 = 0.0;
                    y[j-1] += temp1 * a[(j-1)*(lda)+(j-1)];
                    for(i=j+1;i<=n;i++) {
                        y[i-1] += temp1*a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[i-1];
                    }
                    y[j-1] += alpha *temp2;
                }
            } else {
                /* non-unit incr. */
                jx = kx;
                jy = ky;
                for(j=1;j<=n;j++) {
                    temp1 = alpha * x[jx-1];
                    temp2 = 0.0;
                    y[jy-1] += temp1 * a[(j-1)*(lda)+(j-1)];
                    ix = jx;
                    iy = jy;
                    for(i=j+1;i<=n;i++) {
                        ix += incx;
                        iy += incy;
                        y[iy-1] += temp1 * a[(j-1)*(lda)+(i-1)];
                        temp2 += a[(j-1)*(lda)+(i-1)] * x[ix-1];
                    }
                    y[jy-1] += alpha*temp2;
                    jx += incx;
                    jy += incy;
                }
            }
        }
        return;
}    
}
}
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)(const char *    uplo,
                      int *     n__,
                      float *  alpha__,
                      float *  x,
                      int *     incx__,
                      float *  y,
                      int *     incy__,
                      float *  a,
                      int *     lda__)
{
    int kx,ky,ix,iy,jx,jy,j,i;
    float temp1,temp2;
    const char ch=std::toupper(*uplo);
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    float alpha = *alpha__;
    
    if(n<=0 || std::abs(alpha)<PLUMED_GMX_FLOAT_MIN || incx==0 || incy==0 ||
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
                if( std::abs(x[j-1])>PLUMED_GMX_FLOAT_MIN  || std::abs(y[j-1])>PLUMED_GMX_FLOAT_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=1;i<=j;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>PLUMED_GMX_FLOAT_MIN || std::abs(y[jy-1])>PLUMED_GMX_FLOAT_MIN ) {
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
                if( std::abs(x[j-1])>PLUMED_GMX_FLOAT_MIN  || std::abs(y[j-1])>PLUMED_GMX_FLOAT_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=j;i<=n;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>PLUMED_GMX_FLOAT_MIN || std::abs(y[jy-1])>PLUMED_GMX_FLOAT_MIN ) {
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
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(ssyr2k,SSYR2K)(const char *uplo, 
                        const char *trans,
                        int *n__,
                        int *k__,
                        float *alpha__,
                        float *a,
                        int *lda__,
                        float *b,
                        int *ldb__,
                        float *beta__,
                        float *c,
                        int *ldc__)
{
  char ch1,ch2;
  int i,j,l;
  float temp1,temp2;

  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  float alpha = *alpha__;
  float beta  = *beta__;
  
  ch1 = std::toupper(*uplo);
  ch2 = std::toupper(*trans);

  if(n==0 || ( ( std::abs(alpha)<PLUMED_GMX_FLOAT_MIN || k==0 ) && std::abs(beta-1.0)<PLUMED_GMX_FLOAT_EPS))
    return;

  if(std::abs(alpha)<PLUMED_GMX_FLOAT_MIN) {
    if(ch1=='U') {
      if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    } else {
      /* lower */
      if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    }
    return;
  }

  if(ch2=='N') {
    if(ch1=='U') {
      for(j=1;j<=n;j++) {
	if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	  for(i=1;i<=j;i++)
	     c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS)
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( std::abs(a[(l-1)*(lda)+(j-1)])>PLUMED_GMX_FLOAT_MIN ||
	      std::abs(b[(l-1)*(ldb)+(j-1)])>PLUMED_GMX_FLOAT_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=1;i<=j;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    } else {
      /* lower */
      for(j=1;j<=n;j++) {
	if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(std::abs(beta-1.0)>PLUMED_GMX_FLOAT_EPS)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( std::abs(a[(l-1)*(lda)+(j-1)])>PLUMED_GMX_FLOAT_MIN ||
	      std::abs(b[(l-1)*(ldb)+(j-1)])>PLUMED_GMX_FLOAT_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=j;i<=n;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    }
  } else {
    /* transpose */
    if(ch1=='U') {
      for(j=1;j<=n;j++) 
	for(i=1;i<=j;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    } else {
      /* lower */
      for(j=1;j<=n;j++) 
	for(i=j;i<=n;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(std::abs(beta)<PLUMED_GMX_FLOAT_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    }
  }
  return;
}
}
}
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(strmm,STRMM)(const char *side, 
                      const char *uplo, 
                      const char *transa, 
                      const char *diag, 
                      int *m__, 
                      int *n__, 
                      float *alpha__, 
                      float *a, 
                      int *lda__, 
                      float *b, 
                      int *ldb__)
{
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    
    int m = *m__;
    int n = *n__;
    int lda = *lda__;
    int ldb = *ldb__;
    float alpha = *alpha__;
    
    /* Local variables */
    int i__, j, k;
    float temp;
    int lside;
    int upper;
    int nounit;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = (*side=='L' || *side=='l');

    nounit = (*diag=='N' || *diag=='n');
    upper = (*uplo=='U' || *uplo=='u');

    if (n == 0) {
	return;
    }
    if (std::abs(alpha)<PLUMED_GMX_FLOAT_MIN) {
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
	    }
	}
	return;
    }
    if (lside) {
	if (*transa=='N' || *transa=='n') {
	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (k = 1; k <= i__2; ++k) {
			if ( std::abs(b[k + j * b_dim1])>PLUMED_GMX_FLOAT_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = m; k >= 1; --k) {
			if (std::abs(b[k + j * b_dim1])>PLUMED_GMX_FLOAT_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (*transa=='N' || *transa=='n') {

	    if (upper) {
		for (j = n; j >= 1; --j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if ( std::abs(a[k + j * a_dim1])>PLUMED_GMX_FLOAT_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__2 = n;
		    for (k = j + 1; k <= i__2; ++k) {
			if ( std::abs(a[k + j * a_dim1])>PLUMED_GMX_FLOAT_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if ( std::abs(a[j + k * a_dim1])>PLUMED_GMX_FLOAT_MIN ) {
			    temp = alpha * a[j + k * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if ( std::abs(temp-1.0)>PLUMED_GMX_FLOAT_EPS) {
			i__2 = m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    } else {
		for (k = n; k >= 1; --k) {
		    i__1 = n;
		    for (j = k + 1; j <= i__1; ++j) {
			if ( std::abs(a[j + k * a_dim1])>PLUMED_GMX_FLOAT_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if ( std::abs(temp-1.0)>PLUMED_GMX_FLOAT_EPS) {
			i__1 = m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
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
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void 
PLUMED_BLAS_F77_FUNC(strmv,STRMV)(const char *uplo, 
                      const char *trans,
                      const char *diag, 
                      int *n__, 
                      float *a, 
                      int *lda__, 
                      float *x, 
                      int *incx__)
{
    int a_dim1, a_offset, i__1, i__2;

    int i__, j, ix, jx, kx;
    float temp;
    int nounit;
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    if (n == 0) {
	return;
    }

    nounit = (*diag=='n' || *diag=='N');

    if (incx <= 0) {
	kx = 1 - (n - 1) * incx;
    } else {
	kx = 1;
    }

    if (*trans=='N' || *trans=='n') {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (std::abs(x[j])>PLUMED_GMX_FLOAT_MIN) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (std::abs(x[jx])>PLUMED_GMX_FLOAT_MIN) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix += incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += incx;
		}
	    }
	} else {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    if (std::abs(x[j])>PLUMED_GMX_FLOAT_MIN) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		kx += (n - 1) * incx;
		jx = kx;
		for (j = n; j >= 1; --j) {
		    if (std::abs(x[jx])>PLUMED_GMX_FLOAT_MIN) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix -= incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= incx;
		}
	    }
	}
    } else {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx + (n - 1) * incx;
		for (j = n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx -= incx;
		}
	    }
	} else {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx += incx;
		}
	    }
	}
    }

    return;

}


}
}
#include <cctype>
#include <cmath>

#include "real.h"
#include "blas.h"

namespace PLMD{
namespace blas{
void
PLUMED_BLAS_F77_FUNC(strsm,STRSM)(const char * side,
                      const char * uplo,
                      const char * transa,
                      const char * diag,
                      int *  m__,
                      int *  n__,
                      float *alpha__,
                      float *a,
                      int *  lda__,
                      float *b,
                      int *  ldb__)
{
  const char xside  = std::toupper(*side);
  const char xuplo  = std::toupper(*uplo);
  const char xtrans = std::toupper(*transa);
  const char xdiag  = std::toupper(*diag);
  int i,j,k;
  float temp;

  int m = *m__;
  int n = *n__;
  int lda = *lda__;
  int ldb = *ldb__;
  float alpha = *alpha__;
  
  if(n<=0)
    return;

  
  if(std::abs(alpha)<PLUMED_GMX_FLOAT_MIN) { 
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS) {
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  }
	  for(k=m-1;k>=0;k--) {
	    if( std::abs(b[j*(ldb)+k])>PLUMED_GMX_FLOAT_MIN) {
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<m;k++) {
	    if( std::abs(b[j*(ldb)+k])>PLUMED_GMX_FLOAT_MIN) {
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<j;k++) {
	    if( std::abs(a[j*(lda)+k])>PLUMED_GMX_FLOAT_MIN) {
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
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=j+1;k<n;k++) {
	    if( std::abs(a[j*(lda)+k])>PLUMED_GMX_FLOAT_MIN ) {
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
	    if( std::abs(a[k*(lda)+j])>PLUMED_GMX_FLOAT_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS)
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
	    if( std::abs(a[k*(lda)+j])>PLUMED_GMX_FLOAT_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(std::abs(alpha-1.0)>PLUMED_GMX_FLOAT_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      }      
    }
  }    
}
}
}
#endif

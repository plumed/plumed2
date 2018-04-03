/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib LAPACK library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference LAPACK implementation.

The reference LAPACK implementation is available from http://www.netlib.org/lapack 

LAPACK does not come with a formal named "license", but a general statement saying:

"The reference LAPACK is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified LAPACK files as the original netlib versions, so do what you want with them.

However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)
#include <cctype>
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)(const char *uplo, 
	const char *compq, 
	int *n,
	double *d__, 
	double *e, 
	double *u, 
	int *ldu,
	double *vt, 
	int *ldvt,
	double *q,
	int *iq,
	double *work, 
	int *iwork, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    int i__, j, k;
    double p, r__;
    int z__, ic, ii, kk;
    double cs;
    int is, iu;
    double sn;
    int nm1;
    double eps;
    int ivt, difl, difr, ierr, perm, mlvl, sqre;
    int poles, iuplo, nsize, start;
    int givcol;
    int icompq;
    double orgnrm;
    int givnum, givptr, qstart, smlsiz, wstart, smlszp;
    double zero = 0.0;
    double one = 1.0;
    int c_0 = 0;
    int c_1 = 1;

    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --q;
    --iq;
    --work;
    --iwork;

    k = iu = z__ = ic = is = ivt = difl = difr = perm = 0;
    poles = givnum = givptr = givcol = 0;
    
    smlsiz = DBDSDC_SMALLSIZE;
    *info = 0;

    iuplo = (*uplo=='U' || *uplo=='u') ? 1 : 2;

    switch(*compq) {
    case 'n':
    case 'N':
      icompq = 0;
      break;
    case 'p':
    case 'P':
      icompq = 1;
      break;
    case 'i':
    case 'I':
      icompq = 2;
      break;
    default:
      return;
    }

    if (*n <= 0) 
	return;
    
    if (*n == 1) {
	if (icompq == 1) {
	  q[1] = (d__[1]>0) ? 1.0 : -1.0;
	  q[smlsiz * *n + 1] = 1.;
	} else if (icompq == 2) {
	  u[u_dim1 + 1] = (d__[1]>0) ? 1.0 : -1.0;
	  vt[vt_dim1 + 1] = 1.;
	}
	d__[1] = std::abs(d__[1]);
	return;
    }
    nm1 = *n - 1;
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &d__[1], &c_1, &q[1], &c_1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &e[1], &c_1, &q[*n + 1], &c_1);
    }
    if (iuplo == 2) {
	qstart = 5;
	wstart = (*n << 1) - 1;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (icompq == 1) {
		q[i__ + (*n << 1)] = cs;
		q[i__ + *n * 3] = sn;
	    } else if (icompq == 2) {
		work[i__] = cs;
		work[nm1 + i__] = -sn;
	    }
	}
    }
    if (icompq == 0) {
      PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U",&c_0,n,&c_0,&c_0,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
	      &u[u_offset], ldu, &u[u_offset], ldu, &work[wstart], info);
	goto L40;
    }
    if (*n <= smlsiz) {
	if (icompq == 2) {
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U",&c_0,n,n,n,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
		    &u[u_offset],ldu,&u[u_offset],ldu,&work[wstart],info);
	} else if (icompq == 1) {
	    iu = 1;
	    ivt = iu + *n;
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &q[iu + (qstart - 1) * *n], n);
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &q[ivt + (qstart - 1) * *n], n);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &c_0, n, n, n, &c_0, &d__[1], &e[1], 
		    &q[ivt + (qstart - 1) * *n], n, &q[iu + (qstart - 1) * *n], 
		    n, &q[iu + (qstart - 1) * *n], n, &work[wstart], info);
	}
	goto L40;
    }

    if (icompq == 2) {
	PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
    }

    orgnrm = PLUMED_BLAS_F77_FUNC(dlanst,DLANST)("M", n, &d__[1], &e[1]);
    if ( std::abs(orgnrm)<PLUMED_GMX_DOUBLE_MIN) {
	return;
    }
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &orgnrm, &one, n, &c_1, &d__[1], n, &ierr);
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &orgnrm, &one, &nm1, &c_1, &e[1], &nm1, &ierr);

    eps = PLUMED_GMX_DOUBLE_EPS;

    mlvl = (int) (std::log((double) (*n) / (double) (smlsiz + 1)) / 
		  std::log(2.)) + 1;
    smlszp = smlsiz + 1;

    if (icompq == 1) {
	iu = 1;
	ivt = smlsiz + 1;
	difl = ivt + smlszp;
	difr = difl + mlvl;
	z__ = difr + (mlvl << 1);
	ic = z__ + mlvl;
	is = ic + 1;
	poles = is + 1;
	givnum = poles + (mlvl << 1);

	k = 1;
	givptr = 2;
	perm = 3;
	givcol = perm + mlvl;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(d__[i__]) < eps) 
	    d__[i__] = (d__[i__]>0) ? eps : -eps;
    }

    start = 1;
    sqre = 0;

    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__]) < eps || i__ == nm1) {
	    if (i__ < nm1) {
		nsize = i__ - start + 1;
	    } else if (std::abs(e[i__]) >= eps) {
		nsize = *n - start + 1;
	    } else {
		nsize = i__ - start + 1;
		if (icompq == 2) {
		    u[*n + *n * u_dim1] = (d__[*n]>0) ? 1.0 : -1.0; 
		    vt[*n + *n * vt_dim1] = 1.;
		} else if (icompq == 1) {
		    q[*n + (qstart - 1) * *n] = (d__[*n]>0) ? 1.0 : -1.0; 
		    q[*n + (smlsiz + qstart - 1) * *n] = 1.;
		}
		d__[*n] = std::abs(d__[*n]);
	    }
	    if (icompq == 2) {
		PLUMED_BLAS_F77_FUNC(dlasd0,DLASD0)(&nsize, &sqre, &d__[start], &e[start], 
			&u[start + start * u_dim1], ldu, 
			&vt[start + start * vt_dim1], 
			ldvt, &smlsiz, &iwork[1], &work[wstart], info);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlasda,DLASDA)(&icompq, &smlsiz, &nsize, &sqre, &d__[start], 
			&e[start], &q[start + (iu + qstart - 2) * *n], n, 
			&q[start + (ivt + qstart - 2) * *n], &iq[start + k * *n],
			&q[start + (difl + qstart - 2) * *n], 
			&q[start + (difr + qstart - 2) * *n], 
			&q[start + (z__ + qstart - 2) * *n], 
			&q[start + (poles + qstart - 2) * *n], 
			&iq[start + givptr * *n], &iq[start + givcol * *n], n, 
			&iq[start + perm * *n], 
			&q[start + (givnum + qstart - 2) * *n], 
			&q[start + (ic + qstart - 2) * *n], 
			&q[start + (is + qstart - 2) * *n], &work[wstart], 
			&iwork[1], info);
		if (*info != 0) {
		    return;
		}
	    }
	    start = i__ + 1;
	}
    }
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &one, &orgnrm, n, &c_1, &d__[1], n, &ierr);
L40:
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	kk = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] > p) {
		kk = j;
		p = d__[j];
	    }
	}
	if (kk != i__) {
	    d__[kk] = d__[i__];
	    d__[i__] = p;
	    if (icompq == 1) {
		iq[i__] = kk;
	    } else if (icompq == 2) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &u[i__ * u_dim1 + 1],&c_1,&u[kk*u_dim1+1],&c_1);
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &vt[i__ + vt_dim1], ldvt, &vt[kk + vt_dim1], ldvt);
	    }
	} else if (icompq == 1) {
	    iq[i__] = i__;
	}
    }
    if (icompq == 1) {
	if (iuplo == 1) {
	    iq[*n] = 1;
	} else {
	    iq[*n] = 0;
	}
    }
    if (iuplo == 2 && icompq == 2) {
	PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "B", n, n, &work[1], &work[*n], &u[u_offset], ldu);
    }

    return;
}
}
}
#include <cctype>
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dbdsqr,DBDSQR)(const char *uplo,
                        int *n,
                        int *ncvt,
                        int *nru, 
                        int *ncc, 
                        double *d__,
                        double *e,
                        double *vt, 
                        int *ldvt,
                        double *u, 
                        int *ldu,
                        double *c__, 
                        int *ldc,
                        double *work,
                        int *info)
{
    const char xuplo = std::toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    double r__1, r__2, r__3, r__4;
    double c_b15 = -.125;

    int c__1 = 1;
    double c_b49 = 1.f;
    double c_b72 = -1.f;

    double f, g, h__;
    int i__, j, m;
    double r__, cs;
    int ll;
    double sn, mu;
    int nm1, nm12, nm13, lll;
    double eps, sll, tol, abse;
    int idir;
    double abss;
    int oldm;
    double cosl;
    int isub, iter;
    double unfl, sinl, cosr, smin, smax, sinr;
    double oldcs;
    int oldll;
    double shift, sigmn, oldsn = 0.;
    int maxit;
    double sminl;
    double sigmx;
    int lower;
    double sminoa;
    double thresh;
    int rotate;
    double tolmul;
    int itmp1,itmp2;
    
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    
    lower = (xuplo == 'L');
    if ( (xuplo!='U') && !lower) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ncvt < 0) {
	*info = -3;
    } else if (*nru < 0) {
	*info = -4;
    } else if (*ncc < 0) {
	*info = -5;
    } else if ( ((*ncvt == 0) && (*ldvt < 1)) || ((*ncvt > 0) && (*ldvt < itmp1)) ) {
	*info = -9;
    } else if (*ldu < itmp2) {
	*info = -11;
    } else if ( ((*ncc == 0) && (*ldc < 1)) || ((*ncc > 0) && (*ldc < itmp1))) {
	*info = -13;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	goto L160;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

    if (! rotate) {
	PLUMED_BLAS_F77_FUNC(dlasq1,DLASQ1)(n, &d__[1], &e[1], &work[1], info);
	return;
    }

    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;

    eps = PLUMED_GMX_DOUBLE_EPS;
    unfl = PLUMED_GMX_DOUBLE_MIN/PLUMED_GMX_DOUBLE_EPS;

    if (lower) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    work[i__] = cs;
	    work[nm1 + i__] = sn;
	}

	if (*nru > 0) {
	    PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, n, &work[1], &work[*n], &u[u_offset], 
		    ldu);
	}
	if (*ncc > 0) {
	    PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", n, ncc, &work[1], &work[*n], &c__[c_offset],
		     ldc);
	}
    }

    r__3 = 100.f, r__4 = std::pow(PLUMED_GMX_DOUBLE_EPS,c_b15);
    r__1 = 10.f, r__2 = (r__3<r__4) ? r__3 : r__4;
    tolmul = (r__1>r__2) ? r__1 : r__2;
    tol = tolmul * eps;
    smax = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = d__[i__], std::abs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = e[i__], std::abs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    sminl = 0.f;
    if (tol >= 0.f) {
	sminoa = std::abs(d__[1]);
	if (sminoa == 0.f) {
	    goto L50;
	}
	mu = sminoa;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mu = (r__2 = d__[i__], std::abs(r__2)) * (mu / (mu + (r__1 = e[i__ - 
		    1], std::abs(r__1))));
	    sminoa = (sminoa<mu) ? sminoa : mu;
	    if (sminoa == 0.f) {
		goto L50;
	    }
	}
L50:
	sminoa /=  std::sqrt((double) (*n));
	r__1 = tol * sminoa, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    } else {
	r__1 = std::abs(tol) * smax, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    }
    maxit = *n * 6 * *n;
    iter = 0;
    oldll = -1;
    oldm = -1;
    m = *n;

L60:

    if (m <= 1) {
	goto L160;
    }
    if (iter > maxit) {
	goto L200;
    }

    if (tol < 0.f && (r__1 = d__[m], std::abs(r__1)) <= thresh) {
	d__[m] = 0.f;
    }
    smax = (r__1 = d__[m], std::abs(r__1));
    smin = smax;
    i__1 = m - 1;
    for (lll = 1; lll <= i__1; ++lll) {
	ll = m - lll;
	abss = (r__1 = d__[ll], std::abs(r__1));
	abse = (r__1 = e[ll], std::abs(r__1));
	if (tol < 0.f && abss <= thresh) {
	    d__[ll] = 0.f;
	}
	if (abse <= thresh) {
	    goto L80;
	}
	smin = (smin<abss) ? smin : abss;
	r__1 = (smax>abss) ? smax : abss;
	smax = (r__1>abse) ? r__1 : abse;
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.f;
    if (ll == m - 1) {
	--m;
	goto L60;
    }
L90:
    ++ll;
    if (ll == m - 1) {
	PLUMED_BLAS_F77_FUNC(dlasv2,DLASV2)(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
	d__[m - 1] = sigmx;
	e[m - 1] = 0.f;
	d__[m] = sigmn;
	if (*ncvt > 0) {
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
	}
	if (*nru > 0) {
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
	}
	if (*ncc > 0) {
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
	}
	m += -2;
	goto L60;
    }
    if (ll > oldm || m < oldll) {
	if ((r__1 = d__[ll], std::abs(r__1)) >= (r__2 = d__[m], std::abs(r__2))) {
	    idir = 1;
	} else {
	    idir = 2;
	}
    }
    if (idir == 1) {

        if( (std::abs(e[m-1]) <= std::abs(tol) * std::abs(d__[m])) ||
            (tol<0.0 && std::abs(e[m-1])<=thresh)) {
	    e[m - 1] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[ll], std::abs(r__1));
	    sminl = mu;
	    i__1 = m - 1;
	    for (lll = ll; lll <= i__1; ++lll) {
		if ((r__1 = e[lll], std::abs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		mu = (r__2 = d__[lll + 1], std::abs(r__2)) * (mu / (mu + (r__1 = 
			e[lll], std::abs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    } else {
        if( (std::abs(e[ll]) <= std::abs(tol)*std::abs(d__[ll])) ||
            (tol<0.0 && std::abs(e[ll])<=thresh)) {
	    e[ll] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[m], std::abs(r__1));
	    sminl = mu;
	    i__1 = ll;
	    for (lll = m - 1; lll >= i__1; --lll) {
		if ((r__1 = e[lll], std::abs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		mu = (r__2 = d__[lll], std::abs(r__2)) * (mu / (mu + (r__1 = e[
			lll], std::abs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    }
    oldll = ll;
    oldm = m;

    r__1 = eps, r__2 = tol * .01f;
    if (tol >= 0.f && *n * tol * (sminl / smax) <= ((r__1>r__2) ? r__1 : r__2)) {
	shift = 0.f;
    } else {
	if (idir == 1) {
	    sll = (r__1 = d__[ll], std::abs(r__1));
	    PLUMED_BLAS_F77_FUNC(dlas2,DLAS2)(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
	} else {
	    sll = (r__1 = d__[m], std::abs(r__1));
	    PLUMED_BLAS_F77_FUNC(dlas2,DLAS2)(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
	}
	if (sll > 0.f) {
	    r__1 = shift / sll;
	    if (r__1 * r__1 < eps) {
		shift = 0.f;
	    }
	}
    }
    iter = iter + m - ll;
    if (shift == 0.f) {
	if (idir == 1) {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		r__1 = d__[i__] * cs;
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&r__1, &e[i__], &cs, &sn, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ + 1] * sn;
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll + 1] = cs;
		work[i__ - ll + 1 + nm1] = sn;
		work[i__ - ll + 1 + nm12] = oldcs;
		work[i__ - ll + 1 + nm13] = oldsn;
	    }
	    h__ = d__[m] * cs;
	    d__[m] = h__ * oldcs;
	    e[m - 1] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], std::abs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		r__1 = d__[i__] * cs;
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&r__1, &e[i__ - 1], &cs, &sn, &r__);
		if (i__ < m) {
		    e[i__] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ - 1] * sn;
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll] = cs;
		work[i__ - ll + nm1] = -sn;
		work[i__ - ll + nm12] = oldcs;
		work[i__ - ll + nm13] = -oldsn;
	    }
	    h__ = d__[ll] * cs;
	    d__[ll] = h__ * oldcs;
	    e[ll] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[ll], std::abs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	}
    } else {

	if (idir == 1) {
	    f = ((r__1 = d__[ll], std::abs(r__1)) - shift) * ( ((d__[ll] > 0) ? c_b49 : -c_b49) + shift / d__[ll]);
	    g = e[ll];
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__];
		e[i__] = cosr * e[i__] - sinr * d__[i__];
		g = sinr * d__[i__ + 1];
		d__[i__ + 1] = cosr * d__[i__ + 1];
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__] + sinl * d__[i__ + 1];
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
		if (i__ < m - 1) {
		    g = sinl * e[i__ + 1];
		    e[i__ + 1] = cosl * e[i__ + 1];
		}
		work[i__ - ll + 1] = cosr;
		work[i__ - ll + 1 + nm1] = sinr;
		work[i__ - ll + 1 + nm12] = cosl;
		work[i__ - ll + 1 + nm13] = sinl;
	    }
	    e[m - 1] = f;

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], std::abs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {

	    f = ((r__1 = d__[m], std::abs(r__1)) - shift) * ( ((d__[m] > 0) ? c_b49 : -c_b49) + shift / d__[m]);
	    g = e[m - 1];
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ < m) {
		    e[i__] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__ - 1];
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
		g = sinr * d__[i__ - 1];
		d__[i__ - 1] = cosr * d__[i__ - 1];
		PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
		if (i__ > ll + 1) {
		    g = sinl * e[i__ - 2];
		    e[i__ - 2] = cosl * e[i__ - 2];
		}
		work[i__ - ll] = cosr;
		work[i__ - ll + nm1] = -sinr;
		work[i__ - ll + nm12] = cosl;
		work[i__ - ll + nm13] = -sinl;
	    }
	    e[ll] = f;

	    if ((r__1 = e[ll], std::abs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	}
    }

    goto L60;

L160:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] < 0.f) {
	    d__[i__] = -d__[i__];

	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
	    }
	}
    }

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = 1;
	smin = d__[1];
	i__2 = *n + 1 - i__;
	for (j = 2; j <= i__2; ++j) {
	    if (d__[j] <= smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != *n + 1 - i__) {
	    d__[isub] = d__[*n + 1 - i__];
	    d__[*n + 1 - i__] = smin;
	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
	    }
	    if (*ncc > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
	    }
	}
    }
    goto L220;

L200:
    *info = 0;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    ++(*info);
	}
    }
L220:
    return;

}


}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgebd2,DGEBD2)(int *m,
	int *n,
	double *a,
	int *lda,
	double *d,
	double *e,
	double *tauq,
	double *taup,
	double *work,
	int *info)
{
  int i,i1,i2,i3;
    
    *info = 0;

  if(*m>=*n) {
    /* reduce to upper bidiag. form */
    for(i=0;i<*n;i++) {
      i1 = *m - i;
      i2 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      i3 = 1;
      PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*n-1)) {

	i1 = *n - i -1;
	i2 = ( (i+2) < (*n-1)) ? (i+2) : (*n-1); 
	PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i1,&(a[(i+1)*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));

	e[i] = a[(i+1)*(*lda)+i];
	a[(i+1)*(*lda)+i] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("R",&i1,&i2,&(a[(i+1)*(*lda)+i]),lda,&(taup[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i+1)*(*lda)+i] = e[i];
      } else
	taup[i] = 0.0;
    }
  } else {
    /* reduce to lower bidiag. form */
    for(i=0;i<*m;i++) {
      i1 = *n - i;
      i2 = ( (i+1) < (*n-1)) ? (i+1) : (*n-1);
      PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;

      i2 = *m - i - 1;
      i3 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("R",&i2,&i1,&(a[i*(*lda)+i]),lda,&(taup[i]),&(a[(i)*(*lda)+i3]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*m-1)) {

	i1 = *m - i - 1;
	i2 = ( (i+2) < (*m-1)) ? (i+2) : (*m-1);
	i3 = 1;
	PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i1,&(a[(i)*(*lda)+i+1]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));

	e[i] = a[(i)*(*lda)+i+1];
	a[(i)*(*lda)+i+1] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	i3 = 1;
	PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[(i)*(*lda)+i+1]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i)*(*lda)+i+1] = e[i];
      } else
	tauq[i] = 0.0;
    }
  }
  return;
}
}
}
#include "lapack.h"
#include "blas/blas.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(int *m, 
	int *n, 
	double *a, 
	int *lda, 
	double *d__, 
	double *e,
	double *tauq, 
	double *taup,
	double *work, 
	int *lwork,
	int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i_1, i_2, i_3, i_4;

    /* Local variables */
    int i_, j, nx,nb;
    double ws;
    int nbmin, iinfo, minmn;
    int ldwrkx, ldwrky;
    double one = 1.0;
    double minusone = -1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;

    nb = DGEBRD_BLOCKSIZE;
    *info = 0;
    if (*lwork==-1) {
      work[1] = (double) ( (*m + *n) * nb);
      return;
    }
    minmn = (*m < *n) ? *m : *n;
    if (minmn == 0) {
      work[1] = 1.;
      return;
    }

    ws = (*m > *n) ? *m : *n;
    ldwrkx = *m;
    ldwrky = *n;

    if (nb > 1 && nb < minmn) {
	nx = DGEBRD_CROSSOVER;
	if (nx < minmn) {
	    ws = (double) ((*m + *n) * nb);
	    if ((double) (*lwork) < ws) {
	      nbmin = DGEBRD_MINBLOCKSIZE;
		if (*lwork >= (*m + *n) * nbmin) {
		    nb = *lwork / (*m + *n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }

    i_1 = minmn - nx;
    i_2 = nb;
    for (i_ = 1; i_2 < 0 ? i_ >= i_1 : i_ <= i_1; i_ += i_2) {

	i_3 = *m - i_ + 1;
	i_4 = *n - i_ + 1;
	PLUMED_BLAS_F77_FUNC(dlabrd,DLABRD)(&i_3, &i_4, &nb, &a[i_ + i_ * a_dim1], lda, &d__[i_], 
		&e[i_], &tauq[i_], &taup[i_], &work[1], &ldwrkx, 
		&work[ldwrkx * nb + 1], &ldwrky);

	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "T", &i_3, &i_4, &nb, &minusone, 
	       &a[i_ + nb + i_ * a_dim1], lda, &work[ldwrkx * nb + nb + 1],
	       &ldwrky, &one, &a[i_ + nb + (i_ + nb) * a_dim1], lda);
	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", &i_3, &i_4, &nb, &minusone, &work[nb + 1], &ldwrkx,
	       &a[i_ + (i_ + nb) * a_dim1], lda, &one, 
	       &a[i_ + nb + (i_ + nb) * a_dim1], lda);

	if (*m >= *n) {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + (j + 1) * a_dim1] = e[j];
	    }
	} else {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + 1 + j * a_dim1] = e[j];
	    }
	}
    }

    i_2 = *m - i_ + 1;
    i_1 = *n - i_ + 1;
    PLUMED_BLAS_F77_FUNC(dgebd2,DGEBD2)(&i_2, &i_1, &a[i_ + i_ * a_dim1], lda, &d__[i_], &e[i_], &
	    tauq[i_], &taup[i_], &work[1], &iinfo);
    work[1] = ws;
    return;

}
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dgelq2,DGELQ2)(int *m, 
                        int *n, 
                        double *a,
                        int *lda, 
                        double *tau, 
                        double *work, 
                        int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    int i__, k;
    double aii;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    
    i__4 = (*m > 1) ? *m : 1;
    
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < i__4) {
	*info = -4;
    }
    if (*info != 0) {
	return;
    }

    
    k = (*m < *n ) ? *m : *n;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__ + 1;
	i__3 = i__ + 1;
    i__4 = (i__3 < *n) ? i__3 : *n;
	PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + i__4 * a_dim1],
                            lda, &tau[i__]);
	if (i__ < *m) {
	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.f;
	    i__2 = *m - i__;
	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("R", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, 
                              &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    a[i__ + i__ * a_dim1] = aii;
	}
    }
    return;
}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgelqf,DGELQF)(int *m,
	int *n, 
	double *a, 
	int *lda, 
	double *tau,
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, k, ib, nb, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    nb = DGELQF_BLOCKSIZE;
    lwkopt = *m * nb;
    work[1] = (double) lwkopt;

    if (*lwork==-1) {
	return;
    }

    k =(*m < *n) ? *m : *n;
    if (k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {
	nx = DGELQF_CROSSOVER;
	if (nx < k) {
	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DGELQF_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__3 = k - i__ + 1;
	    ib = (i__3 < nb) ? i__3 : nb;

	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dgelq2,DGELQ2)(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
	    if (i__ + ib <= *m) {

		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__3 = *m - i__ - ib + 1;
		i__4 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }
	}
    } else {
	i__ = 1;
    }

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(dgelq2,DGELQ2)(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
    }

    work[1] = (double) iws;
    return;

}
}
}
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgeqr2,DGEQR2)(int *m,
	int *n,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *info)
{
  int k = (*m < *n) ? *m : *n;
  int i,i1,i2,i3;
  double aii;

  *info = 0;
  
  for(i=0;i<k;i++) {
    i1 = *m - i;
    i2 = ( (i+1) < (*m-1) ) ? (i+1) : (*m-1);
    i3 = 1;
    PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tau[i]));
    if(i<(*n-1)) {
      aii = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tau[i]),
	     &(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = aii;
    }
  }
  return;
}
}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dgeqrf,DGEQRF)(int *m, 
	int *n, 
	double *a, 
	int *lda, 
	double *tau,
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, k, ib, nb, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    nb = DGEQRF_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (double) lwkopt;
        if (*lwork==-1)
	return;
    

    k = (*m < *n) ? *m : *n;
    if (k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {
	
      nx = DGEQRF_CROSSOVER;
	if (nx < k) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DGEQRF_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {
	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

	    i__3 = k - i__ + 1;
	    ib = (i__3 < nb) ? i__3 : nb;

	    i__3 = *m - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dgeqr2,DGEQR2)(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
	    if (i__ + ib <= *n) {

		i__3 = *m - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork);
	    }
	}
    } else {
	i__ = 1;
    }

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(dgeqr2,DGEQR2)(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
    }

    work[1] = (double) iws;
    return;

} 

}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dgesdd,DGESDD)(const char *jobz, 
	int *m, 
	int *n, 
	double *a, 
	int *lda, 
	double *s,
	double *u, 
	int *ldu, 
	double *vt, 
	int *ldvt, 
	double *work,
	int *lwork, 
	int *iwork, 
	int *info)
{
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int ie, iu;
    double dum[1], eps;
    int ivt, iscl;
    double anrm;
    int idum[1], ierr, itau;
    int minmn, wrkbl, itaup, itauq, mnthr;
    int nwork;
    int wntqn;
    int bdspac;
    double bignum;
    int ldwrku, maxwrk, ldwkvt;
    double smlnum,minval, safemin;
    int lquery;
    int c__0 = 0;
    int c__1 = 1;
    double zero = 0.0;
    double one = 1.0;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;

    *info = 0;
    minmn = (*m < *n) ? *m : *n;
    mnthr = (int) (minmn * 11. / 6.);
    wntqn = (*jobz=='o' || *jobz=='O');

    maxwrk = 1;
    lquery = *lwork == -1;

    if (*info == 0 && *m > 0 && *n > 0) {
	if (*m >= *n) {

	    if (wntqn) {
		bdspac = *n * 7;
	    } else {
		bdspac = *n * 3 * *n + (*n << 2);
	    }
	    if (*m >= mnthr) {
		if (wntqn) {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = bdspac + *n;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = *n + (*m << 5);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *n * *n;
		}
	    } else {

		wrkbl = *n * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	} else {

	    if (wntqn) {
		bdspac = *m * 7;
	    } else {
		bdspac = *m * 3 * *m + (*m*4);
	    }
	    if (*n >= mnthr) {
		if (wntqn) {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = bdspac + *m;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = *m + (*n*32);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;

		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *m * *m;
		}
	    } else {
		wrkbl = *m * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	}
	work[1] = (double) maxwrk;
    }

    
    if( lquery != 0)
    {
        return;
    }
    

    if (*m == 0 || *n == 0) {
	if (*lwork >= 1) {
	    work[1] = 1.;
	}
	return;
    }
    eps = PLUMED_GMX_DOUBLE_EPS;
    minval = PLUMED_GMX_DOUBLE_MIN;
    safemin = minval / eps;
    smlnum =  std::sqrt(safemin) / eps;


    bignum = 1. / smlnum;


    anrm = PLUMED_BLAS_F77_FUNC(dlange,DLANGE)("M", m, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G",&c__0,&c__0,&anrm,&smlnum,m,n,&a[a_offset],lda,&ierr);
    } else if (anrm > bignum) {
	iscl = 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G",&c__0,&c__0,&anrm,&bignum,m,n,&a[a_offset],lda,&ierr);
    }

    if (*m >= *n) {
	if (*m >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgeqrf,DGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = 1;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *n;

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {
		iu = 1;

		ldwrku = *n;
		itau = iu + ldwrku * *n;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgeqrf,DGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("L", m, n, &a[a_offset], lda, &u[u_offset], ldu);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dorgqr,DORGQR)(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = itau;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);

		PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", m, n, n, &one, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &zero, &a[a_offset], lda);

		PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("F", m, n, &a[a_offset], lda, &u[u_offset], ldu);

	    }

	} else {
	    ie = 1;
	    itauq = ie + *n;
	    itaup = itauq + *n;
	    nwork = itaup + *n;

	    i__1 = *lwork - nwork + 1;
	    PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {

		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("F", m, m, &zero, &zero, &u[u_offset], ldu);
		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *m - *n;
		i__2 = *m - *n;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("F", &i__1, &i__2, &zero, &one, &u[*n + 1 + (*n + 
			1) * u_dim1], ldu);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset],ldvt,&work[nwork],&i__1,&ierr);
	    }

	}

    } else {

	if (*n >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgelqf,DGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = 1;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *m;

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {

		ivt = 1;

		ldwkvt = *m;
		itau = ivt + ldwkvt * *m;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgelqf,DGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dorglq,DORGLQ)(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = itau;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr);

		PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", m, n, m, &one, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &zero, &a[a_offset], lda);

		PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

	    }

	} else {

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    nwork = itaup + *m;

	    i__1 = *lwork - nwork + 1;
	    PLUMED_BLAS_F77_FUNC(dgebrd,DGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("F", n, n, &zero, &zero, &vt[vt_offset], ldvt);
		PLUMED_BLAS_F77_FUNC(dbdsdc,DBDSDC)("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *n - *m;
		i__2 = *n - *m;
		PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("F", &i__1, &i__2, &zero, &one, &vt[*m + 1 + (*m + 
			1) * vt_dim1], ldvt);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);
	    }

	}

    }

    if (iscl == 1) {
	if (anrm > bignum) {
	    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
	if (anrm < smlnum) {
	    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
    }

    work[1] = (double) maxwrk;

    return;

}


}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgetf2,DGETF2)(int *m,
	int *n,
	double *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int j,jp,k,t1,t2,t3;
  double minusone;
  double tmp;

  minusone = -1.0;

  if(*m<=0 || *n<=0)
    return;

  k = (*m < *n) ? *m : *n;
  for(j=1;j<=k;j++) {
    t1 = *m-j+1;
    t2 = 1;
    jp = j - 1 + PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)(&t1,&(a[(j-1)*(*lda)+(j-1)]),&t2);
    ipiv[j-1] = jp;
    if( std::abs(a[(j-1)*(*lda)+(jp-1)])>PLUMED_GMX_DOUBLE_MIN ) {
      if(jp != j)
	PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n,&(a[ j-1 ]),lda,&(a[ jp-1 ]),lda);
      
      if(j<*m) {
	t1 = *m-j;
	t2 = 1;
	tmp = 1.0/a[(j-1)*(*lda)+(j-1)];
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&t1,&tmp,&(a[(j-1)*(*lda)+(j)]),&t2);
      }
    } else {
      *info = j;
    }

    if(j<k) {
      t1 = *m-j;
      t2 = *n-j;
      t3 = 1;
      PLUMED_BLAS_F77_FUNC(dger,DGER)(&t1,&t2,&minusone,&(a[(j-1)*(*lda)+(j)]),&t3,
	    &(a[(j)*(*lda)+(j-1)]),lda, &(a[(j)*(*lda)+(j)]),lda);
    }
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgetrf,DGETRF)(int *m,
	int *n,
	double *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int mindim,jb;
  int i,j,k,l;
  int iinfo;
  double minusone = -1.0;
  double one = 1.0;

  if(*m<=0 || *n<=0)
    return;

  *info = 0;

  mindim = (*m < *n) ? *m : *n;

  if(DGETRF_BLOCKSIZE>=mindim) {

    /* unblocked code */
    PLUMED_BLAS_F77_FUNC(dgetf2,DGETF2)(m,n,a,lda,ipiv,info);

  } else {

    /* blocked case */

    for(j=1;j<=mindim;j+=DGETRF_BLOCKSIZE) {
      jb = ( DGETRF_BLOCKSIZE < (mindim-j+1)) ? DGETRF_BLOCKSIZE : (mindim-j+1);
      /* factor diag. and subdiag blocks and test for singularity */
      k = *m-j+1;
      PLUMED_BLAS_F77_FUNC(dgetf2,DGETF2)(&k,&jb,&(a[(j-1)*(*lda)+(j-1)]),lda,&(ipiv[j-1]),&iinfo);
      
      if(*info==0 && iinfo>0)
	*info = iinfo + j - 1;

      /* adjust pivot indices */
      k = (*m < (j+jb-1)) ? *m : (j+jb-1);
      for(i=j;i<=k;i++)
	ipiv[i-1] += j - 1;

      /* Apply to columns 1 throughj j-1 */
      k = j - 1;
      i = j + jb - 1;
      l = 1;
      PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(&k,a,lda,&j,&i,ipiv,&l);
      if((j+jb)<=*n) {
	/* Apply to cols. j+jb through n */
	k = *n-j-jb+1;
	i = j+jb-1;
	l = 1;
	PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(&k,&(a[(j+jb-1)*(*lda)+0]),lda,&j,&i,ipiv,&l);
	/* Compute block row of U */
	k = *n-j-jb+1;
	PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Left","Lower","No transpose","Unit",&jb,&k,&one,
	       &(a[(j-1)*(*lda)+(j-1)]),lda,&(a[(j+jb-1)*(*lda)+(j-1)]),lda);

	if((j+jb)<=*m) {
	  /* Update trailing submatrix */
	  k = *m-j-jb+1;
	  i = *n-j-jb+1;
	  PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose","No transpose",&k,&i,&jb,&minusone,
		 &(a[(j-1)*(*lda)+(j+jb-1)]),lda,
		 &(a[(j+jb-1)*(*lda)+(j-1)]),lda,&one,
		 &(a[(j+jb-1)*(*lda)+(j+jb-1)]),lda);
	}

      }
    }
  }
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgetri,DGETRI)(int *n, 
	double *a, 
	int *lda, 
	int *ipiv, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, jb, nb, jj, jp, nn, iws;
    int nbmin;
    int ldwork;
    int lwkopt;
    int c__1 = 1;
    double c_b20 = -1.;
    double c_b22 = 1.;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;

    *info = 0;
    nb = DGETRI_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (double) lwkopt;

    if (*n < 0) {
	*info = -1;
    } else if (*lda < (*n)) {
	*info = -3;
    } else if (*lwork < (*n) && *lwork!=-1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (*lwork == -1) {
	return;
    }

    if (*n == 0) {
	return;
    }

    PLUMED_BLAS_F77_FUNC(dtrtri,DTRTRI)("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if (*info > 0) {
	return;
    }

    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
	i__1 = ldwork * nb;
	iws = (i__1>1) ? i__1 : 1;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DGETRI_MINBLOCKSIZE;
	}
    } else {
	iws = *n;
    }

    if (nb < nbmin || nb >= *n) {

	for (j = *n; j >= 1; --j) {

	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		work[i__] = a[i__ + j * a_dim1];
		a[i__ + j * a_dim1] = 0.;
	    }

	    if (j < *n) {
		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 
			+ 1], lda, &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 
			+ 1], &c__1);
	    }
	}
    } else {

	nn = (*n - 1) / nb * nb + 1;
	i__1 = -nb;
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
	    i__2 = nb, i__3 = *n - j + 1;
	    jb = (i__2<i__3) ? i__2 : i__3;

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= i__2; ++jj) {
		i__3 = *n;
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
		    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
		    a[i__ + jj * a_dim1] = 0.;
		}
	    }

	    if (j + jb <= *n) {
		i__2 = *n - j - jb + 1;
		PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &
			ldwork, &c_b22, &a[j * a_dim1 + 1], lda);
	    }
	    PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda);
	}
    }

    for (j = *n - 1; j >= 1; --j) {
	jp = ipiv[j];
	if (jp != j) {
	    PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
	}
    }

    work[1] = (double) iws;
    return;

}


}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dgetrs,DGETRS)(const char *trans, 
	int *n, 
	int *nrhs, 
	double *a, 
	int *lda, 
	int *ipiv,
	double *b, 
	int *ldb, 
	int *info)
{
    int a_dim1, a_offset, b_dim1, b_offset;
    int notran;
    int c__1 = 1;
    int c_n1 = -1;
    double one = 1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    *info = 0;
    notran = (*trans=='N' || *trans=='n');

    if (*n <= 0 || *nrhs <= 0) 
	return;

    if (notran) {
	PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
	PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Left", "Lower", "No transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
    } else {
	PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
	PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Left", "Lower", "Transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlabrd,DLABRD)(int *m, 
	int *n, 
	int *nb,
	double *a, 
	int *lda, 
	double *d__,
	double *e,
	double *tauq, 
	double *taup,
	double *x,
	int *ldx,
	double *y,
	int *ldy)
{
    int a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset;
    int i__1, i__2, i__3;
    double one = 1.0;
    double minusone = -1.0;
    double zero = 0.0;
    int c__1 = 1;
    int i__;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }

    if (*m >= *n) {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + a_dim1], lda,
		     &y[i__ + y_dim1], ldy, &one, &a[i__ + i__ * a_dim1], &c__1);
	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + x_dim1], ldx,
		   &a[i__*a_dim1+1],&c__1,&one,&a[i__+i__*a_dim1],&c__1);

	    i__2 = *m - i__ + 1;
	    i__3 = i__ + 1;
	    if(*m<i__3)
	      i__3 = *m;
	    PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__3 + i__ * a_dim1], 
		    &c__1, &tauq[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *n) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__ + 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + (i__ + 1) * 
			a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &
			y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + a_dim1], 
			lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &x[i__ + x_dim1], 
			ldx, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, 
			&y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &a[i__ + a_dim1], lda, &one, &a[i__ + (
			i__ + 1) * a_dim1], lda);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &one, &a[
			i__ + (i__ + 1) * a_dim1], lda);

		i__2 = *n - i__;
		i__3 = i__ + 2;
		if(*n<i__3)
		  i__3 = *n;
		PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + (i__ + 1) * a_dim1], 
			&a[i__ + i__3 * a_dim1], lda, &taup[i__]);
		e[i__] = a[i__ + (i__ + 1) * a_dim1];
		a[i__ + (i__ + 1) * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ 
			+ 1) * a_dim1], lda, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &zero, &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__, &one, &y[i__ + 1 + y_dim1], 
			ldy, &a[i__ + (i__ + 1) * a_dim1], lda, &zero, &x[
			i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &one, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			zero, &x[i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
	    }
	}
    } else {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + y_dim1], ldy,
		     &a[i__ + a_dim1], lda, &one, &a[i__ + i__ * a_dim1],lda);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &minusone, &a[i__ * a_dim1 + 1], 
		    lda, &x[i__ + x_dim1], ldx, &one,&a[i__+i__*a_dim1],lda);

	    i__2 = *n - i__ + 1;
	    i__3 = i__ + 1;
	    if(*n<i__3)
	      i__3 = *n;
	    PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + i__ * a_dim1], 
		    &a[i__ + i__3 * a_dim1], lda, &taup[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose",&i__2,&i__3,&one,&a[i__+1+i__*a_dim1], 
		       lda, &a[i__ + i__ * a_dim1], lda, &zero, 
		       &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &y[i__ + y_dim1], 
			ldy, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ * 
			x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ * a_dim1 + 
			1], lda, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ *
			 x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &y[i__ + y_dim1], ldy, &one, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &a[i__ * a_dim1 + 1], &c__1, &one, &a[
			i__ + 1 + i__ * a_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ + 2;
		if(*m<i__3)
		  i__3 = *m;
		PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + 1 + i__ * a_dim1], 
			&a[i__3 + i__ * a_dim1], &c__1, &tauq[i__]);
		e[i__] = a[i__ + 1 + i__ * a_dim1];
		a[i__ + 1 + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ + 
			1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, 
			&zero, &y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__, &one, &x[i__ + 1 + x_dim1], 
			ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__, &i__2, &minusone, &a[(i__ + 1) * a_dim1 
			+ 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, &y[i__ 
			+ 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
	    }
	}
    }
    return;
} 

}
}
#include <cctype>
#include "lapack.h"

/* LAPACK */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)(const char *uplo,
	int *m,
	int *n,
	double *a,
	int *lda,
	double *b,
	int *ldb)
{
  int i,j,minjm;
  const char ch=std::toupper(*uplo);

  if(ch=='U') {
    for(j=0;j<*n;j++) {
      minjm = (j < (*m-1)) ? j : (*m-1);
      for(i=0;i<=minjm;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else if(ch=='L') {
    for(j=0;j<*n;j++) {
      for(i=j;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }    
  }
}
}
}
#include <cmath>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)(double *a, 
       double *b,
       double *c__, 
       double *rt1, 
       double *rt2)
{
    double d__1;
    double ab, df, tb, sm, rt, adf, acmn, acmx;


    sm = *a + *c__;
    df = *a - *c__;
    adf = std::abs(df);
    tb = *b + *b;
    ab = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf *  std::sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab *  std::sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab *  std::sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
    }
    return;

}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaebz,DLAEBZ)(int *ijob,
	int *nitmax,
	int *n, 
	int *mmax,
	int *minp,
	int *nbmin,
	double *abstol, 
	double *reltol, 
	double *pivmin, 
	double *d__,
	double *e,
	double *e2, 
	int *nval,
	double *ab, 
	double *c__, 
	int *mout, 
	int *nab,
	double *work,
	int *iwork, 
	int *info)
{
    int nab_dim1, nab_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    double d__1, d__2, d__3, d__4;

    int j, kf, ji, kl, jp, jit;
    double tmp1, tmp2;
    int itmp1, itmp2, kfnew, klnew;

    nab_dim1 = *mmax;
    nab_offset = 1 + nab_dim1;
    nab -= nab_offset;
    ab_dim1 = *mmax;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --d__;
    --e;
    --e2;
    --nval;
    --c__;
    --work;
    --iwork;

    *info = 0;
    if (*ijob < 1 || *ijob > 3) {
	*info = -1;
	return;
    }

    if (*ijob == 1) {

	*mout = 0;

	i__1 = *minp;
	for (ji = 1; ji <= i__1; ++ji) {
	    for (jp = 1; jp <= 2; ++jp) {
		tmp1 = d__[1] - ab[ji + jp * ab_dim1];
		if (std::abs(tmp1) < *pivmin) {
		    tmp1 = -(*pivmin);
		}
		nab[ji + jp * nab_dim1] = 0;
		if (tmp1 <= 0.) {
		    nab[ji + jp * nab_dim1] = 1;
		}

		i__2 = *n;
		for (j = 2; j <= i__2; ++j) {
		    tmp1 = d__[j] - e2[j - 1] / tmp1 - ab[ji + jp * ab_dim1];
		    if (std::abs(tmp1) < *pivmin) {
			tmp1 = -(*pivmin);
		    }
		    if (tmp1 <= 0.) {
			++nab[ji + jp * nab_dim1];
		    }
		}
	    }
	    *mout = *mout + nab[ji + (nab_dim1 << 1)] - nab[ji + nab_dim1];
	}
	return;
    }

    kf = 1;
    kl = *minp;

    if (*ijob == 2) {
	i__1 = *minp;
	for (ji = 1; ji <= i__1; ++ji) {
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
	}
    }

    i__1 = *nitmax;
    for (jit = 1; jit <= i__1; ++jit) {

	if (kl - kf + 1 >= *nbmin && *nbmin > 0) {

	    i__2 = kl;
	    for (ji = kf; ji <= i__2; ++ji) {

		work[ji] = d__[1] - c__[ji];
		iwork[ji] = 0;
		if (work[ji] <= *pivmin) {
		    iwork[ji] = 1;
		    d__1 = work[ji], d__2 = -(*pivmin);
		    work[ji] = (d__1<d__2) ? d__1 : d__2;
		}

		i__3 = *n;
		for (j = 2; j <= i__3; ++j) {
		    work[ji] = d__[j] - e2[j - 1] / work[ji] - c__[ji];
		    if (work[ji] <= *pivmin) {
			++iwork[ji];
			d__1 = work[ji], d__2 = -(*pivmin);
			work[ji] = (d__1<d__2) ? d__1 : d__2;
		    }
		}
	    }

	    if (*ijob <= 2) {

		klnew = kl;
		i__2 = kl;
		for (ji = kf; ji <= i__2; ++ji) {

		  i__5 = nab[ji + nab_dim1];
		  i__6 = iwork[ji];
		  i__3 = nab[ji + (nab_dim1 << 1)];
		  i__4 = (i__5>i__6) ? i__5 : i__6;
		    iwork[ji] = (i__3<i__4) ? i__3 : i__4;

		    if (iwork[ji] == nab[ji + (nab_dim1 << 1)]) {

			ab[ji + (ab_dim1 << 1)] = c__[ji];

		    } else if (iwork[ji] == nab[ji + nab_dim1]) {

			ab[ji + ab_dim1] = c__[ji];
		    } else {
			++klnew;
			if (klnew <= *mmax) {

			    ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 
				    1)];
			    nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 
				    << 1)];
			    ab[klnew + ab_dim1] = c__[ji];
			    nab[klnew + nab_dim1] = iwork[ji];
			    ab[ji + (ab_dim1 << 1)] = c__[ji];
			    nab[ji + (nab_dim1 << 1)] = iwork[ji];
			} else {
			    *info = *mmax + 1;
			}
		    }
		}
		if (*info != 0) {
		    return;
		}
		kl = klnew;
	    } else {

		i__2 = kl;
		for (ji = kf; ji <= i__2; ++ji) {
		    if (iwork[ji] <= nval[ji]) {
			ab[ji + ab_dim1] = c__[ji];
			nab[ji + nab_dim1] = iwork[ji];
		    }
		    if (iwork[ji] >= nval[ji]) {
			ab[ji + (ab_dim1 << 1)] = c__[ji];
			nab[ji + (nab_dim1 << 1)] = iwork[ji];
		    }
		}
	    }

	} else {

	    klnew = kl;
	    i__2 = kl;
	    for (ji = kf; ji <= i__2; ++ji) {

		tmp1 = c__[ji];
		tmp2 = d__[1] - tmp1;
		itmp1 = 0;
		if (tmp2 <= *pivmin) {
		    itmp1 = 1;
		    d__1 = tmp2, d__2 = -(*pivmin);
		    tmp2 = (d__1<d__2) ? d__1 : d__2;
		}

		i__3 = *n;
		for (j = 2; j <= i__3; ++j) {
		    tmp2 = d__[j] - e2[j - 1] / tmp2 - tmp1;
		    if (tmp2 <= *pivmin) {
			++itmp1;
			d__1 = tmp2, d__2 = -(*pivmin);
			tmp2 = (d__1<d__2) ? d__1 : d__2;
		    }
		}

		if (*ijob <= 2) {

		    i__5 = nab[ji + nab_dim1];
		    i__3 = nab[ji + (nab_dim1 << 1)];
		    i__4 = (i__5>itmp1) ? i__5 : itmp1;
		    itmp1 = (i__3<i__4) ? i__3 : i__4;

		    if (itmp1 == nab[ji + (nab_dim1 << 1)]) {

			ab[ji + (ab_dim1 << 1)] = tmp1;

		    } else if (itmp1 == nab[ji + nab_dim1]) {

			ab[ji + ab_dim1] = tmp1;
		    } else if (klnew < *mmax) {

			++klnew;
			ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 1)];
			nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 << 
				1)];
			ab[klnew + ab_dim1] = tmp1;
			nab[klnew + nab_dim1] = itmp1;
			ab[ji + (ab_dim1 << 1)] = tmp1;
			nab[ji + (nab_dim1 << 1)] = itmp1;
		    } else {
			*info = *mmax + 1;
			return;
		    }
		} else {

		    if (itmp1 <= nval[ji]) {
			ab[ji + ab_dim1] = tmp1;
			nab[ji + nab_dim1] = itmp1;
		    }
		    if (itmp1 >= nval[ji]) {
			ab[ji + (ab_dim1 << 1)] = tmp1;
			nab[ji + (nab_dim1 << 1)] = itmp1;
		    }
		}
	    }
	    kl = klnew;

	}

	kfnew = kf;
	i__2 = kl;
	for (ji = kf; ji <= i__2; ++ji) {
	    tmp1 = std::abs(ab[ji + (ab_dim1 << 1)] - ab[ji + ab_dim1]);
	    d__3 = std::abs(ab[ji + (ab_dim1 << 1)]);
	    d__4 = std::abs(ab[ji + ab_dim1]);
	    tmp2 = (d__3>d__4) ? d__3 : d__4;
	    d__1 = (*abstol>*pivmin) ? *abstol : *pivmin;
	    d__2 = *reltol * tmp2;
	    if (tmp1 < ((d__1>d__2) ? d__1 : d__2) || nab[ji + nab_dim1] >= nab[ji + (
		    nab_dim1 << 1)]) {

		if (ji > kfnew) {
		    tmp1 = ab[ji + ab_dim1];
		    tmp2 = ab[ji + (ab_dim1 << 1)];
		    itmp1 = nab[ji + nab_dim1];
		    itmp2 = nab[ji + (nab_dim1 << 1)];
		    ab[ji + ab_dim1] = ab[kfnew + ab_dim1];
		    ab[ji + (ab_dim1 << 1)] = ab[kfnew + (ab_dim1 << 1)];
		    nab[ji + nab_dim1] = nab[kfnew + nab_dim1];
		    nab[ji + (nab_dim1 << 1)] = nab[kfnew + (nab_dim1 << 1)];
		    ab[kfnew + ab_dim1] = tmp1;
		    ab[kfnew + (ab_dim1 << 1)] = tmp2;
		    nab[kfnew + nab_dim1] = itmp1;
		    nab[kfnew + (nab_dim1 << 1)] = itmp2;
		    if (*ijob == 3) {
			itmp1 = nval[ji];
			nval[ji] = nval[kfnew];
			nval[kfnew] = itmp1;
		    }
		}
		++kfnew;
	    }
	}
	kf = kfnew;

	i__2 = kl;
	for (ji = kf; ji <= i__2; ++ji) {
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
	}

	if (kf > kl) {
	    break;
	}
    }

    i__1 = kl + 1 - kf;
    if(i__1>0)
      *info = i__1;

    *mout = kl;

    return;

}


}
}
#include <cmath>

#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaed6,DLAED6)(int *kniter, 
                        int *orgati, 
                        double *rho, 
                        double *d__,
                        double *z__, 
                        double *finit, 
                        double *tau, 
                        int *info)
{
    int i__1;
    double r__1, r__2, r__3, r__4;

    double a, b, c__, f;
    int i__;
    double fc, df, ddf, eta, eps, base;
    int iter;
    double temp, temp1, temp2, temp3, temp4;
    int scale;
    int niter;
    double small1, small2, sminv1, sminv2, dscale[3], sclfac;
    double zscale[3], erretm;
    double safemin;
    double sclinv = 0;
    
    --z__;
    --d__;

    *info = 0;

    niter = 1;
    *tau = 0.f;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (d__[3] - d__[2]) / 2.f;
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
	} else {
	    temp = (d__[1] - d__[2]) / 2.f;
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
	}
        r__1 = std::abs(a), r__2 = std::abs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = std::abs(c__);
        temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    *tau = b / a;
	} else if (a <= 0.f) {
	    *tau = (a -  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1)))) / (
		    c__ * 2.f);
	} else {
	    *tau = b * 2.f / (a +  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1))));
	}

	temp = *rho + z__[1] / (d__[1] - *tau) + z__[2] / (d__[2] - *tau) + 
		z__[3] / (d__[3] - *tau);
	if (std::abs(*finit) <= std::abs(temp)) {
	    *tau = 0.f;
	}
    }

    eps = PLUMED_GMX_DOUBLE_EPS;
    base = 2;
    safemin = PLUMED_GMX_DOUBLE_MIN*(1.0+PLUMED_GMX_DOUBLE_EPS);
    i__1 = static_cast<int>(std::log(safemin) / std::log(base) / 3.f);
    small1 = std::pow(base, static_cast<double>(i__1));
    sminv1 = 1.f / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;

    if (*orgati) {
	r__3 = (r__1 = d__[2] - *tau, std::abs(r__1)), r__4 = (r__2 = d__[3] - *
		tau, std::abs(r__2));
        temp = (r__3<r__4) ? r__3 : r__4;
    } else {
	r__3 = (r__1 = d__[1] - *tau, std::abs(r__1)), r__4 = (r__2 = d__[2] - *
		tau, std::abs(r__2));
	temp = (r__3<r__4) ? r__3 : r__4;
    }
    scale = 0;
    if (temp <= small1) {
	scale = 1;
	if (temp <= small2) {

	    sclfac = sminv2;
	    sclinv = small2;
	} else {

	    sclfac = sminv1;
	    sclinv = small1;

	}

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__] * sclfac;
	    zscale[i__ - 1] = z__[i__] * sclfac;
	}
	*tau *= sclfac;
    } else {

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__];
	    zscale[i__ - 1] = z__[i__];
	}
    }
    fc = 0.f;
    df = 0.f;
    ddf = 0.f;
    for (i__ = 1; i__ <= 3; ++i__) {
	temp = 1.f / (dscale[i__ - 1] - *tau);
	temp1 = zscale[i__ - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i__ - 1];
	df += temp2;
	ddf += temp3;
    }
    f = *finit + *tau * fc;

    if (std::abs(f) <= 0.f) {
	goto L60;
    }
    iter = niter + 1;
    for (niter = iter; niter <= 20; ++niter) {
	if (*orgati) {
	    temp1 = dscale[1] - *tau;
	    temp2 = dscale[2] - *tau;
	} else {
	    temp1 = dscale[0] - *tau;
	    temp2 = dscale[1] - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	r__1 = std::abs(a), r__2 = std::abs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = std::abs(c__);
	temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    eta = b / a;
	} else if (a <= 0.f) {
	    eta = (a -  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1)))) / ( c__ * 2.f);
	} else {
	    eta = b * 2.f / (a +  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs( r__1))));
	}
	if (f * eta >= 0.f) {
	    eta = -f / df;
	}
	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0.f && temp >= dscale[2]) {
		eta = (dscale[2] - *tau) / 2.f;
	    }

	    if (eta < 0.f && temp <= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	} else {
	    if (eta > 0.f && temp >= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	    if (eta < 0.f && temp <= dscale[0]) {
		eta = (dscale[0] - *tau) / 2.f;
	    }
	}
	*tau += eta;
	fc = 0.f;
	erretm = 0.f;
	df = 0.f;
	ddf = 0.f;
	for (i__ = 1; i__ <= 3; ++i__) {
	    temp = 1.f / (dscale[i__ - 1] - *tau);
	    temp1 = zscale[i__ - 1] * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    temp4 = temp1 / dscale[i__ - 1];
	    fc += temp4;
	    erretm += std::abs(temp4);
	    df += temp2;
	    ddf += temp3;
	}
	f = *finit + *tau * fc;
	erretm = (std::abs(*finit) + std::abs(*tau) * erretm) * 8.f + std::abs(*tau) * df;
	if (std::abs(f) <= eps * erretm) {
	    goto L60;
	}
    }
    *info = 1;
L60:
    if (scale) {
	*tau *= sclinv;
    }
    return;
} 


}
}
#include <cmath>
#include "real.h"

#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaev2,DLAEV2)(double *   a, 
	double *   b, 
	double *   c__, 
	double *   rt1, 
	double *   rt2, 
	double *   cs1, 
	double *   sn1)
{
    double d__1;

    double ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    int sgn1, sgn2;
    double acmn, acmx;

    sm = *a + *c__;
    df = *a - *c__;
    adf = std::abs(df);
    tb = *b + *b;
    ab = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf *  std::sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab *  std::sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab *  std::sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }
    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = std::abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. /  std::sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (std::abs(ab)<PLUMED_GMX_DOUBLE_MIN) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. /  std::sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;

}


}
}
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlagtf,DLAGTF)(int *n, 
	double *a, 
	double *lambda, 
	double *b, 
	double *c__, 
	double *tol, 
	double *d__, 
	int *in, 
	int *info)
{
    int i__1;

    int k;
    double tl, eps, piv1, piv2, temp, mult, scale1, scale2;

    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (*n < 0) {
	*info = -1;
	return;
    }

    if (*n == 0) 
	return;
    
    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (std::abs(a[1])<PLUMED_GMX_DOUBLE_MIN) {
	    in[1] = 1;
	}
	return;
    }

    eps = PLUMED_GMX_DOUBLE_EPS;

    tl = (*tol>eps) ? *tol : eps;
    scale1 = std::abs(a[1]) + std::abs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = std::abs(c__[k]) + std::abs(a[k + 1]);
	if (k < *n - 1) {
	    scale2 += std::abs(b[k + 1]);
	}
	if (std::abs(a[k])<PLUMED_GMX_DOUBLE_MIN) {
	    piv1 = 0.;
	} else {
	    piv1 = std::abs(a[k]) / scale1;
	}
	if (std::abs(c__[k])<PLUMED_GMX_DOUBLE_MIN) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = std::abs(c__[k]) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c__[k] /= a[k];
		a[k + 1] -= c__[k] * b[k];
		if (k < *n - 1) {
		    d__[k] = 0.;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c__[k];
		a[k] = c__[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < *n - 1) {
		    d__[k] = b[k + 1];
		    b[k + 1] = -mult * d__[k];
		}
		b[k] = temp;
		c__[k] = mult;
	    }
	}
	if (((piv1>piv2) ? piv1 : piv2) <= tl && in[*n] == 0) {
	    in[*n] = k;
	}
    }
    if (std::abs(a[*n]) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return;

}


}
}
#include <stdlib.h>
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlagts,DLAGTS)(int *job, 
	int *n, 
	double *a, 
	double *b, 
	double *c__, 
	double *d__, 
	int *in, 
	double *y, 
	double *tol, 
	int *info)
{
    int i__1;
    double d__1, d__2, d__4, d__5;

    int k;
    double ak, eps, temp, pert, absak, sfmin;
    double bignum,minval;
    --y;
    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (abs(*job) > 2 || *job == 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }

    if (*n == 0) {
	return;
    }
    eps = PLUMED_GMX_DOUBLE_EPS;
    minval = PLUMED_GMX_DOUBLE_MIN;
    sfmin = minval / eps;

    bignum = 1. / sfmin;

    if (*job < 0) {
	if (*tol <= 0.) {
	    *tol = std::abs(a[1]);
	    if (*n > 1) {
		d__1 = *tol;
		d__2 = std::abs(a[2]);
		d__1 = (d__1>d__2) ? d__1 : d__2;
		d__2 = std::abs(b[1]);
		*tol = (d__1>d__2) ? d__1 : d__2;
	    }
	    i__1 = *n;
	    for (k = 3; k <= i__1; ++k) {
	      d__4 = *tol;
	      d__5 = std::abs(a[k]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = std::abs(b[k - 1]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = std::abs(d__[k - 2]);
	      *tol = (d__4>d__5) ? d__4 : d__5;
	    }
	    *tol *= eps;
	    if (std::abs(*tol)<PLUMED_GMX_DOUBLE_MIN) {
		*tol = eps;
	    }
	}
    }

    if (1 == abs(*job)) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    if (in[k - 1] == 0) {
		y[k] -= c__[k - 1] * y[k - 1];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
	if (*job == 1) {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k == *n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_DOUBLE_MIN || std::abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    for (k = *n; k >= 1; --k) {
		if (k + 2 <= *n) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k + 1 == *n) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;
L40:
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_DOUBLE_MIN || std::abs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L40;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L40;
		    }
		}
		y[k] = temp / ak;
	    }
	}
    } else {

	if (*job == 2) {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_DOUBLE_MIN || std::abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;

L70:
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_DOUBLE_MIN || std::abs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L70;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L70;
		    }
		}
		y[k] = temp / ak;
	    }
	}

	for (k = *n; k >= 2; --k) {
	    if (in[k - 1] == 0) {
		y[k - 1] -= c__[k - 1] * y[k];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
    }

    return;
}


}
}
#include "lapack.h"


/* LAPACK */


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)(int *n1,
                        int *n2,
                        double *a,
                        int *dtrd1,
                        int *dtrd2,
                        int *index)
{
  int n1sv = *n1;
  int n2sv = *n2;
  int i,ind1,ind2;

  if(*dtrd1>0)
    ind1 = 0;
  else
    ind1 = *n1-1;

  if(*dtrd2>0)
    ind2 = *n1;
  else
    ind2 = *n1+*n2-1;

  i = 0;
  
  while(n1sv>0 && n2sv>0) {
    if(a[ind1]<=a[ind2]) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
      n1sv--;
    } else {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
      n2sv--;
    }
  }

  if(n1sv==0) {
    for(n1sv=1;n1sv<=n2sv;n1sv++) {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
    } 
  } else {
    for(n2sv=1;n2sv<=n1sv;n2sv++) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
    } 
  }
  return;
}
}
}
#include <cctype>
#include <cmath>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
double
PLUMED_BLAS_F77_FUNC(dlange,DLANGE)(const char *norm,
	int *m,
	int *n,
	double *a,
	int *lda,
	double *work)
{
  const char ch=std::toupper(*norm);
  double dtemp,sum,max,val,scale;
  int i,j;

  switch(ch) {
  case 'M':
    max = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++) {
	dtemp = std::abs(a[j*(*lda)+i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;

  case 'O':
  case '1':
    max = 0.0;
    for(j=0;j<*n;j++) {
      sum = 0.0;
      for(i=0;i<*m;i++) 
	sum += std::abs(a[j*(*lda)+i]);
      if(sum>max)
	max = sum;
    }
    val = max;
    break;

  case 'I':
    for(i=0;i<*m;i++)
      work[i] = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++)
	work[i] += std::abs(a[j*(*lda)+i]);
    max = 0;
    for(i=0;i<*m;i++)
      if(work[i]>max)
	max=work[i];
    val = max;
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = 1;
    for(j=0;j<*n;j++) 
      PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(m,&(a[j*(*lda)+0]),&i,&scale,&sum);
    val = scale* std::sqrt(sum);
    break;

  default:
    val = 0.0;
    break;
  }
  return val;
}
}
}
#include <cctype>
#include <cmath>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
double
PLUMED_BLAS_F77_FUNC(dlanst,DLANST)(const char *norm,
	int *n,
	double *d,
	double *e)
{
  const char ch=std::toupper(*norm);
  double dtemp,max,val,scale,sum;
  int i,j;


  if(*n<=0)
    return 0.0;
  
  switch(ch) {
  case 'M':
    max = std::abs(d[*n-1]);
      for(i=0;i<(*n-1);i++) {
	dtemp = std::abs(d[i]);
	if(dtemp>max)
	  max = dtemp;
	dtemp = std::abs(e[i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;
    
  case 'O':
  case '1':
  case 'I':

    if(*n==1)
      val = std::abs(d[0]);
    else {
      max = std::abs(d[0]) + std::abs(e[0]);
      dtemp = std::abs(e[*n-2]) + std::abs(d[*n-1]);
      if(dtemp>max)
	max = dtemp;
      for(i=1;i<(*n-1);i++) {
	dtemp = std::abs(d[i]) + std::abs(e[i]) + std::abs(e[i-1]);
	if(dtemp>max)
	  max = dtemp;
      }
      val = max;
    }
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = *n-1;
    j = 1;
    if(*n>1) {
      PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(&i,e,&j,&scale,&sum);
      sum *= 2;
    }
    PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(n,d,&j,&scale,&sum);
    val = scale *  std::sqrt(sum);
    break;
    
  default:
    val = 0.0;
    break;
  }
  return val;
}
}
}
#include <cmath>


#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
double 
PLUMED_BLAS_F77_FUNC(dlansy,DLANSY)(const char *norm, const char *uplo, int *n, double *a, int 
	*lda, double *work)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    double ret_val, d__1, d__2, d__3;
    int c__1 = 1;

    /* Local variables */
    int i__, j;
    double sum, absa, scale;
    double value =0.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;

    if (*n == 0) {
	value = 0.;
    } else if (*norm=='M' || *norm=='m') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = std::abs(a[i__ + j * a_dim1]);
		  value = (d__2>d__3) ? d__2 : d__3;
		}
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = std::abs(a[i__ + j * a_dim1]);
		    value =  (d__2>d__3) ? d__2 : d__3;
		}
	    }
	}
    } else if (*norm=='I' || *norm=='i' || *norm=='O' || *norm=='o' || *norm=='1') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    absa = std::abs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		work[j] = sum + std::abs(a[j + j * a_dim1]);
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = value, d__2 = work[i__];
		value =  (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.;
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = work[j] + std::abs(a[j + j * a_dim1]);
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    absa = std::abs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		if(sum>value)
		  value = sum;
	    }
	}
    } else if (*norm=='F' || *norm=='f' || *norm=='E' || *norm=='e') {

	scale = 0.;
	sum = 1.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(n, &a[a_offset], &i__1, &scale, &sum);
	value = scale *  std::sqrt(sum);
    }

    ret_val = value;
    return ret_val;
}


}
}
#include <cmath>
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
double
PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(double * x, double * y)
{
  double xabs,yabs;
  double w,z;

  xabs = std::abs(*x);
  yabs = std::abs(*y);
  
  if(xabs>yabs) {
    w = xabs;
    z = yabs;
  } else {
    w = yabs;
    z = xabs;
  }

  if( std::abs(z)<PLUMED_GMX_DOUBLE_MIN) 
    return w;
  else {
    z = z/w;
    return w* std::sqrt(1.0+z*z);
  }
}
  
}
}
#include <cmath>

#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;

void PLUMED_BLAS_F77_FUNC(dlar1vx,DLAR1VX)(int *n, 
	      int *b1, 
	      int *bn,
	      double *sigma, 
	      double *d__, 
	      double *l, 
	      double *ld, 
	      double *lld, 
	      double *eval, 
	      double *gersch, 
	      double *z__, 
	      double *ztz, 
	      double *mingma, 
	      int *r__, 
	      int *isuppz, 
	      double *work)
{
    int i__1;

    int i__, j;
    double s;
    int r1, r2;
    int to;
    double eps, tmp;
    int indp, inds, from;
    double dplus;
    int sawnan;
    int indumn;
    double dminus;

    --work;
    --isuppz;
    --z__;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    eps = PLUMED_GMX_DOUBLE_EPS;
    if (*r__ == 0) {

	r1 = *b1;
	r2 = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r1 = i__;
		goto L20;
	    }
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r2 = i__;
		goto L40;
	    }
	}
    } else {
	r1 = *r__;
	r2 = *r__;
    }

L40:
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    sawnan = 0;

    if (*b1 == 1) {
	work[inds] = 0.;
    } else {
	work[inds] = lld[*b1 - 1];
    }
    s = work[inds] - *sigma;
    i__1 = r2 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	work[i__] = ld[i__] / dplus;
	work[inds + i__] = s * work[i__] * l[i__];
	s = work[inds + i__] - *sigma;
    }

    if (std::isnan(s)) {

	sawnan = 1;
	j = *b1 + 1;
L60:
    if (!std::isnan(work[inds + j])) {
	    ++j;
	    goto L60;
	}
	work[inds + j] = lld[j];
	s = work[inds + j] - *sigma;
	i__1 = r2 - 1;
	for (i__ = j + 1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    work[i__] = ld[i__] / dplus;
	    if (std::abs(work[i__])<PLUMED_GMX_DOUBLE_MIN) {
		work[inds + i__] = lld[i__];
	    } else {
		work[inds + i__] = s * work[i__] * l[i__];
	    }
	    s = work[inds + i__] - *sigma;
	}
    }

    work[indp + *bn - 1] = d__[*bn] - *sigma;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
	dminus = lld[i__] + work[indp + i__];
	tmp = d__[i__] / dminus;
	work[indumn + i__] = l[i__] * tmp;
	work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
    }
    tmp = work[indp + r1 - 1];
    if (std::isnan(tmp)) {

	sawnan = 1;
	j = *bn - 3;
L90:
    if (!std::isnan(work[indp + j])) {
	    --j;
	    goto L90;
	}
	work[indp + j] = d__[j + 1] - *sigma;
	i__1 = r1;
	for (i__ = j; i__ >= i__1; --i__) {
	    dminus = lld[i__] + work[indp + i__];
	    tmp = d__[i__] / dminus;
	    work[indumn + i__] = l[i__] * tmp;
	    if (std::abs(tmp)<PLUMED_GMX_DOUBLE_MIN) {
		work[indp + i__ - 1] = d__[i__] - *sigma;
	    } else {
		work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
	    }
	}
    }

    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (std::abs(*mingma)<PLUMED_GMX_DOUBLE_MIN) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = work[inds + i__] + work[indp + i__];
	if (std::abs(tmp)<PLUMED_GMX_DOUBLE_MIN) {
	    tmp = eps * work[inds + i__];
	}
	if (std::abs(tmp) < std::abs(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
    }

    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;
    if (! sawnan) {
	from = *r__ - 1;
	i__1 = *r__ - 32;
	to = (i__1>(*b1)) ? i__1 : (*b1);
L120:
	if (from >= *b1) {
	    i__1 = to;
	    for (i__ = from; i__ >= i__1; --i__) {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (std::abs(z__[to]) <= eps && std::abs(z__[to + 1]) <= eps) {
		isuppz[1] = to + 2;
	    } else {
		from = to - 1;
		i__1 = to - 32;
		to = (i__1>*b1) ? i__1 : *b1;
		goto L120;
	    }
	}
	from = *r__ + 1;
	i__1 = *r__ + 32;
	to = (i__1<*bn) ? i__1 : *bn;
L140:
	if (from <= *bn) {
	    i__1 = to;
	    for (i__ = from; i__ <= i__1; ++i__) {
		z__[i__] = -(work[indumn + i__ - 1] * z__[i__ - 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (std::abs(z__[to]) <= eps && std::abs(z__[to - 1]) <= eps) {
		isuppz[2] = to - 2;
	    } else {
		from = to + 1;
		i__1 = to + 32;
		to = (i__1<*bn) ? i__1 : *bn;
		goto L140;
	    }
	}
    } else {
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (std::abs(z__[i__ + 1])<PLUMED_GMX_DOUBLE_MIN) {
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
	    } else {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
	    }
	    if (std::abs(z__[i__]) <= eps && std::abs(z__[i__ + 1]) <= eps) {
		isuppz[1] = i__ + 2;
		goto L170;
	    }
	    *ztz += z__[i__] * z__[i__];
	}
L170:
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (std::abs(z__[i__])<PLUMED_GMX_DOUBLE_MIN) {
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
	    } else {
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    }
	    if (std::abs(z__[i__]) <= eps && std::abs(z__[i__ + 1]) <= eps) {
		isuppz[2] = i__ - 1;
		break;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
	}
    }

    return;

}


}
}
#include <cctype>
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarf,DLARF)(const char *side,
       int *m,
       int *n,
       double *v,
       int *incv,
       double *tau,
       double *c,
       int *ldc,
       double *work)
{
  const char ch=std::toupper(*side);
  double one = 1.0;
  double zero = 0.0;
  double minustau = -(*tau);
  int i1 = 1;


  if(ch=='L') {
    if(std::abs(*tau)>PLUMED_GMX_DOUBLE_MIN) {
      PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("T",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      PLUMED_BLAS_F77_FUNC(dger,DGER)(m,n,&minustau,v,incv,work,&i1,c,ldc);
    }
  } else {
    if(std::abs(*tau)>PLUMED_GMX_DOUBLE_MIN) {
      PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      PLUMED_BLAS_F77_FUNC(dger,DGER)(m,n,&minustau,work,&i1,v,incv,c,ldc);
    }
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)(const char *side, 
	const char *trans, 
	const char *direct, 
	const char *storev, 
	int *m, 
	int *n, 
	int *k, 
	double *v, 
	int *ldv, 
	double *t, 
	int *ldt, 
	double *c__,
	int *ldc, 
	double *work, 
	int *ldwork)
{
    int c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    int i__, j;
    char transt[1];
    int c__1 = 1;
    double one = 1.0;
    double minusone = -1.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }
    if (*trans=='N' || *trans=='n') {
      *(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }
    
    if (*storev=='C' || *storev=='c') {

	if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &one, &work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &one, &c__[(*k + 1) * c_dim1 + 1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}
	    }

	} else {

	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc)
			    ;
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    one, &work[work_offset], ldwork);
		}
		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    one, &c__[c_offset], ldc)
			    ;
		}
		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}
	    }
	}

    } else  if (*storev=='r' || *storev=='R') {
      if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}
		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &one, &work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &c__[(*k + 1) * c_dim1 
			    + 1], ldc);
		}
		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    }

	} else {

	    if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &one, &c__[c_offset], ldc);
		}

		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
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

  if(std::abs(xnorm)<PLUMED_GMX_DOUBLE_MIN) {
    *tau = 0.0;
  } else {

    t = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(alpha,&xnorm);

    if(*alpha<0)
      beta = t;
    else
      beta = -t;

    minval = PLUMED_GMX_DOUBLE_MIN;
    
    safmin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS) / PLUMED_GMX_DOUBLE_EPS;

        
    if(std::abs(beta)<safmin) {

      knt = 0;
      rsafmn = 1.0 / safmin;
      
      while(std::abs(beta)<safmin) {
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
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)(const char *direct, 
	const char *storev, 
	int *n, 
	int *k, 
	double *v, 
	int *ldv, 
	double *tau, 
	double *t, 
	int *ldt)
{
    /* System generated locals */
    int t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    int i__, j;
    double vii;
    int c__1 = 1;
    double zero = 0.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    if (*n == 0) {
	return;
    }

    if (*direct=='F' || *direct=='f') {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (std::abs(tau[i__])<PLUMED_GMX_DOUBLE_MIN) {

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		vii = v[i__ + i__ * v_dim1];
		v[i__ + i__ * v_dim1] = 1.;
		if (*storev=='C' || *storev=='c') {

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    d__1 = -tau[i__];
		    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1],
			     ldv, &v[i__ + i__ * v_dim1], &c__1, &zero, &t[
			    i__ * t_dim1 + 1], &c__1);
		} else {

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    d__1 = -tau[i__];
		    PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &d__1, &v[i__ * 
			    v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &
			    zero, &t[i__ * t_dim1 + 1], &c__1);
		}
		v[i__ + i__ * v_dim1] = vii;


		i__2 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (std::abs(tau[i__])<PLUMED_GMX_DOUBLE_MIN) {

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		if (i__ < *k) {
		    if (*storev=='C' || *storev=='c') {
			vii = v[*n - *k + i__ + i__ * v_dim1];
			v[*n - *k + i__ + i__ * v_dim1] = 1.;

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			d__1 = -tau[i__];
			PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1) 
				* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &
				c__1, &zero, &t[i__ + 1 + i__ * t_dim1], &
				c__1);
			v[*n - *k + i__ + i__ * v_dim1] = vii;
		    } else {
			vii = v[i__ + (*n - *k + i__) * v_dim1];
			v[i__ + (*n - *k + i__) * v_dim1] = 1.;

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			d__1 = -tau[i__];
			PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("No transpose", &i__1, &i__2, &d__1, &v[i__ + 
				1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &
				zero, &t[i__ + 1 + i__ * t_dim1], &c__1);
			v[i__ + (*n - *k + i__) * v_dim1] = vii;
		    }

		    i__1 = *k - i__;
		    PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)("Lower", "No transpose", "Non-unit", &i__1, &t[i__ 
			    + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
			     t_dim1], &c__1)
			    ;
		}
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    }
    return;


}
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarnv,DLARNV)(int *idist, 
	int *iseed, 
	int *n, 
	double *x)
{
    int i__1, i__2, i__3;

    int i__;
    double u[128];
    int il, iv, il2;

    --x;
    --iseed;

    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
	i__2 = 64, i__3 = *n - iv + 1;
	il = (i__2<i__3) ? i__2 : i__3;
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

	PLUMED_BLAS_F77_FUNC(dlaruv,DLARUV)(&iseed[1], &il2, u);

	if (*idist == 1) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1];
	    }
	} else if (*idist == 2) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
	    }
	} else if (*idist == 3) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
                x[iv + i__ - 1] =  std::sqrt(std::log(u[(i__ << 1) - 2]) * -2.) * 
		  std::cos(u[(i__ << 1) - 1] * (double)6.2831853071795864769252867663);
	    }
	}
    }
    return;

}
}
}
#include <cmath>

#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarrbx,DLARRBX)(int *n, 
	 double *d__, 
	 double *l, 
	 double *ld, 
	 double *lld, 
	 int *ifirst, 
	 int *ilast, 
	 double *rtol1, 
	 double *rtol2, 
	 int *offset, 
	 double *w, 
	 double *wgap, 
	 double *werr, 
	 double *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    double d__1, d__2;

    int i__, j, k, p;
    double s;
    int i1, i2, ii, kk;
    double fac, gap, mid;
    int cnt;
    double tmp, left;
    int nint, prev, next, nleft;
    double right, width, dplus;
    int nright, olnint;
    k = 0;
    right = 0.0;

    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    *info = 0;
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    i1 = *ifirst;
    i2 = *ifirst;
    prev = 0;
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	iwork[k - 1] = 1;
	i2 = i__;
    }

    i__ = i1;
    nint = 0;
L30:
    if (i__ <= i2) {
	ii = i__ - *offset;
	if (iwork[(i__ << 1) - 1] == 1) {
	    fac = 1.;
	    left = w[ii] - werr[ii];


L40:
	    if (i__ > i1 && left <= right) {
		left = right;
		cnt = i__ - 1;
	    } else {
		s = -left;
		cnt = 0;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    s = s * lld[j] / dplus - left;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
        if (std::isnan(s)) {

		    cnt = 0;
		    s = -left;
		    i__1 = *n - 1;
		    for (j = 1; j <= i__1; ++j) {
			dplus = d__[j] + s;
			if (dplus < 0.) {
			    ++cnt;
			}
			tmp = lld[j] / dplus;
			if (std::abs(tmp)<PLUMED_GMX_DOUBLE_MIN) {
			    s = lld[j] - left;
			} else {
			    s = s * tmp - left;
			}
		    }
		    dplus = d__[*n] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		if (cnt > i__ - 1) {
		    left -= werr[ii] * fac;
		    fac *= 2.;
		    goto L40;
		}
	    }
	    nleft = cnt + 1;
	    i1 = (i1<nleft) ? i1 : nleft;
	    fac = 1.;
	    right = w[ii] + werr[ii];
L60:
	    s = -right;
	    cnt = 0;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		dplus = d__[j] + s;
		s = s * lld[j] / dplus - right;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    if (std::isnan(s)) {

		cnt = 0;
		s = -right;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		    tmp = lld[j] / dplus;
		    if (std::abs(tmp)<PLUMED_GMX_DOUBLE_MIN) {
			s = lld[j] - right;
		    } else {
			s = s * tmp - right;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt < i__) {
		right += werr[ii] * fac;
		fac *= 2.;
		goto L60;
	    }
	    cnt = (cnt<i2) ? cnt : i2;
	    ++nint;
	    k = nleft << 1;
	    work[k - 1] = left;
	    work[k] = right;
	    i__ = cnt + 1;
	    iwork[k - 1] = i__;
	    iwork[k] = cnt;
	    if (prev != nleft - 1) {
		work[k - 2] = left;
	    }
	    prev = nleft;
	} else {
	    right = work[i__ * 2];

	    ++iwork[k - 1];
	    prev = i__;
	    ++i__;
	}
	goto L30;
    }
    if (i__ <= *n && iwork[(i__ << 1) - 1] != -1) {
	work[(i__ << 1) - 1] = work[prev * 2];
    }

L80:
    prev = i1 - 1;
    olnint = nint;
    i__ = i1;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
	k = i__ << 1;
	left = work[k - 1];
	right = work[k];
	next = iwork[k - 1];
	nright = iwork[k];
	mid = (left + right) * .5;
	width = right - mid;
	d__1 = std::abs(left);
	d__2 = std::abs(right);
	tmp = (d__1>d__2) ? d__1 : d__2;

	gap = 0.;
	if (i__ == nright) {
	    if (prev > 0 && next <= *n) {
		d__1 = left - work[k - 2], d__2 = work[k + 1] - right;
		gap = (d__1<d__2) ? d__1 : d__2;
	    } else if (prev > 0) {
		gap = left - work[k - 2];
	    } else if (next <= *n) {
		gap = work[k + 1] - right;
	    }
	}
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
	if (width < ((d__1>d__2) ? d__1 : d__2)) {
	    --nint;
	    iwork[k - 1] = 0;
	    kk = k;
	    i__2 = nright;
	    for (j = i__ + 1; j <= i__2; ++j) {
		kk += 2;
		iwork[kk - 1] = 0;
		work[kk - 1] = left;
		work[kk] = right;
		wgap[j - 1 - *offset] = 0.;
	    }
	    if (i1 == i__) {
		i1 = next;
	    } else {
		iwork[(prev << 1) - 1] = next;
	    }
	    i__ = next;
	    continue;
	}
	prev = i__;

	s = -mid;
	cnt = 0;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    dplus = d__[j] + s;
	    s = s * lld[j] / dplus - mid;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	dplus = d__[*n] + s;
	if (dplus < 0.) {
	    ++cnt;
	}
	if (std::isnan(s)) {
	    cnt = 0;
	    s = -mid;
	    i__2 = *n - 1;
	    for (j = 1; j <= i__2; ++j) {
		dplus = d__[j] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		tmp = lld[j] / dplus;
		if (std::abs(tmp)<PLUMED_GMX_DOUBLE_MIN) {
		    s = lld[j] - mid;
		} else {
		    s = s * tmp - mid;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	i__2 = i__ - 1, i__3 = (nright<cnt) ? nright : cnt;
	cnt = (i__2>i__3) ? i__2 : i__3;
	if (cnt == i__ - 1) {
	    work[k - 1] = mid;
	} else if (cnt == nright) {
	    work[k] = mid;
	} else {
	    iwork[k] = cnt;
	    ++cnt;
	    iwork[k - 1] = cnt;
	    kk = cnt << 1;
	    iwork[kk - 1] = next;
	    iwork[kk] = nright;
	    work[k] = mid;
	    work[kk - 1] = mid;
	    work[kk] = right;
	    prev = cnt;
	    if (cnt - 1 > i__) {
		work[kk - 2] = mid;
	    }
	    if (cnt > *ifirst && cnt <= *ilast) {
		++nint;
	    } else if (cnt <= *ifirst) {
		i1 = cnt;
	    }
	}
	i__ = next;
    }
    if (nint > 0) {
	goto L80;
    }
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
	if (iwork[k - 1] != -1) {
	    w[ii] = (work[k - 1] + work[k]) * .5;
	    werr[ii] = work[k] - w[ii];
	    if (i__ != *ilast) {
		wgap[ii] = work[k + 1] - work[k];
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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarrex,DLARREX)(const char *range,
	 int *n, 
	 double *vl, 
	 double *vu, 
	 int *il, 
	 int *iu, 
	 double *d__, 
	 double *e, 
	 double *tol, 
	 int *nsplit, 
	 int *isplit, 
	 int *m, 
	 double *w, 
	 int *iblock, 
	 int *indexw, 
	 double *gersch, 
	 double *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    double d__1, d__2;
    int c__1 = 1;
    int c__0 = 0;

    int i__, j, k;
    double s, gl;
    int in;
    double gu;
    int cnt;
    double eps, tau, nrm, tmp, vvl, vvu, offd;
    int iend, jblk, till, itmp;
    double rtol, delta, sigma;
    int iinfo;
    double width;
    int ibegin;
    int irange;
    double sgndef;
    int maxcnt;
    --iwork;
    --work;
    --gersch;
    --indexw;
    --iblock;
    --w;
    --isplit;
    --e;
    --d__;

    sigma = 0;
    irange = 0;
    sgndef = 0;
    maxcnt = 0;

    *info = 0;

    if (*range=='A' || *range=='a')
	irange = 1;
    else if (*range=='V' || *range=='v')
	irange = 2;
    else if (*range=='I' || *range=='i')
	irange = 3;
    

    *m = 0;
    eps = PLUMED_GMX_DOUBLE_EPS;

    *nsplit = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__]) <= *tol) {
	    isplit[*nsplit] = i__;
	    ++(*nsplit);
	}
    }
    isplit[*nsplit] = *n;

    ibegin = 1;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];
	if (ibegin == iend) {
	    ++(*m);
	    w[*m] = d__[ibegin];
	    iblock[*m] = jblk;
	    indexw[*m] = 1;
	    e[iend] = 0.;
	    ibegin = iend + 1;
	    goto L170;
	}
	in = iend - ibegin + 1;

	gl = d__[ibegin] - std::abs(e[ibegin]);
	gu = d__[ibegin] + std::abs(e[ibegin]);
	gersch[(ibegin << 1) - 1] = gl;
	gersch[ibegin * 2] = gu;
	gersch[(iend << 1) - 1] = d__[iend] - std::abs(e[iend - 1]);
	gersch[iend * 2] = d__[iend] + std::abs(e[iend - 1]);
	d__1 = gersch[(iend << 1) - 1];
	gl = (d__1<gl) ? d__1 : gl;
	d__1 = gersch[iend * 2];
	gu = (d__1>gu) ? d__1 : gu;
	i__2 = iend - 1;
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
	    offd = std::abs(e[i__ - 1]) + std::abs(e[i__]);
	    gersch[(i__ << 1) - 1] = d__[i__] - offd;
	    d__1 = gersch[(i__ << 1) - 1];
	    gl = (d__1<gl) ? d__1 : gl;
	    gersch[i__ * 2] = d__[i__] + offd;
	    d__1 = gersch[i__ * 2];
	    gu = (d__1>gu) ? d__1 : gu;
	}
	d__1 = std::abs(gl), d__2 = std::abs(gu);
	nrm = (d__1>d__2) ? d__1 : d__2;

	width = gu - gl;
	i__2 = iend - 1;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    work[i__] = e[i__] * e[i__];
	}
	for (j = 1; j <= 2; ++j) {
	    if (j == 1) {
		tau = gl + width * .25;
	    } else {
		tau = gu - width * .25;
	    }
	    tmp = d__[ibegin] - tau;
	    if (tmp < 0.) {
		cnt = 1;
	    } else {
		cnt = 0;
	    }
	    i__2 = iend;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		tmp = d__[i__] - tau - work[i__ - 1] / tmp;
		if (tmp < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt == 0) {
		gl = tau;
	    } else if (cnt == in) {
		gu = tau;
	    }
	    if (j == 1) {
		maxcnt = cnt;
		sigma = gl;
		sgndef = 1.;
	    } else {
		if (in - cnt > maxcnt) {
		    sigma = gu;
		    sgndef = -1.;
		}
	    }
	}

	work[in * 3] = 1.;
	delta = eps;
	tau = sgndef * nrm;
L60:
	sigma -= delta * tau;
	work[1] = d__[ibegin] - sigma;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(in << 1) + i__] = 1. / work[i__];
	    tmp = e[j] * work[(in << 1) + i__];
	    work[i__ + 1] = d__[j + 1] - sigma - tmp * e[j];
	    work[in + i__] = tmp;
	    ++j;
	}
	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
        if (tmp < 0. || std::abs(work[(in << 1) + i__])<PLUMED_GMX_DOUBLE_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L60;
	    }
	}

	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[in * 3 + i__] = work[i__] * work[in + i__];
	    work[(in << 2) + i__] = work[in * 3 + i__] * work[in + i__];
	}
	if (sgndef > 0.) {
	    cnt = 1;
	    work[1] = (gl + gu) / 2. - sigma;
	    work[in + 1] = 0.;
	    work[(in << 1) + 1] = (gu - gl) / 2.;
	} else {
	    cnt = in;
	    work[in] = (gl + gu) / 2. - sigma;
	    work[in * 2] = 0.;
	    work[in * 3] = (gu - gl) / 2.;
	}
	rtol = eps * 4.;
	PLUMED_BLAS_F77_FUNC(dlarrbx,DLARRBX)(&in, &d__[ibegin], &e[ibegin], &work[in * 3 + 1], &work[(in <<
		 2) + 1], &cnt, &cnt, &rtol, &rtol, &c__0, &work[1], &work[in 
		+ 1], &work[(in << 1) + 1], &work[in * 5 + 1], &iwork[1], &
		iinfo);
	if (sgndef > 0.) {
	    tau = work[1] - work[(in << 1) + 1];
	} else {
	    tau = work[in] + work[in * 3];
	}

	work[in * 3] = 1.;
	delta = eps * 2.;
L100:
	tau *= 1. - delta;

	s = -tau;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__] = d__[j] + s;
	    work[(in << 1) + i__] = 1. / work[i__];
	    work[in + i__] = e[j] * d__[j] * work[(in << 1) + i__];
	    s = s * work[in + i__] * e[j] - tau;
	    ++j;
	}
	work[in] = d__[iend] + s;

	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
	    if (tmp < 0. || std::abs(work[(in << 1) + i__])<PLUMED_GMX_DOUBLE_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L100;
	    }
	}

	sigma += tau;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	e[iend] = sigma;
	tmp = (double) in * 4. * eps * (std::abs(sigma) + std::abs(tau));
	i__2 = iend;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    gersch[(i__ << 1) - 1] = gersch[(i__ << 1) - 1] - sigma - tmp;
	    gersch[i__ * 2] = gersch[i__ * 2] - sigma + tmp;
	}

	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(i__ << 1) - 1] = std::abs(d__[j]);
	    work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
	    ++j;
	}
	work[(in << 1) - 1] = std::abs(d__[iend]);

	PLUMED_BLAS_F77_FUNC(dlasq2,DLASQ2)(&in, &work[1], info);
	if (*info != 0) {
	    return;
	}

	if (sgndef > 0.) {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = work[in - i__ + 1];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	} else {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = -work[i__];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	}
	ibegin = iend + 1;
L170:
	;
    }
    if (irange == 2) {
	*m = 0;
	ibegin = 1;
	i__1 = *nsplit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iend = isplit[i__];
	    vvl = *vl - e[iend];
	    vvu = *vu - e[iend];
	    i__2 = iend;
	    for (j = ibegin; j <= i__2; ++j) {
		if (vvl <= w[j] && w[j] <= vvu) {
		    ++(*m);
		    w[*m] = w[j];
		    iblock[*m] = i__;
		    indexw[*m] = j - ibegin + 1;
		}
	    }
	    ibegin = iend + 1;
	}
    } else if (irange == 3) {
	*m = *iu - *il + 1;
	if (*nsplit == 1) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = w[*il + i__ - 1];
		indexw[i__] = *il + i__ - 1;
	    }
	} else {
	    ibegin = 1;
	    i__1 = *nsplit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iend = isplit[i__];
		i__2 = iend;
		for (j = ibegin; j <= i__2; ++j) {
		    work[j] = w[j] + e[iend];
		}
		ibegin = iend + 1;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = i__;
		iwork[*n + i__] = iblock[i__];
	    }
	    PLUMED_BLAS_F77_FUNC(dlasrt2,DLASRT2)("I", n, &work[1], &iwork[1], &iinfo);
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		itmp = iwork[*il + i__ - 1];
		work[i__] = w[itmp];
		iblock[i__] = iwork[*n + itmp];
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[*n + i__] = iwork[*il + i__ - 1];
		iwork[i__] = i__;
	    }
	    PLUMED_BLAS_F77_FUNC(ilasrt2,ILASRT2)("I", m, &iblock[1], &iwork[1], &iinfo);
	    j = 1;
	    itmp = iblock[j];
	    cnt = iwork[*n + iwork[j]];
	    if (itmp == 1) {
		ibegin = 1;
	    } else {
		ibegin = isplit[itmp - 1] + 1;
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = work[iwork[i__]];
		if (iblock[i__] != itmp || i__ == *m) {
		    if (iblock[i__] == itmp) {
			till = *m;
		    } else {
			till = i__ - 1;
		    }
		    i__2 = till - j + 1;
		    PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("I", &i__2, &w[j], &iinfo);
		    cnt = cnt - ibegin + 1;
		    i__2 = till;
		    for (k = j; k <= i__2; ++k) {
			indexw[k] = cnt + k - j;
		    }
		    j = i__;
		    itmp = iblock[j];
		    cnt = iwork[*n + iwork[j]];
		    ibegin = isplit[itmp - 1] + 1;
		    if (i__ == *m && till < *m) {
			indexw[*m] = cnt - ibegin + 1;
		    }
		} else {
		    i__2 = cnt, i__3 = iwork[*n + iwork[i__]];
		    cnt = (i__2<i__3) ? i__2 : i__3;
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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarrfx,DLARRFX)(int *n, 
	double *d__, 
	double *l, 
	double *ld, 
	double *lld, 
	int *ifirst, 
	int *ilast, 
	double *w, 
	double *sigma, 
	double *dplus, 
	double *lplus, 
	double *work,
	int *info)
{
    int i1 = 1;
    int i__1;
    double d__2, d__3;

    int i__;
    double s, eps, tmp, dmax1, dmax2, delta;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;
    *info = 0;
    eps = PLUMED_GMX_DOUBLE_EPS;
    *sigma = w[*ifirst];
    delta = eps * 2.;

L10:
    s = -(*sigma);
    dplus[1] = d__[1] + s;
    dmax1 = std::abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	s = s * lplus[i__] * l[i__] - *sigma;
	dplus[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax1, d__3 = std::abs(dplus[i__ + 1]);
	dmax1 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax1)) {
	*sigma -= std::abs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    s = -tmp;
    work[1] = d__[1] + s;
    dmax2 = std::abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	s = s * work[*n + i__] * l[i__] - tmp;
	work[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax2, d__3 = std::abs(work[i__ + 1]);
	dmax2 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax2)) {
	tmp += std::abs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }
    if (dmax2 < dmax1) {
	*sigma = tmp;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &work[1], &i1, &dplus[1], &i1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &work[*n + 1], &i1, &lplus[1], &i1);
    }

    return;
}
}
}
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlarrvx,DLARRVX)(int *n, 
	double *d__, 
	double *l, 
	int *isplit,
	int *m, 
	double *w,
	int *iblock, 
	int *indexw, 
	double *gersch, 
	double *tol, 
	double *z__, 
	int *ldz, 
	int *isuppz, 
	double *work, 
	int *iwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    double d__1, d__2;
    double c_b5 = 0.;
    int c__1 = 1;
    int c__2 = 2;

    int i__, j, k, p, q;
    int im, in;
    double gap, eps, tmp;
    int zto;
    double ztz;
    int iend, jblk;
    int wend, iter, temp[1], ktot;
    int itmp1, itmp2;
    int indld;
    double sigma;
    int ndone, iinfo, iindr;
    double resid;
    int nomgs;
    int nclus;
    int zfrom, iindc1, iindc2;
    double lambda;
    int ibegin;
    int indgap, indlld;
    double mingma;
    int oldien, oldncl, wbegin;
    double relgap;
    int oldcls;
    int ndepth, inderr, iindwk;
    int newcls, oldfst;
    double minrgp=0.0;
    int indwrk, oldlst;
    double reltol;
    int newfrs, newftt, parity;
    double mgstol, nrminv, rqcorr;
    int newlst, newsiz;


    --d__;
    --l;
    --isplit;
    --w;
    --iblock;
    --indexw;
    --gersch;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    inderr = *n;
    indld = *n << 1;
    indlld = *n * 3;
    indgap = *n << 2;
    indwrk = *n * 5 + 1;

    iindr = *n;
    iindc1 = *n << 1;
    iindc2 = *n * 3;
    iindwk = (*n << 2) + 1;

    eps = PLUMED_GMX_DOUBLE_EPS;

    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("Full", n, m, &c_b5, &c_b5, &z__[z_offset], ldz);
    mgstol = eps * 100.;

    ibegin = 1;
    wbegin = 1;
    i__1 = iblock[*m];
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];

	wend = wbegin - 1;
L171:
	if (wend < *m) {
	    if (iblock[wend + 1] == jblk) {
		++wend;
		goto L171;
	    }
	}
	if (wend < wbegin) {
	    ibegin = iend + 1;
	    continue;
	}

	if (ibegin == iend) {
	    z__[ibegin + wbegin * z_dim1] = 1.;
	    isuppz[(wbegin << 1) - 1] = ibegin;
	    isuppz[wbegin * 2] = ibegin;
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	    continue;
	}
	oldien = ibegin - 1;
	in = iend - oldien;
	d__1 = .001, d__2 = 1. / (double) in;
	reltol = (d__1<d__2) ? d__1 : d__2;
	im = wend - wbegin + 1;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&im, &w[wbegin], &c__1, &work[1], &c__1);
	i__2 = im - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[inderr + i__] = eps * std::abs(work[i__]);
	    work[indgap + i__] = work[i__ + 1] - work[i__];
	}
	work[inderr + im] = eps * std::abs(work[im]);
	d__2 = std::abs(work[im]);
	work[indgap + im] = (d__2>eps) ? d__2 : eps;
	ndone = 0;

	ndepth = 0;
	parity = 1;
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;

L40:
	if (ndone < im) {
	    oldncl = nclus;
	    nclus = 0;
	    parity = 1 - parity;
	    if (parity == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
	    i__2 = oldncl;
	    for (i__ = 1; i__ <= i__2; ++i__) {

		j = oldcls + (i__ << 1);
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
		    j = wbegin + oldfst - 1;
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
		    i__3 = in - 1;
		    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
		    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz);
		}
		k = ibegin;
		i__3 = in - 1;
		for (j = 1; j <= i__3; ++j) {
		    tmp = d__[k] * l[k];
		    work[indld + j] = tmp;
		    work[indlld + j] = tmp * l[k];
		    ++k;
		}
		if (ndepth > 0) {

		    p = indexw[wbegin - 1 + oldfst];
		    q = indexw[wbegin - 1 + oldlst];
		    d__1 = eps * 4.;
		    i__3 = p - oldfst;
		    PLUMED_BLAS_F77_FUNC(dlarrbx,DLARRBX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 1], &
			    work[indlld + 1], &p, &q, &reltol, &d__1, &i__3, &
			    work[1], &work[indgap + 1], &work[inderr + 1], &
			    work[indwrk + in], &iwork[iindwk], &iinfo);
		}
		newfrs = oldfst;
		i__3 = oldlst;
		for (j = oldfst; j <= i__3; ++j) {
		    if (j == oldlst || work[indgap + j] >= 
			reltol * std::abs(work[j])) {
			newlst = j;
		    } else {

			relgap = work[indgap + j] / std::abs(work[j]);
			if (j == newfrs) {
			    minrgp = relgap;
			} else {
			    minrgp = (minrgp<relgap) ? minrgp : relgap;
			}
			continue;
		    }
		    newsiz = newlst - newfrs + 1;
		    newftt = wbegin + newfrs - 1;
		    nomgs = newsiz == 1 || newsiz > 1 || minrgp < mgstol;
		    if (newsiz > 1 && nomgs) {

			PLUMED_BLAS_F77_FUNC(dlarrfx,DLARRFX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				1], &work[indlld + 1], &newfrs, &newlst, &
				work[1], &sigma, &z__[ibegin + newftt * 
				z_dim1], &z__[ibegin + (newftt + 1) * z_dim1],
				 &work[indwrk], info);
			if (*info == 0) {
			    tmp = eps * std::abs(sigma);
			    i__4 = newlst;
			    for (k = newfrs; k <= i__4; ++k) {
				work[k] -= sigma;
				d__1 = work[indgap + k];
				work[indgap + k] = (d__1>tmp) ? d__1 : tmp;
				work[inderr + k] += tmp;
			    }
			    ++nclus;
			    k = newcls + (nclus << 1);
			    iwork[k - 1] = newfrs;
			    iwork[k] = newlst;
			} else {
			    *info = 0;
			    if (minrgp < mgstol) {
				work[indwrk] = d__[ibegin];
				i__4 = in - 1;
				for (k = 1; k <= i__4; ++k) {
				    work[indwrk + k] = d__[ibegin + k] + work[
					    indlld + k];
				}
				i__4 = newsiz;
				for (k = 1; k <= i__4; ++k) {
				    iwork[iindwk + k - 1] = 1;
				}
				i__4 = newlst;
				for (k = newfrs; k <= i__4; ++k) {
				    isuppz[2*(oldien + k) - 1] = 1;
				    isuppz[(oldien + k) * 2] = in;
				}
				temp[0] = in;
				PLUMED_BLAS_F77_FUNC(dstein,DSTEIN)(&in, &work[indwrk], &work[indld + 1], 
					&newsiz, &work[newfrs], &iwork[iindwk]
					, temp, &z__[ibegin + newftt * z_dim1]
					, ldz, &work[indwrk + in], &iwork[
					iindwk + in], &iwork[iindwk + (in*2)], &iinfo);
				if (iinfo != 0) {
				    *info = 2;
				    return;
				}
				ndone += newsiz;
			    }
			}
		    } else {
			ktot = newftt;
			i__4 = newlst;
			for (k = newfrs; k <= i__4; ++k) {
			    iter = 0;
L90:
			    lambda = work[k];

			    PLUMED_BLAS_F77_FUNC(dlar1vx,DLAR1VX)(&in, &c__1, &in, &lambda, &d__[ibegin], &
				    l[ibegin], &work[indld + 1], &work[indlld 
				    + 1], &w[wbegin + k - 1], &gersch[(oldien 
				    << 1) + 1], &z__[ibegin + ktot * z_dim1], 
				    &ztz, &mingma, &iwork[iindr + ktot], &
				    isuppz[(ktot << 1) - 1], &work[indwrk]);
			    tmp = 1. / ztz;
			    nrminv =  std::sqrt(tmp);
			    resid = std::abs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
				d__1 = work[indgap + k - 1], d__2 = work[
					indgap + k];
				gap = (d__1<d__2) ? d__1 : d__2;
			    }
			    ++iter;
			    if (resid > *tol * gap && std::abs(rqcorr) > eps * 4. *
				     std::abs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    goto L90;
				}
			    }
			    iwork[ktot] = 1;
			    if (newsiz == 1) {
				++ndone;
			    }
			    zfrom = isuppz[(ktot << 1) - 1];
			    zto = isuppz[ktot * 2];
			    i__5 = zto - zfrom + 1;
			    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__5, &nrminv, &z__[ibegin + zfrom - 1 + 
				    ktot * z_dim1], &c__1);
			    ++ktot;
			}
			if (newsiz > 1) {
			    itmp1 = isuppz[(newftt << 1) - 1];
			    itmp2 = isuppz[newftt * 2];
			    ktot = oldien + newlst;
			    i__4 = ktot;
			    for (p = newftt + 1; p <= i__4; ++p) {
				i__5 = p - 1;
				for (q = newftt; q <= i__5; ++q) {
				    tmp = -PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&in, &z__[ibegin + p * 
					    z_dim1], &c__1, &z__[ibegin + q * 
					    z_dim1], &c__1);
				    PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&in, &tmp, &z__[ibegin + q * 
					    z_dim1], &c__1, &z__[ibegin + p * 
					    z_dim1], &c__1);
				}
				tmp = 1. / PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(&in, &z__[ibegin + p * 
					z_dim1], &c__1);
				PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&in, &tmp, &z__[ibegin + p * z_dim1], &
					c__1);
				i__5 = itmp1, i__6 = isuppz[(p << 1) - 1];
				itmp1 = (i__5<i__6) ? i__5 : i__6;
				i__5 = itmp2, i__6 = isuppz[p * 2];
				itmp2 = (i__5>i__6) ? i__5 : i__6;
			    }
			    i__4 = ktot;
			    for (p = newftt; p <= i__4; ++p) {
				isuppz[(p << 1) - 1] = itmp1;
				isuppz[p * 2] = itmp2;
			    }
			    ndone += newsiz;
			}
		    }
		    newfrs = j + 1;
		}
	    }
	    ++ndepth;
	    goto L40;
	}
	j = wbegin << 1;
	i__2 = wend;
	for (i__ = wbegin; i__ <= i__2; ++i__) {
	    isuppz[j - 1] += oldien;
	    isuppz[j] += oldien;
	    j += 2;

	}
	ibegin = iend + 1;
	wbegin = wend + 1;
    }

    return;

} 
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(double *f,
	double *g,
	double *cs,
	double *sn,
	double *r)
{
  double minval,safemin, safemin2, safemx2, eps;
  double f1,g1,f1a,g1a,scale;
  int i,n,count;

  eps = PLUMED_GMX_DOUBLE_EPS;
  minval = PLUMED_GMX_DOUBLE_MIN;
  safemin = minval*(1.0+eps);
  n = static_cast<int>(0.5*std::log( safemin/eps ) / std::log(2.0));
  safemin2 = std::pow(2.0,static_cast<double>(n));

  safemx2 = 1.0 / safemin2;

  if(std::abs(*g)<PLUMED_GMX_DOUBLE_MIN) {
    *cs = 1.0;
    *sn = 0.0;
    *r = *f;
  } else if (std::abs(*f)<PLUMED_GMX_DOUBLE_MIN) {
    *cs = 0.0;
    *sn = 1.0;
    *r = *g;
  } else {
    f1 = *f;
    g1 = *g;
    f1a = std::abs(f1);
    g1a = std::abs(g1);
    scale = (f1a > g1a) ? f1a : g1a;
    if(scale >= safemx2) {
      count = 0;
      while(scale >= safemx2) {
	count++;
	f1 *= safemin2;
	g1 *= safemin2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemx2;
    } else if (scale<=safemin2) {
      count = 0;
      while(scale <= safemin2) {
	count++;
	f1 *= safemx2;
	g1 *= safemx2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemin2;
    } else {
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
    }
    if(std::abs(*f)>std::abs(*g) && *cs<0.0) {
      *cs *= -1.0;
      *sn *= -1.0;
      *r  *= -1.0;
    }
  }
  return;
}
      
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaruv,DLARUV)(int *iseed, int *n, double *x)
{
  const int
    mm[512] = {
      494,2637,255,2008,1253,
      3344,4084,1739,3143,3468,688,1657,1238,3166,1292,3422,1270,2016,
      154,2862,697,1706,491,931,1444,444,3577,3944,2184,1661,3482,657,
      3023,3618,1267,1828,164,3798,3087,2400,2870,3876,1905,1593,1797,
      1234,3460,328,2861,1950,617,2070,3331,769,1558,2412,2800,189,287,
      2045,1227,2838,209,2770,3654,3993,192,2253,3491,2889,2857,2094,
      1818,688,1407,634,3231,815,3524,1914,516,164,303,2144,3480,119,
      3357,837,2826,2332,2089,3780,1700,3712,150,2000,3375,1621,3090,
      3765,1149,3146,33,3082,2741,359,3316,1749,185,2784,2202,2199,1364,
      1244,2020,3160,2785,2772,1217,1822,1245,2252,3904,2774,997,2573,
      1148,545,322,789,1440,752,2859,123,1848,643,2405,2638,2344,46,
      3814,913,3649,339,3808,822,2832,3078,3633,2970,637,2249,2081,4019,
      1478,242,481,2075,4058,622,3376,812,234,641,4005,1122,3135,2640,
      2302,40,1832,2247,2034,2637,1287,1691,496,1597,2394,2584,1843,336,
      1472,2407,433,2096,1761,2810,566,442,41,1238,1086,603,840,3168,
      1499,1084,3438,2408,1589,2391,288,26,512,1456,171,1677,2657,2270,
      2587,2961,1970,1817,676,1410,3723,2803,3185,184,663,499,3784,1631,
      1925,3912,1398,1349,1441,2224,2411,1907,3192,2786,382,37,759,2948,
      1862,3802,2423,2051,2295,1332,1832,2405,3638,3661,327,3660,716,
      1842,3987,1368,1848,2366,2508,3754,1766,3572,2893,307,1297,3966,
      758,2598,3406,2922,1038,2934,2091,2451,1580,1958,2055,1507,1078,
      3273,17,854,2916,3971,2889,3831,2621,1541,893,736,3992,787,2125,
      2364,2460,257,1574,3912,1216,3248,3401,2124,2762,149,2245,166,466,
      4018,1399,190,2879,153,2320,18,712,2159,2318,2091,3443,1510,449,
      1956,2201,3137,3399,1321,2271,3667,2703,629,2365,2431,1113,3922,
      2554,184,2099,3228,4012,1921,3452,3901,572,3309,3171,817,3039,
      1696,1256,3715,2077,3019,1497,1101,717,51,981,1978,1813,3881,76,
      3846,3694,1682,124,1660,3997,479,1141,886,3514,1301,3604,1888,
      1836,1990,2058,692,1194,20,3285,2046,2107,3508,3525,3801,2549,
      1145,2253,305,3301,1065,3133,2913,3285,1241,1197,3729,2501,1673,
      541,2753,949,2361,1165,4081,2725,3305,3069,3617,3733,409,2157,
      1361,3973,1865,2525,1409,3445,3577,77,3761,2149,1449,3005,225,85,
      3673,3117,3089,1349,2057,413,65,1845,697,3085,3441,1573,3689,2941,
      929,533,2841,4077,721,2821,2249,2397,2817,245,1913,1997,3121,997,
      1833,2877,1633,981,2009,941,2449,197,2441,285,1473,2741,3129,909,
      2801,421,4073,2813,2337,1429,1177,1901,81,1669,2633,2269,129,1141,
      249,3917,2481,3941,2217,2749,3041,1877,345,2861,1809,3141,2825,
      157,2881,3637,1465,2829,2161,3365,361,2685,3745,2325,3609,3821,
      3537,517,3017,2141,1537 
    };

    int i__1;

    int i__, i1, i2, i3, i4, it1, it2, it3, it4;


    --iseed;
    --x;

    it1 = it2 = it3 = it4 = 0;

    i1 = iseed[1];
    i2 = iseed[2];
    i3 = iseed[3];
    i4 = iseed[4];

    i__1 = (*n<128) ? *n : 128;
    for (i__ = 1; i__ <= i__1; ++i__) {

	it4 = i4 * mm[i__ + 383];
	it3 = it4 / 4096;
	it4 -= it3 << 12;
	it3 = it3 + i3 * mm[i__ + 383] + i4 * mm[i__ + 255];
	it2 = it3 / 4096;
	it3 -= it2 << 12;
	it2 = it2 + i2 * mm[i__ + 383] + i3 * mm[i__ + 255] + 
	  i4 * mm[i__ + 127];
	it1 = it2 / 4096;
	it2 -= it1 << 12;
	it1 = it1 + i1 * mm[i__ + 383] + i2 * mm[i__ + 255] + 
	  i3 * mm[i__ +	127] + i4 * mm[i__ - 1];
	it1 %= 4096;

	x[i__] = ((double) it1 + ((double) it2 + ((double) it3 + (
		double) it4 * 2.44140625e-4) * 2.44140625e-4) * 
		2.44140625e-4) * 2.44140625e-4;
    }

    iseed[1] = it1;
    iseed[2] = it2;
    iseed[3] = it3;
    iseed[4] = it4;
    return;

} 
}
}
#include <cmath>
#include "real.h"

#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlas2,DLAS2)(double *f,
       double *g,
       double *h,
       double *ssmin,
       double *ssmax)
{
  double fa = std::abs(*f);
  double ga = std::abs(*g);
  double ha = std::abs(*h);
  double fhmin,fhmax,tmax,tmin,tmp1,tmp2;
  double as,at,au,c;

  fhmin = (fa<ha) ? fa : ha;
  fhmax = (fa>ha) ? fa : ha;
  
  if(std::abs(fhmin)<PLUMED_GMX_DOUBLE_MIN) {
    *ssmin = 0.0;
    if(std::abs(fhmax)<PLUMED_GMX_DOUBLE_MIN) 
      *ssmax = ga;
    else {
      tmax = (fhmax>ga) ? fhmax : ga;
      tmin = (fhmax<ga) ? fhmax : ga;
      tmp1 = tmin / tmax;
      tmp1 = tmp1 * tmp1;
      *ssmax = tmax* std::sqrt(1.0 + tmp1);
    }
  } else {
    if(ga<fhmax) {
      as = 1.0 + fhmin / fhmax;
      at = (fhmax-fhmin) / fhmax;
      au = (ga/fhmax);
      au = au * au;
      c = 2.0 / (  std::sqrt(as*as+au) + std::sqrt(at*at+au) );
      *ssmin = fhmin * c;
      *ssmax = fhmax / c;
    } else {
      au = fhmax / ga;
      if(std::abs(au)<PLUMED_GMX_DOUBLE_MIN) {
	*ssmin = (fhmin*fhmax)/ga;
	*ssmax = ga;
      } else {
	as = 1.0 + fhmin / fhmax;
	at = (fhmax-fhmin)/fhmax;
	tmp1 = as*au;
	tmp2 = at*au;
	c = 1.0 / (  std::sqrt(1.0+tmp1*tmp1) + std::sqrt(1.0+tmp2*tmp2));
	*ssmin = (fhmin*c)*au;
	*ssmin = *ssmin + *ssmin;
	*ssmax = ga / (c+c);
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

#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)(const char *type,
                        int *kl,
                        int *ku,
                        double *cfrom,
                        double *cto,
                        int *m,
                        int *n,
                        double *a,
                        int *lda,
                        int *info)
{
  const char ch=std::toupper(*type);
  int i,j,k,l,k1,k2,k3,k4;
  int done=0;
  double minval,smlnum,bignum;
  double cfromc, ctoc, cfrom1, cto1, mul;

  if(*n<=0 || *m<=0)
    return;

  minval = PLUMED_GMX_DOUBLE_MIN;
  smlnum = minval / PLUMED_GMX_DOUBLE_EPS;
  bignum = 1.0 / smlnum;

  cfromc = *cfrom;
  ctoc   = *cto;

  while(!done) {
    
    cfrom1 = cfromc * smlnum;
    cto1   = ctoc / bignum;

    if(std::abs(cfrom1)>std::abs(ctoc) && std::abs(ctoc)>PLUMED_GMX_DOUBLE_MIN) {
      mul = smlnum;
      done = 0;
      cfromc = cfrom1;
    } else if(std::abs(cto1)>std::abs(cfromc)) {
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
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd0,DLASD0)(int *n, 
	int *sqre, 
	double *d__, 
	double *e, 
	double *u, 
	int *ldu, 
	double *vt, 
	int *ldvt,
	int *smlsiz, 
	int *iwork,
	double *work, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int i__, j, m, i1, ic, lf, nd, ll, nl, nr, im1, ncc, nlf, nrf, 
	    iwk, lvl, ndb1, nlp1, nrp1;
    double beta;
    int idxq, nlvl;
    double alpha;
    int inode, ndiml, idxqc, ndimr, itemp, sqrei;
    int c__0 = 0;


    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --iwork;
    --work;

    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -2;
    }

    m = *n + *sqre;

    if (*ldu < *n) {
	*info = -6;
    } else if (*ldvt < m) {
	*info = -8;
    } else if (*smlsiz < 3) {
	*info = -9;
    }
    if (*info != 0) {
	return;
    }

    if (*n <= *smlsiz) {
	PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], 
		ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info);
	return;
    }

    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;
    PLUMED_BLAS_F77_FUNC(dlasdt,DLASDT)(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

    ndb1 = (nd + 1) / 2;
    ncc = 0;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {

	i1 = i__ - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nrp1 = nr + 1;
	nlf = ic - nl;
	nrf = ic + 1;
	sqrei = 1;
	PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[
		nlf + nlf * vt_dim1], ldvt, &u[nlf + nlf * u_dim1], ldu, &u[
		nlf + nlf * u_dim1], ldu, &work[1], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + nlf - 2;
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    iwork[itemp + j] = j;
	}
	if (i__ == nd) {
	    sqrei = *sqre;
	} else {
	    sqrei = 1;
	}
	nrp1 = nr + sqrei;
	PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[
		nrf + nrf * vt_dim1], ldvt, &u[nrf + nrf * u_dim1], ldu, &u[
		nrf + nrf * u_dim1], ldu, &work[1], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + ic;
	i__2 = nr;
	for (j = 1; j <= i__2; ++j) {
	    iwork[itemp + j - 1] = j;
	}
    }

    for (lvl = nlvl; lvl >= 1; --lvl) {

	if (lvl == 1) {
	    lf = 1;
	    ll = 1;
	} else {
	    i__1 = lvl - 1;
	    lf = (1 << i__1);
	    ll = (lf << 1) - 1;
	}
	i__1 = ll;
	for (i__ = lf; i__ <= i__1; ++i__) {
	    im1 = i__ - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    if (*sqre == 0 && i__ == ll) {
		sqrei = *sqre;
	    } else {
		sqrei = 1;
	    }
	    idxqc = idxq + nlf - 1;
	    alpha = d__[ic];
	    beta = e[ic];
	    PLUMED_BLAS_F77_FUNC(dlasd1,DLASD1)(&nl, &nr, &sqrei, &d__[nlf], &alpha, &beta, &u[nlf + nlf *
		     u_dim1], ldu, &vt[nlf + nlf * vt_dim1], ldvt, &iwork[
		    idxqc], &iwork[iwk], &work[1], info);
	    if (*info != 0) {
		return;
	    }
	}
    }

    return;

}
}
}
#include <cmath>
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

    d__1 = std::abs(*alpha);
    d__2 = std::abs(*beta);
    orgnrm = (d__1>d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(d__[i__]) > orgnrm) {
	    orgnrm = std::abs(d__[i__]);
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
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd2,DLASD2)(int *nl, 
                        int *nr, 
                        int *sqre, 
                        int *k, 
                        double *d__, 
                        double *z__, 
                        double *alpha, 
                        double *beta, 
                        double *u, 
                        int *ldu, 
                        double *vt, 
                        int *ldvt, 
                        double *dsigma, 
                        double *u2, 
                        int *ldu2, 
                        double *vt2, 
                        int *ldvt2, 
                        int *idxp, 
                        int *idx, 
                        int *idxc, 
                        int *idxq, 
                        int *coltyp, 
                        int *info)
{
    int u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset;
    int vt2_dim1, vt2_offset, i__1;
    double d__1, d__2;

    double c__;
    int i__, j, m, n;
    double s;
    int k2;
    double z1;
    int ct, jp;
    double eps, tau, tol;
    int psm[4], nlp1, nlp2, idxi, idxj;
    int ctot[4], idxjp;
    int jprev = 0;
    double hlftol;
    double zero = 0.0;
    int c__1 = 1;


    --d__;
    --z__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --dsigma;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxp;
    --idx;
    --idxc;
    --idxq;
    --coltyp;

    *info = 0;

    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
    z__[1] = z1;
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
    }

    i__1 = nlp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	coltyp[i__] = 1;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	coltyp[i__] = 2;
    }

    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	u2[i__ + u2_dim1] = z__[idxq[i__]];
	idxc[i__] = coltyp[idxq[i__]];
    }

    PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = u2[idxi + u2_dim1];
	coltyp[i__] = idxc[idxi];
    }

    eps = PLUMED_GMX_DOUBLE_EPS;
    d__1 = std::abs(*alpha), d__2 = std::abs(*beta);
    tol = (d__1 > d__2) ? d__1 : d__2;
    d__2 = std::abs(d__[n]);
    tol = eps * 8. * ((d__2 > tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (std::abs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    coltyp[j] = 4;
	    if (j == n) {
		goto L120;
	    }
	} else {
	    jprev = j;
	    goto L90;
	}
    }
L90:
    j = jprev;
L100:
    ++j;
    if (j > n) {
	goto L110;
    }
    if (std::abs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
	coltyp[j] = 4;
    } else {

	if (std::abs(d__[j] - d__[jprev]) <= tol) {

            s = z__[jprev];
	    c__ = z__[j];

	    tau = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&c__, &s);
	    c__ /= tau;
	    s = -s / tau;
	    z__[j] = tau;
	    z__[jprev] = 0.;

	    idxjp = idxq[idx[jprev] + 1];
	    idxj = idxq[idx[j] + 1];
	    if (idxjp <= nlp1) {
		--idxjp;
	    }
	    if (idxj <= nlp1) {
		--idxj;
	    }
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], &
		    c__1, &c__, &s);
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, &
		    c__, &s);
	    if (coltyp[j] != coltyp[jprev]) {
		coltyp[j] = 3;
	    }
	    coltyp[jprev] = 4;
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    u2[*k + u2_dim1] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L100;
L110:

    ++(*k);
    u2[*k + u2_dim1] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L120:

    for (j = 1; j <= 4; ++j) {
	ctot[j - 1] = 0;
    }
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	ct = coltyp[j];
	++ctot[ct - 1];
    }

    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	ct = coltyp[jp];
	idxc[psm[ct - 1]] = j;
	++psm[ct - 1];
    }

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	idxj = idxq[idx[idxp[idxc[j]]] + 1];
	if (idxj <= nlp1) {
	    --idxj;
	}
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
    }

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (std::abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    c__ = 1.;
	    s = 0.;
	    z__[1] = tol;
	} else {
	    c__ = z1 / z__[1];
	    s = z__[m] / z__[1];
	}
    } else {
	if (std::abs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);

    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &n, &c__1, &zero, &zero, &u2[u2_offset], ldu2);
    u2[nlp1 + u2_dim1] = 1.;
    if (m > n) {
	i__1 = nlp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
	    vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
	    vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
	}
    } else {
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
    }
    if (m > n) {
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
    }

    if (n > *k) {
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1)
		 * u_dim1 + 1], ldu);
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + 
		vt_dim1], ldvt);
    }
    for (j = 1; j <= 4; ++j) {
	coltyp[j] = ctot[j - 1];
    }

    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd3,DLASD3)(int *nl, 
	int *nr,
	int *sqre, 
	int *k, 
	double *d__, 
	double *q, 
	int *ldq, 
	double *dsigma, 
	double *u, 
	int *ldu, 
	double *u2, 
	int *ldu2, 
	double *vt, 
	int *ldvt, 
	double *vt2, 
	int *ldvt2, 
	int *idxc, 
	int *ctot, 
	double *z__, 
	int *info)
{
    int q_dim1, q_offset, u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, 
	    vt_offset, vt2_dim1, vt2_offset, i__1, i__2;
    double d__2;

    int i__, j, m, n, jc;
    double rho;
    int nlp1, nlp2, nrp1;
    double temp;
    int ctemp;
    int ktemp;
    int c__1 = 1;
    int c__0 = 0;
    double zero = 0.0;
    double one = 1.0;

    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dsigma;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxc;
    --ctot;
    --z__;

    /* Function Body */
    *info = 0;

    if (*nl < 1) {
	*info = -1;
    } else if (*nr < 1) {
	*info = -2;
    } else if (*sqre != 1 && *sqre != 0) {
	*info = -3;
    }

    n = *nl + *nr + 1;
    m = n + *sqre;
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    if (*k == 1) {
	d__[1] = std::abs(z__[1]);
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
	if (z__[1] > 0.) {
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
	} else {
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		u[i__ + u_dim1] = -u2[i__ + u2_dim1];
	    }
	}
	return;
    }

    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &z__[1], &c__1, &q[q_offset], &c__1);

    rho = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(k, &z__[1], &c__1);
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;


    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	PLUMED_BLAS_F77_FUNC(dlasd4,DLASD4)(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j],
		 &vt[j * vt_dim1 + 1], info);

	if (*info != 0) {
	    return;
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
	}
	i__2 = *k - 1;
	for (j = i__; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
	}
	d__2 =  std::sqrt(std::abs(z__[i__]));
	z__[i__] = (q[i__ + q_dim1] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * 
		vt_dim1 + 1];
	u[i__ * u_dim1 + 1] = -1.;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ 
		    * vt_dim1];
	    u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
	}
	temp = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(k, &u[i__ * u_dim1 + 1], &c__1);
	q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
	}
    }

    if (*k == 2) {
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", &n, k, k, &one, &u2[u2_offset], ldu2, &q[q_offset],
		 ldq, &zero, &u[u_offset], ldu);
	goto L100;
    }
    if (ctot[1] > 0) {
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[1], &one, &u2[(u2_dim1 << 1) + 1], 
		ldu2, &q[q_dim1 + 2], ldq, &zero, &u[u_dim1 + 1], ldu);
	if (ctot[3] > 0) {
	    ktemp = ctot[1] + 2 + ctot[2];
	    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1]
		    , ldu2, &q[ktemp + q_dim1], ldq, &one, &u[u_dim1 + 1], 
		    ldu);
	}
    } else if (ctot[3] > 0) {
	ktemp = ctot[1] + 2 + ctot[2];
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1], 
		ldu2, &q[ktemp + q_dim1], ldq, &zero, &u[u_dim1 + 1], ldu);
    } else {
	PLUMED_BLAS_F77_FUNC(dlacpy,DLACPY)("F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu);
    }
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
    ktemp = ctot[1] + 2;
    ctemp = ctot[2] + ctot[3];
    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", nr, k, &ctemp, &one, &u2[nlp2 + ktemp * u2_dim1], ldu2,
	     &q[ktemp + q_dim1], ldq, &zero, &u[nlp2 + u_dim1], ldu);

L100:
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(k, &vt[i__ * vt_dim1 + 1], &c__1);
	q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
	}
    }

    if (*k == 2) {
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", k, &m, k, &one, &q[q_offset], ldq, &vt2[vt2_offset]
		, ldvt2, &zero, &vt[vt_offset], ldvt);
	return;
    }
    ktemp = ctot[1] + 1;
    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", k, &nlp1, &ktemp, &one, &q[q_dim1 + 1], ldq, &vt2[
	    vt2_dim1 + 1], ldvt2, &zero, &vt[vt_dim1 + 1], ldvt);
    ktemp = ctot[1] + 2 + ctot[2];
    if (ktemp <= *ldvt2) {
	PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", k, &nlp1, &ctot[3], &one, &q[ktemp * q_dim1 + 1], 
		ldq, &vt2[ktemp + vt2_dim1], ldvt2, &one, &vt[vt_dim1 + 1], 
		ldvt);
    }

    ktemp = ctot[1] + 1;
    nrp1 = *nr + *sqre;
    if (ktemp > 1) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
	}
    }
    ctemp = ctot[2] + 1 + ctot[3];
    PLUMED_BLAS_F77_FUNC(dgemm,DGEMM)("N", "N", k, &nrp1, &ctemp, &one, &q[ktemp * q_dim1 + 1], ldq, &
	    vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &zero, &vt[nlp2 * vt_dim1 + 
	    1], ldvt);

    return;


}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd4,DLASD4)(int *n, 
	int *i__, 
	double *d__, 
	double *z__, 
	double *delta, 
	double *rho, 
	double *sigma, 
	double *work, 
	int *info)
{
    int i__1;
    double d__1;

    double a, b, c__;
    int j;
    double w, dd[3];
    int ii;
    double dw, zz[3];
    int ip1;
    double eta, phi, eps, tau, psi;
    int iim1, iip1;
    double dphi, dpsi;
    int iter;
    double temp, prew, sg2lb, sg2ub, temp1, temp2, dtiim, delsq, 
	    dtiip;
    int niter;
    double dtisq;
    int swtch;
    double dtnsq;
    double delsq2, dtnsq1;
    int swtch3;
    int orgati;
    double erretm, dtipsq, rhoinv;

    --work;
    --delta;
    --z__;
    --d__;

    *info = 0;
    if (*n == 1) {

	*sigma =  std::sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
	delta[1] = 1.;
	work[1] = 1.;
	return;
    }
    if (*n == 2) {
	PLUMED_BLAS_F77_FUNC(dlasd5,DLASD5)(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
	return;
    }

    eps = PLUMED_GMX_DOUBLE_EPS;
    rhoinv = 1. / *rho;

    if (*i__ == *n) {

	ii = *n - 1;
	niter = 1;

	temp = *rho / 2.;

	temp1 = temp / (d__[*n] +  std::sqrt(d__[*n] * d__[*n] + temp));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] = d__[j] + d__[*n] + temp1;
	    delta[j] = d__[j] - d__[*n] - temp1;
	}

	psi = 0.;
	i__1 = *n - 2;
	for (j = 1; j <= i__1; ++j) {
	    psi += z__[j] * z__[j] / (delta[j] * work[j]);
	}

	c__ = rhoinv + psi;
	w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) + z__[*n] * z__[*
		n] / (delta[*n] * work[*n]);

	if (w <= 0.) {
	    temp1 =  std::sqrt(d__[*n] * d__[*n] + *rho);
	    temp = z__[*n - 1] * z__[*n - 1] / ((d__[*n - 1] + temp1) * (d__[*
		    n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) + z__[*n] * 
		    z__[*n] / *rho;

	    if (c__ <= temp) {
		tau = *rho;
	    } else {
		delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
		a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*
			n];
		b = z__[*n] * z__[*n] * delsq;
		if (a < 0.) {
		    tau = b * 2. / ( std::sqrt(a * a + b * 4. * c__) - a);
		} else {
		    tau = (a +  std::sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
		}
	    }

	} else {
	    delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
	    a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
	    b = z__[*n] * z__[*n] * delsq;

	    if (a < 0.) {
		tau = b * 2. / ( std::sqrt(a * a + b * 4. * c__) - a);
	    } else {
		tau = (a +  std::sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
	    }

	}

	eta = tau / (d__[*n] +  std::sqrt(d__[*n] * d__[*n] + tau));

	*sigma = d__[*n] + eta;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    delta[j] = d__[j] - d__[*i__] - eta;
	    work[j] = d__[j] + d__[*i__] + eta;
	}

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = ii;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (delta[j] * work[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	temp = z__[*n] / (delta[*n] * work[*n]);
	phi = z__[*n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (dpsi 
		+ dphi);

	w = rhoinv + phi + psi;

	if (std::abs(w) <= eps * erretm) {
	    goto L240;
	}

	++niter;
	dtnsq1 = work[*n - 1] * delta[*n - 1];
	dtnsq = work[*n] * delta[*n];
	c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
	a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
	b = dtnsq * dtnsq1 * w;
	if (c__ < 0.) {
	    c__ = std::abs(c__);
	}
	if ( std::abs(c__)<PLUMED_GMX_DOUBLE_MIN) {
	    eta = *rho - *sigma * *sigma;
	} else if (a >= 0.) {
	    eta = (a +  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__  * 2.);
	} else {
	  eta = b * 2. / (a -  std::sqrt(std::abs(a * a - b * 4. * c__)));
	}

	if (w * eta > 0.) {
	    eta = -w / (dpsi + dphi);
	}
	temp = eta - dtnsq;
	if (temp > *rho) {
	    eta = *rho + dtnsq;
	}

	tau += eta;
	eta /= *sigma +  std::sqrt(eta + *sigma * *sigma);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    delta[j] -= eta;
	    work[j] += eta;
	}

	*sigma += eta;

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = ii;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	temp = z__[*n] / (work[*n] * delta[*n]);
	phi = z__[*n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (dpsi 
		+ dphi);

	w = rhoinv + phi + psi;

	iter = niter + 1;

	for (niter = iter; niter <= 20; ++niter) {

	    if (std::abs(w) <= eps * erretm) {
		goto L240;
	    }
	    dtnsq1 = work[*n - 1] * delta[*n - 1];
	    dtnsq = work[*n] * delta[*n];
	    c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
	    a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
	    b = dtnsq1 * dtnsq * w;
	    if (a >= 0.) {
		eta = (a +  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    } else {
	      eta = b * 2. / (a -  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    }

	    if (w * eta > 0.) {
		eta = -w / (dpsi + dphi);
	    }
	    temp = eta - dtnsq;
	    if (temp <= 0.) {
		eta /= 2.;
	    }

	    tau += eta;
	    eta /= *sigma +  std::sqrt(eta + *sigma * *sigma);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		delta[j] -= eta;
		work[j] += eta;
	    }

	    *sigma += eta;

	    dpsi = 0.;
	    psi = 0.;
	    erretm = 0.;
	    i__1 = ii;
	    for (j = 1; j <= i__1; ++j) {
		temp = z__[j] / (work[j] * delta[j]);
		psi += z__[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = std::abs(erretm);

	    temp = z__[*n] / (work[*n] * delta[*n]);
	    phi = z__[*n] * temp;
	    dphi = temp * temp;
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (
		    dpsi + dphi);

	    w = rhoinv + phi + psi;
	}

	*info = 1;
	goto L240;

    } else {

	niter = 1;
	ip1 = *i__ + 1;

	delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
	delsq2 = delsq / 2.;
	temp = delsq2 / (d__[*i__] +  std::sqrt(d__[*i__] * d__[*i__] + delsq2));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] = d__[j] + d__[*i__] + temp;
	    delta[j] = d__[j] - d__[*i__] - temp;
	}

	psi = 0.;
	i__1 = *i__ - 1;
	for (j = 1; j <= i__1; ++j) {
	    psi += z__[j] * z__[j] / (work[j] * delta[j]);
	}

	phi = 0.;
	i__1 = *i__ + 2;
	for (j = *n; j >= i__1; --j) {
	    phi += z__[j] * z__[j] / (work[j] * delta[j]);
	}
	c__ = rhoinv + psi + phi;
	w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) + z__[
		ip1] * z__[ip1] / (work[ip1] * delta[ip1]);

	if (w > 0.) {

	    orgati = 1;
	    sg2lb = 0.;
	    sg2ub = delsq2;
	    a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
	    b = z__[*i__] * z__[*i__] * delsq;
	    if (a > 0.) {
		tau = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    } else {
		tau = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    }
	    eta = tau / (d__[*i__] +  std::sqrt(d__[*i__] * d__[*i__] + tau));
	} else {

	    orgati = 0;
	    sg2lb = -delsq2;
	    sg2ub = 0.;
	    a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
	    b = z__[ip1] * z__[ip1] * delsq;
	    if (a < 0.) {
		tau = b * 2. / (a -  std::sqrt(std::abs(a * a + b * 4. * c__)));
	    } else {
		tau = -(a +  std::sqrt(std::abs(a * a + b * 4. * c__))) /	(c__ * 2.);
	    }
	    eta = tau / (d__[ip1] +  std::sqrt(std::abs(d__[ip1] * d__[ip1] + tau)));
	}

	if (orgati) {
	    ii = *i__;
	    *sigma = d__[*i__] + eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] = d__[j] + d__[*i__] + eta;
		delta[j] = d__[j] - d__[*i__] - eta;
	    }
	} else {
	    ii = *i__ + 1;
	    *sigma = d__[ip1] + eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] = d__[j] + d__[ip1] + eta;
		delta[j] = d__[j] - d__[ip1] - eta;
	    }
	}
	iim1 = ii - 1;
	iip1 = ii + 1;

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = iim1;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	dphi = 0.;
	phi = 0.;
	i__1 = iip1;
	for (j = *n; j >= i__1; --j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    phi += z__[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}

	w = rhoinv + phi + psi;

	swtch3 = 0;
	if (orgati) {
	    if (w < 0.) {
		swtch3 = 1;
	    }
	} else {
	    if (w > 0.) {
		swtch3 = 1;
	    }
	}
	if (ii == 1 || ii == *n) {
	    swtch3 = 0;
	}

	temp = z__[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z__[ii] * temp;
	w += temp;
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. + 
		std::abs(tau) * dw;

	if (std::abs(w) <= eps * erretm) {
	    goto L240;
	}

	if (w <= 0.) {
	    sg2lb = (sg2lb > tau) ? sg2lb : tau;
	} else {
	    sg2ub = (sg2ub < tau) ? sg2ub : tau;
	}

	++niter;
	if (! swtch3) {
	    dtipsq = work[ip1] * delta[ip1];
	    dtisq = work[*i__] * delta[*i__];
	    if (orgati) {
		d__1 = z__[*i__] / dtisq;
		c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
	    } else {
		d__1 = z__[ip1] / dtipsq;
		c__ = w - dtisq * dw - delsq * (d__1 * d__1);
	    }
	    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
	    b = dtipsq * dtisq * w;
	    if ( std::abs(c__)<PLUMED_GMX_DOUBLE_MIN) {
		if ( std::abs(a)<PLUMED_GMX_DOUBLE_MIN) {
		    if (orgati) {
			a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + 
				dphi);
		    } else {
			a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				dphi);
		    }
		}
		eta = b / a;
	    } else if (a <= 0.) {
		eta = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    } else {
		eta = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    }
	} else {

	    dtiim = work[iim1] * delta[iim1];
	    dtiip = work[iip1] * delta[iip1];
	    temp = rhoinv + psi + phi;
	    if (orgati) {
		temp1 = z__[iim1] / dtiim;
		temp1 *= temp1;
		c__ = temp - dtiip * (dpsi + dphi) - (d__[iim1] - d__[iip1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
		zz[0] = z__[iim1] * z__[iim1];
		if (dpsi < temp1) {
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
		}
	    } else {
		temp1 = z__[iip1] / dtiip;
		temp1 *= temp1;
		c__ = temp - dtiim * (dpsi + dphi) - (d__[iip1] - d__[iim1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
		if (dphi < temp1) {
		    zz[0] = dtiim * dtiim * dpsi;
		} else {
		    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
		}
		zz[2] = z__[iip1] * z__[iip1];
	    }
	    zz[1] = z__[ii] * z__[ii];
	    dd[0] = dtiim;
	    dd[1] = delta[ii] * work[ii];
	    dd[2] = dtiip;
	    PLUMED_BLAS_F77_FUNC(dlaed6,DLAED6)(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
	    if (*info != 0) {
		goto L240;
	    }
	}

	if (w * eta >= 0.) {
	    eta = -w / dw;
	}
	if (orgati) {
	    temp1 = work[*i__] * delta[*i__];
	    temp = eta - temp1;
	} else {
	    temp1 = work[ip1] * delta[ip1];
	    temp = eta - temp1;
	}
	if (temp > sg2ub || temp < sg2lb) {
	    if (w < 0.) {
		eta = (sg2ub - tau) / 2.;
	    } else {
		eta = (sg2lb - tau) / 2.;
	    }
	}

	tau += eta;
	eta /= *sigma +  std::sqrt(*sigma * *sigma + eta);

	prew = w;

	*sigma += eta;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] += eta;
	    delta[j] -= eta;
	}

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = iim1;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	dphi = 0.;
	phi = 0.;
	i__1 = iip1;
	for (j = *n; j >= i__1; --j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    phi += z__[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}

	temp = z__[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z__[ii] * temp;
	w = rhoinv + phi + psi + temp;
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. + 
		std::abs(tau) * dw;

	if (w <= 0.) {
	    sg2lb = (sg2lb > tau) ? sg2lb : tau;
	} else {
	    sg2ub = (sg2ub < tau) ? sg2ub : tau;
	}

	swtch = 0;
	if (orgati) {
	    if (-w > std::abs(prew) / 10.) {
		swtch = 1;
	    }
	} else {
	    if (w > std::abs(prew) / 10.) {
		swtch = 1;
	    }
	}

	iter = niter + 1;

	for (niter = iter; niter <= 20; ++niter) {

	    if (std::abs(w) <= eps * erretm) {
		goto L240;
	    }

	    if (! swtch3) {
		dtipsq = work[ip1] * delta[ip1];
		dtisq = work[*i__] * delta[*i__];
		if (! swtch) {
		    if (orgati) {
			d__1 = z__[*i__] / dtisq;
			c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
		    } else {
			d__1 = z__[ip1] / dtipsq;
			c__ = w - dtisq * dw - delsq * (d__1 * d__1);
		    }
		} else {
		    temp = z__[ii] / (work[ii] * delta[ii]);
		    if (orgati) {
			dpsi += temp * temp;
		    } else {
			dphi += temp * temp;
		    }
		    c__ = w - dtisq * dpsi - dtipsq * dphi;
		}
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
		b = dtipsq * dtisq * w;
		if (std::abs(c__)<PLUMED_GMX_DOUBLE_MIN) {
		    if (std::abs(a)<PLUMED_GMX_DOUBLE_MIN) {
			if (! swtch) {
			    if (orgati) {
				a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * 
					(dpsi + dphi);
			    } else {
				a = z__[ip1] * z__[ip1] + dtisq * dtisq * (
					dpsi + dphi);
			    }
			} else {
			    a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
			}
		    }
		    eta = b / a;
		} else if (a <= 0.) {
		  eta = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
		} else {
		  eta = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
		}
	    } else {

		dtiim = work[iim1] * delta[iim1];
		dtiip = work[iip1] * delta[iip1];
		temp = rhoinv + psi + phi;
		if (swtch) {
		    c__ = temp - dtiim * dpsi - dtiip * dphi;
		    zz[0] = dtiim * dtiim * dpsi;
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    if (orgati) {
			temp1 = z__[iim1] / dtiim;
			temp1 *= temp1;
			temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[
				iip1]) * temp1;
			c__ = temp - dtiip * (dpsi + dphi) - temp2;
			zz[0] = z__[iim1] * z__[iim1];
			if (dpsi < temp1) {
			    zz[2] = dtiip * dtiip * dphi;
			} else {
			    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
			}
		    } else {
			temp1 = z__[iip1] / dtiip;
			temp1 *= temp1;
			temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[
				iip1]) * temp1;
			c__ = temp - dtiim * (dpsi + dphi) - temp2;
			if (dphi < temp1) {
			    zz[0] = dtiim * dtiim * dpsi;
			} else {
			    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
			}
			zz[2] = z__[iip1] * z__[iip1];
		    }
		}
		dd[0] = dtiim;
		dd[1] = delta[ii] * work[ii];
		dd[2] = dtiip;
		PLUMED_BLAS_F77_FUNC(dlaed6,DLAED6)(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
		if (*info != 0) {
		    goto L240;
		}
	    }

	    if (w * eta >= 0.) {
		eta = -w / dw;
	    }
	    if (orgati) {
		temp1 = work[*i__] * delta[*i__];
		temp = eta - temp1;
	    } else {
		temp1 = work[ip1] * delta[ip1];
		temp = eta - temp1;
	    }
	    if (temp > sg2ub || temp < sg2lb) {
		if (w < 0.) {
		    eta = (sg2ub - tau) / 2.;
		} else {
		    eta = (sg2lb - tau) / 2.;
		}
	    }

	    tau += eta;
	    eta /= *sigma +  std::sqrt(*sigma * *sigma + eta);

	    *sigma += eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] += eta;
		delta[j] -= eta;
	    }

	    prew = w;

	    dpsi = 0.;
	    psi = 0.;
	    erretm = 0.;
	    i__1 = iim1;
	    for (j = 1; j <= i__1; ++j) {
		temp = z__[j] / (work[j] * delta[j]);
		psi += z__[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = std::abs(erretm);

	    dphi = 0.;
	    phi = 0.;
	    i__1 = iip1;
	    for (j = *n; j >= i__1; --j) {
		temp = z__[j] / (work[j] * delta[j]);
		phi += z__[j] * temp;
		dphi += temp * temp;
		erretm += phi;
	    }

	    temp = z__[ii] / (work[ii] * delta[ii]);
	    dw = dpsi + dphi + temp * temp;
	    temp = z__[ii] * temp;
	    w = rhoinv + phi + psi + temp;
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. 
		    + std::abs(tau) * dw;
	    if (w * prew > 0. && std::abs(w) > std::abs(prew) / 10.) {
		swtch = ! swtch;
	    }

	    if (w <= 0.) {
		sg2lb = (sg2lb > tau) ? sg2lb : tau;
	    } else {
		sg2ub = (sg2ub < tau) ? sg2ub : tau;
	    }
	}

	*info = 1;

    }

L240:
    return;

} 
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd5,DLASD5)(int *i__, 
	double *d__, 
	double *z__, 
	double *delta, 
	double *rho, 
	double *dsigma, 
	double *work)
{
    double b, c__, w, del, tau, delsq;

    --work;
    --delta;
    --z__;
    --d__;

    del = d__[2] - d__[1];
    delsq = del * (d__[2] + d__[1]);
    if (*i__ == 1) {
	w = *rho * 4. * (z__[2] * z__[2] / (d__[1] + d__[2] * 3.) - z__[1] * 
		z__[1] / (d__[1] * 3. + d__[2])) / del + 1.;
	if (w > 0.) {
	    b = delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[1] * z__[1] * delsq;

	    tau = c__ * 2. / (b +  std::sqrt(std::abs(b * b - c__ * 4.)));

	    tau /= d__[1] +  std::sqrt(d__[1] * d__[1] + tau);
	    *dsigma = d__[1] + tau;
	    delta[1] = -tau;
	    delta[2] = del - tau;
	    work[1] = d__[1] * 2. + tau;
	    work[2] = d__[1] + tau + d__[2];
	} else {
	    b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[2] * z__[2] * delsq;

	    if (b > 0.) {
		tau = c__ * -2. / (b +  std::sqrt(b * b + c__ * 4.));
	    } else {
		tau = (b -  std::sqrt(b * b + c__ * 4.)) / 2.;
	    }

	    tau /= d__[2] +  std::sqrt(std::abs(d__[2] * d__[2] + tau));
	    *dsigma = d__[2] + tau;
	    delta[1] = -(del + tau);
	    delta[2] = -tau;
	    work[1] = d__[1] + tau + d__[2];
	    work[2] = d__[2] * 2. + tau;
	}
    } else {

	b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	c__ = *rho * z__[2] * z__[2] * delsq;

	if (b > 0.) {
	    tau = (b +  std::sqrt(b * b + c__ * 4.)) / 2.;
	} else {
	    tau = c__ * 2. / (-b +  std::sqrt(b * b + c__ * 4.));
	}
	tau /= d__[2] +  std::sqrt(d__[2] * d__[2] + tau);
	*dsigma = d__[2] + tau;
	delta[1] = -(del + tau);
	delta[2] = -tau;
	work[1] = d__[1] + tau + d__[2];
	work[2] = d__[2] * 2. + tau;
    }
    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd6,DLASD6)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	double *d__, 
	double *vf, 
	double *vl, 
	double *alpha, 
	double *beta, 
	int *idxq, 
	int *perm, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	double *givnum,
	int *ldgnum, 
	double *poles, 
	double *difl, 
	double *difr, 
	double *z__, 
	int *k, 
	double *c__, 
	double *s, 
	double *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, i__1;
    double d__1, d__2;

    int i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    int isigma;
    double orgnrm;
    int c__0 = 0;
    double one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;

    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;

    d__1 = std::abs(*alpha); 
    d__2 = std::abs(*beta);
    orgnrm = (d__1 > d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = std::abs(d__[i__]);
	if (d__1 > orgnrm)
	    orgnrm = d__1;
    }
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    PLUMED_BLAS_F77_FUNC(dlasd7,DLASD7)(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], &
	    work[ivfw], &vl[1], &work[ivlw], alpha, beta, &work[isigma], &
	    iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[
	    givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, 
	    info);

    PLUMED_BLAS_F77_FUNC(dlasd8,DLASD8)(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], 
	    ldgnum, &work[isigma], &work[iw], info);

    if (*icompq == 1) {
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
	PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }

    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = *k;
    n2 = n - *k;
    PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}


}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd7,DLASD7)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	int *k, 
	double *d__, 
	double *z__, 
	double *zw, 
	double *vf, 
	double *vfw,
	double *vl, 
	double *vlw,
	double *alpha, 
	double *beta,
	double *dsigma, 
	int *idx, 
	int *idxp,
	int *idxq, 
	int *perm, 
	int *givptr,
	int *givcol, 
	int *ldgcol, 
	double *givnum,
	int *ldgnum, 
	double *c__, 
	double *s, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    double d__1, d__2;

    int i__, j, m, n, k2;
    double z1;
    int jp;
    double eps, tau, tol;
    int nlp1, nlp2, idxi, idxj;
    int idxjp;
    int jprev = 0;
    double hlftol;
    int c__1 = 1;

    --d__;
    --z__;
    --zw;
    --vf;
    --vfw;
    --vl;
    --vlw;
    --dsigma;
    --idx;
    --idxp;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*icompq == 1) {
	*givptr = 0;
    }

    z1 = *alpha * vl[nlp1];
    vl[nlp1] = 0.;
    tau = vf[nlp1];
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vl[i__];
	vl[i__] = 0.;
	vf[i__ + 1] = vf[i__];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }
    vf[1] = tau;

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vf[i__];
	vf[i__] = 0.;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	zw[i__] = z__[idxq[i__]];
	vfw[i__] = vf[idxq[i__]];
	vlw[i__] = vl[idxq[i__]];
    }

    PLUMED_BLAS_F77_FUNC(dlamrg,DLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = zw[idxi];
	vf[i__] = vfw[idxi];
	vl[i__] = vlw[idxi];
    }

    eps = PLUMED_GMX_DOUBLE_EPS;

    d__1 = std::abs(*alpha);
    d__2 = std::abs(*beta);
    tol = (d__1>d__2) ? d__1 : d__2;
    d__2 = std::abs(d__[n]);
    tol = eps * 64. * ((d__2>tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (std::abs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    if (j == n) {
		goto L100;
	    }
	} else {
	    jprev = j;
	    goto L70;
	}
    }
L70:
    j = jprev;
L80:
    ++j;
    if (j > n) {
	goto L90;
    }
    if (std::abs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
    } else {

	if (std::abs(d__[j] - d__[jprev]) <= tol) {

	    *s = z__[jprev];
	    *c__ = z__[j];

	    tau = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(c__, s);
	    z__[j] = tau;
	    z__[jprev] = 0.;
	    *c__ /= tau;
	    *s = -(*s) / tau;


	    if (*icompq == 1) {
		++(*givptr);
		idxjp = idxq[idx[jprev] + 1];
		idxj = idxq[idx[j] + 1];
		if (idxjp <= nlp1) {
		    --idxjp;
		}
		if (idxj <= nlp1) {
		    --idxj;
		}
		givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
		givcol[*givptr + givcol_dim1] = idxj;
		givnum[*givptr + (givnum_dim1 << 1)] = *c__;
		givnum[*givptr + givnum_dim1] = *s;
	    }
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
	    PLUMED_BLAS_F77_FUNC(drot,DROT)(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    zw[*k] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L80;
L90:

    ++(*k);
    zw[*k] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L100:

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	vfw[j] = vf[jp];
	vlw[j] = vl[jp];
    }
    if (*icompq == 1) {
	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    jp = idxp[j];
	    perm[j] = idxq[idx[jp] + 1];
	    if (perm[j] <= nlp1) {
		--perm[j];
	    }
	}
    }
    i__1 = n - *k;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (std::abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    *c__ = 1.;
	    *s = 0.;
	    z__[1] = tol;
	} else {
	    *c__ = z1 / z__[1];
	    *s = -z__[m] / z__[1];
	}
	PLUMED_BLAS_F77_FUNC(drot,DROT)(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
	PLUMED_BLAS_F77_FUNC(drot,DROT)(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
    } else {
	if (std::abs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &zw[2], &c__1, &z__[2], &c__1);
    i__1 = n - 1;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
    i__1 = n - 1;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &vlw[2], &c__1, &vl[2], &c__1);

    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasd8,DLASD8)(int *icompq, 
	int *k, 
	double *d__, 
     	double *z__, 
	double *vf, 
	double *vl, 
	double *difl, 
	double *difr, 
	int *lddifr, 
	double *dsigma, 
	double *work, 
	int *info)
{
    int difr_dim1, difr_offset, i__1, i__2;
    double d__2;

    int i__, j;
    double dj, rho;
    int iwk1, iwk2, iwk3;
    double temp;
    int iwk2i, iwk3i;
    double diflj, difrj, dsigj;
    double dsigjp;
    int c__1 = 1;
    int c__0 = 0;
    double one = 1.;

    /* avoid warnings on high gcc optimization levels */
    difrj = dsigjp = 0;

     --d__;
    --z__;
    --vf;
    --vl;
    --difl;
    difr_dim1 = *lddifr;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    --dsigma;
    --work;

    *info = 0;

    if (*k == 1) {
	d__[1] = std::abs(z__[1]);
	difl[1] = d__[1];
	if (*icompq == 1) {
	    difl[2] = 1.;
	    difr[(difr_dim1 << 1) + 1] = 1.;
	}
	return;
    }

    iwk1 = 1;
    iwk2 = iwk1 + *k;
    iwk3 = iwk2 + *k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;

    rho = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(k, &z__[1], &c__1);
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;

    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", k, &c__1, &one, &one, &work[iwk3], k);

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	PLUMED_BLAS_F77_FUNC(dlasd4,DLASD4)(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[
		iwk2], info);

	if (*info != 0) {
	    return;
	}
	work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
	difl[j] = -work[j];
	difr[j + difr_dim1] = -work[j + 1];
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__2 =  std::sqrt(std::abs(work[iwk3i + i__]));
	z__[i__] = (z__[i__] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	diflj = difl[j];
	dj = d__[j];
	dsigj = -dsigma[j];
	if (j < *k) {
	    difrj = -difr[j + difr_dim1];
	    dsigjp = -dsigma[j + 1];
	}
	work[j] = -z__[j] / diflj / (dsigma[j] + dj);
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
        work[i__] = z__[i__] / (dsigma[i__] + dsigj - diflj) / ( dsigma[i__] + dj);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    work[i__] = z__[i__] / (dsigma[i__] + dsigjp - difrj) / (dsigma[i__] + dj);
	}
	temp = PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(k, &work[1], &c__1);
	work[iwk2i + j] = PLUMED_BLAS_F77_FUNC(ddot,DDOT)(k, &work[1], &c__1, &vf[1], &c__1) / temp;
	work[iwk3i + j] = PLUMED_BLAS_F77_FUNC(ddot,DDOT)(k, &work[1], &c__1, &vl[1], &c__1) / temp;
	if (*icompq == 1) {
	    difr[j + (difr_dim1 << 1)] = temp;
	}
    }

    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &work[iwk2], &c__1, &vf[1], &c__1);
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(k, &work[iwk3], &c__1, &vl[1], &c__1);

    return;

} 
}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasda,DLASDA)(int *icompq, 
	int *smlsiz, 
	int *n, 
	int *sqre, 
	double *d__, 
	double *e, 
	double *u, 
	int *ldu, 
	double *vt, 
	int *k, 
	double *difl, 
	double *difr, 
	double *z__, 
	double *poles, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	int *perm, 
	double *givnum, 
	double *c__, 
	double *s, 
	double *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, 
	    difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, 
	    z_dim1, z_offset, i__1, i__2;

    int i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, 
	    nlf, nrf, vfi, iwk, vli, lvl, nru, ndb1, nlp1, lvl2, nrp1;
    double beta;
    int idxq, nlvl;
    double alpha;
    int inode, ndiml, ndimr, idxqi, itemp;
    int sqrei;
    int nwork1, nwork2;
    int smlszp;
    int c__0 = 0;
    double zero = 0.0;
    double one = 1.;
    int c__1 = 1;
    --d__;
    --e;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --work;
    --iwork;
    *info = 0;

    m = *n + *sqre;

    if (*n <= *smlsiz) {
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		    vt_offset], ldu, &u[u_offset], ldu, &u[u_offset], ldu, &
		    work[1], info);
	} else {
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldu, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], 
		    info);
	}
	return;
    }

    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;

    ncc = 0;
    nru = 0;

    smlszp = *smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;

    PLUMED_BLAS_F77_FUNC(dlasdt,DLASDT)(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
	i1 = i__ - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nlf = ic - nl;
	nrf = ic + 1;
	idxqi = idxq + nlf - 2;
	vfi = vf + nlf - 1;
	vli = vl + nlf - 1;
	sqrei = 1;
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nlp1, &nlp1, &zero, &one, &work[nwork1], &smlszp);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &
		    work[nwork1], &smlszp, &work[nwork2], &nl, &work[nwork2], 
		    &nl, &work[nwork2], info);
	    itemp = nwork1 + nl * smlszp;
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nl, &nl, &zero, &one, &u[nlf + u_dim1], ldu);
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nlp1, &nlp1, &zero, &one, &vt[nlf + vt_dim1], 
		    ldu);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &
		    vt[nlf + vt_dim1], ldu, &u[nlf + u_dim1], ldu, &u[nlf + 
		    u_dim1], ldu, &work[nwork1], info);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
	if (i__ == nd && *sqre == 0) {
	    sqrei = 0;
	} else {
	    sqrei = 1;
	}
	idxqi += nlp1;
	vfi += nlp1;
	vli += nlp1;
	nrp1 = nr + sqrei;
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nrp1, &nrp1, &zero, &one, &work[nwork1], &smlszp);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &
		    work[nwork1], &smlszp, &work[nwork2], &nr, &work[nwork2], 
		    &nr, &work[nwork2], info);
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nr, &nr, &zero, &one, &u[nrf + u_dim1], ldu);
	    PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("A", &nrp1, &nrp1, &zero, &one, &vt[nrf + vt_dim1], 
		    ldu);
	    PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &
		    vt[nrf + vt_dim1], ldu, &u[nrf + u_dim1], ldu, &u[nrf + 
		    u_dim1], ldu, &work[nwork1], info);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nr;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
    }

    j = (1 << nlvl);

    for (lvl = nlvl; lvl >= 1; --lvl) {
	lvl2 = (lvl << 1) - 1;

	if (lvl == 1) {
	    lf = 1;
	    ll = 1;
	} else {
	    lf = (1 << (lvl-1));
	    ll = (lf << 1) - 1;
	}
	i__1 = ll;
	for (i__ = lf; i__ <= i__1; ++i__) {
	    im1 = i__ - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    if (i__ == ll) {
		sqrei = *sqre;
	    } else {
		sqrei = 1;
	    }
	    vfi = vf + nlf - 1;
	    vli = vl + nlf - 1;
	    idxqi = idxq + nlf - 1;
	    alpha = d__[ic];
	    beta = e[ic];
	    if (*icompq == 0) {
		PLUMED_BLAS_F77_FUNC(dlasd6,DLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[
			perm_offset], &givptr[1], &givcol[givcol_offset], 
			ldgcol, &givnum[givnum_offset], ldu, &poles[
			poles_offset], &difl[difl_offset], &difr[difr_offset],
			 &z__[z_offset], &k[1], &c__[1], &s[1], &work[nwork1],
			 &iwork[iwk], info);
	    } else {
		--j;
		PLUMED_BLAS_F77_FUNC(dlasd6,DLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[nlf + 
			lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * 
			givcol_dim1], ldgcol, &givnum[nlf + lvl2 * 
			givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &
			difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * 
			difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[j], 
			&s[j], &work[nwork1], &iwork[iwk], info);
	    }
	    if (*info != 0) {
		return;
	    }
	}
    }

    return;

}


}
}
#include <cctype>

#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasdq,DLASDQ)(const char *uplo,
                        int *sqre,
                        int *n,
                        int *ncvt,
                        int *nru,
                        int *ncc,
                        double *d__,
                        double *e, 
                        double *vt, 
                        int *ldvt, 
                        double *u,
                        int *ldu, 
                        double *c__,
                        int *ldc,
                        double *work, 
                        int *info)
{
    const char xuplo=std::toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    int c__1 = 1;
    int itmp1,itmp2;
    int i__, j;
    double r__, cs, sn;
    int np1, isub;
    double smin;
    int sqre1;
    int iuplo;
    int rotate;

    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    iuplo = 0;
    if (xuplo == 'U') {
	iuplo = 1;
    }
    if (xuplo == 'L') {
	iuplo = 2;
    }
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    if (iuplo == 0) {
	*info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ncvt < 0) {
	*info = -4;
    } else if (*nru < 0) {
	*info = -5;
    } else if (*ncc < 0) {
	*info = -6;
    } else if ((*ncvt == 0 && *ldvt < 1) || (*ncvt > 0 && *ldvt < itmp1)) {
	*info = -10;
    } else if (*ldu < itmp2) {
	*info = -12;
    } else if ((*ncc == 0 && *ldc < 1) || (*ncc > 0 && *ldc < itmp1)) {
	*info = -14;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    np1 = *n + 1;
    sqre1 = *sqre;

    if (iuplo == 1 && sqre1 == 1) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}
	PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	d__[*n] = r__;
	e[*n] = 0.f;
	if (rotate) {
	    work[*n] = cs;
	    work[*n + *n] = sn;
	}
	iuplo = 2;
	sqre1 = 0;

	if (*ncvt > 0) {
	    PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &np1, ncvt, &work[1], &work[np1], &vt[
		    vt_offset], ldvt);
	}
    }
    if (iuplo == 2) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}

	if (sqre1 == 1) {
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	    d__[*n] = r__;
	    if (rotate) {
		work[*n] = cs;
		work[*n + *n] = sn;
	    }
	}
	if (*nru > 0) {
	    if (sqre1 == 0) {
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, n, &work[1], &work[np1], &u[
			u_offset], ldu);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &np1, &work[1], &work[np1], &u[
			u_offset], ldu);
	    }
	}
	if (*ncc > 0) {
	    if (sqre1 == 0) {
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", n, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("L", "V", "F", &np1, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    }
	}
    }

    PLUMED_BLAS_F77_FUNC(dbdsqr,DBDSQR)("U", n, ncvt, nru, ncc, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[
	    u_offset], ldu, &c__[c_offset], ldc, &work[1], info);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = i__;
	smin = d__[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (d__[j] < smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != i__) {
	    d__[isub] = d__[i__];
	    d__[i__] = smin;
	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[i__ + vt_dim1], 
			ldvt);
	    }
	    if (*nru > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[i__ * u_dim1 + 1]
			, &c__1);
	    }
	    if (*ncc > 0) {
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[i__ + c_dim1], ldc)
			;
	    }
	}
    }

    return;
}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlasdt,DLASDT)(int *n,
	int *lvl,
	int *nd,
	int *inode,
	int *ndiml,
	int *ndimr,
	int *msub)
{
  int maxn = (*n > 1) ? *n : 1;
  double temp;
  int i,il,ir,llst,nlvl,ncrnt;

  temp = std::log( ((double) maxn) / ((double)(*msub+1))) / std::log(2.0);
  
  *lvl = 1 + (int) temp;

  i = *n / 2;
  inode[0] = i + 1;
  ndiml[0] = i;
  ndimr[0] = *n - i - 1;
  il = -1;
  ir = 0;
  llst = 1;

  for(nlvl=1;nlvl<*lvl;nlvl++) {
    for(i=0;i<llst;i++) {
      il += 2;
      ir += 2;
      ncrnt = llst + i - 1;
      ndiml[il] = ndiml[ncrnt] / 2;
      ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
      inode[il] = inode[ncrnt] - ndimr[il] - 1;
      ndiml[ir] = ndimr[ncrnt] / 2;
      ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
      inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
    }
    llst *= 2;
  }
  *nd = llst*2 - 1;
  return;
}
}
}
#include <cctype>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaset,DLASET)(const char *uplo,
	int *m,
	int *n,
	double *alpha,
	double *beta,
	double *a,
	int *lda)
{
  int i,j,k;
  const char ch=std::toupper(*uplo);

  if(ch=='U') {
    for(j=1;j<*n;j++) {
      k = (j < *m) ? j : *m;
      for(i=0;i<k;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else if(ch=='L') {
    k = (*m < *n) ? *m : *n;
    for(j=0;j<k;j++) {
      for(i=j+1;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }    
  }

  k = (*m < *n) ? *m : *n;
  for(i=0;i<k;i++)
    a[i*(*lda)+i] = *beta;
}
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlasq1,DLASQ1)(int *n,
	double *d,
	double *e,
	double *work,
	int *info)
{
  double sigmx = 0.0;
  int i,j,k,iinfo;
  double minval,safemin;
  double dtemp,scale;
  double eps;

  eps = PLUMED_GMX_DOUBLE_EPS;
  minval = PLUMED_GMX_DOUBLE_MIN;
  safemin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS);
  *info = 0;

  if(*n<0) {
    *info = -2;
    return;
  }
  
  for(i=0;i<*n-1;i++) {
    d[i] = std::abs(d[i]);
    dtemp = std::abs(e[i]);
    if(dtemp>sigmx)
      sigmx=dtemp;
  }
  d[*n-1] = std::abs(d[*n-1]);
  
  if(std::abs(sigmx)<PLUMED_GMX_DOUBLE_MIN) {
    PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("D",n,d,&iinfo);
    return;
  }

  for(i=0;i<*n;i++) {
    if(d[i]>sigmx)
      sigmx=d[i];
  }

  /* Copy d and e into work (z format) and scale.
   * Squaring input data makes scaling by a power of the
   * radix pointless.
   */
  scale =  std::sqrt(eps/safemin);
  i = 1;
  j = 2;
  PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n,d,&i,work,&j);
  k = *n-1;
  PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&k,e,&i,work+1,&j);
  i = 0;
  j = 2*(*n)-1;
  k = 1;
  PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G",&i,&i,&sigmx,&scale,&j,&k,work,&j,&iinfo);


  /* Compute q and e elements */
  for(i=0;i<2*(*n)-1;i++)
    work[i] = work[i]*work[i];

  work[2*(*n)-1] = 0.0;

  PLUMED_BLAS_F77_FUNC(dlasq2,DLASQ2)(n,work,info);

  j = 0;
  k = 1;
  if(*info==0) {
    for(i=0;i<*n;i++)
      d[i]= std::sqrt(work[i]);
    PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G",&j,&j,&scale,&sigmx,n,&k,d,n,&iinfo);
  }
  return;
}
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#ifdef _MSC_VER
#pragma warning(disable: 4723) /*division by zero - is used on purpose here*/
#endif

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasq2,DLASQ2)(int *n, 
                        double *z__, 
                        int *info)
{
    int i__1, i__2, i__3;
    double d__1, d__2;

    double d__, e;
    int k;
    double s, t;
    int i0, i4, n0, pp;
    double dee, eps, tol;
    int ipn4;
    double tol2;
    int ieee;
    int nbig;
    double dmin__, emin, emax;
    int kmin, ndiv, iter;
    double qmin, temp, qmax, zmax;
    int splt, nfail;
    double desig, trace, sigma;
    int iinfo;
    double deemin;
    int iwhila, iwhilb;
    double oldemn, safmin, minval;
    double posinf,neginf,negzro,newzro;
    double zero = 0.0;
    double one = 1.0;

    --z__;

    *info = 0;
    eps = PLUMED_GMX_DOUBLE_EPS;
    minval = PLUMED_GMX_DOUBLE_MIN;
    safmin = minval*(1.0+eps);

    tol = eps * 100.;

    d__1 = tol;
    tol2 = d__1 * d__1;

    if (*n < 0) {
	*info = -1;
	return;
    } else if (*n == 0) {
	return;
    } else if (*n == 1) {

	if (z__[1] < 0.) {
	    *info = -201;
	}
	return;
    } else if (*n == 2) {

	if (z__[2] < 0. || z__[3] < 0.) {
	    *info = -2;
	    return;
	} else if (z__[3] > z__[1]) {
	    d__ = z__[3];
	    z__[3] = z__[1];
	    z__[1] = d__;
	}
	z__[5] = z__[1] + z__[2] + z__[3];
	if (z__[2] > z__[3] * tol2) {
	    t = (z__[1] - z__[3] + z__[2]) * .5;
	    s = z__[3] * (z__[2] / t);
	    if (s <= t) {
		s = z__[3] * (z__[2] / (t * ( std::sqrt(s / t + 1.) + 1.)));
	    } else {
		s = z__[3] * (z__[2] / (t +  std::sqrt(t) * std::sqrt(t + s)));
	    }
	    t = z__[1] + (s + z__[2]);
	    z__[3] *= z__[1] / t;
	    z__[1] = t;
	}
	z__[2] = z__[3];
	z__[6] = z__[2] + z__[1];
	return;
    }
    z__[*n * 2] = 0.;
    emin = z__[2];
    qmax = 0.;
    zmax = 0.;
    d__ = 0.;
    e = 0.;

    i__1 = 2*(*n - 1);
    for (k = 1; k <= i__1; k += 2) {
	if (z__[k] < 0.) {
	    *info = -(k + 200);
	    return;
	} else if (z__[k + 1] < 0.) {
	    *info = -(k + 201);
	    return;
	}
	d__ += z__[k];
	e += z__[k + 1];
	d__1 = qmax, d__2 = z__[k];
	qmax = (d__1>d__2) ? d__1 : d__2;
	d__1 = emin, d__2 = z__[k + 1];
	emin = (d__1<d__2) ? d__1 : d__2;
	d__1 = (qmax>zmax) ? qmax : zmax;
	d__2 = z__[k + 1];
	zmax = (d__1>d__2) ? d__1 : d__2;
    }
    if (z__[(*n << 1) - 1] < 0.) {
	*info = -((*n << 1) + 199);
	return;
    }
    d__ += z__[(*n << 1) - 1];
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
    qmax = (d__1>d__2) ? d__1 : d__2;

    if (std::abs(e)<PLUMED_GMX_DOUBLE_MIN) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    z__[k] = z__[(k << 1) - 1];
	}
	PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("D", n, &z__[1], &iinfo);
	z__[(*n << 1) - 1] = d__;
	return;
    }

    trace = d__ + e;

    if (std::abs(trace)<PLUMED_GMX_DOUBLE_MIN) {
	z__[(*n << 1) - 1] = 0.;
	return;
    }

    ieee = 1;
    posinf = one/zero;
    if(posinf<=1.0)
      ieee = 0;
    neginf = -one/zero;
    if(neginf>=0.0)
      ieee = 0;
    negzro = one/(neginf+one);
    if(std::abs(negzro)>PLUMED_GMX_DOUBLE_MIN)
      ieee = 0;
    neginf = one/negzro;
    if(neginf>=0)
      ieee = 0;
    newzro = negzro + zero;
    if(std::abs(newzro-zero)>PLUMED_GMX_DOUBLE_MIN)
      ieee = 0;
    posinf = one /newzro;
    if(posinf<=one)
      ieee = 0;
    neginf = neginf*posinf;
    if(neginf>=zero)
      ieee = 0;
    posinf = posinf*posinf;
    if(posinf<=1.0)
      ieee = 0;

    for (k = *n << 1; k >= 2; k += -2) {
	z__[k * 2] = 0.;
	z__[(k << 1) - 1] = z__[k];
	z__[(k << 1) - 2] = 0.;
	z__[(k << 1) - 3] = z__[k - 1];
    }

    i0 = 1;
    n0 = *n;

    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
	ipn4 = 4*(i0 + n0);
	i__1 = 2*(i0 + n0 - 1);
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
	    temp = z__[i4 - 3];
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
	    z__[ipn4 - i4 - 3] = temp;
	    temp = z__[i4 - 1];
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
	    z__[ipn4 - i4 - 5] = temp;
	}
    }

    pp = 0;

    for (k = 1; k <= 2; ++k) {

	d__ = z__[(n0 << 2) + pp - 3];
	i__1 = (i0 << 2) + pp;
	for (i4 = 4*(n0 - 1) + pp; i4 >= i__1; i4 += -4) {
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		d__ = z__[i4 - 3];
	    } else {
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
	    }
	}

	emin = z__[(i0 << 2) + pp + 1];
	d__ = z__[(i0 << 2) + pp - 3];
	i__1 = 4*(n0 - 1) + pp;
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		z__[i4 - (pp << 1) - 2] = d__;
		z__[i4 - (pp << 1)] = 0.;
		d__ = z__[i4 + 1];
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
	    }
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	z__[(n0 << 2) - pp - 2] = d__;


	qmax = z__[(i0 << 2) - pp - 2];
	i__1 = (n0 << 2) - pp - 2;
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
	    d__1 = qmax, d__2 = z__[i4];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	}

	pp = 1 - pp;
    }

    iter = 2;
    nfail = 0;
    ndiv = 2*(n0 - i0);

    i__1 = *n + 1;
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
	if (n0 < 1) {
	    goto L170;
	}

	desig = 0.;
	if (n0 == *n) {
	    sigma = 0.;
	} else {
	    sigma = -z__[(n0 << 2) - 1];
	}
	if (sigma < 0.) {
	    *info = 1;
	    return;
	}

	emax = 0.;
	if (n0 > i0) {
	    emin = std::abs(z__[(n0 << 2) - 5]);
	} else {
	    emin = 0.;
	}
	qmin = z__[(n0 << 2) - 3];
	qmax = qmin;
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
	    if (z__[i4 - 5] <= 0.) {
		goto L100;
	    }
	    if (qmin >= emax * 4.) {
		d__1 = qmin, d__2 = z__[i4 - 3];
		qmin = (d__1<d__2) ? d__1 : d__2;
		d__1 = emax, d__2 = z__[i4 - 5];
		emax = (d__1>d__2) ? d__1 : d__2;
	    }
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	    d__1 = emin, d__2 = z__[i4 - 5];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	i4 = 4;

L100:
	i0 = i4 / 4;
	pp = 0;

	if (n0 - i0 > 1) {
	    dee = z__[(i0 << 2) - 3];
	    deemin = dee;
	    kmin = i0;
	    i__2 = (n0 << 2) - 3;
	    for (i4 = (i0 << 2) - 3; i4 <= i__2; i4 += 4) {
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
		if (dee <= deemin) {
		    deemin = dee;
		    kmin = (i4 + 3) / 4;
		}
	    }
	    if (2*(kmin - i0) < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
		ipn4 = 4*(i0 + n0);
		pp = 2;
		i__2 = 2*(i0 + n0 - 1);
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
		    temp = z__[i4 - 3];
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
		    z__[ipn4 - i4 - 3] = temp;
		    temp = z__[i4 - 2];
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
		    z__[ipn4 - i4 - 2] = temp;
		    temp = z__[i4 - 1];
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
		    z__[ipn4 - i4 - 5] = temp;
		    temp = z__[i4];
		    z__[i4] = z__[ipn4 - i4 - 4];
		    z__[ipn4 - i4 - 4] = temp;
		}
	    }
	}


	d__1 = 0., d__2 = qmin -  std::sqrt(qmin) * 2. * std::sqrt(emax);
	dmin__ = -((d__1>d__2) ? d__1 : d__2);

	nbig = (n0 - i0 + 1) * 30;
	i__2 = nbig;
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
	    if (i0 > n0) {
		goto L150;
	    }

	    PLUMED_BLAS_F77_FUNC(dlasq3,DLASQ3)(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee);

	    pp = 1 - pp;

	    if (pp == 0 && n0 - i0 >= 3) {
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
		    splt = i0 - 1;
		    qmax = z__[(i0 << 2) - 3];
		    emin = z__[(i0 << 2) - 1];
		    oldemn = z__[i0 * 4];
		    i__3 = 4*(n0 - 3);
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
			    z__[i4 - 1] = -sigma;
			    splt = i4 / 4;
			    qmax = 0.;
			    emin = z__[i4 + 3];
			    oldemn = z__[i4 + 4];
			} else {
			    d__1 = qmax, d__2 = z__[i4 + 1];
			    qmax = (d__1>d__2) ? d__1 : d__2;
			    d__1 = emin, d__2 = z__[i4 - 1];
			    emin = (d__1<d__2) ? d__1 : d__2;
			    d__1 = oldemn, d__2 = z__[i4];
			    oldemn = (d__1<d__2) ? d__1 : d__2;
			}
		    }
		    z__[(n0 << 2) - 1] = emin;
		    z__[n0 * 4] = oldemn;
		    i0 = splt + 1;
		}
	    }
	}

	*info = 2;
	return;

L150:
	;
    }

    *info = 3;
    return;


L170:

    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	z__[k] = z__[(k << 2) - 3];
    }

    PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("D", n, &z__[1], &iinfo);

    e = 0.;
    for (k = *n; k >= 1; --k) {
	e += z__[k];
    }


    z__[(*n << 1) + 1] = trace;
    z__[(*n << 1) + 2] = e;
    z__[(*n << 1) + 3] = (double) iter;
    i__1 = *n;
    z__[(*n << 1) + 4] = (double) ndiv / (double) (i__1 * i__1);
    z__[(*n << 1) + 5] = nfail * 100. / (double) iter;

    return;

}



}
}
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlasq3,DLASQ3)(int *i0, 
                        int *n0, 
                        double *z__, 
                        int *pp, 
                        double *dmin__, 
                        double *sigma,
                        double *desig,
                        double *qmax, 
                        int *nfail, 
                        int *iter, 
                        int *ndiv, 
	int *ieee)
{

    int ttype = 0;
    double dmin1 = 0.;
    double dmin2 = 0.;
    double dn = 0.;
    double dn1 = 0.;
    double dn2 = 0.;
    double tau = 0.;

    int i__1;
    double d__1, d__2;
    double s, t;
    int j4, nn;
    double eps, tol;
    int n0in, ipn4;
    double tol2, temp;
    --z__;

    n0in = *n0;
    eps = PLUMED_GMX_DOUBLE_EPS;
    tol = eps * 100.;
    d__1 = tol;
    tol2 = d__1 * d__1;


L10:

    if (*n0 < *i0) {
	return;
    }
    if (*n0 == *i0) {
	goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0 == *i0 + 1) {
	goto L40;
    }

    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 
	    4] > tol2 * z__[nn - 7]) {
	goto L30;
    }

L20:

    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;

L30:

    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[
	    nn - 11]) {
	goto L50;
    }

L40:

    if (z__[nn - 3] > z__[nn - 7]) {
	s = z__[nn - 3];
	z__[nn - 3] = z__[nn - 7];
	z__[nn - 7] = s;
    }
    if (z__[nn - 5] > z__[nn - 3] * tol2) {
	t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
	s = z__[nn - 3] * (z__[nn - 5] / t);
	if (s <= t) {
	    s = z__[nn - 3] * (z__[nn - 5] / (t * ( std::sqrt(s / t + 1.) + 1.)));
	} else {
	    s = z__[nn - 3] * (z__[nn - 5] / (t +  std::sqrt(t) * std::sqrt(t + s)));
	}
	t = z__[nn - 7] + (s + z__[nn - 5]);
	z__[nn - 3] *= z__[nn - 7] / t;
	z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;

L50:
    if (*pp == 2) {
	*pp = 0;
    }

    if (*dmin__ <= 0. || *n0 < n0in) {
	if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
	    ipn4 = 4*(*i0 + *n0);
	    i__1 = 2*(*i0 + *n0 - 1);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		temp = z__[j4 - 3];
		z__[j4 - 3] = z__[ipn4 - j4 - 3];
		z__[ipn4 - j4 - 3] = temp;
		temp = z__[j4 - 2];
		z__[j4 - 2] = z__[ipn4 - j4 - 2];
		z__[ipn4 - j4 - 2] = temp;
		temp = z__[j4 - 1];
		z__[j4 - 1] = z__[ipn4 - j4 - 5];
		z__[ipn4 - j4 - 5] = temp;
		temp = z__[j4];
		z__[j4] = z__[ipn4 - j4 - 4];
		z__[ipn4 - j4 - 4] = temp;
	    }
	    if (*n0 - *i0 <= 4) {
		z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
		z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
	    }
	    d__1 = dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
	    dmin2 = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1]
		    , d__1 = ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
	    z__[(*n0 << 2) + *pp - 1] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 =
		     ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) - *pp + 4];
	    z__[(*n0 << 2) - *pp] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = *qmax;
	    d__2 = z__[(*i0 << 2) + *pp - 3];
	    d__1 = (d__1>d__2) ? d__1 : d__2;
	    d__2 = z__[(*i0 << 2) + *pp + 1];
	    *qmax = ((d__1>d__2) ? d__1 : d__2);
	    *dmin__ = -0.;
	}
    }


    PLUMED_BLAS_F77_FUNC(dlasq4,DLASQ4)(i0, n0, &z__[1], pp, &n0in, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, &tau, &ttype);

L70:

    PLUMED_BLAS_F77_FUNC(dlasq5,DLASQ5)(i0, n0, &z__[1], pp, &tau, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, ieee);

    *ndiv += *n0 - *i0 + 2;
    ++(*iter);

    if (*dmin__ >= 0. && dmin1 > 0.) {

	goto L90;

    } else if (*dmin__ < 0. && dmin1 > 0. && z__[4*(*n0 - 1) - *pp] < tol *
	     (*sigma + dn1) && std::abs(dn) < tol * *sigma) {

	z__[4*(*n0 - 1) - *pp + 2] = 0.;
	*dmin__ = 0.;
	goto L90;
    } else if (*dmin__ < 0.) {

	++(*nfail);
	if (ttype < -22) {

	    tau = 0.;
	} else if (dmin1 > 0.) {

	    tau = (tau + *dmin__) * (1. - eps * 2.);
	    ttype += -11;
	} else {

	    tau *= .25;
	    ttype += -12;
	}
	goto L70;
    }
    else {
        
        goto L80;
    }

L80:
    PLUMED_BLAS_F77_FUNC(dlasq6,DLASQ6)(i0, n0, &z__[1], pp, dmin__, &dmin1, &dmin2, &dn, &dn1, &dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    tau = 0.;

L90:
    if (tau < *sigma) {
	*desig += tau;
	t = *sigma + *desig;
	*desig -= t - *sigma;
    } else {
	t = *sigma + tau;
	*desig = *sigma - (t - tau) + *desig;
    }
    *sigma = t;

    return;
}
}
}
#include <cmath>
#include "real.h"

#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasq4,DLASQ4)(int *i0, 
	int *n0, 
	double *z__, 
	int *pp, 
	int *n0in, 
	double *dmin__, 
	double *dmin1, 
	double *dmin2, 
	double *dn, 
	double *dn1, 
	double *dn2, 
	double *tau, 
	int *ttype)
{
    double g = 0.;
    int i__1;
    double d__1, d__2;

    double s, a2, b1, b2;
    int i4, nn, np;
    double gam, gap1, gap2;


    if (*dmin__ <= 0.) {
	*tau = -(*dmin__);
	*ttype = -1;
	return;
    }

    s = 0.0;

    nn = (*n0 << 2) + *pp;
    if (*n0in == *n0) {

	if ( std::abs(*dmin__ - *dn)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin__ + *dn) ||
         std::abs(*dmin__ - *dn1)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin__ + *dn1)) {

	    b1 =  std::sqrt(z__[nn - 3]) * std::sqrt(z__[nn - 5]);
	    b2 =  std::sqrt(z__[nn - 7]) * std::sqrt(z__[nn - 9]);
	    a2 = z__[nn - 7] + z__[nn - 5];

        if ( std::abs(*dmin__ - *dn)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin__ + *dn) &&
             std::abs(*dmin1 - *dn1)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin1 + *dn1)) {

            gap2 = *dmin2 - a2 - *dmin2 * .25;
		if (gap2 > 0. && gap2 > b2) {
		    gap1 = a2 - *dn - b2 / gap2 * b2;
		} else {
		    gap1 = a2 - *dn - (b1 + b2);
		}
		if (gap1 > 0. && gap1 > b1) {
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -2;
		} else {
		    s = 0.;
		    if (*dn > b1) {
			s = *dn - b1;
		    }
		    if (a2 > b1 + b2) {
			d__1 = s, d__2 = a2 - (b1 + b2);
			s = (d__1<d__2) ? d__1 : d__2;
		    }
		    d__1 = s, d__2 = *dmin__ * .333;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -3;
		}
	    } else {


		*ttype = -4;
		s = *dmin__ * .25;
		if (std::abs(*dmin__ - *dn)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin__ + *dn)) {
		    gam = *dn;
		    a2 = 0.;
		    if (z__[nn - 5] > z__[nn - 7]) {
			return;
		    }
		    b2 = z__[nn - 5] / z__[nn - 7];
		    np = nn - 9;
		} else {
		    np = nn - (*pp << 1);
		    gam = *dn1;
		    if (z__[np - 4] > z__[np - 2]) {
			return;
		    }
		    a2 = z__[np - 4] / z__[np - 2];
		    if (z__[nn - 9] > z__[nn - 11]) {
			return;
		    }
		    b2 = z__[nn - 9] / z__[nn - 11];
		    np = nn - 13;
		}


		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = np; i4 >= i__1; i4 += -4) {
		    if (std::abs(b2)<PLUMED_GMX_DOUBLE_MIN) {
			goto L20;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L20;
		    }
		}
L20:
		a2 *= 1.05;


		if (a2 < .563) {
		    s = gam * (1. -  std::sqrt(a2)) / (a2 + 1.);
		}
	    }
	} else if (std::abs(*dmin__ - *dn2)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin__ + *dn2)) {

	    *ttype = -5;
	    s = *dmin__ * .25;

	    np = nn - (*pp << 1);
	    b1 = z__[np - 2];
	    b2 = z__[np - 6];
	    gam = *dn2;
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
		return;
	    }
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);


	    if (*n0 - *i0 > 2) {
		b2 = z__[nn - 13] / z__[nn - 15];
		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
		    if (std::abs(b2)<PLUMED_GMX_DOUBLE_MIN) {
			goto L40;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L40;
		    }
		}
L40:
		a2 *= 1.05;
	    }

	    if (a2 < .563) {
		s = gam * (1. -  std::sqrt(a2)) / (a2 + 1.);
	    }
	} else {

	    if (*ttype == -6) {
		g += (1. - g) * .333;
	    } else if (*ttype == -18) {
		g = .083250000000000005;
	    } else {
		g = .25;
	    }
	    s = g * *dmin__;
	    *ttype = -6;
	}

    } else if (*n0in == *n0 + 1) {

        if ( std::abs(*dmin1 - *dn1)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin1 + *dn1) &&
             std::abs(*dmin2 - *dn2)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin2 + *dn2)) {

	    *ttype = -7;
	    s = *dmin1 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (std::abs(b2)<PLUMED_GMX_DOUBLE_MIN) {
		goto L60;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		a2 = b1;
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (((a2>b1) ? a2 : b1) * 100. < b2) {
		    goto L60;
		}
	    }
L60:
	    b2 =  std::sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
	    gap2 = *dmin2 * .5 - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
		*ttype = -8;
	    }
	} else {

	    s = *dmin1 * .25;
	    if (std::abs(*dmin1 - *dn1)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin1 + *dn1)) {
		s = *dmin1 * .5;
	    }
	    *ttype = -9;
	}

    } else if (*n0in == *n0 + 2) {

	if (std::abs(*dmin2 - *dn2)<PLUMED_GMX_DOUBLE_EPS*std::abs(*dmin2 + *dn2) &&
        z__[nn - 5] * 2. < z__[nn - 7]) {
	    *ttype = -10;
	    s = *dmin2 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (std::abs(b2)<PLUMED_GMX_DOUBLE_MIN) {
		goto L80;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (b1 * 100. < b2) {
		    goto L80;
		}
	    }
L80:
	    b2 =  std::sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
	    gap2 = z__[nn - 7] + z__[nn - 9] -  std::sqrt(z__[nn - 11]) * std::sqrt(z__[
		    nn - 9]) - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    s = *dmin2 * .25;
	    *ttype = -11;
	}
    } else if (*n0in > *n0 + 2) {

	s = 0.;
	*ttype = -12;
    }

    *tau = s;
    return;

}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlasq5,DLASQ5)(int *i0, 
	int *n0,
	double *z__, 
	int *pp, 
	double *tau,
	double *dmin__, 
	double *dmin1, 
	double *dmin2, 
	double *dn,
	double *dnm1, 
	double *dnm2,
	int *ieee)
{
    int i__1;
    double d__1, d__2;

    double d__;
    int    j4, j4p2;
    double emin, temp;

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4] - *tau;
    *dmin__ = d__;
    *dmin1 = -z__[j4];

    if (*ieee) {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		temp = z__[j4 + 1] / z__[j4 - 2];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4] = z__[j4 - 1] * temp;
		d__1 = z__[j4];
                if(d__1<emin)
                  emin = d__1;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		temp = z__[j4 + 2] / z__[j4 - 3];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4 - 1] = z__[j4] * temp;
		d__1 = z__[j4 - 1];
                if(d__1<emin)
                  emin = d__1;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
        if(*dn<*dmin__)
          *dmin__ = *dn;

    } else {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		    d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		    d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4 - 1];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	if (*dnm2 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	}
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	if (*dnm1 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	}
        if(*dn<*dmin__)
          *dmin__ = *dn;

    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;

}

}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasq6,DLASQ6)(int *i0, 
	int *n0, 
	double *z__, 
	int *pp, 
	double *dmin__, 
	double *dmin1, 
	double *dmin2,
	double *dn, 
	double *dnm1, 
	double *dnm2)
{
    int i__1;
    double d__1, d__2;

    /* Local variables */
    double d__;
    int j4, j4p2;
    double emin, temp;
    const double safemin = PLUMED_GMX_DOUBLE_MIN*(1.0+PLUMED_GMX_DOUBLE_EPS);

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;

    if (*pp == 0) {
	i__1 = 4*(*n0 - 3);
	for (j4 = *i0*4; j4 <= i__1; j4 += 4) {
	    z__[j4 - 2] = d__ + z__[j4 - 1];
	    if (std::abs(z__[j4 - 2])<PLUMED_GMX_DOUBLE_MIN) {
		z__[j4] = 0.;
		d__ = z__[j4 + 1];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 1] < z__[j4 - 2] && safemin * z__[j4 
		    - 2] < z__[j4 + 1]) {
		temp = z__[j4 + 1] / z__[j4 - 2];
		z__[j4] = z__[j4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
	    }
	    if(d__<*dmin__)
	      *dmin__ = d__;

	    d__1 = emin, d__2 = z__[j4];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
    } else {
	i__1 = 4*(*n0 - 3);
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
	    z__[j4 - 3] = d__ + z__[j4];
	    if (std::abs(z__[j4 - 3])<PLUMED_GMX_DOUBLE_MIN) {
		z__[j4 - 1] = 0.;
		d__ = z__[j4 + 2];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 2] < z__[j4 - 3] && safemin * z__[j4 
		    - 3] < z__[j4 + 2]) {
		temp = z__[j4 + 2] / z__[j4 - 3];
		z__[j4 - 1] = z__[j4] * temp;
		d__ *= temp;
	    } else {
		z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
	    }
	    if(d__<*dmin__)
	      *dmin__ = d__;
	    d__1 = emin, d__2 = z__[j4 - 1];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
    }

    *dnm2 = d__;
    *dmin2 = *dmin__;
    j4 = 4*(*n0 - 2) - *pp;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm2 + z__[j4p2];
    if (std::abs(z__[j4 - 2])<PLUMED_GMX_DOUBLE_MIN) {
	z__[j4] = 0.;
	*dnm1 = z__[j4p2 + 2];
	*dmin__ = *dnm1;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dnm1 = *dnm2 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    if(*dnm1<*dmin__)
      *dmin__ = *dnm1;

    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (std::abs(z__[j4 - 2])<PLUMED_GMX_DOUBLE_MIN) {
	z__[j4] = 0.;
	*dn = z__[j4p2 + 2];
	*dmin__ = *dn;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dn = *dnm1 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    if(*dn<*dmin__)
      *dmin__ = *dn;

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;


} 
}
}
#include <cmath>

#include "real.h"
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
    int i__, j;
    double temp;
    double ctemp, stemp;

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */

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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_DOUBLE_EPS || std::abs(stemp)>PLUMED_GMX_DOUBLE_MIN) {
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
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)(const char *id, 
	int *n, 
	double *d__, 
	int *info)
{
    int i__1, i__2;

    int i__, j;
    double d1, d2, d3;
    int dir;
    double tmp;
    int endd;
    int stack[64];
    double dmnmx;
    int start;
    int stkpnt;

    --d__;

    *info = 0;
    dir = -1;
    if (*id=='D' || *id=='d') 
	dir = 0;
    else if (*id=='I' || *id=='i') 
	dir = 1;
   
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }
    if (*n <= 1) {
	return;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {


	if (dir == 0) {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L30;
		    }
		}
L30:
		;
	    }

	} else {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L50;
		    }
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return;

}
}
}
#include "lapack.h"
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;

void PLUMED_BLAS_F77_FUNC(dlasrt2,DLASRT2)(const char *id, 
	      int *n, 
	      double *d__, 
	      int * key, 
	      int *info)
{
    int i__1, i__2;

    int i__, j;
    double d1, d2, d3;
    int dir;
    double tmp;
    int endd;
    int stack[64];
    double dmnmx;
    int start;
    int tmpkey, stkpnt;

    --key;
    --d__;

    *info = 0;
    dir = -1;
    if (*id=='D' || *id=='d')
	dir = 0;
    else if (*id=='I' || *id=='i')
	dir = 1;
    
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
        return;
    }

    if (*n <= 1) {
	return;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start > 0) {

	if (dir == 0) {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			break;
		    }
		}
	    }

	} else {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			break;
		    }
		}
	    }

	}

    } else if (endd - start > 20) {

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }

    return;
}
}
}
#include <cmath>
#include "real.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlassq,DLASSQ)(int *n,
                        double *x,
                        int *incx,
                        double *scale,
                        double *sumsq)
{
  int ix;
  double absxi,t;

  if(*n>0) {
    for(ix=0;ix<=(*n-1)*(*incx);ix+=*incx) {
      if(std::abs(x[ix])>PLUMED_GMX_DOUBLE_MIN) {
	absxi = std::abs(x[ix]);
	if(*scale<absxi) {
	  t = *scale/absxi;
	  t = t*t;
	  *sumsq = 1.0 + (*sumsq)*t;
	  *scale = absxi;
	} else {
	  t = absxi/(*scale);
	  *sumsq += t*t;
	}
      }
    }
  }
  return;
}
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dlasv2,DLASV2)(double *f, 
                        double *g, 
                        double *h__, 
                        double *ssmin, 
                        double *ssmax, 
                        double *snr, 
                        double *csr, 
                        double *snl, 
                        double *csl)
{
    double d__1;

    double a, d__, l, m, r__, s, t, fa, ga, ha, ft, gt, ht, mm, tt,
	     clt, crt, slt, srt;
    int pmax;
    double temp;
    int swap;
    double tsign=1.0;
    int gasmal;

    ft = *f;
    fa = std::abs(ft);
    ht = *h__;
    ha = std::abs(*h__);

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

    }
    gt = *g;
    ga = std::abs(gt);
    if (std::abs(ga)<PLUMED_GMX_DOUBLE_MIN) {

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = 1;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < PLUMED_GMX_DOUBLE_EPS) {

		gasmal = 0;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

	    d__ = fa - ha;
	    if ( std::abs( fa - d__ )<PLUMED_GMX_DOUBLE_EPS*std::abs( fa + d__ )) {
		l = 1.;
	    } else {
		l = d__ / fa;
	    }

	    m = gt / ft;
	    t = 2. - l;

	    mm = m * m;
	    tt = t * t;
	    s =  std::sqrt(tt + mm);

	    if ( std::abs(l)<PLUMED_GMX_DOUBLE_MIN) {
		r__ = std::abs(m);
	    } else {
		r__ =  std::sqrt(l * l + mm);
	    }
	    a = (s + r__) * .5;

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if ( std::abs(mm)<PLUMED_GMX_DOUBLE_MIN) {

		if (std::abs(l)<PLUMED_GMX_DOUBLE_MIN) {
		    t = ( (ft>0) ? 2.0 : -2.0) * ( (gt>0) ? 1.0 : -1.0);
		} else {
		    t = gt / ( (ft>0) ? d__ : -d__) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
	    }
	    l =  std::sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

    if (pmax == 1) {
	tsign = ( (*csr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*f>0) ? 1.0 : -1.0);
    }
    if (pmax == 2) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*g>0) ? 1.0 : -1.0);
    }
    if (pmax == 3) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*snl>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    }
    if(tsign<0)
      *ssmax *= -1.0;
    d__1 = tsign * ( (*f>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    if(d__1<0)
      *ssmin *= -1.0;
    return;

}
}
}
#include "lapack.h"

/* LAPACK */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlaswp,DLASWP)(int *n,
	double *a,
	int *lda,
	int *k1,
	int *k2,
	int *ipiv,
	int *incx)
{
  int ix0,i1,i2,inc,n32;
  int ix,i,j,ip,k;
  double temp;

  if(*incx>0) {
    ix0 = *k1 - 1;
    i1 = *k1 - 1;
    i2 = *k2;
    inc = 1;
  } else if(*incx<0) {
    ix0 = *incx * (1- *k2);
    i1 = *k2 - 1;
    i2 = *k1;
    inc = -1;
  } else
    return;

  n32 = *n / 32;
  
  n32 *= 32;


  if(n32!=0) {
    for(j=0;j<n32;j+=32) {
      ix = ix0;
      for(i=i1;i<i2;i+=inc,ix+=*incx) {
	ip = ipiv[ix] - 1;
	if(ip != i) {
	  for(k=j;k<j+32;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	  }
	}
      }
    }
  }
  if(n32!=*n) {
    ix = ix0;
    for(i=i1;i<i2;i+=inc,ix+=*incx) {
      ip = ipiv[ix] - 1;
      if(ip != i) {
	for(k=n32;k<*n;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	}
      }
    }
  }
  return;
}
}
}
#include <cctype>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dlatrd,DLATRD)(const char *  uplo,
       int  *   n,
       int  *   nb,
       double * a,
       int *    lda,
       double * e,
       double * tau,
       double * w,
       int *    ldw)
{
  int i,iw;
  int ti1,ti2,ti3;
  double one,zero,minusone,alpha;
  const char ch=std::toupper(*uplo);

  one=1.0;
  minusone=-1.0;
  zero=0.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n;i>=(*n-*nb+1);i--) {
      iw = i -*n + *nb;
      
      if(i<*n) {
	ti1 = *n-i;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&i,&ti1,&minusone, &(a[ i*(*lda) + 0]),lda,&(w[iw*(*ldw)+(i-1)]),
	       ldw,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&i,&ti1,&minusone, &(w[ iw*(*ldw) + 0]),ldw,&(a[i*(*lda)+(i-1)]),
	       lda,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
      }

      if(i>1) {
	/*  Generate elementary reflector H(i) to annihilate
	 *              A(1:i-2,i) 
	 */
	ti1 = i-1;
	ti2 = 1;

	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&ti1,&(a[(i-1)*(*lda)+(i-2)]),&(a[(i-1)*(*lda)+0]),&ti2,&(tau[i-2]));
      
	e[i-2] = a[(i-1)*(*lda)+(i-2)];
	a[(i-1)*(*lda)+(i-2)] = 1.0;

	/* Compute W(1:i-1,i) */
	ti1 = i-1;
	ti2 = 1;

	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)("U",&ti1,&one,a,lda,&(a[(i-1)*(*lda)+0]),&ti2,&zero,
	       &(w[(iw-1)*(*ldw)+0]),&ti2);
	if(i<*n) {
	  ti1 = i-1;
	  ti2 = *n-i;
	  ti3 = 1;
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("T",&ti1,&ti2,&one,&(w[iw*(*ldw)+0]),ldw,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone,&(a[i*(*lda)+0]),lda,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("T",&ti1,&ti2,&one,&(a[i*(*lda)+0]),lda,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone,&(w[iw*(*ldw)+0]),ldw,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	}
      
	ti1 = i-1;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&ti1,&(tau[i-2]),&(w[(iw-1)*(*ldw)+0]),&ti2);
      
	alpha = -0.5*tau[i-2]*PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&ti1,&(w[(iw-1)*(*ldw)+0]),&ti2,
				    &(a[(i-1)*(*lda)+0]),&ti2);
      
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+0]),&ti2,&(w[(iw-1)*(*ldw)+0]),&ti2);

      }
    }
  } else {
    /* lower */
    for(i=1;i<=*nb;i++) {

      ti1 = *n-i+1;
      ti2 = i-1;
      ti3 = 1;
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone, &(a[ i-1 ]),lda,&(w[ i-1 ]),
	       ldw,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone, &(w[ i-1 ]),ldw,&(a[ i-1 ]),
	       lda,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);

      if(i<*n) {
	ti1 = *n - i;
	ti2 = (*n < i+2 ) ? *n : (i+2);
	ti3 = 1;
	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+(ti2-1)]),&ti3,&(tau[i-1]));
	e[i-1] = a[(i-1)*(*lda)+(i)];
	a[(i-1)*(*lda)+(i)] = 1.0;
	
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)("L",&ti1,&one,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),&ti2,
	       &zero,&(w[(i-1)*(*ldw)+i]),&ti2);
	ti1 = *n - i;
	ti2 = i-1;
	ti3 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("T",&ti1,&ti2,&one,&(w[ i ]),ldw,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone,&(a[ i ]),lda,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("T",&ti1,&ti2,&one,&(a[ i ]),lda,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(dgemv,DGEMV)("N",&ti1,&ti2,&minusone,&(w[ i ]),ldw,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);

	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&ti1,&(tau[i-1]),&(w[(i-1)*(*ldw)+i]),&ti3);
	alpha = -0.5*tau[i-1]*PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&ti1,&(w[(i-1)*(*ldw)+i]),&ti3,
				   &(a[(i-1)*(*lda)+i]),&ti3);
	
	PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti3,&(w[(i-1)*(*ldw)+i]),&ti3);
      }
    }
  }
  return;
}
	


  
}
}
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dorg2r,DORG2R)(int *m, 
                        int *n,
                        int *k, 
                        double *a, 
                        int *lda,
                        double *tau,
                        double *work,
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    double r__1;
    int c__1 = 1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;

    if (*n <= 0) {
        return;
    }

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a[l + j * a_dim1] = 0.0;
	}
	a[j + j * a_dim1] = 1.0;
    }
    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    a[i__ + i__ * a_dim1] = 1.0;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__;
	    PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("L", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, 
                              &tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
	}
	if (i__ < *m) {
	    i__1 = *m - i__;
	    r__1 = -tau[i__];
	    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__1, &r__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[l + i__ * a_dim1] = 0.0;
	}
    }
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dorgbr,DORGBR)(const char *vect,
	int *m,
	int *n,
	int *k,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *lwork,
	int *info)
{
  int wantq,iinfo,j,i,i1,wrksz;
  int mn = (*m < *n) ? *m : *n;

  wantq = (*vect=='Q' || *vect=='q');

  *info = 0;
  wrksz = mn*DORGBR_BLOCKSIZE;
  if(*lwork==-1) {
    work[0] = wrksz;
    return;
  }
  
  if(*m==0 || *n==0)
    return;

  if(wantq) {
    if(*m>=*k)
      PLUMED_BLAS_F77_FUNC(dorgqr,DORGQR)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      for(j=*m;j>=2;j--) {
	a[(j-1)*(*lda)+0] = 0.0;
	for(i=j+1;i<=*m;i++)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-2)*(*lda)+(i-1)]; 
      }
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      if(*m>1) {
	i1 = *m-1;
	PLUMED_BLAS_F77_FUNC(dorgqr,DORGQR)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  } else {
    if(*k<*n)
      PLUMED_BLAS_F77_FUNC(dorglq,DORGLQ)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      for(j=2;j<=*n;j++) {
	for(i=j-1;i>=2;i--)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-1)*(*lda)+(i-2)]; 
	a[(j-1)*(*lda)+0] = 0.0;
      }
      if(*n>1) {
	i1 = *n-1;
	PLUMED_BLAS_F77_FUNC(dorglq,DORGLQ)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  }
  work[0] = wrksz;
  return;
}
 
}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dorgl2,DORGL2)(int *m,
                        int *n, 
                        int *k, 
                        double *a, 
                        int *lda, 
                        double *tau, 
                        double *work, 
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    double r__1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    i__ = (*m > 1) ? *m : 1;
    
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < i__) {
	*info = -5;
    }
    if (*info != 0) {
	return;
    }
    if (*m <= 0) {
	return;
    }

    if (*k < *m) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (l = *k + 1; l <= i__2; ++l) {
		a[l + j * a_dim1] = 0.0;
	    }
	    if (j > *k && j <= *m) {
		a[j + j * a_dim1] = 1.0;
	    }
	}
    }

    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.0;
		i__1 = *m - i__;
		i__2 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarf,DLARF)("R", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, 
               &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    }
	    i__1 = *n - i__;
	    r__1 = -tau[i__];
	    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__1, &r__1, &a[i__ + (i__ + 1) * a_dim1], lda);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[i__ + l * a_dim1] = 0.0;
	}
    }
    return;

}



}
}
#include "lapack.h"

#define DORGLQ_BLOCKSIZE    32
#define DORGLQ_MINBLOCKSIZE 2
#define DORGLQ_CROSSOVER    128


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dorglq,DORGLQ)(int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;

    int ldwork, lwkopt;
    int lquery;
    
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    ki = 0;
    nb = DORGLQ_BLOCKSIZE;
    lwkopt = (*m) * nb;
    work[1] = (double) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*m) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

	nx = DORGLQ_CROSSOVER;
	if (nx < *k) {

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DORGLQ_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = kk;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }
    if (kk < *m) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	PLUMED_BLAS_F77_FUNC(dorgl2,DORGL2)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *m) {

		i__2 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ - ib + 1;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)("Right", "Transpose", "Forward", "Rowwise", &i__2, &
			i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }

	    i__2 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dorgl2,DORGL2)(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = i__ + ib - 1;
		for (l = i__; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (double) iws;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dorgqr,DORGQR)(int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;
    int lquery;
 
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    ki = 0;
    *info = 0;
    nb = DORGQR_BLOCKSIZE;
    lwkopt = (*n) * nb;
    work[1] = (double) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*n <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

	nx = DORGQR_CROSSOVER;
	if (nx < *k) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DORGQR_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = *n;
	for (j = kk + 1; j <= i__1; ++j) {
	    i__2 = kk;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	PLUMED_BLAS_F77_FUNC(dorg2r,DORG2R)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *n) {

		i__2 = *m - i__ + 1;
		PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ + 1;
		i__3 = *n - i__ - ib + 1;
		PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
			1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &
			work[ib + 1], &ldwork);
	    }

	    i__2 = *m - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dorg2r,DORG2R)(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ + ib - 1;
	    for (j = i__; j <= i__2; ++j) {
		i__3 = i__ - 1;
		for (l = 1; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (double) iws;
    return;

} 
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dorm2l,DORM2L)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	double *a,
	int *lda, 
	double *tau,
	double *c__,
	int *ldc, 
	double *work, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    int c__1 = 1;

    int i__, i1, i2, i3, mi, ni, nq;
    double aii;
    int left;
    int notran;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (*info != 0) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	return;
    }

    if ((left && notran) || (! left && ! notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
    } else {
	mi = *m;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - *k + i__;
	} else {

	    ni = *n - *k + i__;
	}

	aii = a[nq - *k + i__ + i__ * a_dim1];
	a[nq - *k + i__ + i__ * a_dim1] = 1.;
	PLUMED_BLAS_F77_FUNC(dlarf,DLARF)(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[
		c_offset], ldc, &work[1]);
	a[nq - *k + i__ + i__ * a_dim1] = aii;
    }
    return;
}
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dorm2r,DORM2R)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc, 
	double *work, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    int i__, i1, i2, i3, ic, jc, mi, ni;
    double aii;
    int left;
    int notran;
    int c__1 = 1;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');

    ic = jc = 0;

    if (*m <= 0 || *n <= 0 || *k <= 0) {
	return;
    }

    if ((left && !notran) || (!left && notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

	    ni = *n - i__ + 1;
	    jc = i__;
	}


	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	PLUMED_BLAS_F77_FUNC(dlarf,DLARF)(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1]);
	a[i__ + i__ * a_dim1] = aii;
    }
    return;

} 
}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dormbr,DORMBR)(const char *vect, 
	const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1;
 

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int notran;
    int applyq;
    char transt[1];
    int lwkopt;
    int lquery;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    applyq = (*vect=='Q' || *vect=='q');
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMQR_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    work[1] = 1.;
    if (*m == 0 || *n == 0) {
	return;
    }

    if (applyq) {

	if (nq >= *k) {

	    PLUMED_BLAS_F77_FUNC(dormqr,DORMQR)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    PLUMED_BLAS_F77_FUNC(dormqr,DORMQR)(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
	}
    } else {

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}
	if (nq > *k) {

	    PLUMED_BLAS_F77_FUNC(dormlq,DORMLQ)(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    PLUMED_BLAS_F77_FUNC(dormlq,DORMLQ)(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo);
	}
    }
    work[1] = (double) lwkopt;
    return;


}


}
}
#include <cctype>
#include "lapack.h"
#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dorml2,DORML2)(const char *side,
	const char *trans,
	int *m,
	int *n,
	int *k,
	double *a,
	int *lda,
	double *tau,
	double *c,
	int *ldc,
	double *work,
    int *info)
{
  const char xside=std::toupper(*side);
  const char xtrans=std::toupper(*trans);
  int i,i1,i2,i3,ni,mi,ic,jc;
  double aii;

  if(*m<=0 || *n<=0 || *k<=0)
    return;

  ic = jc = 0;

  if((xside=='L' && xtrans=='N') || (xside!='L' && xtrans!='N')) {
    i1 = 0;
    i2 = *k;
    i3 = 1;
  } else {
    i1 = *k-1;
    i2 = -1;
    i3 = -1;
  }
  
  if(xside=='L') {
    ni = *n;
    jc = 0;
  } else {
    mi = *m;
    ic = 0;
  }

  for(i=i1;i!=i2;i+=i3) {
    if(xside=='L') {
      mi = *m - i;
      ic = i;
    } else {
      ni = *n - i;
      jc = i;
    }
    aii = a[i*(*lda)+i];
    a[i*(*lda)+i] = 1.0;
    PLUMED_BLAS_F77_FUNC(dlarf,DLARF)(side,&mi,&ni,&(a[i*(*lda)+i]),lda,tau+i,
	   &(c[jc*(*ldc)+ic]),ldc,work);
    a[i*(*lda)+i] = aii;
  }
  return;
}
	     
}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dormlq,DORMLQ)(const char *side, 
	const char *trans,
	int *m, 
	int *n, 
	int *k,
	double *a,
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, 
	    i__5;
  

    int i__;
    double t[4160]	/* was [65][64] */;
    int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork;
    char transt[1];
    int lwkopt;
    int lquery;
    int ldt = 65;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    ic = jc = 0;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMLQ_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
       	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMLQ_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {


	PLUMED_BLAS_F77_FUNC(dorml2,DORML2)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && notran) || (!left && !notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], 
		    lda, &tau[i__], t, &ldt);
	    if (left) {

		mi = *m - i__ + 1;
		ic = i__;
	    } else {

		ni = *n - i__ + 1;
		jc = i__;
	    }

	    PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ 
		    + i__ * a_dim1], lda, t, &ldt, &c__[ic + jc * c_dim1], 
		    ldc, &work[1], &ldwork);
	}
    }
    work[1] = (double) lwkopt;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dormql,DORMQL)(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;
    int c__65 = 65;

    int i__;
    double t[4160];
    int i1, i2, i3, ib, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork, lwkopt;
    int lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMQL_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMQL_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {

	PLUMED_BLAS_F77_FUNC(dorm2l,DORM2L)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && notran) || (! left && ! notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	} else {
	    mi = *m;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - *k + i__ + ib - 1;
	    PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], t, &c__65);
	    if (left) {

		mi = *m - *k + i__ + ib - 1;
	    } else {

		ni = *n - *k + i__ + ib - 1;
	    }

	    PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, t, &c__65, &c__[c_offset], ldc, &
		    work[1], &ldwork);
	}
    }
    work[1] = (double) lwkopt;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(dormqr,DORMQR)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc, 
	double *work, 
	int *lwork, 
	int *info)
{
   int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;

    int i__;
    double t[4160];
    int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork, lwkopt;
    int lquery;
    int ldt = 65;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

     ic = jc = 0;
     nb = DORMQR_BLOCKSIZE;
     lwkopt = nw * nb;
     work[1] = (double) lwkopt;

    if (*info != 0) {
	return;
    } else if (lquery) {
      return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMQR_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {

	PLUMED_BLAS_F77_FUNC(dorm2r,DORM2R)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && !notran) || (!left && notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dlarft,DLARFT)("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], t, &ldt);
	    if (left) {

		mi = *m - i__ + 1;
		ic = i__;
	    } else {
		ni = *n - i__ + 1;
		jc = i__;
	    }

	    PLUMED_BLAS_F77_FUNC(dlarfb,DLARFB)(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, t, &ldt, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork);
	}
    }
    work[1] = (double) lwkopt;
    return;


}


}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dormtr,DORMTR)(const char *side, 
	const char *uplo,
	const char *trans, 
	int *m, 
	int *n,
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc,
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__2;

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int upper;
    int lwkopt;
    int lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    upper = (*uplo=='U' || *uplo=='u');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }


    nb = DORMQL_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
	i__2 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || nq == 1) {
	work[1] = 1.;
	return;
    }

    if (left) {
	mi = *m - 1;
	ni = *n;
    } else {
	mi = *m;
	ni = *n - 1;
    }

    if (upper) {
	i__2 = nq - 1;
	PLUMED_BLAS_F77_FUNC(dormql,DORMQL)(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &
		tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo);
    } else {
	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	i__2 = nq - 1;
	PLUMED_BLAS_F77_FUNC(dormqr,DORMQR)(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &
		c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
    }
    work[1] = (double) lwkopt;
    return;

}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dstebz,DSTEBZ)(const char *range, 
                        const char *order,
                        int *n,
                        double *vl, 
                        double *vu, 
                        int *il,
                        int *iu,
                        double *abstol, 
                        double *d__,
                        double *e,
                        int *m, 
                        int *nsplit, 
                        double *w,
                        int *iblock,
                        int *isplit,
                        double *work, 
                        int *iwork, 
                        int *info)
{
    int i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4, d__5;
    int c__1 = 1;
    int c__3 = 3;
    int c__2 = 2;
    int c__0 = 0;

    int j, ib, jb, ie, je, nb;
    double gl;
    int im, in;
    double gu;
    int iw;
    double wl, wu;
    int nwl;
    double ulp, wlu, wul;
    int nwu;
    double tmp1, tmp2;
    int iend, ioff, iout, itmp1, jdisc;
    int iinfo;
    double atoli;
    int iwoff;
    double bnorm;
    int itmax;
    double wkill, rtoli, tnorm;
    int ibegin;
    int irange, idiscl;
    int idumma[1];
    int idiscu, iorder;
    int ncnvrg;
    double pivmin;
    int toofew;
    const double safemn = PLUMED_GMX_DOUBLE_MIN*(1.0+PLUMED_GMX_DOUBLE_EPS);

    --iwork;
    --work;
    --isplit;
    --iblock;
    --w;
    --e;
    --d__;

    *info = 0;

    if (*range=='A' || *range=='a') {
	irange = 1;
    } else if (*range=='V' || *range=='v') {
	irange = 2;
    } else if (*range=='I' || *range=='i') {
	irange = 3;
    } else {
	irange = 0;
    }

    if (*order=='B' || *order=='b') {
	iorder = 2;
    } else if (*order=='E' || *order=='e') {
	iorder = 1;
    } else {
	iorder = 0;
    }

    if (irange <= 0) {
	*info = -1;
    } else if (iorder <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (irange == 2) {
	if (*vl >= *vu) {
	    *info = -5;
	}
    } else if (irange == 3 && (*il < 1 || *il > (*n))) {
	*info = -6;
    } else if (irange == 3 && (*iu < ((*n<*il) ? *n : *il) || *iu > *n)) {
	*info = -7;
    }

    if (*info != 0) {
	return;
    }

    *info = 0;
    ncnvrg = 0;
    toofew = 0;

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (irange == 3 && *il == 1 && *iu == *n) {
	irange = 1;
    }

    ulp = 2*PLUMED_GMX_DOUBLE_EPS;
    rtoli = ulp * 2.;
    nb = DSTEBZ_BLOCKSIZE;
    // cppcheck-suppress knownConditionTrueFalse
    if (nb <= 1) {
	nb = 0;
    }

    if (*n == 1) {
	*nsplit = 1;
	isplit[1] = 1;
	if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
	    *m = 0;
	} else {
	    w[1] = d__[1];
	    iblock[1] = 1;
	    *m = 1;
	}
	return;
    }

    *nsplit = 1;
    work[*n] = 0.;
    pivmin = 1.;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	d__1 = e[j - 1];
	tmp1 = d__1 * d__1;
	d__2 = ulp;
	if (std::abs(d__[j] * d__[j - 1]) * (d__2 * d__2) + safemn 
		> tmp1) {
	    isplit[*nsplit] = j - 1;
	    ++(*nsplit);
	    work[j - 1] = 0.;
	} else {
	    work[j - 1] = tmp1;
	    pivmin = (pivmin>tmp1) ? pivmin : tmp1;
	}
    }
    isplit[*nsplit] = *n;
    pivmin *= safemn;

    if (irange == 3) {

	gu = d__[1];
	gl = d__[1];
	tmp1 = 0.;

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    tmp2 =  std::sqrt(work[j]);
	    d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
	    gu = (d__1>d__2) ? d__1 : d__2;
	    d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
	    gl = (d__1<d__2) ? d__1 : d__2;
	    tmp1 = tmp2;
	}

	d__1 = gu, d__2 = d__[*n] + tmp1;
	gu = (d__1>d__2) ? d__1 : d__2;
	d__1 = gl, d__2 = d__[*n] - tmp1;
	gl = (d__1<d__2) ? d__1 : d__2;
	d__1 = std::abs(gl);
	d__2 = std::abs(gu);
	tnorm = (d__1>d__2) ? d__1 : d__2;
	gl = gl - tnorm * 2. * ulp * *n - pivmin * 4.;
	gu = gu + tnorm * 2. * ulp * *n + pivmin * 2.;

	itmax = (int) ((std::log(tnorm + pivmin) - std::log(pivmin)) / std::log(2.)) + 2;
	if (*abstol <= 0.) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = *abstol;
	}

	work[*n + 1] = gl;
	work[*n + 2] = gl;
	work[*n + 3] = gu;
	work[*n + 4] = gu;
	work[*n + 5] = gl;
	work[*n + 6] = gu;
	iwork[1] = -1;
	iwork[2] = -1;
	iwork[3] = *n + 1;
	iwork[4] = *n + 1;
	iwork[5] = *il - 1;
	iwork[6] = *iu;

	PLUMED_BLAS_F77_FUNC(dlaebz,DLAEBZ)(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, 
		&d__[1], &e[1], &work[1], &iwork[5], &work[*n + 1], &work[*n 
		+ 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

	if (iwork[6] == *iu) {
	    wl = work[*n + 1];
	    wlu = work[*n + 3];
	    nwl = iwork[1];
	    wu = work[*n + 4];
	    wul = work[*n + 2];
	    nwu = iwork[4];
	} else {
	    wl = work[*n + 2];
	    wlu = work[*n + 4];
	    nwl = iwork[2];
	    wu = work[*n + 3];
	    wul = work[*n + 1];
	    nwu = iwork[3];
	}

	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
	    *info = 4;
	    return;
	}
    } else {


      /* avoid warnings for high gcc optimization */
      wlu = wul = 1.0;

	d__3 = std::abs(d__[1]) + std::abs(e[1]);
	d__4 = std::abs(d__[*n]) + std::abs(e[*n - 1]);
	tnorm = (d__3>d__4) ? d__3 : d__4;

	i__1 = *n - 1;
	for (j = 2; j <= i__1; ++j) {
	    d__4 = tnorm;
	    d__5 = std::abs(d__[j]) + std::abs(e[j - 1]) + std::abs(e[j]);
	    tnorm = (d__4>d__5) ? d__4 : d__5;
	}

	if (*abstol <= 0.) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = *abstol;
	}

	if (irange == 2) {
	    wl = *vl;
	    wu = *vu;
	} else {
	    wl = 0.;
	    wu = 0.;
	}
    }

    *m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;

    i__1 = *nsplit;
    for (jb = 1; jb <= i__1; ++jb) {
	ioff = iend;
	ibegin = ioff + 1;
	iend = isplit[jb];
	in = iend - ioff;

	if (in == 1) {

	    if (irange == 1 || wl >= d__[ibegin] - pivmin) {
		++nwl;
	    }
	    if (irange == 1 || wu >= d__[ibegin] - pivmin) {
		++nwu;
	    }
	    if (irange == 1 || ((wl < d__[ibegin] - pivmin) && (wu >= d__[ibegin] - pivmin))) {
		++(*m);
		w[*m] = d__[ibegin];
		iblock[*m] = jb;
	    }
	} else {

	    gu = d__[ibegin];
	    gl = d__[ibegin];
	    tmp1 = 0.;

	    i__2 = iend - 1;
	    for (j = ibegin; j <= i__2; ++j) {
		tmp2 = std::abs(e[j]);
		d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
		gu = (d__1>d__2) ? d__1 : d__2;
		d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
		gl = (d__1<d__2) ? d__1 : d__2;
		tmp1 = tmp2;
	    }

	    d__1 = gu, d__2 = d__[iend] + tmp1;
	    gu = (d__1>d__2) ? d__1 : d__2;
	    d__1 = gl, d__2 = d__[iend] - tmp1;
	    gl = (d__1<d__2) ? d__1 : d__2;
	    d__1 = std::abs(gl);
	    d__2 = std::abs(gu);
	    bnorm = (d__1>d__2) ? d__1 : d__2;
	    gl = gl - bnorm * 2. * ulp * in - pivmin * 2.;
	    gu = gu + bnorm * 2. * ulp * in + pivmin * 2.;

	    if (*abstol <= 0.) {
		d__1 = std::abs(gl);
		d__2 = std::abs(gu);
		atoli = ulp * ((d__1>d__2) ? d__1 : d__2);
	    } else {
		atoli = *abstol;
	    }

	    if (irange > 1) {
		if (gu < wl) {
		    nwl += in;
		    nwu += in;
		}
		gl = (gl>wl) ? gl : wl;
		gu = (gu<wu) ? gu : wu;
		if (gl >= gu) {
		}
		continue;
	    }

	    work[*n + 1] = gl;
	    work[*n + in + 1] = gu;
	    PLUMED_BLAS_F77_FUNC(dlaebz,DLAEBZ)(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);

	    nwl += iwork[1];
	    nwu += iwork[in + 1];
	    iwoff = *m - iwork[1];

	    itmax = (int) ((std::log(gu - gl + pivmin) - std::log(pivmin)) / std::log(2.)
		    ) + 2;
	    PLUMED_BLAS_F77_FUNC(dlaebz,DLAEBZ)(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);

	    i__2 = iout;
	    for (j = 1; j <= i__2; ++j) {
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

		if (j > iout - iinfo) {
		    ncnvrg = 1;
		    ib = -jb;
		} else {
		    ib = jb;
		}
		i__3 = iwork[j + in] + iwoff;
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
		    w[je] = tmp1;
		    iblock[je] = ib;
		}
	    }

	    *m += im;
	}
    }

    if (irange == 3) {
	im = 0;
	idiscl = *il - 1 - nwl;
	idiscu = nwu - *iu;

	if (idiscl > 0 || idiscu > 0) {
	    i__1 = *m;
	    for (je = 1; je <= i__1; ++je) {
		if (w[je] <= wlu && idiscl > 0) {
		    --idiscl;
		} else if (w[je] >= wul && idiscu > 0) {
		    --idiscu;
		} else {
		    ++im;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    *m = im;
	}
	if (idiscl > 0 || idiscu > 0) {

	    if (idiscl > 0) {
		wkill = wu;
		i__1 = idiscl;
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
		    iw = 0;
		    i__2 = *m;
		    for (je = 1; je <= i__2; ++je) {
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}
		    }
		    iblock[iw] = 0;
		}
	    }
	    if (idiscu > 0) {

		wkill = wl;
		i__1 = idiscu;
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
		    iw = 0;
		    i__2 = *m;
		    for (je = 1; je <= i__2; ++je) {
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}
		    }
		    iblock[iw] = 0;
		}
	    }
	    im = 0;
	    i__1 = *m;
	    for (je = 1; je <= i__1; ++je) {
		if (iblock[je] != 0) {
		    ++im;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    *m = im;
	}
	if (idiscl < 0 || idiscu < 0) {
	    toofew = 1;
	}
    }

    if (iorder == 1 && *nsplit > 1) {
	i__1 = *m - 1;
	for (je = 1; je <= i__1; ++je) {
	    ie = 0;
	    tmp1 = w[je];
	    i__2 = *m;
	    for (j = je + 1; j <= i__2; ++j) {
		if (w[j] < tmp1) {
		    ie = j;
		    tmp1 = w[j];
		}
	    }

	    if (ie != 0) {
		itmp1 = iblock[ie];
		w[ie] = w[je];
		iblock[ie] = iblock[je];
		w[je] = tmp1;
		iblock[je] = itmp1;
	    }
	}
    }

    *info = 0;
    if (ncnvrg) {
	++(*info);
    }
    if (toofew) {
	*info += 2;
    }
    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dstegr,DSTEGR)(const char *jobz, 
	const char *range, 
	int *n, 
	double *d__, 
	double *e, 
	double *vl, 
	double *vu, 
	int *il, 
	int *iu, 
	double *abstol, 
	int *m, 
	double *w, 
	double *z__, 
	int *ldz, 
	int *isuppz,
	double *work, 
	int *lwork, 
	int *iwork, 
	int *liwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;
    int c__1 = 1;

    int i__, j;
    int jj;
    double eps, tol, tmp, rmin, rmax;
    int itmp;
    double tnrm;
    double scale;
    int iinfo, iindw;
    int lwmin;
    int wantz;
    int iindbl;
    int valeig,alleig,indeig;
    double safmin,minval;
    double bignum;
    int iindwk, indgrs;
    double thresh;
    int iinspl, indwrk, liwmin, nsplit;
    double smlnum;
    int lquery;


    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    lquery = *lwork == -1 || *liwork == -1;
    lwmin = *n * 17;
    liwmin = *n * 10;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && (*il < 1 || *il > *n)) {
	*info = -8;
    } else if (indeig && (*iu < *il || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -14;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -17;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -19;
    }
    if (*info == 0) {
	work[1] = (double) lwmin;
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    minval = PLUMED_GMX_DOUBLE_MIN;
    safmin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS);
    eps = PLUMED_GMX_DOUBLE_EPS;
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin =  std::sqrt(smlnum);
    d__1 =  std::sqrt(bignum), d__2 = 1. / std::sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;
    scale = 1.;
    tnrm = PLUMED_BLAS_F77_FUNC(dlanst,DLANST)("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin) {
	scale = rmin / tnrm;
    } else if (tnrm > rmax) {
	scale = rmax / tnrm;
    }
    if ( std::abs(scale-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(n, &scale, &d__[1], &c__1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__1, &scale, &e[1], &c__1);
	tnrm *= scale;
    }
    indgrs = 1;
    indwrk = (*n << 1) + 1;

    iinspl = 1;
    iindbl = *n + 1;
    iindw = (*n << 1) + 1;
    iindwk = *n * 3 + 1;

    thresh = eps * tnrm;
    PLUMED_BLAS_F77_FUNC(dlarrex,DLARREX)(range, n, vl, vu, il, iu, &d__[1], &e[1], &thresh, &nsplit, &
	    iwork[iinspl], m, &w[1], &iwork[iindbl], &iwork[iindw], &work[
	    indgrs], &work[indwrk], &iwork[iindwk], &iinfo);
    
    if (iinfo != 0) {
	*info = 1;
	return;
    }

    if (wantz) {
	d__1 = *abstol, d__2 = (double) (*n) * eps;
	tol = (d__1>d__2) ? d__1 : d__2;
	PLUMED_BLAS_F77_FUNC(dlarrvx,DLARRVX)(n, &d__[1], &e[1], &iwork[iinspl], m, &w[1], &iwork[iindbl], &
		iwork[iindw], &work[indgrs], &tol, &z__[z_offset], ldz, &
		isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
	if (iinfo != 0) {
	    *info = 2;
	    return;
	}
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	itmp = iwork[iindbl + j - 1];
	w[j] += e[iwork[iinspl + itmp - 1]];
    } 

    if (std::abs(scale-1.0)>PLUMED_GMX_DOUBLE_EPS) {
	d__1 = 1. / scale;
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(m, &d__1, &w[1], &c__1);
    }
    if (nsplit > 1) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp) {
		    i__ = jj;
		    tmp = w[jj];
		}
	    }
	    if (i__ != 0) {
		w[i__] = w[j];
		w[j] = tmp;
		if (wantz) {
		    PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 
			    + 1], &c__1);
		    itmp = isuppz[(i__ << 1) - 1];
		    isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
		    isuppz[(j << 1) - 1] = itmp;
		    itmp = isuppz[i__ * 2];
		    isuppz[i__ * 2] = isuppz[j * 2];
		    isuppz[j * 2] = itmp;
		}
	    }
	}
    }

    work[1] = (double) lwmin;
    iwork[1] = liwmin;
    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dstein,DSTEIN)(int *n, 
	double *d__, 
	double *e, 
	int *m, 
	double *w, 
	int *iblock,
	int *isplit, 
	double *z__,
	int *ldz, 
	double *work,
	int *iwork, 
	int *ifail,
	int *info)
{
    int z_dim1, z_offset, i__1, i__2, i__3;
    double d__2, d__3, d__4, d__5;

    int i__, j, b1, j1, bn;
    double xj, scl, eps, sep, nrm, tol;
    int its;
    double xjm, ztr, eps1;
    int jblk, nblk;
    int jmax;

    int iseed[4], gpind, iinfo;
    double ortol;
    int indrv1, indrv2, indrv3, indrv4, indrv5;
    int nrmchk;
    int blksiz;
    double onenrm, dtpcrt, pertol;
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
	return;
    }

    if (*n == 0 || *m == 0) {
	return;
    } else if (*n == 1) {
	z__[z_dim1 + 1] = 1.;
	return;
    }

    eps = PLUMED_GMX_DOUBLE_EPS;

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

	onenrm = std::abs(d__[b1]) + std::abs(e[b1]);
	d__3 = onenrm;
	d__4 = std::abs(d__[bn]) + std::abs(e[bn - 1]);
	onenrm = (d__3>d__4) ? d__3 : d__4;
	i__2 = bn - 1;
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
	  d__4 = onenrm;
	  d__5 = std::abs(d__[i__]) + std::abs(e[i__ - 1]) + std::abs(e[i__]);
	    onenrm = (d__4>d__5) ? d__4 : d__5;
	}
	ortol = onenrm * .001;

	dtpcrt =  std::sqrt(.1 / blksiz);

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
		eps1 = std::abs(eps * xj);
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

	    PLUMED_BLAS_F77_FUNC(dlarnv,DLARNV)(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
	    i__3 = blksiz - 1;
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
	    i__3 = blksiz - 1;
	    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

	    tol = 0.;
	    PLUMED_BLAS_F77_FUNC(dlagtf,DLAGTF)(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

	    d__2 = eps;
	    d__3 = std::abs(work[indrv4 + blksiz]);
	    scl = blksiz * onenrm * ((d__2>d__3) ? d__2 : d__3) / PLUMED_BLAS_F77_FUNC(dasum,DASUM)(&blksiz, &work[
		    indrv1 + 1], &c__1);
	    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&blksiz, &scl, &work[indrv1 + 1], &c__1);

	    PLUMED_BLAS_F77_FUNC(dlagts,DLAGTS)(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

	    if (jblk == 1) {
		goto L90;
	    }
	    if (std::abs(xj - xjm) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i__ = gpind; i__ <= i__3; ++i__) {
		    ztr = -PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
			    i__ * z_dim1], &c__1);
		    PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &
			    work[indrv1 + 1], &c__1);
		}
	    }

L90:
	    jmax = PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)(&blksiz, &work[indrv1 + 1], &c__1);
	    nrm = std::abs(work[indrv1 + jmax]);

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
	    scl = 1. / PLUMED_BLAS_F77_FUNC(dnrm2,DNRM2)(&blksiz, &work[indrv1 + 1], &c__1);
	    jmax = PLUMED_BLAS_F77_FUNC(idamax,IDAMAX)(&blksiz, &work[indrv1 + 1], &c__1);
	    if (work[indrv1 + jmax] < 0.) {
		scl = -scl;
	    }
	    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&blksiz, &scl, &work[indrv1 + 1], &c__1);
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
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dsteqr,DSTEQR)(const char *    compz, 
	int *     n, 
	double *  d__, 
	double *  e, 
	double *  z__, 
	int *     ldz, 
	double *  work, 
	int *     info)
{
    double c_b9 = 0.;
    double c_b10 = 1.;
    int c__0 = 0;
    int c__1 = 1;
    int c__2 = 2;
    int z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;

    double b, c__, f, g;
    int i__, j, k, l, m;
    double p, r__, s;
    int l1, ii, mm, lm1, mm1, nm1;
    double rt1, rt2, eps;
    int lsv;
    double tst, eps2;
    int lend, jtot;
    double anorm;
    int lendm1, lendp1;
    int iscale;
    double safmin,minval;
    double safmax;
    int lendsv;
    double ssfmin;
    int nmaxit, icompz;
    double ssfmax;


    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    *info = 0;

    if (*compz=='N' || *compz=='n') {
	icompz = 0;
    } else if (*compz=='V' || *compz=='v') {
	icompz = 1;
    } else if (*compz=='I' || *compz=='i') {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < ((*n>1) ? *n : 1))) {
	*info = -6;
    }
    if (*info != 0) {
	return;
    }


    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    eps = PLUMED_GMX_DOUBLE_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    minval = PLUMED_GMX_DOUBLE_MIN;
    safmin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS);

    safmax = 1. / safmin;
    ssfmax =  std::sqrt(safmax) / 3.;
    ssfmin =  std::sqrt(safmin) / eps2;

    if (icompz == 2) {
	PLUMED_BLAS_F77_FUNC(dlaset,DLASET)("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= i__1; ++m) {
	    tst = std::abs(e[m]);
	    if (std::abs(tst)<PLUMED_GMX_DOUBLE_MIN) {
		goto L30;
	    }
	    if (tst <=  std::sqrt(std::abs(d__[m])) * std::sqrt(std::abs(d__[m + 1])) * eps) {
		e[m] = 0.;
		goto L30;
	    }
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    i__1 = lend - l + 1;
    anorm = PLUMED_BLAS_F77_FUNC(dlanst,DLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (std::abs(anorm)<PLUMED_GMX_DOUBLE_MIN) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    if (std::abs(d__[lend]) < std::abs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= i__1; ++m) {
  	        d__2 = std::abs(e[m]);
		tst = d__2 * d__2;
		if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m+ 1]) + safmin) {
		    goto L60;
		}
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L80;
	}

	if (m == l + 1) {
	    if (icompz > 0) {
		PLUMED_BLAS_F77_FUNC(dlaev2,DLAEV2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
		work[l] = c__;
		work[*n - 1 + l] = s;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
	    }
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l + 1] - p) / (e[l] * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l] / (g + ( (g>0) ? r__ : -r__ ) );

	s = 1.;
	c__ = 1.;
	p = 0.;

	mm1 = m - 1;
	i__1 = l;
	for (i__ = mm1; i__ >= i__1; --i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m - 1) {
		e[i__ + 1] = r__;
	    }
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = -s;
	    }
	}

	if (icompz > 0) {
	    mm = m - l + 1;
	    PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[l] = g;
	goto L40;

L80:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= i__1; --m) {
		d__2 = std::abs(e[m - 1]);
		tst = d__2 * d__2;
		if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m- 1]) + safmin) {
		    goto L110;
		}
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L130;
	}
	if (m == l - 1) {
	    if (icompz > 0) {
		PLUMED_BLAS_F77_FUNC(dlaev2,DLAEV2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
		work[m] = c__;
		work[*n - 1 + m] = s;
		PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz);
	    } else {
		PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
	    }
	    d__[l - 1] = rt1;
	    d__[l] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l - 1] / (g + ( (g>0) ? r__ : -r__ ));

	s = 1.;
	c__ = 1.;
	p = 0.;

	lm1 = l - 1;
	i__1 = lm1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    PLUMED_BLAS_F77_FUNC(dlartg,DLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m) {
		e[i__ - 1] = r__;
	    }
	    g = d__[i__] - p;
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = s;
	    }
	}

	if (icompz > 0) {
	    mm = l - m + 1;
	    PLUMED_BLAS_F77_FUNC(dlasr,DLASR)("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[lm1] = g;
	goto L90;

L130:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__])>PLUMED_GMX_DOUBLE_MIN) {
	    ++(*info);
	}
    }
    goto L190;

L160:
    if (icompz == 0) {

	PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("I", n, &d__[1], info);

    } else {

	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
	    i__ = ii - 1;
	    k = i__;
	    p = d__[i__];
	    i__2 = *n;
	    for (j = ii; j <= i__2; ++j) {
		if (d__[j] < p) {
		    k = j;
		    p = d__[j];
		}
	    }
	    if (k != i__) {
		d__[k] = d__[i__];
		d__[i__] = p;
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

L190:
    return;
}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dsterf,DSTERF)(int *n, 
	double *d__, 
	double *e, 
	int *info)
{
    int i__1;
    double d__1;

    double c__;
    int i__, l, m;
    double p, r__, s;
    int l1;
    double bb, rt1, rt2, eps, rte;
    int lsv;
    double eps2, oldc;
    int lend, jtot;
    double gamma, alpha, sigma, anorm;
      int iscale;
    double oldgam;
    double safmax;
    int lendsv;
    double ssfmin;
    int nmaxit;
    double ssfmax;
    int c__0 = 0;
    int c__1 = 1;
    double c_b32 = 1.;
    const double safmin = PLUMED_GMX_DOUBLE_MIN*(1.0+PLUMED_GMX_DOUBLE_EPS);

    --e;
    --d__;

    *info = 0;

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	return;
    }
    if (*n <= 1) {
	return;
    }

    eps = PLUMED_GMX_DOUBLE_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmax = 1. / safmin;
    ssfmax =  std::sqrt(safmax) / 3.;
    ssfmin =  std::sqrt(safmin) / eps2;

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

    l1 = 1;

L10:
    if (l1 > *n) {
      PLUMED_BLAS_F77_FUNC(dlasrt,DLASRT)("I", n, &d__[1], info);
      return;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
	if (std::abs(e[m]) <=  std::sqrt(std::abs(d__[m])) * 
		 std::sqrt(std::abs(d__[m + 1])) * eps) {
	    e[m] = 0.;
	    goto L30;
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    i__1 = lend - l + 1;
    anorm = PLUMED_BLAS_F77_FUNC(dlanst,DLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
	d__1 = e[i__];
	e[i__] = d__1 * d__1;
    }

    if (std::abs(d__[lend]) < std::abs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

L50:
	if (l != lend) {
	    i__1 = lend - 1;
	    for (m = l; m <= i__1; ++m) {
		if (std::abs(e[m]) <= eps2 * std::abs(d__[m] * d__[m + 1])) {
		    goto L70;
		}
	    }
	}
	m = lend;

L70:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L90;
	}
	if (m == l + 1) {
	    rte =  std::sqrt(e[l]);
	    PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte =  std::sqrt(e[l]);
	sigma = (d__[l + 1] - p) / (rte * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l;
	for (i__ = m - 1; i__ >= i__1; --i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m - 1) {
		e[i__ + 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__ + 1] = oldgam + (alpha - gamma);
	    if (std::abs(c__)>PLUMED_GMX_DOUBLE_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l] = s * p;
	d__[l] = sigma + gamma;
	goto L50;

L90:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

L100:
	i__1 = lend + 1;
	for (m = l; m >= i__1; --m) {
	    if (std::abs(e[m - 1]) <= eps2 * std::abs(d__[m] * d__[m - 1])) {
		goto L120;
	    }
	}
	m = lend;

L120:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L140;
	}

	if (m == l - 1) {
	    rte =  std::sqrt(e[l - 1]);
	    PLUMED_BLAS_F77_FUNC(dlae2,DLAE2)(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l - 1] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte =  std::sqrt(e[l - 1]);
	sigma = (d__[l - 1] - p) / (rte * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(dlapy2,DLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l - 1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m) {
		e[i__ - 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__ + 1];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__] = oldgam + (alpha - gamma);
	    if (std::abs(c__)>PLUMED_GMX_DOUBLE_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l - 1] = s * p;
	d__[l] = sigma + gamma;
	goto L100;

L140:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__])>PLUMED_GMX_DOUBLE_MIN) {
	    ++(*info);
	}
    }
    return;
}


}
}
#include "lapack.h"


/* Normally, DSTEVR is the LAPACK wrapper which calls one
 * of the eigenvalue methods. However, our code includes a
 * version of DSTEGR which is never than LAPACK 3.0 and can
 * handle requests for a subset of eigenvalues/vectors too,
 * and it should not need to call DSTEIN.
 * Just in case somebody has a faster version in their lapack
 * library we still call the driver routine, but in our own
 * case this is just a wrapper to dstegr.
 */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dstevr,DSTEVR)(const char *jobz, 
	const char *range,
	int *n,
	double *d,
	double *e,
	double *vl, 
	double *vu,
	int *il, 
	int *iu, 
	double *abstol,
	int *m,
	double *w, 
	double *z,
	int *ldz,
	int *isuppz, 
	double *work, 
	int *lwork, 
	int *iwork,
	int *liwork, 
	int *info)
{
  PLUMED_BLAS_F77_FUNC(dstegr,DSTEGR)(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w,
	  z, ldz, isuppz, work, lwork, iwork, liwork, info);
  

    return;

}


}
}
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dsyevr,DSYEVR)(const char *jobz, const char *range, const char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    int c__1 = 1;
    int i__, j, nb, jj;
    double eps, tmp1;
    int indd, inde;
    double anrm;
    int imax;
    double rmin, rmax;
    int itmp1, inddd, indee;
    double sigma;
    int iinfo;
    int indwk;
    int lwmin;
    int lower, wantz;
    int alleig, indeig;
    int iscale, indibl, indifl;
    int valeig;
    double safmin,minval;
    double bignum;
    int indtau;
    int indwkn;
    int liwmin;
    int llwrkn, llwork;
    double smlnum;
    int lwkopt;
    int lquery;
    
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    lower = (*uplo=='L' || *uplo=='l');
    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    indibl = 0;
    lquery = *lwork == -1 || *liwork == -1;

    i__1 = 1;
    i__2 = *n * 26;

    if(*n>0) 
      lwmin = *n * 26;
    else
      lwmin = 1;

    if(*n>0) 
      liwmin = *n * 10;
    else
      liwmin = 1;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || (*uplo=='U' || *uplo=='u'))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < ((*n>1) ? *n : 1) ) {
	*info = -6;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -8;
	    }
	} else if (indeig) {
	  if (*il < 1 || *il > ((*n>1) ? *n : 1)) {
		*info = -9;
	    } else if (*iu < ((*n<*il) ? *n : *il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
      if (*ldz < 1 || (wantz && *ldz < *n)) {
	    *info = -15;
	} else if (*lwork < lwmin && ! lquery) {
	    *info = -18;
	} else if (*liwork < liwmin && ! lquery) {
	    *info = -20;
	}
    }

    if (*info == 0) {
      nb = 32;
      /* Computing MAX */
      i__1 = (nb + 1) * *n;
      lwkopt = (i__1>lwmin) ? i__1 : lwmin;
      work[1] = (double) lwkopt;
      iwork[1] = liwmin;
    } else 
      return;

    if (lquery) 
	return;
    
    *m = 0;
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    if (*n == 1) {
	work[1] = 7.;
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = a[a_dim1 + 1];
	} else {
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
		*m = 1;
		w[1] = a[a_dim1 + 1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }
    minval = PLUMED_GMX_DOUBLE_MIN;
    safmin = minval*(1.0+PLUMED_GMX_DOUBLE_EPS);
    eps = PLUMED_GMX_DOUBLE_EPS;

    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin =  std::sqrt(smlnum);

    d__1 =  std::sqrt(bignum), d__2 = 1. / std::sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;

    iscale = 0;
    anrm = PLUMED_BLAS_F77_FUNC(dlansy,DLANSY)("M", uplo, n, &a[a_offset], lda, &work[1]);
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm; 
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&j, &sigma, &a[j * a_dim1 + 1], &c__1);

	    }
	}
    }

    indtau = 1;
    inde = indtau + *n;
    indd = inde + *n;
    indee = indd + *n;
    inddd = indee + *n;
    indifl = inddd + *n;
    indwk = indifl + *n;
    llwork = *lwork - indwk + 1;
    PLUMED_BLAS_F77_FUNC(dsytrd,DSYTRD)(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwk], &llwork, &iinfo);

    i__1 = *n - 1;
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(&i__1, &work[inde], &c__1, &work[indee], &c__1);
    PLUMED_BLAS_F77_FUNC(dcopy,DCOPY)(n, &work[indd], &c__1, &work[inddd], &c__1);

    PLUMED_BLAS_F77_FUNC(dstegr,DSTEGR)(jobz, range, n, &work[inddd], &work[indee], vl, vu, il, iu, 
	    abstol, m, &w[1], &z__[z_offset], ldz, &isuppz[1], 
	    &work[indwk], lwork, &iwork[1], liwork, info);
    if (wantz && *info == 0) {
      indwkn = inde;
      llwrkn = *lwork - indwkn + 1;
      PLUMED_BLAS_F77_FUNC(dormtr,DORMTR)("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
	      , &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
    }

    if (*info != 0) 
      return;

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&imax, &d__1, &w[1], &c__1);
    }

    if (wantz) {
	i__1 = *m - 1;
	
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		PLUMED_BLAS_F77_FUNC(dswap,DSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

    work[1] = (double) lwkopt;
    iwork[1] = liwmin;
    return;

}
}
}
#include <cctype>
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dsytd2,DSYTD2)(const char *    uplo,
	int *     n,
	double *  a,
	int *     lda,
	double *  d,
	double *  e,
	double *  tau,
    int *     info)
{
  double minusone,zero;
  double taui,alpha,tmp;
  int ti1,ti2,ti3,i;
  const char ch=std::toupper(*uplo);

  zero = 0.0;
  minusone = -1.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n-1;i>=1;i--) {

      ti1 = 1;
      PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&i,&(a[i*(*lda)+(i-1)]),&(a[i*(*lda)+0]),&ti1,&taui);
      e[i-1] = a[i*(*lda) + (i-1)];
      if(std::abs(taui)>PLUMED_GMX_DOUBLE_MIN) {
	a[i*(*lda)+(i-1)] = 1.0;
      
	ti1 = 1;
	PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)("U",&i,&taui,a,lda,&(a[i*(*lda)+0]),&ti1,&zero,tau,&ti1);

	tmp = PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&i,tau,&ti1,&(a[i*(*lda)+0]),&ti1);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&i,&alpha,&(a[i*(*lda)+0]),&ti1,tau,&ti1);

	PLUMED_BLAS_F77_FUNC(dsyr2,DSYR2)("U",&i,&minusone,&(a[i*(*lda)+0]),&ti1,tau,&ti1,a,lda);

	a[i*(*lda)+(i-1)] = e[i-1]; 

      }
      d[i] = a[i*(*lda)+i];
      tau[i-1] = taui;
    }
    d[0] = a[0];
    
  } else {
    /* lower */

    for(i=1;i<*n;i++) {

      ti1 = *n - i;
      ti2 = ( *n < i+2) ? *n : i+2;
      ti3 = 1;
      PLUMED_BLAS_F77_FUNC(dlarfg,DLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+ti2-1]),&ti3,&taui);

      e[i-1] = a[(i-1)*(*lda) + (i)];

      if(std::abs(taui)>PLUMED_GMX_DOUBLE_MIN) {
	a[(i-1)*(*lda)+(i)] = 1.0;
      
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(dsymv,DSYMV)(uplo,&ti1,&taui,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),
	       &ti2,&zero,&(tau[i-1]),&ti2);
	
	tmp = PLUMED_BLAS_F77_FUNC(ddot,DDOT)(&ti1,&(tau[i-1]),&ti2,&(a[(i-1)*(*lda)+i]),&ti2);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(daxpy,DAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2);

	PLUMED_BLAS_F77_FUNC(dsyr2,DSYR2)(uplo,&ti1,&minusone,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2,
	       &(a[(i)*(*lda)+i]),lda);

	a[(i-1)*(*lda)+(i)] = e[i-1]; 

      }
      d[i-1] = a[(i-1)*(*lda)+i-1];
      tau[i-1] = taui;
    }
    d[*n-1] = a[(*n-1)*(*lda)+(*n-1)];
 
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dsytrd,DSYTRD)(const char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, double *
	work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, j, nb, kk, nx, iws;
    int nbmin, iinfo;
    int upper;
    int ldwork, lwkopt;
    int lquery;
    double c_b22 = -1.;
    double c_b23 = 1.;


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    upper = (*uplo=='U' || *uplo=='u');
    lquery = (*lwork == -1);

    if (! upper && ! (*uplo=='L' || *uplo=='l')) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ((1>*n) ? 1 : *n)) {
	*info = -4;
    } else if (*lwork < 1 && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {

      nb = DSYTRD_BLOCKSIZE;
      lwkopt = *n * nb;
      work[1] = (double) lwkopt;
    } else
      return;

    if (lquery) 
      return;
  
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    nx = *n;
    if (nb > 1 && nb < *n) {

	nx = DSYTRD_CROSSOVER;
	if (nx < *n) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		i__1 = *lwork / ldwork;
		nb = (i__1>1) ? i__1 : 1;
		nbmin = DSYTRD_MINBLOCKSIZE;
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

	    i__3 = i__ + nb - 1;
	    PLUMED_BLAS_F77_FUNC(dlatrd,DLATRD)(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork);

	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(dsyr2k,DSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 
		    + 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j - 1 + j * a_dim1] = e[j - 1];
		d__[j] = a[j + j * a_dim1];

	    }

	}

	PLUMED_BLAS_F77_FUNC(dsytd2,DSYTD2)(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
    } else {

	i__2 = *n - nx;
	i__1 = nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {


	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(dlatrd,DLATRD)(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
		    tau[i__], &work[1], &ldwork);

	    i__3 = *n - i__ - nb + 1;
	    PLUMED_BLAS_F77_FUNC(dsyr2k,DSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + 
		    i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
		    i__ + nb + (i__ + nb) * a_dim1], lda);


	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j + 1 + j * a_dim1] = e[j];
		d__[j] = a[j + j * a_dim1];

	    }

	}


	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(dsytd2,DSYTD2)(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
		&tau[i__], &iinfo);
    }

    work[1] = (double) lwkopt;
    return;

}


}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dtrti2,DTRTI2)(const char *uplo,
	const char *diag, 
	int *n, 
	double *a,
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__2;

    int j;
    double ajj;
    int upper, nounit;
    int c__1 = 1;


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

    if (upper) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }

	    i__2 = j - 1;
	    PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1);
	    i__2 = j - 1;
	    PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
	}
    } else {

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(dtrmv,DTRMV)("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1);
		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(dscal,DSCAL)(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
	    }
	}
    }
    return;
}
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(dtrtri,DTRTRI)(const char *uplo,
	const char *diag, 
	int *n,
	double *a, 
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__3, i__4, i__5;
    int j, jb, nb, nn;
    double c_b18 = 1.;
    double c_b22 = -1.;

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
	    if (std::abs(a[*info + *info * a_dim1])<PLUMED_GMX_DOUBLE_MIN) {
		return;
	    }
	}
	*info = 0;
    }

    nb = DTRTRI_BLOCKSIZE;
    if (nb <= 1 || nb >= *n) {

	PLUMED_BLAS_F77_FUNC(dtrti2,DTRTI2)(uplo, diag, n, &a[a_offset], lda, info);
    } else {

	if (upper) {

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
		i__4 = nb, i__5 = *n - j + 1;
		jb = (i__4<i__5) ? i__4 : i__5;

		i__4 = j - 1;
		PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda);
		i__4 = j - 1;
		PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], 
			lda);

		PLUMED_BLAS_F77_FUNC(dtrti2,DTRTI2)("Upper", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	} else {

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
		i__1 = nb, i__4 = *n - j + 1;
		jb = (i__1<i__4) ? i__1 : i__4;
		if (j + jb <= *n) {

		    i__1 = *n - j - jb + 1;
		    PLUMED_BLAS_F77_FUNC(dtrmm,DTRMM)("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda);
		    i__1 = *n - j - jb + 1;
		    PLUMED_BLAS_F77_FUNC(dtrsm,DTRSM)("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda);
		}

		PLUMED_BLAS_F77_FUNC(dtrti2,DTRTI2)("Lower", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	}
    }
    return;
}


}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(ilasrt2,ILASRT2)(const char *id, 
	 int *n, 
	 int *d__, 
	 int *key, 
	 int *info)
{
    int i__1, i__2;
    int i__, j, d1, d2, d3, dir, tmp, endd;
    int stack[64], dmnmx, start;
    int tmpkey, stkpnt;

    --key;
    --d__;

    *info = 0;
    dir = -1;
    if (*id=='D' || *id=='d') 
	dir = 0;
    else if (*id=='I' || *id=='i') 
	dir = 1;
    
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }

    if (*n <= 1) {
	return;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start > 0) {

	if (dir == 0) {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			goto L30;
		    }
		}
L30:
		;
	    }

	} else {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			goto L50;
		    }
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }


    return;
}
}
}
#include <cctype>
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)(const char *uplo, 
	const char *compq, 
	int *n,
	float *d__, 
	float *e, 
	float *u, 
	int *ldu,
	float *vt, 
	int *ldvt,
	float *q,
	int *iq,
	float *work, 
	int *iwork, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    int i__, j, k;
    float p, r__;
    int z__, ic, ii, kk;
    float cs;
    int is, iu;
    float sn;
    int nm1;
    float eps;
    int ivt, difl, difr, ierr, perm, mlvl, sqre;
    int poles, iuplo, nsize, start;
    int givcol;
    int icompq;
    float orgnrm;
    int givnum, givptr, qstart, smlsiz, wstart, smlszp;
    float zero = 0.0;
    float one = 1.0;
    int c_0 = 0;
    int c_1 = 1;

    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --q;
    --iq;
    --work;
    --iwork;

    k = iu = z__ = ic = is = ivt = difl = difr = perm = 0;
    poles = givnum = givptr = givcol = 0;
    
    smlsiz = DBDSDC_SMALLSIZE;
    *info = 0;

    iuplo = (*uplo=='U' || *uplo=='u') ? 1 : 2;

    switch(*compq) {
    case 'n':
    case 'N':
      icompq = 0;
      break;
    case 'p':
    case 'P':
      icompq = 1;
      break;
    case 'i':
    case 'I':
      icompq = 2;
      break;
    default:
      return;
    }

    if (*n <= 0) 
	return;
    
    if (*n == 1) {
	if (icompq == 1) {
	  q[1] = (d__[1]>0) ? 1.0 : -1.0;
	  q[smlsiz * *n + 1] = 1.;
	} else if (icompq == 2) {
	  u[u_dim1 + 1] = (d__[1]>0) ? 1.0 : -1.0;
	  vt[vt_dim1 + 1] = 1.;
	}
	d__[1] = std::abs(d__[1]);
	return;
    }
    nm1 = *n - 1;
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &d__[1], &c_1, &q[1], &c_1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &e[1], &c_1, &q[*n + 1], &c_1);
    }
    if (iuplo == 2) {
	qstart = 5;
	wstart = (*n << 1) - 1;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (icompq == 1) {
		q[i__ + (*n << 1)] = cs;
		q[i__ + *n * 3] = sn;
	    } else if (icompq == 2) {
		work[i__] = cs;
		work[nm1 + i__] = -sn;
	    }
	}
    }
    if (icompq == 0) {
      PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U",&c_0,n,&c_0,&c_0,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
	      &u[u_offset], ldu, &u[u_offset], ldu, &work[wstart], info);
	goto L40;
    }
    if (*n <= smlsiz) {
	if (icompq == 2) {
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U",&c_0,n,n,n,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
		    &u[u_offset],ldu,&u[u_offset],ldu,&work[wstart],info);
	} else if (icompq == 1) {
	    iu = 1;
	    ivt = iu + *n;
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &q[iu + (qstart - 1) * *n], n);
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &q[ivt + (qstart - 1) * *n], n);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &c_0, n, n, n, &c_0, &d__[1], &e[1], 
		    &q[ivt + (qstart - 1) * *n], n, &q[iu + (qstart - 1) * *n], 
		    n, &q[iu + (qstart - 1) * *n], n, &work[wstart], info);
	}
	goto L40;
    }

    if (icompq == 2) {
	PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
    }

    orgnrm = PLUMED_BLAS_F77_FUNC(slanst,SLANST)("M", n, &d__[1], &e[1]);
    if ( std::abs(orgnrm)<PLUMED_GMX_FLOAT_MIN) {
	return;
    }
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c_0, &c_0, &orgnrm, &one, n, &c_1, &d__[1], n, &ierr);
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c_0, &c_0, &orgnrm, &one, &nm1, &c_1, &e[1], &nm1, &ierr);

    eps = PLUMED_GMX_FLOAT_EPS;

    mlvl = (int) (std::log((float) (*n) / (float) (smlsiz + 1)) / std::log(2.)) + 1;
    smlszp = smlsiz + 1;

    if (icompq == 1) {
	iu = 1;
	ivt = smlsiz + 1;
	difl = ivt + smlszp;
	difr = difl + mlvl;
	z__ = difr + (mlvl << 1);
	ic = z__ + mlvl;
	is = ic + 1;
	poles = is + 1;
	givnum = poles + (mlvl << 1);

	k = 1;
	givptr = 2;
	perm = 3;
	givcol = perm + mlvl;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(d__[i__]) < eps) 
	    d__[i__] = (d__[i__]>0) ? eps : -eps;
    }

    start = 1;
    sqre = 0;

    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__]) < eps || i__ == nm1) {
	    if (i__ < nm1) {
		nsize = i__ - start + 1;
	    } else if (std::abs(e[i__]) >= eps) {
		nsize = *n - start + 1;
	    } else {
		nsize = i__ - start + 1;
		if (icompq == 2) {
		    u[*n + *n * u_dim1] = (d__[*n]>0) ? 1.0 : -1.0; 
		    vt[*n + *n * vt_dim1] = 1.;
		} else if (icompq == 1) {
		    q[*n + (qstart - 1) * *n] = (d__[*n]>0) ? 1.0 : -1.0; 
		    q[*n + (smlsiz + qstart - 1) * *n] = 1.;
		}
		d__[*n] = std::abs(d__[*n]);
	    }
	    if (icompq == 2) {
		PLUMED_BLAS_F77_FUNC(slasd0,SLASD0)(&nsize, &sqre, &d__[start], &e[start], 
			&u[start + start * u_dim1], ldu, 
			&vt[start + start * vt_dim1], 
			ldvt, &smlsiz, &iwork[1], &work[wstart], info);
	    } else {
		PLUMED_BLAS_F77_FUNC(slasda,SLASDA)(&icompq, &smlsiz, &nsize, &sqre, &d__[start], 
			&e[start], &q[start + (iu + qstart - 2) * *n], n, 
			&q[start + (ivt + qstart - 2) * *n], &iq[start + k * *n],
			&q[start + (difl + qstart - 2) * *n], 
			&q[start + (difr + qstart - 2) * *n], 
			&q[start + (z__ + qstart - 2) * *n], 
			&q[start + (poles + qstart - 2) * *n], 
			&iq[start + givptr * *n], &iq[start + givcol * *n], n, 
			&iq[start + perm * *n], 
			&q[start + (givnum + qstart - 2) * *n], 
			&q[start + (ic + qstart - 2) * *n], 
			&q[start + (is + qstart - 2) * *n], &work[wstart], 
			&iwork[1], info);
		if (*info != 0) {
		    return;
		}
	    }
	    start = i__ + 1;
	}
    }
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c_0, &c_0, &one, &orgnrm, n, &c_1, &d__[1], n, &ierr);
L40:
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	kk = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] > p) {
		kk = j;
		p = d__[j];
	    }
	}
	if (kk != i__) {
	    d__[kk] = d__[i__];
	    d__[i__] = p;
	    if (icompq == 1) {
		iq[i__] = kk;
	    } else if (icompq == 2) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &u[i__ * u_dim1 + 1],&c_1,&u[kk*u_dim1+1],&c_1);
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &vt[i__ + vt_dim1], ldvt, &vt[kk + vt_dim1], ldvt);
	    }
	} else if (icompq == 1) {
	    iq[i__] = i__;
	}
    }
    if (icompq == 1) {
	if (iuplo == 1) {
	    iq[*n] = 1;
	} else {
	    iq[*n] = 0;
	}
    }
    if (iuplo == 2 && icompq == 2) {
	PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "B", n, n, &work[1], &work[*n], &u[u_offset], ldu);
    }

    return;
}
}
}
#include <cctype>
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sbdsqr,SBDSQR)(const char *uplo,
                        int *n,
                        int *ncvt,
                        int *nru, 
                        int *ncc, 
                        float *d__,
                        float *e,
                        float *vt, 
                        int *ldvt,
                        float *u, 
                        int *ldu,
                        float *c__, 
                        int *ldc,
                        float *work,
                        int *info)
{
    const char xuplo = std::toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    float r__1, r__2, r__3, r__4;
    float c_b15 = -.125;

    int c__1 = 1;
    float c_b49 = 1.f;
    float c_b72 = -1.f;

    float f, g, h__;
    int i__, j, m;
    float r__, cs;
    int ll;
    float sn, mu;
    int nm1, nm12, nm13, lll;
    float eps, sll, tol, abse;
    int idir;
    float abss;
    int oldm;
    float cosl;
    int isub, iter;
    float unfl, sinl, cosr, smin, smax, sinr;
    float oldcs;
    int oldll;
    float shift, sigmn, oldsn = 0.;
    int maxit;
    float sminl;
    float sigmx;
    int lower;
    float sminoa;
    float thresh;
    int rotate;
    float tolmul;
    int itmp1,itmp2;
    
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    
    lower = (xuplo == 'L');
    if ( (xuplo!='U') && !lower) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ncvt < 0) {
	*info = -3;
    } else if (*nru < 0) {
	*info = -4;
    } else if (*ncc < 0) {
	*info = -5;
    } else if ( ((*ncvt == 0) && (*ldvt < 1)) || ((*ncvt > 0) && (*ldvt < itmp1)) ) {
        *info = -9;
    } else if (*ldu < itmp2) {
        *info = -11;
    } else if ( ((*ncc == 0) && (*ldc < 1)) || ((*ncc > 0) && (*ldc < itmp1))) {
        *info = -13;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	goto L160;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

    if (! rotate) {
	PLUMED_BLAS_F77_FUNC(slasq1,SLASQ1)(n, &d__[1], &e[1], &work[1], info);
	return;
    }

    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;

    eps = PLUMED_GMX_FLOAT_EPS;
    unfl = PLUMED_GMX_FLOAT_MIN/PLUMED_GMX_FLOAT_EPS;

    if (lower) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    work[i__] = cs;
	    work[nm1 + i__] = sn;
	}

	if (*nru > 0) {
	    PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", nru, n, &work[1], &work[*n], &u[u_offset], 
		    ldu);
	}
	if (*ncc > 0) {
	    PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", n, ncc, &work[1], &work[*n], &c__[c_offset],
		     ldc);
	}
    }

    r__3 = 100.f, r__4 = std::pow(static_cast<float>(PLUMED_GMX_FLOAT_EPS),c_b15);
    r__1 = 10.f, r__2 = (r__3<r__4) ? r__3 : r__4;
    tolmul = (r__1>r__2) ? r__1 : r__2;
    tol = tolmul * eps;
    smax = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = d__[i__], std::abs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = e[i__], std::abs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    sminl = 0.f;
    if (tol >= 0.f) {
	sminoa = std::abs(d__[1]);
	if (sminoa == 0.f) {
	    goto L50;
	}
	mu = sminoa;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mu = (r__2 = d__[i__], std::abs(r__2)) * (mu / (mu + (r__1 = e[i__ - 
		    1], std::abs(r__1))));
	    sminoa = (sminoa<mu) ? sminoa : mu;
	    if (sminoa == 0.f) {
		goto L50;
	    }
	}
L50:
	sminoa /=  std::sqrt((float) (*n));
	r__1 = tol * sminoa, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    } else {
	r__1 = std::abs(tol) * smax, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    }
    maxit = *n * 6 * *n;
    iter = 0;
    oldll = -1;
    oldm = -1;
    m = *n;

L60:

    if (m <= 1) {
	goto L160;
    }
    if (iter > maxit) {
	goto L200;
    }

    if (tol < 0.f && (r__1 = d__[m], std::abs(r__1)) <= thresh) {
	d__[m] = 0.f;
    }
    smax = (r__1 = d__[m], std::abs(r__1));
    smin = smax;
    i__1 = m - 1;
    for (lll = 1; lll <= i__1; ++lll) {
	ll = m - lll;
	abss = (r__1 = d__[ll], std::abs(r__1));
	abse = (r__1 = e[ll], std::abs(r__1));
	if (tol < 0.f && abss <= thresh) {
	    d__[ll] = 0.f;
	}
	if (abse <= thresh) {
	    goto L80;
	}
	smin = (smin<abss) ? smin : abss;
	r__1 = (smax>abss) ? smax : abss;
	smax = (r__1>abse) ? r__1 : abse;
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.f;
    if (ll == m - 1) {
	--m;
	goto L60;
    }
L90:
    ++ll;
    if (ll == m - 1) {
	PLUMED_BLAS_F77_FUNC(slasv2,SLASV2)(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
	d__[m - 1] = sigmx;
	e[m - 1] = 0.f;
	d__[m] = sigmn;
	if (*ncvt > 0) {
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
	}
	if (*nru > 0) {
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
	}
	if (*ncc > 0) {
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
	}
	m += -2;
	goto L60;
    }
    if (ll > oldm || m < oldll) {
	if ((r__1 = d__[ll], std::abs(r__1)) >= (r__2 = d__[m], std::abs(r__2))) {
	    idir = 1;
	} else {
	    idir = 2;
	}
    }
    if (idir == 1) {

        if( (std::abs(e[m-1]) <= std::abs(tol) * std::abs(d__[m])) ||
            (tol<0.0 && std::abs(e[m-1])<=thresh)) {
	    e[m - 1] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[ll], std::abs(r__1));
	    sminl = mu;
	    i__1 = m - 1;
	    for (lll = ll; lll <= i__1; ++lll) {
		if ((r__1 = e[lll], std::abs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		mu = (r__2 = d__[lll + 1], std::abs(r__2)) * (mu / (mu + (r__1 = 
			e[lll], std::abs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    } else {
        if( (std::abs(e[ll]) <= std::abs(tol)*std::abs(d__[ll])) ||
            (tol<0.0 && std::abs(e[ll])<=thresh)) {
	    e[ll] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[m], std::abs(r__1));
	    sminl = mu;
	    i__1 = ll;
	    for (lll = m - 1; lll >= i__1; --lll) {
		if ((r__1 = e[lll], std::abs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		mu = (r__2 = d__[lll], std::abs(r__2)) * (mu / (mu + (r__1 = e[
			lll], std::abs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    }
    oldll = ll;
    oldm = m;

    r__1 = eps, r__2 = tol * .01f;
    if (tol >= 0.f && *n * tol * (sminl / smax) <= ((r__1>r__2) ? r__1 : r__2)) {
	shift = 0.f;
    } else {
	if (idir == 1) {
	    sll = (r__1 = d__[ll], std::abs(r__1));
	    PLUMED_BLAS_F77_FUNC(slas2,SLAS2)(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
	} else {
	    sll = (r__1 = d__[m], std::abs(r__1));
	    PLUMED_BLAS_F77_FUNC(slas2,SLAS2)(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
	}
	if (sll > 0.f) {
	    r__1 = shift / sll;
	    if (r__1 * r__1 < eps) {
		shift = 0.f;
	    }
	}
    }
    iter = iter + m - ll;
    if (shift == 0.f) {
	if (idir == 1) {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		r__1 = d__[i__] * cs;
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&r__1, &e[i__], &cs, &sn, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ + 1] * sn;
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll + 1] = cs;
		work[i__ - ll + 1 + nm1] = sn;
		work[i__ - ll + 1 + nm12] = oldcs;
		work[i__ - ll + 1 + nm13] = oldsn;
	    }
	    h__ = d__[m] * cs;
	    d__[m] = h__ * oldcs;
	    e[m - 1] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], std::abs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		r__1 = d__[i__] * cs;
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&r__1, &e[i__ - 1], &cs, &sn, &r__);
		if (i__ < m) {
		    e[i__] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ - 1] * sn;
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll] = cs;
		work[i__ - ll + nm1] = -sn;
		work[i__ - ll + nm12] = oldcs;
		work[i__ - ll + nm13] = -oldsn;
	    }
	    h__ = d__[ll] * cs;
	    d__[ll] = h__ * oldcs;
	    e[ll] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[ll], std::abs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	}
    } else {

	if (idir == 1) {
	    f = ((r__1 = d__[ll], std::abs(r__1)) - shift) * ( ((d__[ll] > 0) ? c_b49 : -c_b49) + shift / d__[ll]);
	    g = e[ll];
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__];
		e[i__] = cosr * e[i__] - sinr * d__[i__];
		g = sinr * d__[i__ + 1];
		d__[i__ + 1] = cosr * d__[i__ + 1];
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__] + sinl * d__[i__ + 1];
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
		if (i__ < m - 1) {
		    g = sinl * e[i__ + 1];
		    e[i__ + 1] = cosl * e[i__ + 1];
		}
		work[i__ - ll + 1] = cosr;
		work[i__ - ll + 1 + nm1] = sinr;
		work[i__ - ll + 1 + nm12] = cosl;
		work[i__ - ll + 1 + nm13] = sinl;
	    }
	    e[m - 1] = f;

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], std::abs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {

	    f = ((r__1 = d__[m], std::abs(r__1)) - shift) * ( ((d__[m] > 0) ? c_b49 : -c_b49) + shift / d__[m]);
	    g = e[m - 1];
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ < m) {
		    e[i__] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__ - 1];
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
		g = sinr * d__[i__ - 1];
		d__[i__ - 1] = cosr * d__[i__ - 1];
		PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
		if (i__ > ll + 1) {
		    g = sinl * e[i__ - 2];
		    e[i__ - 2] = cosl * e[i__ - 2];
		}
		work[i__ - ll] = cosr;
		work[i__ - ll + nm1] = -sinr;
		work[i__ - ll + nm12] = cosl;
		work[i__ - ll + nm13] = -sinl;
	    }
	    e[ll] = f;

	    if ((r__1 = e[ll], std::abs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	}
    }

    goto L60;

L160:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] < 0.f) {
	    d__[i__] = -d__[i__];

	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
	    }
	}
    }

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = 1;
	smin = d__[1];
	i__2 = *n + 1 - i__;
	for (j = 2; j <= i__2; ++j) {
	    if (d__[j] <= smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != *n + 1 - i__) {
	    d__[isub] = d__[*n + 1 - i__];
	    d__[*n + 1 - i__] = smin;
	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
	    }
	    if (*ncc > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
	    }
	}
    }
    goto L220;

L200:
    *info = 0;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    ++(*info);
	}
    }
L220:
    return;

}


}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgebd2,SGEBD2)(int *m,
	int *n,
	float *a,
	int *lda,
	float *d,
	float *e,
	float *tauq,
	float *taup,
	float *work,
	int *info)
{

  int i,i1,i2,i3;

    *info = 0;

  if(*m>=*n) {
    /* reduce to upper bidiag. form */
    for(i=0;i<*n;i++) {
      i1 = *m - i;
      i2 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      i3 = 1;
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      PLUMED_BLAS_F77_FUNC(slarf,SLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*n-1)) {

	i1 = *n - i -1;
	i2 = ( (i+2) < (*n-1)) ? (i+2) : (*n-1); 
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i1,&(a[(i+1)*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));

	e[i] = a[(i+1)*(*lda)+i];
	a[(i+1)*(*lda)+i] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	PLUMED_BLAS_F77_FUNC(slarf,SLARF)("R",&i1,&i2,&(a[(i+1)*(*lda)+i]),lda,&(taup[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i+1)*(*lda)+i] = e[i];
      } else
	taup[i] = 0.0;
    }
  } else {
    /* reduce to lower bidiag. form */
    for(i=0;i<*m;i++) {
      i1 = *n - i;
      i2 = ( (i+1) < (*n-1)) ? (i+1) : (*n-1);
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;

      i2 = *m - i - 1;
      i3 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      PLUMED_BLAS_F77_FUNC(slarf,SLARF)("R",&i2,&i1,&(a[i*(*lda)+i]),lda,&(taup[i]),&(a[(i)*(*lda)+i3]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*m-1)) {

	i1 = *m - i - 1;
	i2 = ( (i+2) < (*m-1)) ? (i+2) : (*m-1);
	i3 = 1;
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i1,&(a[(i)*(*lda)+i+1]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));

	e[i] = a[(i)*(*lda)+i+1];
	a[(i)*(*lda)+i+1] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	i3 = 1;
	PLUMED_BLAS_F77_FUNC(slarf,SLARF)("L",&i1,&i2,&(a[(i)*(*lda)+i+1]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i)*(*lda)+i+1] = e[i];
      } else
	tauq[i] = 0.0;
    }
  }
  return;
}
}
}
#include "lapack.h"
#include "blas/blas.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(int *m, 
	int *n, 
	float *a, 
	int *lda, 
	float *d__, 
	float *e,
	float *tauq, 
	float *taup,
	float *work, 
	int *lwork,
	int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i_1, i_2, i_3, i_4;

    /* Local variables */
    int i_, j, nx,nb;
    float ws;
    int nbmin, iinfo, minmn;
    int ldwrkx, ldwrky;
    float one = 1.0;
    float minusone = -1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;

    nb = DGEBRD_BLOCKSIZE;
    *info = 0;
    if (*lwork==-1) {
      work[1] = (float) ( (*m + *n) * nb);
      return;
    }
    minmn = (*m < *n) ? *m : *n;
    if (minmn == 0) {
      work[1] = 1.;
      return;
    }

    ws = (*m > *n) ? *m : *n;
    ldwrkx = *m;
    ldwrky = *n;

    if (nb > 1 && nb < minmn) {
	nx = DGEBRD_CROSSOVER;
	if (nx < minmn) {
	    ws = (float) ((*m + *n) * nb);
	    if ((float) (*lwork) < ws) {
	      nbmin = DGEBRD_MINBLOCKSIZE;
		if (*lwork >= (*m + *n) * nbmin) {
		    nb = *lwork / (*m + *n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }

    i_1 = minmn - nx;
    i_2 = nb;
    for (i_ = 1; i_2 < 0 ? i_ >= i_1 : i_ <= i_1; i_ += i_2) {

	i_3 = *m - i_ + 1;
	i_4 = *n - i_ + 1;
	PLUMED_BLAS_F77_FUNC(slabrd,SLABRD)(&i_3, &i_4, &nb, &a[i_ + i_ * a_dim1], lda, &d__[i_], 
		&e[i_], &tauq[i_], &taup[i_], &work[1], &ldwrkx, 
		&work[ldwrkx * nb + 1], &ldwrky);

	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "T", &i_3, &i_4, &nb, &minusone, 
	       &a[i_ + nb + i_ * a_dim1], lda, &work[ldwrkx * nb + nb + 1],
	       &ldwrky, &one, &a[i_ + nb + (i_ + nb) * a_dim1], lda);
	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", &i_3, &i_4, &nb, &minusone, &work[nb + 1], &ldwrkx,
	       &a[i_ + (i_ + nb) * a_dim1], lda, &one, 
	       &a[i_ + nb + (i_ + nb) * a_dim1], lda);

	if (*m >= *n) {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + (j + 1) * a_dim1] = e[j];
	    }
	} else {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + 1 + j * a_dim1] = e[j];
	    }
	}
    }

    i_2 = *m - i_ + 1;
    i_1 = *n - i_ + 1;
    PLUMED_BLAS_F77_FUNC(sgebd2,SGEBD2)(&i_2, &i_1, &a[i_ + i_ * a_dim1], lda, &d__[i_], &e[i_], &
	    tauq[i_], &taup[i_], &work[1], &iinfo);
    work[1] = ws;
    return;

}
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sgelq2,SGELQ2)(int *m, 
                        int *n, 
                        float *a,
                        int *lda, 
                        float *tau, 
                        float *work, 
                        int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    int i__, k;
    float aii;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    
    i__4 = (*m > 1) ? *m : 1;
    
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < i__4) {
	*info = -4;
    }
    if (*info != 0) {
	return;
    }

    
    k = (*m < *n ) ? *m : *n;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__ + 1;
	i__3 = i__ + 1;
    i__4 = (i__3 < *n) ? i__3 : *n;
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + i__4 * a_dim1],
            lda, &tau[i__]);
	if (i__ < *m) {
	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.f;
	    i__2 = *m - i__;
	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(slarf,SLARF)("R", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, 
               &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    a[i__ + i__ * a_dim1] = aii;
	}
    }
    return;
}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgelqf,SGELQF)(int *m,
	int *n, 
	float *a, 
	int *lda, 
	float *tau,
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, k, ib, nb, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    nb = DGELQF_BLOCKSIZE;
    lwkopt = *m * nb;
    work[1] = (float) lwkopt;

    if (*lwork==-1) {
	return;
    }

    k =(*m < *n) ? *m : *n;
    if (k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {
	nx = DGELQF_CROSSOVER;
	if (nx < k) {
	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DGELQF_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__3 = k - i__ + 1;
	    ib = (i__3 < nb) ? i__3 : nb;

	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(sgelq2,SGELQ2)(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
	    if (i__ + ib <= *m) {

		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__3 = *m - i__ - ib + 1;
		i__4 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }
	}
    } else {
	i__ = 1;
    }

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(sgelq2,SGELQ2)(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
    }

    work[1] = (float) iws;
    return;

}
}
}
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgeqr2,SGEQR2)(int *m,
	int *n,
	float *a,
	int *lda,
	float *tau,
	float *work,
	int *info)
{
  int k = (*m < *n) ? *m : *n;
  int i,i1,i2,i3;
  float aii;

  *info = 0;
  
  for(i=0;i<k;i++) {
    i1 = *m - i;
    i2 = ( (i+1) < (*m-1) ) ? (i+1) : (*m-1);
    i3 = 1;
    PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tau[i]));
    if(i<(*n-1)) {
      aii = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      PLUMED_BLAS_F77_FUNC(slarf,SLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tau[i]),
	     &(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = aii;
    }
  }
  return;
}
}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sgeqrf,SGEQRF)(int *m, 
	int *n, 
	float *a, 
	int *lda, 
	float *tau,
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, k, ib, nb, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    nb = DGEQRF_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (float) lwkopt;
        if (*lwork==-1)
	return;
    

    k = (*m < *n) ? *m : *n;
    if (k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {
	
      nx = DGEQRF_CROSSOVER;
	if (nx < k) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DGEQRF_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {
	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

	    i__3 = k - i__ + 1;
	    ib = (i__3 < nb) ? i__3 : nb;

	    i__3 = *m - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(sgeqr2,SGEQR2)(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
	    if (i__ + ib <= *n) {

		i__3 = *m - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork);
	    }
	}
    } else {
	i__ = 1;
    }

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(sgeqr2,SGEQR2)(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
    }

    work[1] = (float) iws;
    return;

} 

}
}
#include <cmath>
#include "real.h"


#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgesdd,SGESDD)(const char *jobz, 
                        int *m, 
                        int *n, 
                        float *a, 
                        int *lda, 
                        float *s,
                        float *u, 
                        int *ldu, 
                        float *vt, 
                        int *ldvt, 
                        float *work,
                        int *lwork, 
                        int *iwork, 
                        int *info)
{
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int ie, iu;
    float dum[1], eps;
    int ivt, iscl;
    float anrm;
    int idum[1], ierr, itau;
    int minmn, wrkbl, itaup, itauq, mnthr;
    int nwork;
    int wntqn;
    int bdspac;
    float bignum;
    int ldwrku, maxwrk, ldwkvt;
    float smlnum,minval, safemin;
    int lquery;
    int c__0 = 0;
    int c__1 = 1;
    float zero = 0.0;
    float one = 1.0;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;

    *info = 0;
    minmn = (*m < *n) ? *m : *n;
    mnthr = (int) (minmn * 11. / 6.);
    wntqn = (*jobz=='o' || *jobz=='O');

    maxwrk = 1;
    lquery = *lwork == -1;

    if (*info == 0 && *m > 0 && *n > 0) {
	if (*m >= *n) {

	    if (wntqn) {
		bdspac = *n * 7;
	    } else {
		bdspac = *n * 3 * *n + (*n << 2);
	    }
	    if (*m >= mnthr) {
		if (wntqn) {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = bdspac + *n;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = *n + (*m << 5);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *n * *n;
		}
	    } else {

		wrkbl = *n * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	} else {

	    if (wntqn) {
		bdspac = *m * 7;
	    } else {
		bdspac = *m * 3 * *m + (*m*4);
	    }
	    if (*n >= mnthr) {
		if (wntqn) {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = bdspac + *m;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = *m + (*n*32);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;

		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *m * *m;
		}
	    } else {
		wrkbl = *m * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	}
	work[1] = (float) maxwrk;
    }
    
    if( lquery != 0)
    {
        return;
    }
    
    if (*m == 0 || *n == 0) {
	if (*lwork >= 1) {
	    work[1] = 1.;
	}
	return;
    }
    eps = PLUMED_GMX_FLOAT_EPS;
    minval = PLUMED_GMX_FLOAT_MIN;
    safemin = minval / eps;
    smlnum =  std::sqrt(safemin) / eps;


    bignum = 1. / smlnum;


    anrm = PLUMED_BLAS_F77_FUNC(slange,SLANGE)("M", m, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G",&c__0,&c__0,&anrm,&smlnum,m,n,&a[a_offset],lda,&ierr);
    } else if (anrm > bignum) {
	iscl = 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G",&c__0,&c__0,&anrm,&bignum,m,n,&a[a_offset],lda,&ierr);
    }

    if (*m >= *n) {
	if (*m >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgeqrf,SGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = 1;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *n;

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {
		iu = 1;

		ldwrku = *n;
		itau = iu + ldwrku * *n;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgeqrf,SGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("L", m, n, &a[a_offset], lda, &u[u_offset], ldu);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sorgqr,SORGQR)(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = itau;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);

		PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", m, n, n, &one, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &zero, &a[a_offset], lda);

		PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("F", m, n, &a[a_offset], lda, &u[u_offset], ldu);

	    }

	} else {
	    ie = 1;
	    itauq = ie + *n;
	    itaup = itauq + *n;
	    nwork = itaup + *n;

	    i__1 = *lwork - nwork + 1;
	    PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {

		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("F", m, m, &zero, &zero, &u[u_offset], ldu);
		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *m - *n;
		i__2 = *m - *n;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("F", &i__1, &i__2, &zero, &one, &u[*n + 1 + (*n + 
			1) * u_dim1], ldu);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset],ldvt,&work[nwork],&i__1,&ierr);
	    }

	}

    } else {

	if (*n >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgelqf,SGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = 1;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *m;

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {

		ivt = 1;

		ldwkvt = *m;
		itau = ivt + ldwkvt * *m;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgelqf,SGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sorglq,SORGLQ)(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = itau;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr);

		PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", m, n, m, &one, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &zero, &a[a_offset], lda);

		PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

	    }

	} else {

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    nwork = itaup + *m;

	    i__1 = *lwork - nwork + 1;
	    PLUMED_BLAS_F77_FUNC(sgebrd,SGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("F", n, n, &zero, &zero, &vt[vt_offset], ldvt);
		PLUMED_BLAS_F77_FUNC(sbdsdc,SBDSDC)("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *n - *m;
		i__2 = *n - *m;
		PLUMED_BLAS_F77_FUNC(slaset,SLASET)("F", &i__1, &i__2, &zero, &one, &vt[*m + 1 + (*m + 
			1) * vt_dim1], ldvt);

		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);
	    }

	}

    }

    if (iscl == 1) {
	if (anrm > bignum) {
	    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
	if (anrm < smlnum) {
	    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
    }

    work[1] = (float) maxwrk;

    return;

}


}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)(int *m,
	int *n,
	float *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int j,jp,k,t1,t2,t3;
  float minusone;
  float tmp;

  minusone = -1.0;

  if(*m<=0 || *n<=0)
    return;

  k = (*m < *n) ? *m : *n;
  for(j=1;j<=k;j++) {
    t1 = *m-j+1;
    t2 = 1;
    jp = j - 1 + PLUMED_BLAS_F77_FUNC(isamax,ISAMAX)(&t1,&(a[(j-1)*(*lda)+(j-1)]),&t2);
    ipiv[j-1] = jp;
    if( std::abs(a[(j-1)*(*lda)+(jp-1)])>PLUMED_GMX_FLOAT_MIN ) {
      if(jp != j)
	PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n,&(a[ j-1 ]),lda,&(a[ jp-1 ]),lda);
      
      if(j<*m) {
	t1 = *m-j;
	t2 = 1;
	tmp = 1.0/a[(j-1)*(*lda)+(j-1)];
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&t1,&tmp,&(a[(j-1)*(*lda)+(j)]),&t2);
      }
    } else {
      *info = j;
    }

    if(j<k) {
      t1 = *m-j;
      t2 = *n-j;
      t3 = 1;
      PLUMED_BLAS_F77_FUNC(sger,SGER)(&t1,&t2,&minusone,&(a[(j-1)*(*lda)+(j)]),&t3,
	    &(a[(j)*(*lda)+(j-1)]),lda, &(a[(j)*(*lda)+(j)]),lda);
    }
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgetrf,SGETRF)(int *m,
	int *n,
	float *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int mindim,jb;
  int i,j,k,l;
  int iinfo;
  float minusone = -1.0;
  float one = 1.0;

  if(*m<=0 || *n<=0)
    return;

  *info = 0;

  mindim = (*m < *n) ? *m : *n;

  if(DGETRF_BLOCKSIZE>=mindim) {

    /* unblocked code */
    PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)(m,n,a,lda,ipiv,info);

  } else {

    /* blocked case */

    for(j=1;j<=mindim;j+=DGETRF_BLOCKSIZE) {
      jb = ( DGETRF_BLOCKSIZE < (mindim-j+1)) ? DGETRF_BLOCKSIZE : (mindim-j+1);
      /* factor diag. and subdiag blocks and test for singularity */
      k = *m-j+1;
      PLUMED_BLAS_F77_FUNC(sgetf2,SGETF2)(&k,&jb,&(a[(j-1)*(*lda)+(j-1)]),lda,&(ipiv[j-1]),&iinfo);
      
      if(*info==0 && iinfo>0)
	*info = iinfo + j - 1;

      /* adjust pivot indices */
      k = (*m < (j+jb-1)) ? *m : (j+jb-1);
      for(i=j;i<=k;i++)
	ipiv[i-1] += j - 1;

      /* Apply to columns 1 throughj j-1 */
      k = j - 1;
      i = j + jb - 1;
      l = 1;
      PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(&k,a,lda,&j,&i,ipiv,&l);
      if((j+jb)<=*n) {
	/* Apply to cols. j+jb through n */
	k = *n-j-jb+1;
	i = j+jb-1;
	l = 1;
	PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(&k,&(a[(j+jb-1)*(*lda)+0]),lda,&j,&i,ipiv,&l);
	/* Compute block row of U */
	k = *n-j-jb+1;
	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left","Lower","No transpose","Unit",&jb,&k,&one,
	       &(a[(j-1)*(*lda)+(j-1)]),lda,&(a[(j+jb-1)*(*lda)+(j-1)]),lda);

	if((j+jb)<=*m) {
	  /* Update trailing submatrix */
	  k = *m-j-jb+1;
	  i = *n-j-jb+1;
	  PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose","No transpose",&k,&i,&jb,&minusone,
		 &(a[(j-1)*(*lda)+(j+jb-1)]),lda,
		 &(a[(j+jb-1)*(*lda)+(j-1)]),lda,&one,
		 &(a[(j+jb-1)*(*lda)+(j+jb-1)]),lda);
	}

      }
    }
  }
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgetri,SGETRI)(int *n, 
	float *a, 
	int *lda, 
	int *ipiv, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, jb, nb, jj, jp, nn, iws;
    int nbmin;
    int ldwork;
    int lwkopt;
    int c__1 = 1;
    float c_b20 = -1.;
    float c_b22 = 1.;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;

    *info = 0;
    nb = DGETRI_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (float) lwkopt;

    if (*n < 0) {
	*info = -1;
    } else if (*lda < (*n)) {
	*info = -3;
    } else if (*lwork < (*n) && *lwork!=-1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (*lwork == -1) {
	return;
    }

    if (*n == 0) {
	return;
    }

    PLUMED_BLAS_F77_FUNC(strtri,STRTRI)("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if (*info > 0) {
	return;
    }

    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
	i__1 = ldwork * nb;
	iws = (i__1>1) ? i__1 : 1;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DGETRI_MINBLOCKSIZE;
	}
    } else {
	iws = *n;
    }

    if (nb < nbmin || nb >= *n) {

	for (j = *n; j >= 1; --j) {

	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		work[i__] = a[i__ + j * a_dim1];
		a[i__ + j * a_dim1] = 0.;
	    }

	    if (j < *n) {
		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 
			+ 1], lda, &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 
			+ 1], &c__1);
	    }
	}
    } else {

	nn = (*n - 1) / nb * nb + 1;
	i__1 = -nb;
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
	    i__2 = nb, i__3 = *n - j + 1;
	    jb = (i__2<i__3) ? i__2 : i__3;

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= i__2; ++jj) {
		i__3 = *n;
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
		    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
		    a[i__ + jj * a_dim1] = 0.;
		}
	    }

	    if (j + jb <= *n) {
		i__2 = *n - j - jb + 1;
		PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &
			ldwork, &c_b22, &a[j * a_dim1 + 1], lda);
	    }
	    PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda);
	}
    }

    for (j = *n - 1; j >= 1; --j) {
	jp = ipiv[j];
	if (jp != j) {
	    PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
	}
    }

    work[1] = (float) iws;
    return;

}


}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sgetrs,SGETRS)(const char *trans, 
	int *n, 
	int *nrhs, 
	float *a, 
	int *lda, 
	int *ipiv,
	float *b, 
	int *ldb, 
	int *info)
{
    int a_dim1, a_offset, b_dim1, b_offset;
    int notran;
    int c__1 = 1;
    int c_n1 = -1;
    float one = 1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    *info = 0;
    notran = (*trans=='N' || *trans=='n');

    if (*n <= 0 || *nrhs <= 0) 
	return;

    if (notran) {
	PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left", "Lower", "No transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
    } else {
	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
	PLUMED_BLAS_F77_FUNC(strsm,STRSM)("Left", "Lower", "Transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slabrd,SLABRD)(int *m, 
	int *n, 
	int *nb,
	float *a, 
	int *lda, 
	float *d__,
	float *e,
	float *tauq, 
	float *taup,
	float *x,
	int *ldx,
	float *y,
	int *ldy)
{
    int a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset;
    int i__1, i__2, i__3;
    float one = 1.0;
    float minusone = -1.0;
    float zero = 0.0;
    int c__1 = 1;
    int i__;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }

    if (*m >= *n) {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + a_dim1], lda,
		     &y[i__ + y_dim1], ldy, &one, &a[i__ + i__ * a_dim1], &c__1);
	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + x_dim1], ldx,
		   &a[i__*a_dim1+1],&c__1,&one,&a[i__+i__*a_dim1],&c__1);

	    i__2 = *m - i__ + 1;
	    i__3 = i__ + 1;
	    if(*m<i__3)
	      i__3 = *m;
	    PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__3 + i__ * a_dim1], 
		    &c__1, &tauq[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *n) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__ + 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + (i__ + 1) * 
			a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &
			y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + a_dim1], 
			lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &x[i__ + x_dim1], 
			ldx, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, 
			&y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &a[i__ + a_dim1], lda, &one, &a[i__ + (
			i__ + 1) * a_dim1], lda);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &one, &a[
			i__ + (i__ + 1) * a_dim1], lda);

		i__2 = *n - i__;
		i__3 = i__ + 2;
		if(*n<i__3)
		  i__3 = *n;
		PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + (i__ + 1) * a_dim1], 
			&a[i__ + i__3 * a_dim1], lda, &taup[i__]);
		e[i__] = a[i__ + (i__ + 1) * a_dim1];
		a[i__ + (i__ + 1) * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ 
			+ 1) * a_dim1], lda, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &zero, &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__, &one, &y[i__ + 1 + y_dim1], 
			ldy, &a[i__ + (i__ + 1) * a_dim1], lda, &zero, &x[
			i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			zero, &x[i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
	    }
	}
    } else {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + y_dim1], ldy,
		     &a[i__ + a_dim1], lda, &one, &a[i__ + i__ * a_dim1],lda);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[i__ * a_dim1 + 1], 
		    lda, &x[i__ + x_dim1], ldx, &one,&a[i__+i__*a_dim1],lda);

	    i__2 = *n - i__ + 1;
	    i__3 = i__ + 1;
	    if(*n<i__3)
	      i__3 = *n;
	    PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], 
		    &a[i__ + i__3 * a_dim1], lda, &taup[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose",&i__2,&i__3,&one,&a[i__+1+i__*a_dim1], 
		       lda, &a[i__ + i__ * a_dim1], lda, &zero, 
		       &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &y[i__ + y_dim1], 
			ldy, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ * 
			x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ * a_dim1 + 
			1], lda, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ *
			 x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &y[i__ + y_dim1], ldy, &one, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &a[i__ * a_dim1 + 1], &c__1, &one, &a[
			i__ + 1 + i__ * a_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ + 2;
		if(*m<i__3)
		  i__3 = *m;
		PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + 1 + i__ * a_dim1], 
			&a[i__3 + i__ * a_dim1], &c__1, &tauq[i__]);
		e[i__] = a[i__ + 1 + i__ * a_dim1];
		a[i__ + 1 + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ + 
			1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, 
			&zero, &y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__, &one, &x[i__ + 1 + x_dim1], 
			ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__, &i__2, &minusone, &a[(i__ + 1) * a_dim1 
			+ 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, &y[i__ 
			+ 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
	    }
	}
    }
    return;
} 

}
}
#include <cctype>
#include "lapack.h"

/* LAPACK */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)(const char *uplo,
	int *m,
	int *n,
	float *a,
	int *lda,
	float *b,
	int *ldb)
{
  int i,j,minjm;
  const char ch=std::toupper(*uplo);

  if(ch=='U') {
    for(j=0;j<*n;j++) {
      minjm = (j < (*m-1)) ? j : (*m-1);
      for(i=0;i<=minjm;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else if(ch=='L') {
    for(j=0;j<*n;j++) {
      for(i=j;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }    
  }
}
}
}
#include <cmath>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slae2,SLAE2)(float *a, 
       float *b,
       float *c__, 
       float *rt1, 
       float *rt2)
{
    float d__1;
    float ab, df, tb, sm, rt, adf, acmn, acmx;


    sm = *a + *c__;
    df = *a - *c__;
    adf = std::abs(df);
    tb = *b + *b;
    ab = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf *  std::sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab *  std::sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab *  std::sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
    }
    return;

}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaebz,SLAEBZ)(int *ijob,
	int *nitmax,
	int *n, 
	int *mmax,
	int *minp,
	int *nbmin,
	float *abstol, 
	float *reltol, 
	float *pivmin, 
	float *d__,
	float *e,
	float *e2, 
	int *nval,
	float *ab, 
	float *c__, 
	int *mout, 
	int *nab,
	float *work,
	int *iwork, 
	int *info)
{
    int nab_dim1, nab_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    float d__1, d__2, d__3, d__4;

    int j, kf, ji, kl, jp, jit;
    float tmp1, tmp2;
    int itmp1, itmp2, kfnew, klnew;

    nab_dim1 = *mmax;
    nab_offset = 1 + nab_dim1;
    nab -= nab_offset;
    ab_dim1 = *mmax;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --d__;
    --e;
    --e2;
    --nval;
    --c__;
    --work;
    --iwork;

    *info = 0;
    if (*ijob < 1 || *ijob > 3) {
	*info = -1;
	return;
    }

    if (*ijob == 1) {

	*mout = 0;

	i__1 = *minp;
	for (ji = 1; ji <= i__1; ++ji) {
	    for (jp = 1; jp <= 2; ++jp) {
		tmp1 = d__[1] - ab[ji + jp * ab_dim1];
		if (std::abs(tmp1) < *pivmin) {
		    tmp1 = -(*pivmin);
		}
		nab[ji + jp * nab_dim1] = 0;
		if (tmp1 <= 0.) {
		    nab[ji + jp * nab_dim1] = 1;
		}

		i__2 = *n;
		for (j = 2; j <= i__2; ++j) {
		    tmp1 = d__[j] - e2[j - 1] / tmp1 - ab[ji + jp * ab_dim1];
		    if (std::abs(tmp1) < *pivmin) {
			tmp1 = -(*pivmin);
		    }
		    if (tmp1 <= 0.) {
			++nab[ji + jp * nab_dim1];
		    }
		}
	    }
	    *mout = *mout + nab[ji + (nab_dim1 << 1)] - nab[ji + nab_dim1];
	}
	return;
    }

    kf = 1;
    kl = *minp;

    if (*ijob == 2) {
	i__1 = *minp;
	for (ji = 1; ji <= i__1; ++ji) {
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
	}
    }

    i__1 = *nitmax;
    for (jit = 1; jit <= i__1; ++jit) {

	if (kl - kf + 1 >= *nbmin && *nbmin > 0) {

	    i__2 = kl;
	    for (ji = kf; ji <= i__2; ++ji) {

		work[ji] = d__[1] - c__[ji];
		iwork[ji] = 0;
		if (work[ji] <= *pivmin) {
		    iwork[ji] = 1;
		    d__1 = work[ji], d__2 = -(*pivmin);
		    work[ji] = (d__1<d__2) ? d__1 : d__2;
		}

		i__3 = *n;
		for (j = 2; j <= i__3; ++j) {
		    work[ji] = d__[j] - e2[j - 1] / work[ji] - c__[ji];
		    if (work[ji] <= *pivmin) {
			++iwork[ji];
			d__1 = work[ji], d__2 = -(*pivmin);
			work[ji] = (d__1<d__2) ? d__1 : d__2;
		    }
		}
	    }

	    if (*ijob <= 2) {

		klnew = kl;
		i__2 = kl;
		for (ji = kf; ji <= i__2; ++ji) {

		  i__5 = nab[ji + nab_dim1];
		  i__6 = iwork[ji];
		  i__3 = nab[ji + (nab_dim1 << 1)];
		  i__4 = (i__5>i__6) ? i__5 : i__6;
		    iwork[ji] = (i__3<i__4) ? i__3 : i__4;

		    if (iwork[ji] == nab[ji + (nab_dim1 << 1)]) {

			ab[ji + (ab_dim1 << 1)] = c__[ji];

		    } else if (iwork[ji] == nab[ji + nab_dim1]) {

			ab[ji + ab_dim1] = c__[ji];
		    } else {
			++klnew;
			if (klnew <= *mmax) {

			    ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 
				    1)];
			    nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 
				    << 1)];
			    ab[klnew + ab_dim1] = c__[ji];
			    nab[klnew + nab_dim1] = iwork[ji];
			    ab[ji + (ab_dim1 << 1)] = c__[ji];
			    nab[ji + (nab_dim1 << 1)] = iwork[ji];
			} else {
			    *info = *mmax + 1;
			}
		    }
		}
		if (*info != 0) {
		    return;
		}
		kl = klnew;
	    } else {

		i__2 = kl;
		for (ji = kf; ji <= i__2; ++ji) {
		    if (iwork[ji] <= nval[ji]) {
			ab[ji + ab_dim1] = c__[ji];
			nab[ji + nab_dim1] = iwork[ji];
		    }
		    if (iwork[ji] >= nval[ji]) {
			ab[ji + (ab_dim1 << 1)] = c__[ji];
			nab[ji + (nab_dim1 << 1)] = iwork[ji];
		    }
		}
	    }

	} else {

	    klnew = kl;
	    i__2 = kl;
	    for (ji = kf; ji <= i__2; ++ji) {

		tmp1 = c__[ji];
		tmp2 = d__[1] - tmp1;
		itmp1 = 0;
		if (tmp2 <= *pivmin) {
		    itmp1 = 1;
		    d__1 = tmp2, d__2 = -(*pivmin);
		    tmp2 = (d__1<d__2) ? d__1 : d__2;
		}

		i__3 = *n;
		for (j = 2; j <= i__3; ++j) {
		    tmp2 = d__[j] - e2[j - 1] / tmp2 - tmp1;
		    if (tmp2 <= *pivmin) {
			++itmp1;
			d__1 = tmp2, d__2 = -(*pivmin);
			tmp2 = (d__1<d__2) ? d__1 : d__2;
		    }
		}

		if (*ijob <= 2) {

		    i__5 = nab[ji + nab_dim1];
		    i__3 = nab[ji + (nab_dim1 << 1)];
		    i__4 = (i__5>itmp1) ? i__5 : itmp1;
		    itmp1 = (i__3<i__4) ? i__3 : i__4;

		    if (itmp1 == nab[ji + (nab_dim1 << 1)]) {

			ab[ji + (ab_dim1 << 1)] = tmp1;

		    } else if (itmp1 == nab[ji + nab_dim1]) {

			ab[ji + ab_dim1] = tmp1;
		    } else if (klnew < *mmax) {

			++klnew;
			ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 1)];
			nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 << 
				1)];
			ab[klnew + ab_dim1] = tmp1;
			nab[klnew + nab_dim1] = itmp1;
			ab[ji + (ab_dim1 << 1)] = tmp1;
			nab[ji + (nab_dim1 << 1)] = itmp1;
		    } else {
			*info = *mmax + 1;
			return;
		    }
		} else {

		    if (itmp1 <= nval[ji]) {
			ab[ji + ab_dim1] = tmp1;
			nab[ji + nab_dim1] = itmp1;
		    }
		    if (itmp1 >= nval[ji]) {
			ab[ji + (ab_dim1 << 1)] = tmp1;
			nab[ji + (nab_dim1 << 1)] = itmp1;
		    }
		}
	    }
	    kl = klnew;

	}

	kfnew = kf;
	i__2 = kl;
	for (ji = kf; ji <= i__2; ++ji) {
	    tmp1 = std::abs(ab[ji + (ab_dim1 << 1)] - ab[ji + ab_dim1]);
	    d__3 = std::abs(ab[ji + (ab_dim1 << 1)]);
	    d__4 = std::abs(ab[ji + ab_dim1]);
	    tmp2 = (d__3>d__4) ? d__3 : d__4;
	    d__1 = (*abstol>*pivmin) ? *abstol : *pivmin;
	    d__2 = *reltol * tmp2;
	    if (tmp1 < ((d__1>d__2) ? d__1 : d__2) || nab[ji + nab_dim1] >= nab[ji + (
		    nab_dim1 << 1)]) {

		if (ji > kfnew) {
		    tmp1 = ab[ji + ab_dim1];
		    tmp2 = ab[ji + (ab_dim1 << 1)];
		    itmp1 = nab[ji + nab_dim1];
		    itmp2 = nab[ji + (nab_dim1 << 1)];
		    ab[ji + ab_dim1] = ab[kfnew + ab_dim1];
		    ab[ji + (ab_dim1 << 1)] = ab[kfnew + (ab_dim1 << 1)];
		    nab[ji + nab_dim1] = nab[kfnew + nab_dim1];
		    nab[ji + (nab_dim1 << 1)] = nab[kfnew + (nab_dim1 << 1)];
		    ab[kfnew + ab_dim1] = tmp1;
		    ab[kfnew + (ab_dim1 << 1)] = tmp2;
		    nab[kfnew + nab_dim1] = itmp1;
		    nab[kfnew + (nab_dim1 << 1)] = itmp2;
		    if (*ijob == 3) {
			itmp1 = nval[ji];
			nval[ji] = nval[kfnew];
			nval[kfnew] = itmp1;
		    }
		}
		++kfnew;
	    }
	}
	kf = kfnew;

	i__2 = kl;
	for (ji = kf; ji <= i__2; ++ji) {
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
	}

	if (kf > kl) {
	    break;
	}
    }

    i__1 = kl + 1 - kf;
    if(i__1>0)
      *info = i__1;

    *mout = kl;

    return;

}


}
}
#include <cmath>

#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaed6,SLAED6)(int *kniter, 
                        int *orgati, 
                        float *rho, 
                        float *d__,
                        float *z__, 
                        float *finit, 
                        float *tau, 
                        int *info)
{
    int i__1;
    float r__1, r__2, r__3, r__4;

    float a, b, c__, f;
    int i__;
    float fc, df, ddf, eta, eps, base;
    int iter;
    float temp, temp1, temp2, temp3, temp4;
    int scale;
    int niter;
    float small1, small2, sminv1, sminv2, dscale[3], sclfac;
    float zscale[3], erretm;
    float safemin;
    float sclinv = 0;
    
    --z__;
    --d__;

    *info = 0;

    niter = 1;
    *tau = 0.f;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (d__[3] - d__[2]) / 2.f;
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
	} else {
	    temp = (d__[1] - d__[2]) / 2.f;
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
	}
        r__1 = std::abs(a), r__2 = std::abs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = std::abs(c__);
        temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    *tau = b / a;
	} else if (a <= 0.f) {
	    *tau = (a -  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1)))) / (
		    c__ * 2.f);
	} else {
	    *tau = b * 2.f / (a +  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1))));
	}

	temp = *rho + z__[1] / (d__[1] - *tau) + z__[2] / (d__[2] - *tau) + 
		z__[3] / (d__[3] - *tau);
	if (std::abs(*finit) <= std::abs(temp)) {
	    *tau = 0.f;
	}
    }

    eps = PLUMED_GMX_FLOAT_EPS;
    base = 2;
    safemin = PLUMED_GMX_FLOAT_MIN*(1.0+PLUMED_GMX_FLOAT_EPS);
    i__1 = static_cast<int>(std::log(safemin) / std::log(base) / 3.f);
    small1 = std::pow(base, static_cast<float>(i__1));
    sminv1 = 1.f / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;

    if (*orgati) {
	r__3 = (r__1 = d__[2] - *tau, std::abs(r__1)), r__4 = (r__2 = d__[3] - *
		tau, std::abs(r__2));
        temp = (r__3<r__4) ? r__3 : r__4;
    } else {
	r__3 = (r__1 = d__[1] - *tau, std::abs(r__1)), r__4 = (r__2 = d__[2] - *
		tau, std::abs(r__2));
	temp = (r__3<r__4) ? r__3 : r__4;
    }
    scale = 0;
    if (temp <= small1) {
	scale = 1;
	if (temp <= small2) {

	    sclfac = sminv2;
	    sclinv = small2;
	} else {

	    sclfac = sminv1;
	    sclinv = small1;

	}

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__] * sclfac;
	    zscale[i__ - 1] = z__[i__] * sclfac;
	}
	*tau *= sclfac;
    } else {

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__];
	    zscale[i__ - 1] = z__[i__];
	}
    }
    fc = 0.f;
    df = 0.f;
    ddf = 0.f;
    for (i__ = 1; i__ <= 3; ++i__) {
	temp = 1.f / (dscale[i__ - 1] - *tau);
	temp1 = zscale[i__ - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i__ - 1];
	df += temp2;
	ddf += temp3;
    }
    f = *finit + *tau * fc;

    if (std::abs(f) <= 0.f) {
	goto L60;
    }
    iter = niter + 1;
    for (niter = iter; niter <= 20; ++niter) {
	if (*orgati) {
	    temp1 = dscale[1] - *tau;
	    temp2 = dscale[2] - *tau;
	} else {
	    temp1 = dscale[0] - *tau;
	    temp2 = dscale[1] - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	r__1 = std::abs(a), r__2 = std::abs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = std::abs(c__);
	temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    eta = b / a;
	} else if (a <= 0.f) {
	    eta = (a -  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs(r__1)))) / ( c__ * 2.f);
	} else {
	    eta = b * 2.f / (a +  std::sqrt((r__1 = a * a - b * 4.f * c__, std::abs( r__1))));
	}
	if (f * eta >= 0.f) {
	    eta = -f / df;
	}
	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0.f && temp >= dscale[2]) {
		eta = (dscale[2] - *tau) / 2.f;
	    }

	    if (eta < 0.f && temp <= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	} else {
	    if (eta > 0.f && temp >= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	    if (eta < 0.f && temp <= dscale[0]) {
		eta = (dscale[0] - *tau) / 2.f;
	    }
	}
	*tau += eta;
	fc = 0.f;
	erretm = 0.f;
	df = 0.f;
	ddf = 0.f;
	for (i__ = 1; i__ <= 3; ++i__) {
	    temp = 1.f / (dscale[i__ - 1] - *tau);
	    temp1 = zscale[i__ - 1] * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    temp4 = temp1 / dscale[i__ - 1];
	    fc += temp4;
	    erretm += std::abs(temp4);
	    df += temp2;
	    ddf += temp3;
	}
	f = *finit + *tau * fc;
	erretm = (std::abs(*finit) + std::abs(*tau) * erretm) * 8.f + std::abs(*tau) * df;
	if (std::abs(f) <= eps * erretm) {
	    goto L60;
	}
    }
    *info = 1;
L60:
    if (scale) {
	*tau *= sclinv;
    }
    return;
} 


}
}
#include <cmath>
#include "real.h"

#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaev2,SLAEV2)(float *   a, 
	float *   b, 
	float *   c__, 
	float *   rt1, 
	float *   rt2, 
	float *   cs1, 
	float *   sn1)
{
    float d__1;

    float ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    int sgn1, sgn2;
    float acmn, acmx;

    sm = *a + *c__;
    df = *a - *c__;
    adf = std::abs(df);
    tb = *b + *b;
    ab = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf *  std::sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab *  std::sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab *  std::sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }
    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = std::abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. /  std::sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (std::abs(ab)<PLUMED_GMX_FLOAT_MIN) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. /  std::sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;

}


}
}
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slagtf,SLAGTF)(int *n, 
	float *a, 
	float *lambda, 
	float *b, 
	float *c__, 
	float *tol, 
	float *d__, 
	int *in, 
	int *info)
{
    int i__1;

    int k;
    float tl, eps, piv1, piv2, temp, mult, scale1, scale2;

    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (*n < 0) {
	*info = -1;
	return;
    }

    if (*n == 0) 
	return;
    
    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (std::abs(a[1])<PLUMED_GMX_FLOAT_MIN) {
	    in[1] = 1;
	}
	return;
    }

    eps = PLUMED_GMX_FLOAT_EPS;

    tl = (*tol>eps) ? *tol : eps;
    scale1 = std::abs(a[1]) + std::abs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = std::abs(c__[k]) + std::abs(a[k + 1]);
	if (k < *n - 1) {
	    scale2 += std::abs(b[k + 1]);
	}
	if (std::abs(a[k])<PLUMED_GMX_FLOAT_MIN) {
	    piv1 = 0.;
	} else {
	    piv1 = std::abs(a[k]) / scale1;
	}
	if (std::abs(c__[k])<PLUMED_GMX_FLOAT_MIN) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = std::abs(c__[k]) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c__[k] /= a[k];
		a[k + 1] -= c__[k] * b[k];
		if (k < *n - 1) {
		    d__[k] = 0.;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c__[k];
		a[k] = c__[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < *n - 1) {
		    d__[k] = b[k + 1];
		    b[k + 1] = -mult * d__[k];
		}
		b[k] = temp;
		c__[k] = mult;
	    }
	}
	if (((piv1>piv2) ? piv1 : piv2) <= tl && in[*n] == 0) {
	    in[*n] = k;
	}
    }
    if (std::abs(a[*n]) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return;

}


}
}
#include <stdlib.h>
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slagts,SLAGTS)(int *job, 
	int *n, 
	float *a, 
	float *b, 
	float *c__, 
	float *d__, 
	int *in, 
	float *y, 
	float *tol, 
	int *info)
{
    int i__1;
    float d__1, d__2, d__4, d__5;

    int k;
    float ak, eps, temp, pert, absak, sfmin;
    float bignum,minval;
    --y;
    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (abs(*job) > 2 || *job == 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }

    if (*n == 0) {
	return;
    }
    eps = PLUMED_GMX_FLOAT_EPS;
    minval = PLUMED_GMX_FLOAT_MIN;
    sfmin = minval / eps;

    bignum = 1. / sfmin;

    if (*job < 0) {
	if (*tol <= 0.) {
	    *tol = std::abs(a[1]);
	    if (*n > 1) {
		d__1 = *tol;
		d__2 = std::abs(a[2]);
		d__1 = (d__1>d__2) ? d__1 : d__2;
		d__2 = std::abs(b[1]);
		*tol = (d__1>d__2) ? d__1 : d__2;
	    }
	    i__1 = *n;
	    for (k = 3; k <= i__1; ++k) {
	      d__4 = *tol;
	      d__5 = std::abs(a[k]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = std::abs(b[k - 1]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = std::abs(d__[k - 2]);
	      *tol = (d__4>d__5) ? d__4 : d__5;
	    }
	    *tol *= eps;
	    if (std::abs(*tol)<PLUMED_GMX_FLOAT_MIN) {
		*tol = eps;
	    }
	}
    }

    if (1 == abs(*job)) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    if (in[k - 1] == 0) {
		y[k] -= c__[k - 1] * y[k - 1];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
	if (*job == 1) {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k == *n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_FLOAT_MIN || std::abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    for (k = *n; k >= 1; --k) {
		if (k + 2 <= *n) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k + 1 == *n) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;
L40:
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_FLOAT_MIN || std::abs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L40;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L40;
		    }
		}
		y[k] = temp / ak;
	    }
	}
    } else {

	if (*job == 2) {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_FLOAT_MIN || std::abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;

L70:
		absak = std::abs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (std::abs(absak)<PLUMED_GMX_FLOAT_MIN || std::abs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L70;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (std::abs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L70;
		    }
		}
		y[k] = temp / ak;
	    }
	}

	for (k = *n; k >= 2; --k) {
	    if (in[k - 1] == 0) {
		y[k - 1] -= c__[k - 1] * y[k];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
    }

    return;
}


}
}
#include "lapack.h"


/* LAPACK */


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)(int *n1,
                        int *n2,
                        float *a,
                        int *dtrd1,
                        int *dtrd2,
                        int *index)
{
  int n1sv = *n1;
  int n2sv = *n2;
  int i,ind1,ind2;

  if(*dtrd1>0)
    ind1 = 0;
  else
    ind1 = *n1-1;

  if(*dtrd2>0)
    ind2 = *n1;
  else
    ind2 = *n1+*n2-1;

  i = 0;
  
  while(n1sv>0 && n2sv>0) {
    if(a[ind1]<=a[ind2]) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
      n1sv--;
    } else {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
      n2sv--;
    }
  }

  if(n1sv==0) {
    for(n1sv=1;n1sv<=n2sv;n1sv++) {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
    } 
  } else {
    for(n2sv=1;n2sv<=n1sv;n2sv++) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
    } 
  }
  return;
}
}
}
#include <cctype>
#include <cmath>

#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
float
PLUMED_BLAS_F77_FUNC(slange,SLANGE)(const char *norm,
	int *m,
	int *n,
	float *a,
	int *lda,
	float *work)
{
  const char ch=std::toupper(*norm);
  float dtemp,sum,max,val,scale;
  int i,j;

  switch(ch) {
  case 'M':
    max = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++) {
	dtemp = std::abs(a[j*(*lda)+i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;

  case 'O':
  case '1':
    max = 0.0;
    for(j=0;j<*n;j++) {
      sum = 0.0;
      for(i=0;i<*m;i++) 
	sum += std::abs(a[j*(*lda)+i]);
      if(sum>max)
	max = sum;
    }
    val = max;
    break;

  case 'I':
    for(i=0;i<*m;i++)
      work[i] = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++)
	work[i] += std::abs(a[j*(*lda)+i]);
    max = 0;
    for(i=0;i<*m;i++)
      if(work[i]>max)
	max=work[i];
    val = max;
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = 1;
    for(j=0;j<*n;j++) 
      PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(m,&(a[j*(*lda)+0]),&i,&scale,&sum);
    val = scale* std::sqrt(sum);
    break;

  default:
    val = 0.0;
    break;
  }
  return val;
}
}
}
#include <cctype>
#include <cmath>

#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
float
PLUMED_BLAS_F77_FUNC(slanst,SLANST)(const char *norm,
	int *n,
	float *d,
	float *e)
{
  const char ch=std::toupper(*norm);
  float dtemp,max,val,scale,sum;
  int i,j;


  if(*n<=0)
    return 0.0;
  
  switch(ch) {
  case 'M':
    max = std::abs(d[*n-1]);
      for(i=0;i<(*n-1);i++) {
	dtemp = std::abs(d[i]);
	if(dtemp>max)
	  max = dtemp;
	dtemp = std::abs(e[i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;
    
  case 'O':
  case '1':
  case 'I':

    if(*n==1)
      val = std::abs(d[0]);
    else {
      max = std::abs(d[0]) + std::abs(e[0]);
      dtemp = std::abs(e[*n-2]) + std::abs(d[*n-1]);
      if(dtemp>max)
	max = dtemp;
      for(i=1;i<(*n-1);i++) {
	dtemp = std::abs(d[i]) + std::abs(e[i]) + std::abs(e[i-1]);
	if(dtemp>max)
	  max = dtemp;
      }
      val = max;
    }
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = *n-1;
    j = 1;
    if(*n>1) {
      PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(&i,e,&j,&scale,&sum);
      sum *= 2;
    }
    PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(n,d,&j,&scale,&sum);
    val = scale *  std::sqrt(sum);
    break;
    
  default:
    val = 0.0;
    break;
  }
  return val;
}
}
}
#include <cmath>


#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
float 
PLUMED_BLAS_F77_FUNC(slansy,SLANSY)(const char *norm, const char *uplo, int *n, float *a, int 
	*lda, float *work)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    float ret_val, d__1, d__2, d__3;
    int c__1 = 1;

    /* Local variables */
    int i__, j;
    float sum, absa, scale;
    float value =0.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;

    if (*n == 0) {
	value = 0.;
    } else if (*norm=='M' || *norm=='m') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = std::abs(a[i__ + j * a_dim1]);
		  value = (d__2>d__3) ? d__2 : d__3;
		}
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = std::abs(a[i__ + j * a_dim1]);
		    value =  (d__2>d__3) ? d__2 : d__3;
		}
	    }
	}
    } else if (*norm=='I' || *norm=='i' || *norm=='O' || *norm=='o' || *norm=='1') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    absa = std::abs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		work[j] = sum + std::abs(a[j + j * a_dim1]);
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = value, d__2 = work[i__];
		value =  (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.;
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = work[j] + std::abs(a[j + j * a_dim1]);
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    absa = std::abs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		if(sum>value)
		  value = sum;
	    }
	}
    } else if (*norm=='F' || *norm=='f' || *norm=='E' || *norm=='e') {

	scale = 0.;
	sum = 1.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(n, &a[a_offset], &i__1, &scale, &sum);
	value = scale *  std::sqrt(sum);
    }

    ret_val = value;
    return ret_val;
}


}
}
#include <cmath>
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
float
PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(float * x, float * y)
{
  float xabs,yabs;
  float w,z;

  xabs = std::abs(*x);
  yabs = std::abs(*y);
  
  if(xabs>yabs) {
    w = xabs;
    z = yabs;
  } else {
    w = yabs;
    z = xabs;
  }

  if( std::abs(z)<PLUMED_GMX_FLOAT_MIN) 
    return w;
  else {
    z = z/w;
    return w* std::sqrt(1.0+z*z);
  }
}
  
}
}
#include <cmath>

#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;

void PLUMED_BLAS_F77_FUNC(slar1vx,SLAR1VX)(int *n, 
	      int *b1, 
	      int *bn,
	      float *sigma, 
	      float *d__, 
	      float *l, 
	      float *ld, 
	      float *lld, 
	      float *eval, 
	      float *gersch, 
	      float *z__, 
	      float *ztz, 
	      float *mingma, 
	      int *r__, 
	      int *isuppz, 
	      float *work)
{
    int i__1;

    int i__, j;
    float s;
    int r1, r2;
    int to;
    float eps, tmp;
    int indp, inds, from;
    float dplus;
    int sawnan;
    int indumn;
    float dminus;

    --work;
    --isuppz;
    --z__;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    eps = PLUMED_GMX_FLOAT_EPS;
    if (*r__ == 0) {

	r1 = *b1;
	r2 = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r1 = i__;
		goto L20;
	    }
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r2 = i__;
		goto L40;
	    }
	}
    } else {
	r1 = *r__;
	r2 = *r__;
    }

L40:
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    sawnan = 0;

    if (*b1 == 1) {
	work[inds] = 0.;
    } else {
	work[inds] = lld[*b1 - 1];
    }
    s = work[inds] - *sigma;
    i__1 = r2 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	work[i__] = ld[i__] / dplus;
	work[inds + i__] = s * work[i__] * l[i__];
	s = work[inds + i__] - *sigma;
    }

    if (std::isnan(s)) {

	sawnan = 1;
	j = *b1 + 1;
L60:
    if (!std::isnan(work[inds + j])) {
	    ++j;
	    goto L60;
	}
	work[inds + j] = lld[j];
	s = work[inds + j] - *sigma;
	i__1 = r2 - 1;
	for (i__ = j + 1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    work[i__] = ld[i__] / dplus;
	    if (std::abs(work[i__])<PLUMED_GMX_FLOAT_MIN) {
		work[inds + i__] = lld[i__];
	    } else {
		work[inds + i__] = s * work[i__] * l[i__];
	    }
	    s = work[inds + i__] - *sigma;
	}
    }

    work[indp + *bn - 1] = d__[*bn] - *sigma;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
	dminus = lld[i__] + work[indp + i__];
	tmp = d__[i__] / dminus;
	work[indumn + i__] = l[i__] * tmp;
	work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
    }
    tmp = work[indp + r1 - 1];
    if (std::isnan(tmp)) {

	sawnan = 1;
	j = *bn - 3;
L90:
    if (!std::isnan(work[indp + j])) {
	    --j;
	    goto L90;
	}
	work[indp + j] = d__[j + 1] - *sigma;
	i__1 = r1;
	for (i__ = j; i__ >= i__1; --i__) {
	    dminus = lld[i__] + work[indp + i__];
	    tmp = d__[i__] / dminus;
	    work[indumn + i__] = l[i__] * tmp;
	    if (std::abs(tmp)<PLUMED_GMX_FLOAT_MIN) {
		work[indp + i__ - 1] = d__[i__] - *sigma;
	    } else {
		work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
	    }
	}
    }

    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (std::abs(*mingma)<PLUMED_GMX_FLOAT_MIN) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = work[inds + i__] + work[indp + i__];
	if (std::abs(tmp)<PLUMED_GMX_FLOAT_MIN) {
	    tmp = eps * work[inds + i__];
	}
	if (std::abs(tmp) < std::abs(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
    }

    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;
    if (! sawnan) {
	from = *r__ - 1;
	i__1 = *r__ - 32;
	to = (i__1>(*b1)) ? i__1 : (*b1);
L120:
	if (from >= *b1) {
	    i__1 = to;
	    for (i__ = from; i__ >= i__1; --i__) {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (std::abs(z__[to]) <= eps && std::abs(z__[to + 1]) <= eps) {
		isuppz[1] = to + 2;
	    } else {
		from = to - 1;
		i__1 = to - 32;
		to = (i__1>*b1) ? i__1 : *b1;
		goto L120;
	    }
	}
	from = *r__ + 1;
	i__1 = *r__ + 32;
	to = (i__1<*bn) ? i__1 : *bn;
L140:
	if (from <= *bn) {
	    i__1 = to;
	    for (i__ = from; i__ <= i__1; ++i__) {
		z__[i__] = -(work[indumn + i__ - 1] * z__[i__ - 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (std::abs(z__[to]) <= eps && std::abs(z__[to - 1]) <= eps) {
		isuppz[2] = to - 2;
	    } else {
		from = to + 1;
		i__1 = to + 32;
		to = (i__1<*bn) ? i__1 : *bn;
		goto L140;
	    }
	}
    } else {
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (std::abs(z__[i__ + 1])<PLUMED_GMX_FLOAT_MIN) {
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
	    } else {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
	    }
	    if (std::abs(z__[i__]) <= eps && std::abs(z__[i__ + 1]) <= eps) {
		isuppz[1] = i__ + 2;
		goto L170;
	    }
	    *ztz += z__[i__] * z__[i__];
	}
L170:
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (std::abs(z__[i__])<PLUMED_GMX_FLOAT_MIN) {
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
	    } else {
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    }
	    if (std::abs(z__[i__]) <= eps && std::abs(z__[i__ + 1]) <= eps) {
		isuppz[2] = i__ - 1;
		break;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
	}
    }

    return;

}


}
}
#include <cctype>
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarf,SLARF)(const char *side,
       int *m,
       int *n,
       float *v,
       int *incv,
       float *tau,
       float *c,
       int *ldc,
       float *work)
{
  const char ch=std::toupper(*side);
  float one = 1.0;
  float zero = 0.0;
  float minustau = -(*tau);
  int i1 = 1;


  if(ch=='L') {
    if(std::abs(*tau)>PLUMED_GMX_FLOAT_MIN) {
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      PLUMED_BLAS_F77_FUNC(sger,SGER)(m,n,&minustau,v,incv,work,&i1,c,ldc);
    }
  } else {
    if(std::abs(*tau)>PLUMED_GMX_FLOAT_MIN) {
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      PLUMED_BLAS_F77_FUNC(sger,SGER)(m,n,&minustau,work,&i1,v,incv,c,ldc);
    }
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)(const char *side, 
	const char *trans, 
	const char *direct, 
	const char *storev, 
	int *m, 
	int *n, 
	int *k, 
	float *v, 
	int *ldv, 
	float *t, 
	int *ldt, 
	float *c__,
	int *ldc, 
	float *work, 
	int *ldwork)
{
    int c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    int i__, j;
    char transt[1];
    int c__1 = 1;
    float one = 1.0;
    float minusone = -1.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }
    if (*trans=='N' || *trans=='n') {
      *(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }
    
    if (*storev=='C' || *storev=='c') {

	if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &one, &work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &one, &c__[(*k + 1) * c_dim1 + 1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}
	    }

	} else {

	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc)
			    ;
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    one, &work[work_offset], ldwork);
		}
		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    one, &c__[c_offset], ldc)
			    ;
		}
		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}
	    }
	}

    } else  if (*storev=='r' || *storev=='R') {
      if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}
		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &one, &work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &c__[(*k + 1) * c_dim1 
			    + 1], ldc);
		}
		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    }

	} else {

	    if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &one, &c__[c_offset], ldc);
		}

		PLUMED_BLAS_F77_FUNC(strmm,STRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(int   *n,
                        float *alpha,
                        float *x,
                        int    *incx,
                        float *tau)
{
  float xnorm,t;
  int    ti1,knt,j;
  float minval,safmin,rsafmn,beta;

  if(*n<=1) {
    *tau = 0;
    return;
  }

  ti1 = *n-1;

  xnorm = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(&ti1,x,incx);

  if(std::abs(xnorm)<PLUMED_GMX_FLOAT_MIN) {
    *tau = 0.0;
  } else {

    t = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(alpha,&xnorm);

    if(*alpha<0)
      beta = t;
    else
      beta = -t;

    minval = PLUMED_GMX_FLOAT_MIN;
    
    safmin = minval*(1.0+PLUMED_GMX_FLOAT_EPS) / PLUMED_GMX_FLOAT_EPS;

        
    if(std::abs(beta)<safmin) {

      knt = 0;
      rsafmn = 1.0 / safmin;
      
      while(std::abs(beta)<safmin) {
	knt++;
	ti1 = *n-1;
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&rsafmn,x,incx);
	beta *= rsafmn;
	*alpha *= rsafmn;
      }
      
      /* safmin <= beta <= 1 now */
      ti1 = *n-1;
      xnorm = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(&ti1,x,incx);
      t = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(alpha,&xnorm);
      
      if(*alpha<0)
	beta = t;
      else
	beta = -t;
      
      *tau = (beta-*alpha)/beta;

      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&t,x,incx);
   
      *alpha = beta;
      for(j=0;j<knt;j++)
	*alpha *= safmin;
    } else {
      *tau = (beta-*alpha)/beta;
      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&t,x,incx);
      *alpha = beta;
    }
  }
   
  return;
}
}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slarft,SLARFT)(const char *direct, 
	const char *storev, 
	int *n, 
	int *k, 
	float *v, 
	int *ldv, 
	float *tau, 
	float *t, 
	int *ldt)
{
    /* System generated locals */
    int t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    float d__1;

    /* Local variables */
    int i__, j;
    float vii;
    int c__1 = 1;
    float zero = 0.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    if (*n == 0) {
	return;
    }

    if (*direct=='F' || *direct=='f') {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (std::abs(tau[i__])<PLUMED_GMX_FLOAT_MIN) {

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		vii = v[i__ + i__ * v_dim1];
		v[i__ + i__ * v_dim1] = 1.;
		if (*storev=='C' || *storev=='c') {

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    d__1 = -tau[i__];
		    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1],
			     ldv, &v[i__ + i__ * v_dim1], &c__1, &zero, &t[
			    i__ * t_dim1 + 1], &c__1);
		} else {

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    d__1 = -tau[i__];
		    PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &d__1, &v[i__ * 
			    v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &
			    zero, &t[i__ * t_dim1 + 1], &c__1);
		}
		v[i__ + i__ * v_dim1] = vii;


		i__2 = i__ - 1;
		PLUMED_BLAS_F77_FUNC(strmv,STRMV)("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (std::abs(tau[i__])<PLUMED_GMX_FLOAT_MIN) {

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		if (i__ < *k) {
		    if (*storev=='C' || *storev=='c') {
			vii = v[*n - *k + i__ + i__ * v_dim1];
			v[*n - *k + i__ + i__ * v_dim1] = 1.;

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			d__1 = -tau[i__];
			PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1) 
				* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &
				c__1, &zero, &t[i__ + 1 + i__ * t_dim1], &
				c__1);
			v[*n - *k + i__ + i__ * v_dim1] = vii;
		    } else {
			vii = v[i__ + (*n - *k + i__) * v_dim1];
			v[i__ + (*n - *k + i__) * v_dim1] = 1.;

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			d__1 = -tau[i__];
			PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("No transpose", &i__1, &i__2, &d__1, &v[i__ + 
				1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &
				zero, &t[i__ + 1 + i__ * t_dim1], &c__1);
			v[i__ + (*n - *k + i__) * v_dim1] = vii;
		    }

		    i__1 = *k - i__;
		    PLUMED_BLAS_F77_FUNC(strmv,STRMV)("Lower", "No transpose", "Non-unit", &i__1, &t[i__ 
			    + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
			     t_dim1], &c__1)
			    ;
		}
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    }
    return;


}
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarnv,SLARNV)(int *idist, 
	int *iseed, 
	int *n, 
	float *x)
{
    int i__1, i__2, i__3;

    int i__;
    float u[128];
    int il, iv, il2;

    --x;
    --iseed;

    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
	i__2 = 64, i__3 = *n - iv + 1;
	il = (i__2<i__3) ? i__2 : i__3;
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

	PLUMED_BLAS_F77_FUNC(slaruv,SLARUV)(&iseed[1], &il2, u);

	if (*idist == 1) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1];
	    }
	} else if (*idist == 2) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
	    }
	} else if (*idist == 3) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
                x[iv + i__ - 1] =  std::sqrt(std::log(u[(i__ << 1) - 2]) * -2.) * 
		  std::cos(u[(i__ << 1) - 1] * (float)6.2831853071795864769252867663);
	    }
	}
    }
    return;

}
}
}
#include <cmath>

#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarrbx,SLARRBX)(int *n, 
	 float *d__, 
	 float *l, 
	 float *ld, 
	 float *lld, 
	 int *ifirst, 
	 int *ilast, 
	 float *rtol1, 
	 float *rtol2, 
	 int *offset, 
	 float *w, 
	 float *wgap, 
	 float *werr, 
	 float *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2;

    int i__, j, k, p;
    float s;
    int i1, i2, ii, kk;
    float fac, gap, mid;
    int cnt;
    float tmp, left;
    int nint, prev, next, nleft;
    float right, width, dplus;
    int nright, olnint;
    k = 0;
    right = 0.0;

    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    *info = 0;
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    i1 = *ifirst;
    i2 = *ifirst;
    prev = 0;
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	iwork[k - 1] = 1;
	i2 = i__;
    }

    i__ = i1;
    nint = 0;
L30:
    if (i__ <= i2) {
	ii = i__ - *offset;
	if (iwork[(i__ << 1) - 1] == 1) {
	    fac = 1.;
	    left = w[ii] - werr[ii];


L40:
	    if (i__ > i1 && left <= right) {
		left = right;
		cnt = i__ - 1;
	    } else {
		s = -left;
		cnt = 0;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    s = s * lld[j] / dplus - left;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		if (std::isnan(s)) {

		    cnt = 0;
		    s = -left;
		    i__1 = *n - 1;
		    for (j = 1; j <= i__1; ++j) {
			dplus = d__[j] + s;
			if (dplus < 0.) {
			    ++cnt;
			}
			tmp = lld[j] / dplus;
			if (std::abs(tmp)<PLUMED_GMX_FLOAT_MIN) {
			    s = lld[j] - left;
			} else {
			    s = s * tmp - left;
			}
		    }
		    dplus = d__[*n] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		if (cnt > i__ - 1) {
		    left -= werr[ii] * fac;
		    fac *= 2.;
		    goto L40;
		}
	    }
	    nleft = cnt + 1;
	    i1 = (i1<nleft) ? i1 : nleft;
	    fac = 1.;
	    right = w[ii] + werr[ii];
L60:
	    s = -right;
	    cnt = 0;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		dplus = d__[j] + s;
		s = s * lld[j] / dplus - right;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    if (std::isnan(s)) {

		cnt = 0;
		s = -right;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		    tmp = lld[j] / dplus;
		    if (std::abs(tmp)<PLUMED_GMX_FLOAT_MIN) {
			s = lld[j] - right;
		    } else {
			s = s * tmp - right;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt < i__) {
		right += werr[ii] * fac;
		fac *= 2.;
		goto L60;
	    }
	    cnt = (cnt<i2) ? cnt : i2;
	    ++nint;
	    k = nleft << 1;
	    work[k - 1] = left;
	    work[k] = right;
	    i__ = cnt + 1;
	    iwork[k - 1] = i__;
	    iwork[k] = cnt;
	    if (prev != nleft - 1) {
		work[k - 2] = left;
	    }
	    prev = nleft;
	} else {
	    right = work[i__ * 2];

	    ++iwork[k - 1];
	    prev = i__;
	    ++i__;
	}
	goto L30;
    }
    if (i__ <= *n && iwork[(i__ << 1) - 1] != -1) {
	work[(i__ << 1) - 1] = work[prev * 2];
    }

L80:
    prev = i1 - 1;
    olnint = nint;
    i__ = i1;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
	k = i__ << 1;
	left = work[k - 1];
	right = work[k];
	next = iwork[k - 1];
	nright = iwork[k];
	mid = (left + right) * .5;
	width = right - mid;
	d__1 = std::abs(left);
	d__2 = std::abs(right);
	tmp = (d__1>d__2) ? d__1 : d__2;

	gap = 0.;
	if (i__ == nright) {
	    if (prev > 0 && next <= *n) {
		d__1 = left - work[k - 2], d__2 = work[k + 1] - right;
		gap = (d__1<d__2) ? d__1 : d__2;
	    } else if (prev > 0) {
		gap = left - work[k - 2];
	    } else if (next <= *n) {
		gap = work[k + 1] - right;
	    }
	}
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
	if (width < ((d__1>d__2) ? d__1 : d__2)) {
	    --nint;
	    iwork[k - 1] = 0;
	    kk = k;
	    i__2 = nright;
	    for (j = i__ + 1; j <= i__2; ++j) {
		kk += 2;
		iwork[kk - 1] = 0;
		work[kk - 1] = left;
		work[kk] = right;
		wgap[j - 1 - *offset] = 0.;
	    }
	    if (i1 == i__) {
		i1 = next;
	    } else {
		iwork[(prev << 1) - 1] = next;
	    }
	    i__ = next;
	    continue;
	}
	prev = i__;

	s = -mid;
	cnt = 0;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    dplus = d__[j] + s;
	    s = s * lld[j] / dplus - mid;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	dplus = d__[*n] + s;
	if (dplus < 0.) {
	    ++cnt;
	}
	if (std::isnan(s)) {
	    cnt = 0;
	    s = -mid;
	    i__2 = *n - 1;
	    for (j = 1; j <= i__2; ++j) {
		dplus = d__[j] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		tmp = lld[j] / dplus;
		if (std::abs(tmp)<PLUMED_GMX_FLOAT_MIN) {
		    s = lld[j] - mid;
		} else {
		    s = s * tmp - mid;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	i__2 = i__ - 1, i__3 = (nright<cnt) ? nright : cnt;
	cnt = (i__2>i__3) ? i__2 : i__3;
	if (cnt == i__ - 1) {
	    work[k - 1] = mid;
	} else if (cnt == nright) {
	    work[k] = mid;
	} else {
	    iwork[k] = cnt;
	    ++cnt;
	    iwork[k - 1] = cnt;
	    kk = cnt << 1;
	    iwork[kk - 1] = next;
	    iwork[kk] = nright;
	    work[k] = mid;
	    work[kk - 1] = mid;
	    work[kk] = right;
	    prev = cnt;
	    if (cnt - 1 > i__) {
		work[kk - 2] = mid;
	    }
	    if (cnt > *ifirst && cnt <= *ilast) {
		++nint;
	    } else if (cnt <= *ifirst) {
		i1 = cnt;
	    }
	}
	i__ = next;
    }
    if (nint > 0) {
	goto L80;
    }
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
	if (iwork[k - 1] != -1) {
	    w[ii] = (work[k - 1] + work[k]) * .5;
	    werr[ii] = work[k] - w[ii];
	    if (i__ != *ilast) {
		wgap[ii] = work[k + 1] - work[k];
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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"



#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarrex,SLARREX)(const char *range,
	 int *n, 
	 float *vl, 
	 float *vu, 
	 int *il, 
	 int *iu, 
	 float *d__, 
	 float *e, 
	 float *tol, 
	 int *nsplit, 
	 int *isplit, 
	 int *m, 
	 float *w, 
	 int *iblock, 
	 int *indexw, 
	 float *gersch, 
	 float *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2;
    int c__1 = 1;
    int c__0 = 0;

    int i__, j, k;
    float s, gl;
    int in;
    float gu;
    int cnt;
    float eps, tau, nrm, tmp, vvl, vvu, offd;
    int iend, jblk, till, itmp;
    float rtol, delta, sigma;
    int iinfo;
    float width;
    int ibegin;
    int irange;
    float sgndef;
    int maxcnt;
    --iwork;
    --work;
    --gersch;
    --indexw;
    --iblock;
    --w;
    --isplit;
    --e;
    --d__;

    sigma = 0;
    irange = 0;
    sgndef = 0;
    maxcnt = 0;

    *info = 0;

    if (*range=='A' || *range=='a')
	irange = 1;
    else if (*range=='V' || *range=='v')
	irange = 2;
    else if (*range=='I' || *range=='i')
	irange = 3;
    

    *m = 0;
    eps = PLUMED_GMX_FLOAT_EPS;

    *nsplit = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__]) <= *tol) {
	    isplit[*nsplit] = i__;
	    ++(*nsplit);
	}
    }
    isplit[*nsplit] = *n;

    ibegin = 1;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];
	if (ibegin == iend) {
	    ++(*m);
	    w[*m] = d__[ibegin];
	    iblock[*m] = jblk;
	    indexw[*m] = 1;
	    e[iend] = 0.;
	    ibegin = iend + 1;
	    goto L170;
	}
	in = iend - ibegin + 1;

	gl = d__[ibegin] - std::abs(e[ibegin]);
	gu = d__[ibegin] + std::abs(e[ibegin]);
	gersch[(ibegin << 1) - 1] = gl;
	gersch[ibegin * 2] = gu;
	gersch[(iend << 1) - 1] = d__[iend] - std::abs(e[iend - 1]);
	gersch[iend * 2] = d__[iend] + std::abs(e[iend - 1]);
	d__1 = gersch[(iend << 1) - 1];
	gl = (d__1<gl) ? d__1 : gl;
	d__1 = gersch[iend * 2];
	gu = (d__1>gu) ? d__1 : gu;
	i__2 = iend - 1;
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
	    offd = std::abs(e[i__ - 1]) + std::abs(e[i__]);
	    gersch[(i__ << 1) - 1] = d__[i__] - offd;
	    d__1 = gersch[(i__ << 1) - 1];
	    gl = (d__1<gl) ? d__1 : gl;
	    gersch[i__ * 2] = d__[i__] + offd;
	    d__1 = gersch[i__ * 2];
	    gu = (d__1>gu) ? d__1 : gu;
	}
	d__1 = std::abs(gl), d__2 = std::abs(gu);
	nrm = (d__1>d__2) ? d__1 : d__2;

	width = gu - gl;
	i__2 = iend - 1;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    work[i__] = e[i__] * e[i__];
	}
	for (j = 1; j <= 2; ++j) {
	    if (j == 1) {
		tau = gl + width * .25;
	    } else {
		tau = gu - width * .25;
	    }
	    tmp = d__[ibegin] - tau;
	    if (tmp < 0.) {
		cnt = 1;
	    } else {
		cnt = 0;
	    }
	    i__2 = iend;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		tmp = d__[i__] - tau - work[i__ - 1] / tmp;
		if (tmp < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt == 0) {
		gl = tau;
	    } else if (cnt == in) {
		gu = tau;
	    }
	    if (j == 1) {
		maxcnt = cnt;
		sigma = gl;
		sgndef = 1.;
	    } else {
		if (in - cnt > maxcnt) {
		    sigma = gu;
		    sgndef = -1.;
		}
	    }
	}

	work[in * 3] = 1.;
	delta = eps;
	tau = sgndef * nrm;
L60:
	sigma -= delta * tau;
	work[1] = d__[ibegin] - sigma;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(in << 1) + i__] = 1. / work[i__];
	    tmp = e[j] * work[(in << 1) + i__];
	    work[i__ + 1] = d__[j + 1] - sigma - tmp * e[j];
	    work[in + i__] = tmp;
	    ++j;
	}
	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
	    if (tmp < 0. || std::abs(work[(in << 1) + i__])<PLUMED_GMX_FLOAT_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L60;
	    }
	}

	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[in * 3 + i__] = work[i__] * work[in + i__];
	    work[(in << 2) + i__] = work[in * 3 + i__] * work[in + i__];
	}
	if (sgndef > 0.) {
	    cnt = 1;
	    work[1] = (gl + gu) / 2. - sigma;
	    work[in + 1] = 0.;
	    work[(in << 1) + 1] = (gu - gl) / 2.;
	} else {
	    cnt = in;
	    work[in] = (gl + gu) / 2. - sigma;
	    work[in * 2] = 0.;
	    work[in * 3] = (gu - gl) / 2.;
	}
	rtol = eps * 4.;
	PLUMED_BLAS_F77_FUNC(slarrbx,SLARRBX)(&in, &d__[ibegin], &e[ibegin], &work[in * 3 + 1], &work[(in <<
		 2) + 1], &cnt, &cnt, &rtol, &rtol, &c__0, &work[1], &work[in 
		+ 1], &work[(in << 1) + 1], &work[in * 5 + 1], &iwork[1], &
		iinfo);
	if (sgndef > 0.) {
	    tau = work[1] - work[(in << 1) + 1];
	} else {
	    tau = work[in] + work[in * 3];
	}

	work[in * 3] = 1.;
	delta = eps * 2.;
L100:
	tau *= 1. - delta;

	s = -tau;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__] = d__[j] + s;
	    work[(in << 1) + i__] = 1. / work[i__];
	    work[in + i__] = e[j] * d__[j] * work[(in << 1) + i__];
	    s = s * work[in + i__] * e[j] - tau;
	    ++j;
	}
	work[in] = d__[iend] + s;

	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
	    if (tmp < 0. || std::abs(work[(in << 1) + i__])<PLUMED_GMX_FLOAT_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L100;
	    }
	}

	sigma += tau;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	e[iend] = sigma;
	tmp = (float) in * 4. * eps * (std::abs(sigma) + std::abs(tau));
	i__2 = iend;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    gersch[(i__ << 1) - 1] = gersch[(i__ << 1) - 1] - sigma - tmp;
	    gersch[i__ * 2] = gersch[i__ * 2] - sigma + tmp;
	}

	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(i__ << 1) - 1] = std::abs(d__[j]);
	    work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
	    ++j;
	}
	work[(in << 1) - 1] = std::abs(d__[iend]);

	PLUMED_BLAS_F77_FUNC(slasq2,SLASQ2)(&in, &work[1], info);
	if (*info != 0) {
	    return;
	}

	if (sgndef > 0.) {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = work[in - i__ + 1];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	} else {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = -work[i__];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	}
	ibegin = iend + 1;
L170:
	;
    }
    if (irange == 2) {
	*m = 0;
	ibegin = 1;
	i__1 = *nsplit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iend = isplit[i__];
	    vvl = *vl - e[iend];
	    vvu = *vu - e[iend];
	    i__2 = iend;
	    for (j = ibegin; j <= i__2; ++j) {
		if (vvl <= w[j] && w[j] <= vvu) {
		    ++(*m);
		    w[*m] = w[j];
		    iblock[*m] = i__;
		    indexw[*m] = j - ibegin + 1;
		}
	    }
	    ibegin = iend + 1;
	}
    } else if (irange == 3) {
	*m = *iu - *il + 1;
	if (*nsplit == 1) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = w[*il + i__ - 1];
		indexw[i__] = *il + i__ - 1;
	    }
	} else {
	    ibegin = 1;
	    i__1 = *nsplit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iend = isplit[i__];
		i__2 = iend;
		for (j = ibegin; j <= i__2; ++j) {
		    work[j] = w[j] + e[iend];
		}
		ibegin = iend + 1;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = i__;
		iwork[*n + i__] = iblock[i__];
	    }
	    PLUMED_BLAS_F77_FUNC(slasrt2,SLASRT2)("I", n, &work[1], &iwork[1], &iinfo);
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		itmp = iwork[*il + i__ - 1];
		work[i__] = w[itmp];
		iblock[i__] = iwork[*n + itmp];
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[*n + i__] = iwork[*il + i__ - 1];
		iwork[i__] = i__;
	    }
	    PLUMED_BLAS_F77_FUNC(ilasrt2,ILASRT2)("I", m, &iblock[1], &iwork[1], &iinfo);
	    j = 1;
	    itmp = iblock[j];
	    cnt = iwork[*n + iwork[j]];
	    if (itmp == 1) {
		ibegin = 1;
	    } else {
		ibegin = isplit[itmp - 1] + 1;
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = work[iwork[i__]];
		if (iblock[i__] != itmp || i__ == *m) {
		    if (iblock[i__] == itmp) {
			till = *m;
		    } else {
			till = i__ - 1;
		    }
		    i__2 = till - j + 1;
		    PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("I", &i__2, &w[j], &iinfo);
		    cnt = cnt - ibegin + 1;
		    i__2 = till;
		    for (k = j; k <= i__2; ++k) {
			indexw[k] = cnt + k - j;
		    }
		    j = i__;
		    itmp = iblock[j];
		    cnt = iwork[*n + iwork[j]];
		    ibegin = isplit[itmp - 1] + 1;
		    if (i__ == *m && till < *m) {
			indexw[*m] = cnt - ibegin + 1;
		    }
		} else {
		    i__2 = cnt, i__3 = iwork[*n + iwork[i__]];
		    cnt = (i__2<i__3) ? i__2 : i__3;
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

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarrfx,SLARRFX)(int *n, 
	float *d__, 
	float *l, 
	float *ld, 
	float *lld, 
	int *ifirst, 
	int *ilast, 
	float *w, 
	float *sigma, 
	float *dplus, 
	float *lplus, 
	float *work,
	int *info)
{
    int i1 = 1;
    int i__1;
    float d__2, d__3;

    int i__;
    float s, eps, tmp, dmax1, dmax2, delta;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;
    *info = 0;
    eps = PLUMED_GMX_FLOAT_EPS;
    *sigma = w[*ifirst];
    delta = eps * 2.;

L10:
    s = -(*sigma);
    dplus[1] = d__[1] + s;
    dmax1 = std::abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	s = s * lplus[i__] * l[i__] - *sigma;
	dplus[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax1, d__3 = std::abs(dplus[i__ + 1]);
	dmax1 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax1)) {
	*sigma -= std::abs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    s = -tmp;
    work[1] = d__[1] + s;
    dmax2 = std::abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	s = s * work[*n + i__] * l[i__] - tmp;
	work[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax2, d__3 = std::abs(work[i__ + 1]);
	dmax2 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax2)) {
	tmp += std::abs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }
    if (dmax2 < dmax1) {
	*sigma = tmp;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &work[1], &i1, &dplus[1], &i1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &work[*n + 1], &i1, &lplus[1], &i1);
    }

    return;
}
}
}
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slarrvx,SLARRVX)(int *n, 
	float *d__, 
	float *l, 
	int *isplit,
	int *m, 
	float *w,
	int *iblock, 
	int *indexw, 
	float *gersch, 
	float *tol, 
	float *z__, 
	int *ldz, 
	int *isuppz, 
	float *work, 
	int *iwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    float d__1, d__2;
    float c_b5 = 0.;
    int c__1 = 1;
    int c__2 = 2;

    int i__, j, k, p, q;
    int im, in;
    float gap, eps, tmp;
    int zto;
    float ztz;
    int iend, jblk;
    int wend, iter, temp[1], ktot;
    int itmp1, itmp2;
    int indld;
    float sigma;
    int ndone, iinfo, iindr;
    float resid;
    int nomgs;
    int nclus;
    int zfrom, iindc1, iindc2;
    float lambda;
    int ibegin;
    int indgap, indlld;
    float mingma;
    int oldien, oldncl, wbegin;
    float relgap;
    int oldcls;
    int ndepth, inderr, iindwk;
    int newcls, oldfst;
    float minrgp=0.0;
    int indwrk, oldlst;
    float reltol;
    int newfrs, newftt, parity;
    float mgstol, nrminv, rqcorr;
    int newlst, newsiz;


    --d__;
    --l;
    --isplit;
    --w;
    --iblock;
    --indexw;
    --gersch;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    inderr = *n;
    indld = *n << 1;
    indlld = *n * 3;
    indgap = *n << 2;
    indwrk = *n * 5 + 1;

    iindr = *n;
    iindc1 = *n << 1;
    iindc2 = *n * 3;
    iindwk = (*n << 2) + 1;

    eps = PLUMED_GMX_FLOAT_EPS;

    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("Full", n, m, &c_b5, &c_b5, &z__[z_offset], ldz);
    mgstol = eps * 100.;

    ibegin = 1;
    wbegin = 1;
    i__1 = iblock[*m];
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];

	wend = wbegin - 1;
L171:
	if (wend < *m) {
	    if (iblock[wend + 1] == jblk) {
		++wend;
		goto L171;
	    }
	}
	if (wend < wbegin) {
	    ibegin = iend + 1;
	    continue;
	}

	if (ibegin == iend) {
	    z__[ibegin + wbegin * z_dim1] = 1.;
	    isuppz[(wbegin << 1) - 1] = ibegin;
	    isuppz[wbegin * 2] = ibegin;
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	    continue;
	}
	oldien = ibegin - 1;
	in = iend - oldien;
	d__1 = .001, d__2 = 1. / (float) in;
	reltol = (d__1<d__2) ? d__1 : d__2;
	im = wend - wbegin + 1;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&im, &w[wbegin], &c__1, &work[1], &c__1);
	i__2 = im - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[inderr + i__] = eps * std::abs(work[i__]);
	    work[indgap + i__] = work[i__ + 1] - work[i__];
	}
	work[inderr + im] = eps * std::abs(work[im]);
	d__2 = std::abs(work[im]);
	work[indgap + im] = (d__2>eps) ? d__2 : eps;
	ndone = 0;

	ndepth = 0;
	parity = 1;
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;

L40:
	if (ndone < im) {
	    oldncl = nclus;
	    nclus = 0;
	    parity = 1 - parity;
	    if (parity == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
	    i__2 = oldncl;
	    for (i__ = 1; i__ <= i__2; ++i__) {

		j = oldcls + (i__ << 1);
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
		    j = wbegin + oldfst - 1;
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
		    i__3 = in - 1;
		    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
		    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz);
		}
		k = ibegin;
		i__3 = in - 1;
		for (j = 1; j <= i__3; ++j) {
		    tmp = d__[k] * l[k];
		    work[indld + j] = tmp;
		    work[indlld + j] = tmp * l[k];
		    ++k;
		}
		if (ndepth > 0) {

		    p = indexw[wbegin - 1 + oldfst];
		    q = indexw[wbegin - 1 + oldlst];
		    d__1 = eps * 4.;
		    i__3 = p - oldfst;
		    PLUMED_BLAS_F77_FUNC(slarrbx,SLARRBX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 1], &
			    work[indlld + 1], &p, &q, &reltol, &d__1, &i__3, &
			    work[1], &work[indgap + 1], &work[inderr + 1], &
			    work[indwrk + in], &iwork[iindwk], &iinfo);
		}
		newfrs = oldfst;
		i__3 = oldlst;
		for (j = oldfst; j <= i__3; ++j) {
		    if (j == oldlst || work[indgap + j] >= 
			reltol * std::abs(work[j])) {
			newlst = j;
		    } else {

			relgap = work[indgap + j] / std::abs(work[j]);
			if (j == newfrs) {
			    minrgp = relgap;
			} else {
			    minrgp = (minrgp<relgap) ? minrgp : relgap;
			}
			continue;
		    }
		    newsiz = newlst - newfrs + 1;
		    newftt = wbegin + newfrs - 1;
		    nomgs = newsiz == 1 || newsiz > 1 || minrgp < mgstol;
		    if (newsiz > 1 && nomgs) {

			PLUMED_BLAS_F77_FUNC(slarrfx,SLARRFX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				1], &work[indlld + 1], &newfrs, &newlst, &
				work[1], &sigma, &z__[ibegin + newftt * 
				z_dim1], &z__[ibegin + (newftt + 1) * z_dim1],
				 &work[indwrk], info);
			if (*info == 0) {
			    tmp = eps * std::abs(sigma);
			    i__4 = newlst;
			    for (k = newfrs; k <= i__4; ++k) {
				work[k] -= sigma;
				d__1 = work[indgap + k];
				work[indgap + k] = (d__1>tmp) ? d__1 : tmp;
				work[inderr + k] += tmp;
			    }
			    ++nclus;
			    k = newcls + (nclus << 1);
			    iwork[k - 1] = newfrs;
			    iwork[k] = newlst;
			} else {
			    *info = 0;
			    if (minrgp < mgstol) {

				work[indwrk] = d__[ibegin];
				i__4 = in - 1;
				for (k = 1; k <= i__4; ++k) {
				    work[indwrk + k] = d__[ibegin + k] + work[
					    indlld + k];
				}
				i__4 = newsiz;
				for (k = 1; k <= i__4; ++k) {
				    iwork[iindwk + k - 1] = 1;
				}
				i__4 = newlst;
				for (k = newfrs; k <= i__4; ++k) {
				    isuppz[2*(oldien + k) - 1] = 1;
				    isuppz[(oldien + k) * 2] = in;
				}
				temp[0] = in;
				PLUMED_BLAS_F77_FUNC(sstein,SSTEIN)(&in, &work[indwrk], &work[indld + 1], 
					&newsiz, &work[newfrs], &iwork[iindwk]
					, temp, &z__[ibegin + newftt * z_dim1]
					, ldz, &work[indwrk + in], &iwork[
					iindwk + in], &iwork[iindwk + (in*2)], &iinfo);
				if (iinfo != 0) {
				    *info = 2;
				    return;
				}
				ndone += newsiz;
			    }
			}
		    } else {
			ktot = newftt;
			i__4 = newlst;
			for (k = newfrs; k <= i__4; ++k) {
			    iter = 0;
L90:
			    lambda = work[k];

			    PLUMED_BLAS_F77_FUNC(slar1vx,SLAR1VX)(&in, &c__1, &in, &lambda, &d__[ibegin], &
				    l[ibegin], &work[indld + 1], &work[indlld 
				    + 1], &w[wbegin + k - 1], &gersch[(oldien 
				    << 1) + 1], &z__[ibegin + ktot * z_dim1], 
				    &ztz, &mingma, &iwork[iindr + ktot], &
				    isuppz[(ktot << 1) - 1], &work[indwrk]);
			    tmp = 1. / ztz;
			    nrminv =  std::sqrt(tmp);
			    resid = std::abs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
				d__1 = work[indgap + k - 1], d__2 = work[
					indgap + k];
				gap = (d__1<d__2) ? d__1 : d__2;
			    }
			    ++iter;
			    if (resid > *tol * gap && std::abs(rqcorr) > eps * 4. *
				     std::abs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    goto L90;
				}
			    }
			    iwork[ktot] = 1;
			    if (newsiz == 1) {
				++ndone;
			    }
			    zfrom = isuppz[(ktot << 1) - 1];
			    zto = isuppz[ktot * 2];
			    i__5 = zto - zfrom + 1;
			    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__5, &nrminv, &z__[ibegin + zfrom - 1 + 
				    ktot * z_dim1], &c__1);
			    ++ktot;
			}
			if (newsiz > 1) {
			    itmp1 = isuppz[(newftt << 1) - 1];
			    itmp2 = isuppz[newftt * 2];
			    ktot = oldien + newlst;
			    i__4 = ktot;
			    for (p = newftt + 1; p <= i__4; ++p) {
				i__5 = p - 1;
				for (q = newftt; q <= i__5; ++q) {
				    tmp = -PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&in, &z__[ibegin + p * 
					    z_dim1], &c__1, &z__[ibegin + q * 
					    z_dim1], &c__1);
				    PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&in, &tmp, &z__[ibegin + q * 
					    z_dim1], &c__1, &z__[ibegin + p * 
					    z_dim1], &c__1);
				}
				tmp = 1. / PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(&in, &z__[ibegin + p * 
					z_dim1], &c__1);
				PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&in, &tmp, &z__[ibegin + p * z_dim1], &
					c__1);
				i__5 = itmp1, i__6 = isuppz[(p << 1) - 1];
				itmp1 = (i__5<i__6) ? i__5 : i__6;
				i__5 = itmp2, i__6 = isuppz[p * 2];
				itmp2 = (i__5>i__6) ? i__5 : i__6;
			    }
			    i__4 = ktot;
			    for (p = newftt; p <= i__4; ++p) {
				isuppz[(p << 1) - 1] = itmp1;
				isuppz[p * 2] = itmp2;
			    }
			    ndone += newsiz;
			}
		    }
		    newfrs = j + 1;
		}
	    }
	    ++ndepth;
	    goto L40;
	}
	j = wbegin << 1;
	i__2 = wend;
	for (i__ = wbegin; i__ <= i__2; ++i__) {
	    isuppz[j - 1] += oldien;
	    isuppz[j] += oldien;
	    j += 2;

	}
	ibegin = iend + 1;
	wbegin = wend + 1;
    }

    return;

} 
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(float *f,
	float *g,
	float *cs,
	float *sn,
	float *r)
{
  float minval,safemin, safemin2, safemx2, eps;
  float f1,g1,f1a,g1a,scale;
  int i,n,count;

  eps = PLUMED_GMX_FLOAT_EPS;
  minval = PLUMED_GMX_FLOAT_MIN;
  safemin = minval*(1.0+eps);
  n = static_cast<int>(0.5*std::log( safemin/eps ) / std::log(2.0));
  safemin2 = std::pow(static_cast<float>(2.0),static_cast<float>(n));

  safemx2 = 1.0 / safemin2;

  if(std::abs(*g)<PLUMED_GMX_FLOAT_MIN) {
    *cs = 1.0;
    *sn = 0.0;
    *r = *f;
  } else if (std::abs(*f)<PLUMED_GMX_FLOAT_MIN) {
    *cs = 0.0;
    *sn = 1.0;
    *r = *g;
  } else {
    f1 = *f;
    g1 = *g;
    f1a = std::abs(f1);
    g1a = std::abs(g1);
    scale = (f1a > g1a) ? f1a : g1a;
    if(scale >= safemx2) {
      count = 0;
      while(scale >= safemx2) {
	count++;
	f1 *= safemin2;
	g1 *= safemin2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemx2;
    } else if (scale<=safemin2) {
      count = 0;
      while(scale <= safemin2) {
	count++;
	f1 *= safemx2;
	g1 *= safemx2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemin2;
    } else {
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
    }
    if(std::abs(*f)>std::abs(*g) && *cs<0.0) {
      *cs *= -1.0;
      *sn *= -1.0;
      *r  *= -1.0;
    }
  }
  return;
}
      
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaruv,SLARUV)(int *iseed, int *n, float *x)
{
  const int
    mm[512] = {
      494,2637,255,2008,1253,
      3344,4084,1739,3143,3468,688,1657,1238,3166,1292,3422,1270,2016,
      154,2862,697,1706,491,931,1444,444,3577,3944,2184,1661,3482,657,
      3023,3618,1267,1828,164,3798,3087,2400,2870,3876,1905,1593,1797,
      1234,3460,328,2861,1950,617,2070,3331,769,1558,2412,2800,189,287,
      2045,1227,2838,209,2770,3654,3993,192,2253,3491,2889,2857,2094,
      1818,688,1407,634,3231,815,3524,1914,516,164,303,2144,3480,119,
      3357,837,2826,2332,2089,3780,1700,3712,150,2000,3375,1621,3090,
      3765,1149,3146,33,3082,2741,359,3316,1749,185,2784,2202,2199,1364,
      1244,2020,3160,2785,2772,1217,1822,1245,2252,3904,2774,997,2573,
      1148,545,322,789,1440,752,2859,123,1848,643,2405,2638,2344,46,
      3814,913,3649,339,3808,822,2832,3078,3633,2970,637,2249,2081,4019,
      1478,242,481,2075,4058,622,3376,812,234,641,4005,1122,3135,2640,
      2302,40,1832,2247,2034,2637,1287,1691,496,1597,2394,2584,1843,336,
      1472,2407,433,2096,1761,2810,566,442,41,1238,1086,603,840,3168,
      1499,1084,3438,2408,1589,2391,288,26,512,1456,171,1677,2657,2270,
      2587,2961,1970,1817,676,1410,3723,2803,3185,184,663,499,3784,1631,
      1925,3912,1398,1349,1441,2224,2411,1907,3192,2786,382,37,759,2948,
      1862,3802,2423,2051,2295,1332,1832,2405,3638,3661,327,3660,716,
      1842,3987,1368,1848,2366,2508,3754,1766,3572,2893,307,1297,3966,
      758,2598,3406,2922,1038,2934,2091,2451,1580,1958,2055,1507,1078,
      3273,17,854,2916,3971,2889,3831,2621,1541,893,736,3992,787,2125,
      2364,2460,257,1574,3912,1216,3248,3401,2124,2762,149,2245,166,466,
      4018,1399,190,2879,153,2320,18,712,2159,2318,2091,3443,1510,449,
      1956,2201,3137,3399,1321,2271,3667,2703,629,2365,2431,1113,3922,
      2554,184,2099,3228,4012,1921,3452,3901,572,3309,3171,817,3039,
      1696,1256,3715,2077,3019,1497,1101,717,51,981,1978,1813,3881,76,
      3846,3694,1682,124,1660,3997,479,1141,886,3514,1301,3604,1888,
      1836,1990,2058,692,1194,20,3285,2046,2107,3508,3525,3801,2549,
      1145,2253,305,3301,1065,3133,2913,3285,1241,1197,3729,2501,1673,
      541,2753,949,2361,1165,4081,2725,3305,3069,3617,3733,409,2157,
      1361,3973,1865,2525,1409,3445,3577,77,3761,2149,1449,3005,225,85,
      3673,3117,3089,1349,2057,413,65,1845,697,3085,3441,1573,3689,2941,
      929,533,2841,4077,721,2821,2249,2397,2817,245,1913,1997,3121,997,
      1833,2877,1633,981,2009,941,2449,197,2441,285,1473,2741,3129,909,
      2801,421,4073,2813,2337,1429,1177,1901,81,1669,2633,2269,129,1141,
      249,3917,2481,3941,2217,2749,3041,1877,345,2861,1809,3141,2825,
      157,2881,3637,1465,2829,2161,3365,361,2685,3745,2325,3609,3821,
      3537,517,3017,2141,1537 
    };

    int i__1;

    int i__, i1, i2, i3, i4, it1, it2, it3, it4;


    --iseed;
    --x;

    it1 = it2 = it3 = it4 = 0;

    i1 = iseed[1];
    i2 = iseed[2];
    i3 = iseed[3];
    i4 = iseed[4];

    i__1 = (*n<128) ? *n : 128;
    for (i__ = 1; i__ <= i__1; ++i__) {

	it4 = i4 * mm[i__ + 383];
	it3 = it4 / 4096;
	it4 -= it3 << 12;
	it3 = it3 + i3 * mm[i__ + 383] + i4 * mm[i__ + 255];
	it2 = it3 / 4096;
	it3 -= it2 << 12;
	it2 = it2 + i2 * mm[i__ + 383] + i3 * mm[i__ + 255] + 
	  i4 * mm[i__ + 127];
	it1 = it2 / 4096;
	it2 -= it1 << 12;
	it1 = it1 + i1 * mm[i__ + 383] + i2 * mm[i__ + 255] + 
	  i3 * mm[i__ +	127] + i4 * mm[i__ - 1];
	it1 %= 4096;

	x[i__] = ((float) it1 + ((float) it2 + ((float) it3 + (
		float) it4 * 2.44140625e-4) * 2.44140625e-4) * 
		2.44140625e-4) * 2.44140625e-4;
    }

    iseed[1] = it1;
    iseed[2] = it2;
    iseed[3] = it3;
    iseed[4] = it4;
    return;

} 
}
}
#include <cmath>
#include "real.h"

#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slas2,SLAS2)(float *f,
       float *g,
       float *h,
       float *ssmin,
       float *ssmax)
{
  float fa = std::abs(*f);
  float ga = std::abs(*g);
  float ha = std::abs(*h);
  float fhmin,fhmax,tmax,tmin,tmp1,tmp2;
  float as,at,au,c;

  fhmin = (fa<ha) ? fa : ha;
  fhmax = (fa>ha) ? fa : ha;
  
  if(std::abs(fhmin)<PLUMED_GMX_FLOAT_MIN) {
    *ssmin = 0.0;
    if(std::abs(fhmax)<PLUMED_GMX_FLOAT_MIN) 
      *ssmax = ga;
    else {
      tmax = (fhmax>ga) ? fhmax : ga;
      tmin = (fhmax<ga) ? fhmax : ga;
      tmp1 = tmin / tmax;
      tmp1 = tmp1 * tmp1;
      *ssmax = tmax* std::sqrt(1.0 + tmp1);
    }
  } else {
    if(ga<fhmax) {
      as = 1.0 + fhmin / fhmax;
      at = (fhmax-fhmin) / fhmax;
      au = (ga/fhmax);
      au = au * au;
      c = 2.0 / (  std::sqrt(as*as+au) + std::sqrt(at*at+au) );
      *ssmin = fhmin * c;
      *ssmax = fhmax / c;
    } else {
      au = fhmax / ga;
      if(std::abs(au)<PLUMED_GMX_FLOAT_MIN) {
	*ssmin = (fhmin*fhmax)/ga;
	*ssmax = ga;
      } else {
	as = 1.0 + fhmin / fhmax;
	at = (fhmax-fhmin)/fhmax;
	tmp1 = as*au;
	tmp2 = at*au;
	c = 1.0 / (  std::sqrt(1.0+tmp1*tmp1) + std::sqrt(1.0+tmp2*tmp2));
	*ssmin = (fhmin*c)*au;
	*ssmin = *ssmin + *ssmin;
	*ssmax = ga / (c+c);
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
  const char ch=std::toupper(*type);
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

    if(std::abs(cfrom1)>std::abs(ctoc) && std::abs(ctoc)>PLUMED_GMX_FLOAT_MIN) {
      mul = smlnum;
      done = 0;
      cfromc = cfrom1;
    } else if(std::abs(cto1)>std::abs(cfromc)) {
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
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd0,SLASD0)(int *n, 
	int *sqre, 
	float *d__, 
	float *e, 
	float *u, 
	int *ldu, 
	float *vt, 
	int *ldvt,
	int *smlsiz, 
	int *iwork,
	float *work, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int i__, j, m, i1, ic, lf, nd, ll, nl, nr, im1, ncc, nlf, nrf, 
	    iwk, lvl, ndb1, nlp1, nrp1;
    float beta;
    int idxq, nlvl;
    float alpha;
    int inode, ndiml, idxqc, ndimr, itemp, sqrei;
    int c__0 = 0;


    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --iwork;
    --work;

    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -2;
    }

    m = *n + *sqre;

    if (*ldu < *n) {
	*info = -6;
    } else if (*ldvt < m) {
	*info = -8;
    } else if (*smlsiz < 3) {
	*info = -9;
    }
    if (*info != 0) {
	return;
    }

    if (*n <= *smlsiz) {
	PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], 
		ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info);
	return;
    }

    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;
    PLUMED_BLAS_F77_FUNC(slasdt,SLASDT)(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

    ndb1 = (nd + 1) / 2;
    ncc = 0;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {

	i1 = i__ - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nrp1 = nr + 1;
	nlf = ic - nl;
	nrf = ic + 1;
	sqrei = 1;
	PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[
		nlf + nlf * vt_dim1], ldvt, &u[nlf + nlf * u_dim1], ldu, &u[
		nlf + nlf * u_dim1], ldu, &work[1], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + nlf - 2;
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    iwork[itemp + j] = j;
	}
	if (i__ == nd) {
	    sqrei = *sqre;
	} else {
	    sqrei = 1;
	}
	nrp1 = nr + sqrei;
	PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[
		nrf + nrf * vt_dim1], ldvt, &u[nrf + nrf * u_dim1], ldu, &u[
		nrf + nrf * u_dim1], ldu, &work[1], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + ic;
	i__2 = nr;
	for (j = 1; j <= i__2; ++j) {
	    iwork[itemp + j - 1] = j;
	}
    }

    for (lvl = nlvl; lvl >= 1; --lvl) {

	if (lvl == 1) {
	    lf = 1;
	    ll = 1;
	} else {
	    i__1 = lvl - 1;
	    lf = (1 << i__1);
	    ll = (lf << 1) - 1;
	}
	i__1 = ll;
	for (i__ = lf; i__ <= i__1; ++i__) {
	    im1 = i__ - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    if (*sqre == 0 && i__ == ll) {
		sqrei = *sqre;
	    } else {
		sqrei = 1;
	    }
	    idxqc = idxq + nlf - 1;
	    alpha = d__[ic];
	    beta = e[ic];
	    PLUMED_BLAS_F77_FUNC(slasd1,SLASD1)(&nl, &nr, &sqrei, &d__[nlf], &alpha, &beta, &u[nlf + nlf *
		     u_dim1], ldu, &vt[nlf + nlf * vt_dim1], ldvt, &iwork[
		    idxqc], &iwork[iwk], &work[1], info);
	    if (*info != 0) {
		return;
	    }
	}
    }

    return;

}
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd1,SLASD1)(int *nl, 
	int *nr, 
	int *sqre, 
	float *d__, 
	float *alpha, 
	float *beta, 
	float *u, 
	int *ldu, 
	float *vt, 
	int *ldvt, 
	int *idxq, 
	int *iwork, 
	float *work, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    float d__1, d__2;

    int i__, k, m, n, n1, n2, iq, iz, iu2, ldq, idx, ldu2, ivt2, 
	    idxc, idxp, ldvt2;
    int isigma;
    float orgnrm;
    int coltyp;
    int c__0 = 0;
    float one = 1.0;
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

    d__1 = std::abs(*alpha);
    d__2 = std::abs(*beta);
    orgnrm = (d__1>d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(d__[i__]) > orgnrm) {
	    orgnrm = std::abs(d__[i__]);
	}
    }
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    PLUMED_BLAS_F77_FUNC(slasd2,SLASD2)(nl, nr, sqre, &k, &d__[1], &work[iz], alpha, beta, &u[u_offset], 
	    ldu, &vt[vt_offset], ldvt, &work[isigma], &work[iu2], &ldu2, &
	    work[ivt2], &ldvt2, &iwork[idxp], &iwork[idx], &iwork[idxc], &
	    idxq[1], &iwork[coltyp], info);

    ldq = k;
    PLUMED_BLAS_F77_FUNC(slasd3,SLASD3)(nl, nr, sqre, &k, &d__[1], &work[iq], &ldq, &work[isigma], &u[
	    u_offset], ldu, &work[iu2], &ldu2, &vt[vt_offset], ldvt, &work[
	    ivt2], &ldvt2, &iwork[idxc], &iwork[coltyp], &work[iz], info);
    if (*info != 0) {
	return;
    }
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = k;
    n2 = n - k;
    PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd2,SLASD2)(int *nl, 
                        int *nr, 
                        int *sqre, 
                        int *k, 
                        float *d__, 
                        float *z__, 
                        float *alpha, 
                        float *beta, 
                        float *u, 
                        int *ldu, 
                        float *vt, 
                        int *ldvt, 
                        float *dsigma, 
                        float *u2, 
                        int *ldu2, 
                        float *vt2, 
                        int *ldvt2, 
                        int *idxp, 
                        int *idx, 
                        int *idxc, 
                        int *idxq, 
                        int *coltyp, 
                        int *info)
{
    int u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset;
    int vt2_dim1, vt2_offset, i__1;
    float d__1, d__2;

    float c__;
    int i__, j, m, n;
    float s;
    int k2;
    float z1;
    int ct, jp;
    float eps, tau, tol;
    int psm[4], nlp1, nlp2, idxi, idxj;
    int ctot[4], idxjp;
    int jprev = 0;
    float hlftol;
    float zero = 0.0;
    int c__1 = 1;


    --d__;
    --z__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --dsigma;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxp;
    --idx;
    --idxc;
    --idxq;
    --coltyp;

    *info = 0;

    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
    z__[1] = z1;
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
    }

    i__1 = nlp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	coltyp[i__] = 1;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	coltyp[i__] = 2;
    }

    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	u2[i__ + u2_dim1] = z__[idxq[i__]];
	idxc[i__] = coltyp[idxq[i__]];
    }

    PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = u2[idxi + u2_dim1];
	coltyp[i__] = idxc[idxi];
    }

    eps = PLUMED_GMX_FLOAT_EPS;
    d__1 = std::abs(*alpha), d__2 = std::abs(*beta);
    tol = (d__1 > d__2) ? d__1 : d__2;
    d__2 = std::abs(d__[n]);
    tol = eps * 8. * ((d__2 > tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (std::abs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    coltyp[j] = 4;
	    if (j == n) {
		goto L120;
	    }
	} else {
	    jprev = j;
	    goto L90;
	}
    }
L90:
    j = jprev;
L100:
    ++j;
    if (j > n) {
	goto L110;
    }
    if (std::abs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
	coltyp[j] = 4;
    } else {

	if (std::abs(d__[j] - d__[jprev]) <= tol) {

            s = z__[jprev];
	    c__ = z__[j];

	    tau = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&c__, &s);
	    c__ /= tau;
	    s = -s / tau;
	    z__[j] = tau;
	    z__[jprev] = 0.;

	    idxjp = idxq[idx[jprev] + 1];
	    idxj = idxq[idx[j] + 1];
	    if (idxjp <= nlp1) {
		--idxjp;
	    }
	    if (idxj <= nlp1) {
		--idxj;
	    }
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], &
		    c__1, &c__, &s);
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, &
		    c__, &s);
	    if (coltyp[j] != coltyp[jprev]) {
		coltyp[j] = 3;
	    }
	    coltyp[jprev] = 4;
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    u2[*k + u2_dim1] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L100;
L110:

    ++(*k);
    u2[*k + u2_dim1] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L120:

    for (j = 1; j <= 4; ++j) {
	ctot[j - 1] = 0;
    }
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	ct = coltyp[j];
	++ctot[ct - 1];
    }

    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	ct = coltyp[jp];
	idxc[psm[ct - 1]] = j;
	++psm[ct - 1];
    }

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	idxj = idxq[idx[idxp[idxc[j]]] + 1];
	if (idxj <= nlp1) {
	    --idxj;
	}
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
    }

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (std::abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    c__ = 1.;
	    s = 0.;
	    z__[1] = tol;
	} else {
	    c__ = z1 / z__[1];
	    s = z__[m] / z__[1];
	}
    } else {
	if (std::abs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);

    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &n, &c__1, &zero, &zero, &u2[u2_offset], ldu2);
    u2[nlp1 + u2_dim1] = 1.;
    if (m > n) {
	i__1 = nlp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
	    vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
	    vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
	}
    } else {
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
    }
    if (m > n) {
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
    }

    if (n > *k) {
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1)
		 * u_dim1 + 1], ldu);
	i__1 = n - *k;
	PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + 
		vt_dim1], ldvt);
    }
    for (j = 1; j <= 4; ++j) {
	coltyp[j] = ctot[j - 1];
    }

    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd3,SLASD3)(int *nl, 
	int *nr,
	int *sqre, 
	int *k, 
	float *d__, 
	float *q, 
	int *ldq, 
	float *dsigma, 
	float *u, 
	int *ldu, 
	float *u2, 
	int *ldu2, 
	float *vt, 
	int *ldvt, 
	float *vt2, 
	int *ldvt2, 
	int *idxc, 
	int *ctot, 
	float *z__, 
	int *info)
{
    int q_dim1, q_offset, u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, 
	    vt_offset, vt2_dim1, vt2_offset, i__1, i__2;
    float d__2;

    int i__, j, m, n, jc;
    float rho;
    int nlp1, nlp2, nrp1;
    float temp;
    int ctemp;
    int ktemp;
    int c__1 = 1;
    int c__0 = 0;
    float zero = 0.0;
    float one = 1.0;

    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dsigma;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxc;
    --ctot;
    --z__;

    /* Function Body */
    *info = 0;

    if (*nl < 1) {
	*info = -1;
    } else if (*nr < 1) {
	*info = -2;
    } else if (*sqre != 1 && *sqre != 0) {
	*info = -3;
    }

    n = *nl + *nr + 1;
    m = n + *sqre;
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    if (*k == 1) {
	d__[1] = std::abs(z__[1]);
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
	if (z__[1] > 0.) {
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
	} else {
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		u[i__ + u_dim1] = -u2[i__ + u2_dim1];
	    }
	}
	return;
    }

    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &z__[1], &c__1, &q[q_offset], &c__1);

    rho = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(k, &z__[1], &c__1);
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;


    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	PLUMED_BLAS_F77_FUNC(slasd4,SLASD4)(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j],
		 &vt[j * vt_dim1 + 1], info);

	if (*info != 0) {
	    return;
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
	}
	i__2 = *k - 1;
	for (j = i__; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
	}
	d__2 =  std::sqrt(std::abs(z__[i__]));
	z__[i__] = (q[i__ + q_dim1] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * 
		vt_dim1 + 1];
	u[i__ * u_dim1 + 1] = -1.;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ 
		    * vt_dim1];
	    u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
	}
	temp = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(k, &u[i__ * u_dim1 + 1], &c__1);
	q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
	}
    }

    if (*k == 2) {
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", &n, k, k, &one, &u2[u2_offset], ldu2, &q[q_offset],
		 ldq, &zero, &u[u_offset], ldu);
	goto L100;
    }
    if (ctot[1] > 0) {
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", nl, k, &ctot[1], &one, &u2[(u2_dim1 << 1) + 1], 
		ldu2, &q[q_dim1 + 2], ldq, &zero, &u[u_dim1 + 1], ldu);
	if (ctot[3] > 0) {
	    ktemp = ctot[1] + 2 + ctot[2];
	    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1]
		    , ldu2, &q[ktemp + q_dim1], ldq, &one, &u[u_dim1 + 1], 
		    ldu);
	}
    } else if (ctot[3] > 0) {
	ktemp = ctot[1] + 2 + ctot[2];
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1], 
		ldu2, &q[ktemp + q_dim1], ldq, &zero, &u[u_dim1 + 1], ldu);
    } else {
	PLUMED_BLAS_F77_FUNC(slacpy,SLACPY)("F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu);
    }
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
    ktemp = ctot[1] + 2;
    ctemp = ctot[2] + ctot[3];
    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", nr, k, &ctemp, &one, &u2[nlp2 + ktemp * u2_dim1], ldu2,
	     &q[ktemp + q_dim1], ldq, &zero, &u[nlp2 + u_dim1], ldu);

L100:
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(k, &vt[i__ * vt_dim1 + 1], &c__1);
	q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
	}
    }

    if (*k == 2) {
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", k, &m, k, &one, &q[q_offset], ldq, &vt2[vt2_offset]
		, ldvt2, &zero, &vt[vt_offset], ldvt);
	return;
    }
    ktemp = ctot[1] + 1;
    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", k, &nlp1, &ktemp, &one, &q[q_dim1 + 1], ldq, &vt2[
	    vt2_dim1 + 1], ldvt2, &zero, &vt[vt_dim1 + 1], ldvt);
    ktemp = ctot[1] + 2 + ctot[2];
    if (ktemp <= *ldvt2) {
	PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", k, &nlp1, &ctot[3], &one, &q[ktemp * q_dim1 + 1], 
		ldq, &vt2[ktemp + vt2_dim1], ldvt2, &one, &vt[vt_dim1 + 1], 
		ldvt);
    }

    ktemp = ctot[1] + 1;
    nrp1 = *nr + *sqre;
    if (ktemp > 1) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
	}
    }
    ctemp = ctot[2] + 1 + ctot[3];
    PLUMED_BLAS_F77_FUNC(sgemm,SGEMM)("N", "N", k, &nrp1, &ctemp, &one, &q[ktemp * q_dim1 + 1], ldq, &
	    vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &zero, &vt[nlp2 * vt_dim1 + 
	    1], ldvt);

    return;


}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd4,SLASD4)(int *n, 
	int *i__, 
	float *d__, 
	float *z__, 
	float *delta, 
	float *rho, 
	float *sigma, 
	float *work, 
	int *info)
{
    int i__1;
    float d__1;

    float a, b, c__;
    int j;
    float w, dd[3];
    int ii;
    float dw, zz[3];
    int ip1;
    float eta, phi, eps, tau, psi;
    int iim1, iip1;
    float dphi, dpsi;
    int iter;
    float temp, prew, sg2lb, sg2ub, temp1, temp2, dtiim, delsq, 
	    dtiip;
    int niter;
    float dtisq;
    int swtch;
    float dtnsq;
    float delsq2, dtnsq1;
    int swtch3;
    int orgati;
    float erretm, dtipsq, rhoinv;

    --work;
    --delta;
    --z__;
    --d__;

    *info = 0;
    if (*n == 1) {

	*sigma =  std::sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
	delta[1] = 1.;
	work[1] = 1.;
	return;
    }
    if (*n == 2) {
	PLUMED_BLAS_F77_FUNC(slasd5,SLASD5)(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
	return;
    }

    eps = PLUMED_GMX_FLOAT_EPS;
    rhoinv = 1. / *rho;

    if (*i__ == *n) {

	ii = *n - 1;
	niter = 1;

	temp = *rho / 2.;

	temp1 = temp / (d__[*n] +  std::sqrt(d__[*n] * d__[*n] + temp));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] = d__[j] + d__[*n] + temp1;
	    delta[j] = d__[j] - d__[*n] - temp1;
	}

	psi = 0.;
	i__1 = *n - 2;
	for (j = 1; j <= i__1; ++j) {
	    psi += z__[j] * z__[j] / (delta[j] * work[j]);
	}

	c__ = rhoinv + psi;
	w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) + z__[*n] * z__[*
		n] / (delta[*n] * work[*n]);

	if (w <= 0.) {
	    temp1 =  std::sqrt(d__[*n] * d__[*n] + *rho);
	    temp = z__[*n - 1] * z__[*n - 1] / ((d__[*n - 1] + temp1) * (d__[*
		    n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) + z__[*n] * 
		    z__[*n] / *rho;

	    if (c__ <= temp) {
		tau = *rho;
	    } else {
		delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
		a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*
			n];
		b = z__[*n] * z__[*n] * delsq;
		if (a < 0.) {
		    tau = b * 2. / ( std::sqrt(a * a + b * 4. * c__) - a);
		} else {
		    tau = (a +  std::sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
		}
	    }

	} else {
	    delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
	    a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
	    b = z__[*n] * z__[*n] * delsq;

	    if (a < 0.) {
		tau = b * 2. / ( std::sqrt(a * a + b * 4. * c__) - a);
	    } else {
		tau = (a +  std::sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
	    }

	}

	eta = tau / (d__[*n] +  std::sqrt(d__[*n] * d__[*n] + tau));

	*sigma = d__[*n] + eta;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    delta[j] = d__[j] - d__[*i__] - eta;
	    work[j] = d__[j] + d__[*i__] + eta;
	}

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = ii;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (delta[j] * work[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	temp = z__[*n] / (delta[*n] * work[*n]);
	phi = z__[*n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (dpsi 
		+ dphi);

	w = rhoinv + phi + psi;

	if (std::abs(w) <= eps * erretm) {
	    goto L240;
	}

	++niter;
	dtnsq1 = work[*n - 1] * delta[*n - 1];
	dtnsq = work[*n] * delta[*n];
	c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
	a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
	b = dtnsq * dtnsq1 * w;
	if (c__ < 0.) {
	    c__ = std::abs(c__);
	}
	if ( std::abs(c__)<PLUMED_GMX_FLOAT_MIN) {
	    eta = *rho - *sigma * *sigma;
	} else if (a >= 0.) {
	    eta = (a +  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__  * 2.);
	} else {
	  eta = b * 2. / (a -  std::sqrt(std::abs(a * a - b * 4. * c__)));
	}

	if (w * eta > 0.) {
	    eta = -w / (dpsi + dphi);
	}
	temp = eta - dtnsq;
	if (temp > *rho) {
	    eta = *rho + dtnsq;
	}

	tau += eta;
	eta /= *sigma +  std::sqrt(eta + *sigma * *sigma);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    delta[j] -= eta;
	    work[j] += eta;
	}

	*sigma += eta;

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = ii;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	temp = z__[*n] / (work[*n] * delta[*n]);
	phi = z__[*n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (dpsi 
		+ dphi);

	w = rhoinv + phi + psi;

	iter = niter + 1;

	for (niter = iter; niter <= 20; ++niter) {

	    if (std::abs(w) <= eps * erretm) {
		goto L240;
	    }
	    dtnsq1 = work[*n - 1] * delta[*n - 1];
	    dtnsq = work[*n] * delta[*n];
	    c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
	    a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
	    b = dtnsq1 * dtnsq * w;
	    if (a >= 0.) {
		eta = (a +  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    } else {
	      eta = b * 2. / (a -  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    }

	    if (w * eta > 0.) {
		eta = -w / (dpsi + dphi);
	    }
	    temp = eta - dtnsq;
	    if (temp <= 0.) {
		eta /= 2.;
	    }

	    tau += eta;
	    eta /= *sigma +  std::sqrt(eta + *sigma * *sigma);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		delta[j] -= eta;
		work[j] += eta;
	    }

	    *sigma += eta;

	    dpsi = 0.;
	    psi = 0.;
	    erretm = 0.;
	    i__1 = ii;
	    for (j = 1; j <= i__1; ++j) {
		temp = z__[j] / (work[j] * delta[j]);
		psi += z__[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = std::abs(erretm);

	    temp = z__[*n] / (work[*n] * delta[*n]);
	    phi = z__[*n] * temp;
	    dphi = temp * temp;
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + std::abs(tau) * (
		    dpsi + dphi);

	    w = rhoinv + phi + psi;
	}

	*info = 1;
	goto L240;

    } else {

	niter = 1;
	ip1 = *i__ + 1;

	delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
	delsq2 = delsq / 2.;
	temp = delsq2 / (d__[*i__] +  std::sqrt(d__[*i__] * d__[*i__] + delsq2));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] = d__[j] + d__[*i__] + temp;
	    delta[j] = d__[j] - d__[*i__] - temp;
	}

	psi = 0.;
	i__1 = *i__ - 1;
	for (j = 1; j <= i__1; ++j) {
	    psi += z__[j] * z__[j] / (work[j] * delta[j]);
	}

	phi = 0.;
	i__1 = *i__ + 2;
	for (j = *n; j >= i__1; --j) {
	    phi += z__[j] * z__[j] / (work[j] * delta[j]);
	}
	c__ = rhoinv + psi + phi;
	w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) + z__[
		ip1] * z__[ip1] / (work[ip1] * delta[ip1]);

	if (w > 0.) {

	    orgati = 1;
	    sg2lb = 0.;
	    sg2ub = delsq2;
	    a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
	    b = z__[*i__] * z__[*i__] * delsq;
	    if (a > 0.) {
		tau = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    } else {
		tau = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    }
	    eta = tau / (d__[*i__] +  std::sqrt(d__[*i__] * d__[*i__] + tau));
	} else {

	    orgati = 0;
	    sg2lb = -delsq2;
	    sg2ub = 0.;
	    a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
	    b = z__[ip1] * z__[ip1] * delsq;
	    if (a < 0.) {
		tau = b * 2. / (a -  std::sqrt(std::abs(a * a + b * 4. * c__)));
	    } else {
		tau = -(a +  std::sqrt(std::abs(a * a + b * 4. * c__))) /	(c__ * 2.);
	    }
	    eta = tau / (d__[ip1] +  std::sqrt(std::abs(d__[ip1] * d__[ip1] + tau)));
	}

	if (orgati) {
	    ii = *i__;
	    *sigma = d__[*i__] + eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] = d__[j] + d__[*i__] + eta;
		delta[j] = d__[j] - d__[*i__] - eta;
	    }
	} else {
	    ii = *i__ + 1;
	    *sigma = d__[ip1] + eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] = d__[j] + d__[ip1] + eta;
		delta[j] = d__[j] - d__[ip1] - eta;
	    }
	}
	iim1 = ii - 1;
	iip1 = ii + 1;

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = iim1;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	dphi = 0.;
	phi = 0.;
	i__1 = iip1;
	for (j = *n; j >= i__1; --j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    phi += z__[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}

	w = rhoinv + phi + psi;

	swtch3 = 0;
	if (orgati) {
	    if (w < 0.) {
		swtch3 = 1;
	    }
	} else {
	    if (w > 0.) {
		swtch3 = 1;
	    }
	}
	if (ii == 1 || ii == *n) {
	    swtch3 = 0;
	}

	temp = z__[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z__[ii] * temp;
	w += temp;
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. + 
		std::abs(tau) * dw;

	if (std::abs(w) <= eps * erretm) {
	    goto L240;
	}

	if (w <= 0.) {
	    sg2lb = (sg2lb > tau) ? sg2lb : tau;
	} else {
	    sg2ub = (sg2ub < tau) ? sg2ub : tau;
	}

	++niter;
	if (! swtch3) {
	    dtipsq = work[ip1] * delta[ip1];
	    dtisq = work[*i__] * delta[*i__];
	    if (orgati) {
		d__1 = z__[*i__] / dtisq;
		c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
	    } else {
		d__1 = z__[ip1] / dtipsq;
		c__ = w - dtisq * dw - delsq * (d__1 * d__1);
	    }
	    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
	    b = dtipsq * dtisq * w;
	    if ( std::abs(c__)<PLUMED_GMX_FLOAT_MIN) {
		if ( std::abs(a)<PLUMED_GMX_FLOAT_MIN) {
		    if (orgati) {
			a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + 
				dphi);
		    } else {
			a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				dphi);
		    }
		}
		eta = b / a;
	    } else if (a <= 0.) {
		eta = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
	    } else {
		eta = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
	    }
	} else {

	    dtiim = work[iim1] * delta[iim1];
	    dtiip = work[iip1] * delta[iip1];
	    temp = rhoinv + psi + phi;
	    if (orgati) {
		temp1 = z__[iim1] / dtiim;
		temp1 *= temp1;
		c__ = temp - dtiip * (dpsi + dphi) - (d__[iim1] - d__[iip1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
		zz[0] = z__[iim1] * z__[iim1];
		if (dpsi < temp1) {
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
		}
	    } else {
		temp1 = z__[iip1] / dtiip;
		temp1 *= temp1;
		c__ = temp - dtiim * (dpsi + dphi) - (d__[iip1] - d__[iim1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
		if (dphi < temp1) {
		    zz[0] = dtiim * dtiim * dpsi;
		} else {
		    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
		}
		zz[2] = z__[iip1] * z__[iip1];
	    }
	    zz[1] = z__[ii] * z__[ii];
	    dd[0] = dtiim;
	    dd[1] = delta[ii] * work[ii];
	    dd[2] = dtiip;
	    PLUMED_BLAS_F77_FUNC(slaed6,SLAED6)(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
	    if (*info != 0) {
		goto L240;
	    }
	}

	if (w * eta >= 0.) {
	    eta = -w / dw;
	}
	if (orgati) {
	    temp1 = work[*i__] * delta[*i__];
	    temp = eta - temp1;
	} else {
	    temp1 = work[ip1] * delta[ip1];
	    temp = eta - temp1;
	}
	if (temp > sg2ub || temp < sg2lb) {
	    if (w < 0.) {
		eta = (sg2ub - tau) / 2.;
	    } else {
		eta = (sg2lb - tau) / 2.;
	    }
	}

	tau += eta;
	eta /= *sigma +  std::sqrt(*sigma * *sigma + eta);

	prew = w;

	*sigma += eta;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work[j] += eta;
	    delta[j] -= eta;
	}

	dpsi = 0.;
	psi = 0.;
	erretm = 0.;
	i__1 = iim1;
	for (j = 1; j <= i__1; ++j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    psi += z__[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = std::abs(erretm);

	dphi = 0.;
	phi = 0.;
	i__1 = iip1;
	for (j = *n; j >= i__1; --j) {
	    temp = z__[j] / (work[j] * delta[j]);
	    phi += z__[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}

	temp = z__[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z__[ii] * temp;
	w = rhoinv + phi + psi + temp;
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. + 
		std::abs(tau) * dw;

	if (w <= 0.) {
	    sg2lb = (sg2lb > tau) ? sg2lb : tau;
	} else {
	    sg2ub = (sg2ub < tau) ? sg2ub : tau;
	}

	swtch = 0;
	if (orgati) {
	    if (-w > std::abs(prew) / 10.) {
		swtch = 1;
	    }
	} else {
	    if (w > std::abs(prew) / 10.) {
		swtch = 1;
	    }
	}

	iter = niter + 1;

	for (niter = iter; niter <= 20; ++niter) {

	    if (std::abs(w) <= eps * erretm) {
		goto L240;
	    }

	    if (! swtch3) {
		dtipsq = work[ip1] * delta[ip1];
		dtisq = work[*i__] * delta[*i__];
		if (! swtch) {
		    if (orgati) {
			d__1 = z__[*i__] / dtisq;
			c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
		    } else {
			d__1 = z__[ip1] / dtipsq;
			c__ = w - dtisq * dw - delsq * (d__1 * d__1);
		    }
		} else {
		    temp = z__[ii] / (work[ii] * delta[ii]);
		    if (orgati) {
			dpsi += temp * temp;
		    } else {
			dphi += temp * temp;
		    }
		    c__ = w - dtisq * dpsi - dtipsq * dphi;
		}
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
		b = dtipsq * dtisq * w;
		if (std::abs(c__)<PLUMED_GMX_FLOAT_MIN) {
		    if (std::abs(a)<PLUMED_GMX_FLOAT_MIN) {
			if (! swtch) {
			    if (orgati) {
				a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * 
					(dpsi + dphi);
			    } else {
				a = z__[ip1] * z__[ip1] + dtisq * dtisq * (
					dpsi + dphi);
			    }
			} else {
			    a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
			}
		    }
		    eta = b / a;
		} else if (a <= 0.) {
		  eta = (a -  std::sqrt(std::abs(a * a - b * 4. * c__))) / (c__ * 2.);
		} else {
		  eta = b * 2. / (a +  std::sqrt(std::abs(a * a - b * 4. * c__)));
		}
	    } else {

		dtiim = work[iim1] * delta[iim1];
		dtiip = work[iip1] * delta[iip1];
		temp = rhoinv + psi + phi;
		if (swtch) {
		    c__ = temp - dtiim * dpsi - dtiip * dphi;
		    zz[0] = dtiim * dtiim * dpsi;
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    if (orgati) {
			temp1 = z__[iim1] / dtiim;
			temp1 *= temp1;
			temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[
				iip1]) * temp1;
			c__ = temp - dtiip * (dpsi + dphi) - temp2;
			zz[0] = z__[iim1] * z__[iim1];
			if (dpsi < temp1) {
			    zz[2] = dtiip * dtiip * dphi;
			} else {
			    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
			}
		    } else {
			temp1 = z__[iip1] / dtiip;
			temp1 *= temp1;
			temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[
				iip1]) * temp1;
			c__ = temp - dtiim * (dpsi + dphi) - temp2;
			if (dphi < temp1) {
			    zz[0] = dtiim * dtiim * dpsi;
			} else {
			    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
			}
			zz[2] = z__[iip1] * z__[iip1];
		    }
		}
		dd[0] = dtiim;
		dd[1] = delta[ii] * work[ii];
		dd[2] = dtiip;
		PLUMED_BLAS_F77_FUNC(slaed6,SLAED6)(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
		if (*info != 0) {
		    goto L240;
		}
	    }

	    if (w * eta >= 0.) {
		eta = -w / dw;
	    }
	    if (orgati) {
		temp1 = work[*i__] * delta[*i__];
		temp = eta - temp1;
	    } else {
		temp1 = work[ip1] * delta[ip1];
		temp = eta - temp1;
	    }
	    if (temp > sg2ub || temp < sg2lb) {
		if (w < 0.) {
		    eta = (sg2ub - tau) / 2.;
		} else {
		    eta = (sg2lb - tau) / 2.;
		}
	    }

	    tau += eta;
	    eta /= *sigma +  std::sqrt(*sigma * *sigma + eta);

	    *sigma += eta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] += eta;
		delta[j] -= eta;
	    }

	    prew = w;

	    dpsi = 0.;
	    psi = 0.;
	    erretm = 0.;
	    i__1 = iim1;
	    for (j = 1; j <= i__1; ++j) {
		temp = z__[j] / (work[j] * delta[j]);
		psi += z__[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = std::abs(erretm);

	    dphi = 0.;
	    phi = 0.;
	    i__1 = iip1;
	    for (j = *n; j >= i__1; --j) {
		temp = z__[j] / (work[j] * delta[j]);
		phi += z__[j] * temp;
		dphi += temp * temp;
		erretm += phi;
	    }

	    temp = z__[ii] / (work[ii] * delta[ii]);
	    dw = dpsi + dphi + temp * temp;
	    temp = z__[ii] * temp;
	    w = rhoinv + phi + psi + temp;
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + std::abs(temp) * 3. 
		    + std::abs(tau) * dw;
	    if (w * prew > 0. && std::abs(w) > std::abs(prew) / 10.) {
		swtch = ! swtch;
	    }

	    if (w <= 0.) {
		sg2lb = (sg2lb > tau) ? sg2lb : tau;
	    } else {
		sg2ub = (sg2ub < tau) ? sg2ub : tau;
	    }
	}

	*info = 1;

    }

L240:
    return;

} 
}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd5,SLASD5)(int *i__, 
	float *d__, 
	float *z__, 
	float *delta, 
	float *rho, 
	float *dsigma, 
	float *work)
{
    float b, c__, w, del, tau, delsq;

    --work;
    --delta;
    --z__;
    --d__;

    del = d__[2] - d__[1];
    delsq = del * (d__[2] + d__[1]);
    if (*i__ == 1) {
	w = *rho * 4. * (z__[2] * z__[2] / (d__[1] + d__[2] * 3.) - z__[1] * 
		z__[1] / (d__[1] * 3. + d__[2])) / del + 1.;
	if (w > 0.) {
	    b = delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[1] * z__[1] * delsq;

	    tau = c__ * 2. / (b +  std::sqrt(std::abs(b * b - c__ * 4.)));

	    tau /= d__[1] +  std::sqrt(d__[1] * d__[1] + tau);
	    *dsigma = d__[1] + tau;
	    delta[1] = -tau;
	    delta[2] = del - tau;
	    work[1] = d__[1] * 2. + tau;
	    work[2] = d__[1] + tau + d__[2];
	} else {
	    b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[2] * z__[2] * delsq;

	    if (b > 0.) {
		tau = c__ * -2. / (b +  std::sqrt(b * b + c__ * 4.));
	    } else {
		tau = (b -  std::sqrt(b * b + c__ * 4.)) / 2.;
	    }

	    tau /= d__[2] +  std::sqrt(std::abs(d__[2] * d__[2] + tau));
	    *dsigma = d__[2] + tau;
	    delta[1] = -(del + tau);
	    delta[2] = -tau;
	    work[1] = d__[1] + tau + d__[2];
	    work[2] = d__[2] * 2. + tau;
	}
    } else {

	b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	c__ = *rho * z__[2] * z__[2] * delsq;

	if (b > 0.) {
	    tau = (b +  std::sqrt(b * b + c__ * 4.)) / 2.;
	} else {
	    tau = c__ * 2. / (-b +  std::sqrt(b * b + c__ * 4.));
	}
	tau /= d__[2] +  std::sqrt(d__[2] * d__[2] + tau);
	*dsigma = d__[2] + tau;
	delta[1] = -(del + tau);
	delta[2] = -tau;
	work[1] = d__[1] + tau + d__[2];
	work[2] = d__[2] * 2. + tau;
    }
    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd6,SLASD6)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	float *d__, 
	float *vf, 
	float *vl, 
	float *alpha, 
	float *beta, 
	int *idxq, 
	int *perm, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	float *givnum,
	int *ldgnum, 
	float *poles, 
	float *difl, 
	float *difr, 
	float *z__, 
	int *k, 
	float *c__, 
	float *s, 
	float *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, i__1;
    float d__1, d__2;

    int i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    int isigma;
    float orgnrm;
    int c__0 = 0;
    float one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;

    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;

    d__1 = std::abs(*alpha); 
    d__2 = std::abs(*beta);
    orgnrm = (d__1 > d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = std::abs(d__[i__]);
	if (d__1 > orgnrm)
	    orgnrm = d__1;
    }
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    PLUMED_BLAS_F77_FUNC(slasd7,SLASD7)(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], &
	    work[ivfw], &vl[1], &work[ivlw], alpha, beta, &work[isigma], &
	    iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[
	    givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, 
	    info);

    PLUMED_BLAS_F77_FUNC(slasd8,SLASD8)(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], 
	    ldgnum, &work[isigma], &work[iw], info);

    if (*icompq == 1) {
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
	PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }

    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = *k;
    n2 = n - *k;
    PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}


}
}
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd7,SLASD7)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	int *k, 
	float *d__, 
	float *z__, 
	float *zw, 
	float *vf, 
	float *vfw,
	float *vl, 
	float *vlw,
	float *alpha, 
	float *beta,
	float *dsigma, 
	int *idx, 
	int *idxp,
	int *idxq, 
	int *perm, 
	int *givptr,
	int *givcol, 
	int *ldgcol, 
	float *givnum,
	int *ldgnum, 
	float *c__, 
	float *s, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    float d__1, d__2;

    int i__, j, m, n, k2;
    float z1;
    int jp;
    float eps, tau, tol;
    int nlp1, nlp2, idxi, idxj;
    int idxjp;
    int jprev = 0;
    float hlftol;
    int c__1 = 1;

    --d__;
    --z__;
    --zw;
    --vf;
    --vfw;
    --vl;
    --vlw;
    --dsigma;
    --idx;
    --idxp;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*icompq == 1) {
	*givptr = 0;
    }

    z1 = *alpha * vl[nlp1];
    vl[nlp1] = 0.;
    tau = vf[nlp1];
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vl[i__];
	vl[i__] = 0.;
	vf[i__ + 1] = vf[i__];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }
    vf[1] = tau;

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vf[i__];
	vf[i__] = 0.;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	zw[i__] = z__[idxq[i__]];
	vfw[i__] = vf[idxq[i__]];
	vlw[i__] = vl[idxq[i__]];
    }

    PLUMED_BLAS_F77_FUNC(slamrg,SLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = zw[idxi];
	vf[i__] = vfw[idxi];
	vl[i__] = vlw[idxi];
    }

    eps = PLUMED_GMX_FLOAT_EPS;

    d__1 = std::abs(*alpha);
    d__2 = std::abs(*beta);
    tol = (d__1>d__2) ? d__1 : d__2;
    d__2 = std::abs(d__[n]);
    tol = eps * 64. * ((d__2>tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (std::abs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    if (j == n) {
		goto L100;
	    }
	} else {
	    jprev = j;
	    goto L70;
	}
    }
L70:
    j = jprev;
L80:
    ++j;
    if (j > n) {
	goto L90;
    }
    if (std::abs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
    } else {

	if (std::abs(d__[j] - d__[jprev]) <= tol) {

	    *s = z__[jprev];
	    *c__ = z__[j];

	    tau = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(c__, s);
	    z__[j] = tau;
	    z__[jprev] = 0.;
	    *c__ /= tau;
	    *s = -(*s) / tau;


	    if (*icompq == 1) {
		++(*givptr);
		idxjp = idxq[idx[jprev] + 1];
		idxj = idxq[idx[j] + 1];
		if (idxjp <= nlp1) {
		    --idxjp;
		}
		if (idxj <= nlp1) {
		    --idxj;
		}
		givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
		givcol[*givptr + givcol_dim1] = idxj;
		givnum[*givptr + (givnum_dim1 << 1)] = *c__;
		givnum[*givptr + givnum_dim1] = *s;
	    }
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
	    PLUMED_BLAS_F77_FUNC(srot,SROT)(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    zw[*k] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L80;
L90:

    ++(*k);
    zw[*k] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L100:

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	vfw[j] = vf[jp];
	vlw[j] = vl[jp];
    }
    if (*icompq == 1) {
	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    jp = idxp[j];
	    perm[j] = idxq[idx[jp] + 1];
	    if (perm[j] <= nlp1) {
		--perm[j];
	    }
	}
    }
    i__1 = n - *k;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (std::abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    *c__ = 1.;
	    *s = 0.;
	    z__[1] = tol;
	} else {
	    *c__ = z1 / z__[1];
	    *s = -z__[m] / z__[1];
	}
	PLUMED_BLAS_F77_FUNC(srot,SROT)(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
	PLUMED_BLAS_F77_FUNC(srot,SROT)(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
    } else {
	if (std::abs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &zw[2], &c__1, &z__[2], &c__1);
    i__1 = n - 1;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
    i__1 = n - 1;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &vlw[2], &c__1, &vl[2], &c__1);

    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasd8,SLASD8)(int *icompq, 
	int *k, 
	float *d__, 
     	float *z__, 
	float *vf, 
	float *vl, 
	float *difl, 
	float *difr, 
	int *lddifr, 
	float *dsigma, 
	float *work, 
	int *info)
{
    int difr_dim1, difr_offset, i__1, i__2;
    float d__2;

    int i__, j;
    float dj, rho;
    int iwk1, iwk2, iwk3;
    float temp;
    int iwk2i, iwk3i;
    float diflj, difrj, dsigj;
    float dsigjp;
    int c__1 = 1;
    int c__0 = 0;
    float one = 1.;

    /* avoid warnings on high gcc optimization levels */
    difrj = dsigjp = 0;

     --d__;
    --z__;
    --vf;
    --vl;
    --difl;
    difr_dim1 = *lddifr;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    --dsigma;
    --work;

    *info = 0;

    if (*k == 1) {
	d__[1] = std::abs(z__[1]);
	difl[1] = d__[1];
	if (*icompq == 1) {
	    difl[2] = 1.;
	    difr[(difr_dim1 << 1) + 1] = 1.;
	}
	return;
    }

    iwk1 = 1;
    iwk2 = iwk1 + *k;
    iwk3 = iwk2 + *k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;

    rho = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(k, &z__[1], &c__1);
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;

    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", k, &c__1, &one, &one, &work[iwk3], k);

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	PLUMED_BLAS_F77_FUNC(slasd4,SLASD4)(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[
		iwk2], info);

	if (*info != 0) {
	    return;
	}
	work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
	difl[j] = -work[j];
	difr[j + difr_dim1] = -work[j + 1];
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__2 =  std::sqrt(std::abs(work[iwk3i + i__]));
	z__[i__] = (z__[i__] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	diflj = difl[j];
	dj = d__[j];
	dsigj = -dsigma[j];
	if (j < *k) {
	    difrj = -difr[j + difr_dim1];
	    dsigjp = -dsigma[j + 1];
	}
	work[j] = -z__[j] / diflj / (dsigma[j] + dj);
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  work[i__] = z__[i__] / (dsigma[i__] + dsigj - diflj) / ( dsigma[i__] + dj);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    work[i__] = z__[i__] / (dsigma[i__] + dsigjp - difrj) / (dsigma[i__] + dj);
	}
	temp = PLUMED_BLAS_F77_FUNC(snrm2,SNRM2)(k, &work[1], &c__1);
	work[iwk2i + j] = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(k, &work[1], &c__1, &vf[1], &c__1) / temp;
	work[iwk3i + j] = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(k, &work[1], &c__1, &vl[1], &c__1) / temp;
	if (*icompq == 1) {
	    difr[j + (difr_dim1 << 1)] = temp;
	}
    }

    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &work[iwk2], &c__1, &vf[1], &c__1);
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(k, &work[iwk3], &c__1, &vl[1], &c__1);

    return;

} 
}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasda,SLASDA)(int *icompq, 
	int *smlsiz, 
	int *n, 
	int *sqre, 
	float *d__, 
	float *e, 
	float *u, 
	int *ldu, 
	float *vt, 
	int *k, 
	float *difl, 
	float *difr, 
	float *z__, 
	float *poles, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	int *perm, 
	float *givnum, 
	float *c__, 
	float *s, 
	float *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, 
	    difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, 
	    z_dim1, z_offset, i__1, i__2;

    int i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, 
	    nlf, nrf, vfi, iwk, vli, lvl, nru, ndb1, nlp1, lvl2, nrp1;
    float beta;
    int idxq, nlvl;
    float alpha;
    int inode, ndiml, ndimr, idxqi, itemp;
    int sqrei;
    int nwork1, nwork2;
    int smlszp;
    int c__0 = 0;
    float zero = 0.0;
    float one = 1.;
    int c__1 = 1;
    --d__;
    --e;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --work;
    --iwork;
    *info = 0;

    m = *n + *sqre;

    if (*n <= *smlsiz) {
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		    vt_offset], ldu, &u[u_offset], ldu, &u[u_offset], ldu, &
		    work[1], info);
	} else {
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldu, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], 
		    info);
	}
	return;
    }

    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;

    ncc = 0;
    nru = 0;

    smlszp = *smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;

    PLUMED_BLAS_F77_FUNC(slasdt,SLASDT)(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
	i1 = i__ - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nlf = ic - nl;
	nrf = ic + 1;
	idxqi = idxq + nlf - 2;
	vfi = vf + nlf - 1;
	vli = vl + nlf - 1;
	sqrei = 1;
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nlp1, &nlp1, &zero, &one, &work[nwork1], &smlszp);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &
		    work[nwork1], &smlszp, &work[nwork2], &nl, &work[nwork2], 
		    &nl, &work[nwork2], info);
	    itemp = nwork1 + nl * smlszp;
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nl, &nl, &zero, &one, &u[nlf + u_dim1], ldu);
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nlp1, &nlp1, &zero, &one, &vt[nlf + vt_dim1], 
		    ldu);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &
		    vt[nlf + vt_dim1], ldu, &u[nlf + u_dim1], ldu, &u[nlf + 
		    u_dim1], ldu, &work[nwork1], info);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
	if (i__ == nd && *sqre == 0) {
	    sqrei = 0;
	} else {
	    sqrei = 1;
	}
	idxqi += nlp1;
	vfi += nlp1;
	vli += nlp1;
	nrp1 = nr + sqrei;
	if (*icompq == 0) {
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nrp1, &nrp1, &zero, &one, &work[nwork1], &smlszp);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &
		    work[nwork1], &smlszp, &work[nwork2], &nr, &work[nwork2], 
		    &nr, &work[nwork2], info);
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nr, &nr, &zero, &one, &u[nrf + u_dim1], ldu);
	    PLUMED_BLAS_F77_FUNC(slaset,SLASET)("A", &nrp1, &nrp1, &zero, &one, &vt[nrf + vt_dim1], 
		    ldu);
	    PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &
		    vt[nrf + vt_dim1], ldu, &u[nrf + u_dim1], ldu, &u[nrf + 
		    u_dim1], ldu, &work[nwork1], info);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
	    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nr;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
    }

    j = (1 << nlvl);

    for (lvl = nlvl; lvl >= 1; --lvl) {
	lvl2 = (lvl << 1) - 1;

	if (lvl == 1) {
	    lf = 1;
	    ll = 1;
	} else {
	    lf = (1 << (lvl-1));
	    ll = (lf << 1) - 1;
	}
	i__1 = ll;
	for (i__ = lf; i__ <= i__1; ++i__) {
	    im1 = i__ - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    if (i__ == ll) {
		sqrei = *sqre;
	    } else {
		sqrei = 1;
	    }
	    vfi = vf + nlf - 1;
	    vli = vl + nlf - 1;
	    idxqi = idxq + nlf - 1;
	    alpha = d__[ic];
	    beta = e[ic];
	    if (*icompq == 0) {
		PLUMED_BLAS_F77_FUNC(slasd6,SLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[
			perm_offset], &givptr[1], &givcol[givcol_offset], 
			ldgcol, &givnum[givnum_offset], ldu, &poles[
			poles_offset], &difl[difl_offset], &difr[difr_offset],
			 &z__[z_offset], &k[1], &c__[1], &s[1], &work[nwork1],
			 &iwork[iwk], info);
	    } else {
		--j;
		PLUMED_BLAS_F77_FUNC(slasd6,SLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[nlf + 
			lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * 
			givcol_dim1], ldgcol, &givnum[nlf + lvl2 * 
			givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &
			difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * 
			difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[j], 
			&s[j], &work[nwork1], &iwork[iwk], info);
	    }
	    if (*info != 0) {
		return;
	    }
	}
    }

    return;

}


}
}
#include <cctype>

#include "blas/blas.h"
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasdq,SLASDQ)(const char *uplo,
                        int *sqre,
                        int *n,
                        int *ncvt,
                        int *nru,
                        int *ncc,
                        float *d__,
                        float *e, 
                        float *vt, 
                        int *ldvt, 
                        float *u,
                        int *ldu, 
                        float *c__,
                        int *ldc,
                        float *work, 
                        int *info)
{
    const char xuplo=std::toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    int c__1 = 1;
    int itmp1,itmp2;
    int i__, j;
    float r__, cs, sn;
    int np1, isub;
    float smin;
    int sqre1;
    int iuplo;
    int rotate;

    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    iuplo = 0;
    if (xuplo == 'U') {
	iuplo = 1;
    }
    if (xuplo == 'L') {
	iuplo = 2;
    }
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    if (iuplo == 0) {
	*info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ncvt < 0) {
	*info = -4;
    } else if (*nru < 0) {
	*info = -5;
    } else if (*ncc < 0) {
	*info = -6;
    } else if ( (*ncvt == 0 && *ldvt < 1) || (*ncvt > 0 && *ldvt < itmp1)) {
	*info = -10;
    } else if (*ldu < itmp2) {
	*info = -12;
    } else if ((*ncc == 0 && *ldc < 1) || (*ncc > 0 && *ldc < itmp1)) {
	*info = -14;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    np1 = *n + 1;
    sqre1 = *sqre;

    if (iuplo == 1 && sqre1 == 1) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}
	PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	d__[*n] = r__;
	e[*n] = 0.f;
	if (rotate) {
	    work[*n] = cs;
	    work[*n + *n] = sn;
	}
	iuplo = 2;
	sqre1 = 0;

	if (*ncvt > 0) {
	    PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &np1, ncvt, &work[1], &work[np1], &vt[
		    vt_offset], ldvt);
	}
    }
    if (iuplo == 2) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}

	if (sqre1 == 1) {
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	    d__[*n] = r__;
	    if (rotate) {
		work[*n] = cs;
		work[*n + *n] = sn;
	    }
	}
	if (*nru > 0) {
	    if (sqre1 == 0) {
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", nru, n, &work[1], &work[np1], &u[
			u_offset], ldu);
	    } else {
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", nru, &np1, &work[1], &work[np1], &u[
			u_offset], ldu);
	    }
	}
	if (*ncc > 0) {
	    if (sqre1 == 0) {
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", n, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    } else {
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("L", "V", "F", &np1, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    }
	}
    }

    PLUMED_BLAS_F77_FUNC(sbdsqr,SBDSQR)("U", n, ncvt, nru, ncc, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[
	    u_offset], ldu, &c__[c_offset], ldc, &work[1], info);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = i__;
	smin = d__[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (d__[j] < smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != i__) {
	    d__[isub] = d__[i__];
	    d__[i__] = smin;
	    if (*ncvt > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[i__ + vt_dim1], 
			ldvt);
	    }
	    if (*nru > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[i__ * u_dim1 + 1]
			, &c__1);
	    }
	    if (*ncc > 0) {
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[i__ + c_dim1], ldc)
			;
	    }
	}
    }

    return;
}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slasdt,SLASDT)(int *n,
	int *lvl,
	int *nd,
	int *inode,
	int *ndiml,
	int *ndimr,
	int *msub)
{
  int maxn = (*n > 1) ? *n : 1;
  float temp;
  int i,il,ir,llst,nlvl,ncrnt;

  temp = std::log( ((float) maxn) / ((float)(*msub+1))) / std::log(2.0);
  
  *lvl = 1 + (int) temp;

  i = *n / 2;
  inode[0] = i + 1;
  ndiml[0] = i;
  ndimr[0] = *n - i - 1;
  il = -1;
  ir = 0;
  llst = 1;

  for(nlvl=1;nlvl<*lvl;nlvl++) {
    for(i=0;i<llst;i++) {
      il += 2;
      ir += 2;
      ncrnt = llst + i - 1;
      ndiml[il] = ndiml[ncrnt] / 2;
      ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
      inode[il] = inode[ncrnt] - ndimr[il] - 1;
      ndiml[ir] = ndimr[ncrnt] / 2;
      ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
      inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
    }
    llst *= 2;
  }
  *nd = llst*2 - 1;
  return;
}
}
}
#include <cctype>
#include "lapack.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaset,SLASET)(const char *uplo,
	int *m,
	int *n,
	float *alpha,
	float *beta,
	float *a,
	int *lda)
{
  int i,j,k;
  const char ch=std::toupper(*uplo);

  if(ch=='U') {
    for(j=1;j<*n;j++) {
      k = (j < *m) ? j : *m;
      for(i=0;i<k;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else if(ch=='L') {
    k = (*m < *n) ? *m : *n;
    for(j=0;j<k;j++) {
      for(i=j+1;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }    
  }

  k = (*m < *n) ? *m : *n;
  for(i=0;i<k;i++)
    a[i*(*lda)+i] = *beta;
}
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slasq1,SLASQ1)(int *n,
	float *d,
	float *e,
	float *work,
	int *info)
{
  float sigmx = 0.0;
  int i,j,k,iinfo;
  float minval,safemin;
  float dtemp,scale;
  float eps;

  eps = PLUMED_GMX_FLOAT_EPS;
  minval = PLUMED_GMX_FLOAT_MIN;
  safemin = minval*(1.0+PLUMED_GMX_FLOAT_EPS);
  *info = 0;

  if(*n<0) {
    *info = -2;
    return;
  }
  
  for(i=0;i<*n-1;i++) {
    d[i] = std::abs(d[i]);
    dtemp = std::abs(e[i]);
    if(dtemp>sigmx)
      sigmx=dtemp;
  }
  d[*n-1] = std::abs(d[*n-1]);
  
  if(std::abs(sigmx)<PLUMED_GMX_FLOAT_MIN) {
    PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("D",n,d,&iinfo);
    return;
  }

  for(i=0;i<*n;i++) {
    if(d[i]>sigmx)
      sigmx=d[i];
  }

  /* Copy d and e into work (z format) and scale.
   * Squaring input data makes scaling by a power of the
   * radix pointless.
   */
  scale =  std::sqrt(eps/safemin);
  i = 1;
  j = 2;
  PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n,d,&i,work,&j);
  k = *n-1;
  PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&k,e,&i,work+1,&j);
  i = 0;
  j = 2*(*n)-1;
  k = 1;
  PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G",&i,&i,&sigmx,&scale,&j,&k,work,&j,&iinfo);


  /* Compute q and e elements */
  for(i=0;i<2*(*n)-1;i++)
    work[i] = work[i]*work[i];

  work[2*(*n)-1] = 0.0;

  PLUMED_BLAS_F77_FUNC(slasq2,SLASQ2)(n,work,info);

  j = 0;
  k = 1;
  if(*info==0) {
    for(i=0;i<*n;i++)
      d[i]= std::sqrt(work[i]);
    PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G",&j,&j,&scale,&sigmx,n,&k,d,n,&iinfo);
  }
  return;
}
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#ifdef _MSC_VER
#pragma warning(disable: 4723) /*division by zero - is used on purpose here*/
#endif

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasq2,SLASQ2)(int *n, 
                        float *z__, 
                        int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2;

    float d__, e;
    int k;
    float s, t;
    int i0, i4, n0, pp;
    float dee, eps, tol;
    int ipn4;
    float tol2;
    int ieee;
    int nbig;
    float dmin__, emin, emax;
    int kmin, ndiv, iter;
    float qmin, temp, qmax, zmax;
    int splt, nfail;
    float desig, trace, sigma;
    int iinfo;
    float deemin;
    int iwhila, iwhilb;
    float oldemn, safmin, minval;
    float posinf,neginf,negzro,newzro;
    float zero = 0.0;
    float one = 1.0;

    --z__;

    *info = 0;
    eps = PLUMED_GMX_FLOAT_EPS;
    minval = PLUMED_GMX_FLOAT_MIN;
    safmin = minval*(1.0+eps);

    tol = eps * 100.;

    d__1 = tol;
    tol2 = d__1 * d__1;

    if (*n < 0) {
	*info = -1;
	return;
    } else if (*n == 0) {
	return;
    } else if (*n == 1) {

	if (z__[1] < 0.) {
	    *info = -201;
	}
	return;
    } else if (*n == 2) {

	if (z__[2] < 0. || z__[3] < 0.) {
	    *info = -2;
	    return;
	} else if (z__[3] > z__[1]) {
	    d__ = z__[3];
	    z__[3] = z__[1];
	    z__[1] = d__;
	}
	z__[5] = z__[1] + z__[2] + z__[3];
	if (z__[2] > z__[3] * tol2) {
	    t = (z__[1] - z__[3] + z__[2]) * .5;
	    s = z__[3] * (z__[2] / t);
	    if (s <= t) {
		s = z__[3] * (z__[2] / (t * ( std::sqrt(s / t + 1.) + 1.)));
	    } else {
		s = z__[3] * (z__[2] / (t +  std::sqrt(t) * std::sqrt(t + s)));
	    }
	    t = z__[1] + (s + z__[2]);
	    z__[3] *= z__[1] / t;
	    z__[1] = t;
	}
	z__[2] = z__[3];
	z__[6] = z__[2] + z__[1];
	return;
    }
    z__[*n * 2] = 0.;
    emin = z__[2];
    qmax = 0.;
    zmax = 0.;
    d__ = 0.;
    e = 0.;

    i__1 = 2*(*n - 1);
    for (k = 1; k <= i__1; k += 2) {
	if (z__[k] < 0.) {
	    *info = -(k + 200);
	    return;
	} else if (z__[k + 1] < 0.) {
	    *info = -(k + 201);
	    return;
	}
	d__ += z__[k];
	e += z__[k + 1];
	d__1 = qmax, d__2 = z__[k];
	qmax = (d__1>d__2) ? d__1 : d__2;
	d__1 = emin, d__2 = z__[k + 1];
	emin = (d__1<d__2) ? d__1 : d__2;
	d__1 = (qmax>zmax) ? qmax : zmax;
	d__2 = z__[k + 1];
	zmax = (d__1>d__2) ? d__1 : d__2;
    }
    if (z__[(*n << 1) - 1] < 0.) {
	*info = -((*n << 1) + 199);
	return;
    }
    d__ += z__[(*n << 1) - 1];
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
    qmax = (d__1>d__2) ? d__1 : d__2;

    if (std::abs(e)<PLUMED_GMX_FLOAT_MIN) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    z__[k] = z__[(k << 1) - 1];
	}
	PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("D", n, &z__[1], &iinfo);
	z__[(*n << 1) - 1] = d__;
	return;
    }

    trace = d__ + e;

    if (std::abs(trace)<PLUMED_GMX_FLOAT_MIN) {
	z__[(*n << 1) - 1] = 0.;
	return;
    }

    ieee = 1;
    posinf = one/zero;
    if(posinf<=1.0)
      ieee = 0;
    neginf = -one/zero;
    if(neginf>=0.0)
      ieee = 0;
    negzro = one/(neginf+one);
    if(std::abs(negzro)>PLUMED_GMX_FLOAT_MIN)
      ieee = 0;
    neginf = one/negzro;
    if(neginf>=0)
      ieee = 0;
    newzro = negzro + zero;
    if(std::abs(newzro-zero)>PLUMED_GMX_FLOAT_MIN)
      ieee = 0;
    posinf = one /newzro;
    if(posinf<=one)
      ieee = 0;
    neginf = neginf*posinf;
    if(neginf>=zero)
      ieee = 0;
    posinf = posinf*posinf;
    if(posinf<=1.0)
      ieee = 0;

    for (k = *n << 1; k >= 2; k += -2) {
	z__[k * 2] = 0.;
	z__[(k << 1) - 1] = z__[k];
	z__[(k << 1) - 2] = 0.;
	z__[(k << 1) - 3] = z__[k - 1];
    }

    i0 = 1;
    n0 = *n;

    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
	ipn4 = 4*(i0 + n0);
	i__1 = 2*(i0 + n0 - 1);
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
	    temp = z__[i4 - 3];
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
	    z__[ipn4 - i4 - 3] = temp;
	    temp = z__[i4 - 1];
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
	    z__[ipn4 - i4 - 5] = temp;
	}
    }

    pp = 0;

    for (k = 1; k <= 2; ++k) {

	d__ = z__[(n0 << 2) + pp - 3];
	i__1 = (i0 << 2) + pp;
	for (i4 = 4*(n0 - 1) + pp; i4 >= i__1; i4 += -4) {
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		d__ = z__[i4 - 3];
	    } else {
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
	    }
	}

	emin = z__[(i0 << 2) + pp + 1];
	d__ = z__[(i0 << 2) + pp - 3];
	i__1 = 4*(n0 - 1) + pp;
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		z__[i4 - (pp << 1) - 2] = d__;
		z__[i4 - (pp << 1)] = 0.;
		d__ = z__[i4 + 1];
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
	    }
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	z__[(n0 << 2) - pp - 2] = d__;


	qmax = z__[(i0 << 2) - pp - 2];
	i__1 = (n0 << 2) - pp - 2;
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
	    d__1 = qmax, d__2 = z__[i4];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	}

	pp = 1 - pp;
    }

    iter = 2;
    nfail = 0;
    ndiv = 2*(n0 - i0);

    i__1 = *n + 1;
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
	if (n0 < 1) {
	    goto L170;
	}

	desig = 0.;
	if (n0 == *n) {
	    sigma = 0.;
	} else {
	    sigma = -z__[(n0 << 2) - 1];
	}
	if (sigma < 0.) {
	    *info = 1;
	    return;
	}

	emax = 0.;
	if (n0 > i0) {
	    emin = std::abs(z__[(n0 << 2) - 5]);
	} else {
	    emin = 0.;
	}
	qmin = z__[(n0 << 2) - 3];
	qmax = qmin;
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
	    if (z__[i4 - 5] <= 0.) {
		goto L100;
	    }
	    if (qmin >= emax * 4.) {
		d__1 = qmin, d__2 = z__[i4 - 3];
		qmin = (d__1<d__2) ? d__1 : d__2;
		d__1 = emax, d__2 = z__[i4 - 5];
		emax = (d__1>d__2) ? d__1 : d__2;
	    }
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	    d__1 = emin, d__2 = z__[i4 - 5];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	i4 = 4;

L100:
	i0 = i4 / 4;
	pp = 0;

	if (n0 - i0 > 1) {
	    dee = z__[(i0 << 2) - 3];
	    deemin = dee;
	    kmin = i0;
	    i__2 = (n0 << 2) - 3;
	    for (i4 = (i0 << 2) - 3; i4 <= i__2; i4 += 4) {
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
		if (dee <= deemin) {
		    deemin = dee;
		    kmin = (i4 + 3) / 4;
		}
	    }
	    if (2*(kmin - i0) < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
		ipn4 = 4*(i0 + n0);
		pp = 2;
		i__2 = 2*(i0 + n0 - 1);
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
		    temp = z__[i4 - 3];
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
		    z__[ipn4 - i4 - 3] = temp;
		    temp = z__[i4 - 2];
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
		    z__[ipn4 - i4 - 2] = temp;
		    temp = z__[i4 - 1];
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
		    z__[ipn4 - i4 - 5] = temp;
		    temp = z__[i4];
		    z__[i4] = z__[ipn4 - i4 - 4];
		    z__[ipn4 - i4 - 4] = temp;
		}
	    }
	}


	d__1 = 0., d__2 = qmin -  std::sqrt(qmin) * 2. * std::sqrt(emax);
	dmin__ = -((d__1>d__2) ? d__1 : d__2);

	nbig = (n0 - i0 + 1) * 30;
	i__2 = nbig;
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
	    if (i0 > n0) {
		goto L150;
	    }

	    PLUMED_BLAS_F77_FUNC(slasq3,SLASQ3)(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee);

	    pp = 1 - pp;

	    if (pp == 0 && n0 - i0 >= 3) {
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
		    splt = i0 - 1;
		    qmax = z__[(i0 << 2) - 3];
		    emin = z__[(i0 << 2) - 1];
		    oldemn = z__[i0 * 4];
		    i__3 = 4*(n0 - 3);
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
			    z__[i4 - 1] = -sigma;
			    splt = i4 / 4;
			    qmax = 0.;
			    emin = z__[i4 + 3];
			    oldemn = z__[i4 + 4];
			} else {
			    d__1 = qmax, d__2 = z__[i4 + 1];
			    qmax = (d__1>d__2) ? d__1 : d__2;
			    d__1 = emin, d__2 = z__[i4 - 1];
			    emin = (d__1<d__2) ? d__1 : d__2;
			    d__1 = oldemn, d__2 = z__[i4];
			    oldemn = (d__1<d__2) ? d__1 : d__2;
			}
		    }
		    z__[(n0 << 2) - 1] = emin;
		    z__[n0 * 4] = oldemn;
		    i0 = splt + 1;
		}
	    }
	}

	*info = 2;
	return;

L150:
	;
    }

    *info = 3;
    return;


L170:

    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	z__[k] = z__[(k << 2) - 3];
    }

    PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("D", n, &z__[1], &iinfo);

    e = 0.;
    for (k = *n; k >= 1; --k) {
	e += z__[k];
    }


    z__[(*n << 1) + 1] = trace;
    z__[(*n << 1) + 2] = e;
    z__[(*n << 1) + 3] = (float) iter;
    i__1 = *n;
    z__[(*n << 1) + 4] = (float) ndiv / (float) (i__1 * i__1);
    z__[(*n << 1) + 5] = nfail * 100. / (float) iter;

    return;

}



}
}
#include <cmath>
#include "real.h"

#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slasq3,SLASQ3)(int *i0, 
                        int *n0, 
                        float *z__, 
                        int *pp, 
                        float *dmin__, 
                        float *sigma,
                        float *desig,
                        float *qmax, 
                        int *nfail, 
                        int *iter, 
                        int *ndiv, 
	int *ieee)
{

    int ttype = 0;
    float dmin1 = 0.;
    float dmin2 = 0.;
    float dn = 0.;
    float dn1 = 0.;
    float dn2 = 0.;
    float tau = 0.;

    int i__1;
    float d__1, d__2;
    float s, t;
    int j4, nn;
    float eps, tol;
    int n0in, ipn4;
    float tol2, temp;
    --z__;

    n0in = *n0;
    eps = PLUMED_GMX_FLOAT_EPS;
    tol = eps * 100.;
    d__1 = tol;
    tol2 = d__1 * d__1;


L10:

    if (*n0 < *i0) {
	return;
    }
    if (*n0 == *i0) {
	goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0 == *i0 + 1) {
	goto L40;
    }

    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 
	    4] > tol2 * z__[nn - 7]) {
	goto L30;
    }

L20:

    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;

L30:

    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[
	    nn - 11]) {
	goto L50;
    }

L40:

    if (z__[nn - 3] > z__[nn - 7]) {
	s = z__[nn - 3];
	z__[nn - 3] = z__[nn - 7];
	z__[nn - 7] = s;
    }
    if (z__[nn - 5] > z__[nn - 3] * tol2) {
	t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
	s = z__[nn - 3] * (z__[nn - 5] / t);
	if (s <= t) {
	    s = z__[nn - 3] * (z__[nn - 5] / (t * ( std::sqrt(s / t + 1.) + 1.)));
	} else {
	    s = z__[nn - 3] * (z__[nn - 5] / (t +  std::sqrt(t) * std::sqrt(t + s)));
	}
	t = z__[nn - 7] + (s + z__[nn - 5]);
	z__[nn - 3] *= z__[nn - 7] / t;
	z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;

L50:
    if (*pp == 2) {
	*pp = 0;
    }

    if (*dmin__ <= 0. || *n0 < n0in) {
	if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
	    ipn4 = 4*(*i0 + *n0);
	    i__1 = 2*(*i0 + *n0 - 1);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		temp = z__[j4 - 3];
		z__[j4 - 3] = z__[ipn4 - j4 - 3];
		z__[ipn4 - j4 - 3] = temp;
		temp = z__[j4 - 2];
		z__[j4 - 2] = z__[ipn4 - j4 - 2];
		z__[ipn4 - j4 - 2] = temp;
		temp = z__[j4 - 1];
		z__[j4 - 1] = z__[ipn4 - j4 - 5];
		z__[ipn4 - j4 - 5] = temp;
		temp = z__[j4];
		z__[j4] = z__[ipn4 - j4 - 4];
		z__[ipn4 - j4 - 4] = temp;
	    }
	    if (*n0 - *i0 <= 4) {
		z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
		z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
	    }
	    d__1 = dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
	    dmin2 = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1]
		    , d__1 = ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
	    z__[(*n0 << 2) + *pp - 1] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 =
		     ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) - *pp + 4];
	    z__[(*n0 << 2) - *pp] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = *qmax;
	    d__2 = z__[(*i0 << 2) + *pp - 3];
	    d__1 = (d__1>d__2) ? d__1 : d__2;
	    d__2 = z__[(*i0 << 2) + *pp + 1];
	    *qmax = ((d__1>d__2) ? d__1 : d__2);
	    *dmin__ = -0.;
	}
    }


    PLUMED_BLAS_F77_FUNC(slasq4,SLASQ4)(i0, n0, &z__[1], pp, &n0in, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, &tau, &ttype);

L70:

    PLUMED_BLAS_F77_FUNC(slasq5,SLASQ5)(i0, n0, &z__[1], pp, &tau, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, ieee);

    *ndiv += *n0 - *i0 + 2;
    ++(*iter);

    if (*dmin__ >= 0. && dmin1 > 0.) {

	goto L90;

    } else if (*dmin__ < 0. && dmin1 > 0. && z__[4*(*n0 - 1) - *pp] < tol *
	     (*sigma + dn1) && std::abs(dn) < tol * *sigma) {

	z__[4*(*n0 - 1) - *pp + 2] = 0.;
	*dmin__ = 0.;
	goto L90;
    } else if (*dmin__ < 0.) {

	++(*nfail);
	if (ttype < -22) {

	    tau = 0.;
	} else if (dmin1 > 0.) {

	    tau = (tau + *dmin__) * (1. - eps * 2.);
	    ttype += -11;
	} else {

	    tau *= .25;
	    ttype += -12;
	}
	goto L70;
    }
    else {
        
        goto L80;
    }

L80:
    PLUMED_BLAS_F77_FUNC(slasq6,SLASQ6)(i0, n0, &z__[1], pp, dmin__, &dmin1, &dmin2, &dn, &dn1, &dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    tau = 0.;

L90:
    if (tau < *sigma) {
	*desig += tau;
	t = *sigma + *desig;
	*desig -= t - *sigma;
    } else {
	t = *sigma + tau;
	*desig = *sigma - (t - tau) + *desig;
    }
    *sigma = t;

    return;
}
}
}
#include <cmath>
#include "real.h"

#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasq4,SLASQ4)(int *i0, 
	int *n0, 
	float *z__, 
	int *pp, 
	int *n0in, 
	float *dmin__, 
	float *dmin1, 
	float *dmin2, 
	float *dn, 
	float *dn1, 
	float *dn2, 
	float *tau, 
	int *ttype)
{
    float g = 0.;
    int i__1;
    float d__1, d__2;

    float s, a2, b1, b2;
    int i4, nn, np;
    float gam, gap1, gap2;


    if (*dmin__ <= 0.) {
	*tau = -(*dmin__);
	*ttype = -1;
	return;
    }

    s = 0.0;

    nn = (*n0 << 2) + *pp;
    if (*n0in == *n0) {

	if ( std::abs(*dmin__ - *dn)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin__ + *dn) ||
         std::abs(*dmin__ - *dn1)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin__ + *dn1)) {

	    b1 =  std::sqrt(z__[nn - 3]) * std::sqrt(z__[nn - 5]);
	    b2 =  std::sqrt(z__[nn - 7]) * std::sqrt(z__[nn - 9]);
	    a2 = z__[nn - 7] + z__[nn - 5];

        if ( std::abs(*dmin__ - *dn)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin__ + *dn) &&
             std::abs(*dmin1 - *dn1)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin1 + *dn1)) {

            gap2 = *dmin2 - a2 - *dmin2 * .25;
		if (gap2 > 0. && gap2 > b2) {
		    gap1 = a2 - *dn - b2 / gap2 * b2;
		} else {
		    gap1 = a2 - *dn - (b1 + b2);
		}
		if (gap1 > 0. && gap1 > b1) {
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -2;
		} else {
		    s = 0.;
		    if (*dn > b1) {
			s = *dn - b1;
		    }
		    if (a2 > b1 + b2) {
			d__1 = s, d__2 = a2 - (b1 + b2);
			s = (d__1<d__2) ? d__1 : d__2;
		    }
		    d__1 = s, d__2 = *dmin__ * .333;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -3;
		}
	    } else {


		*ttype = -4;
		s = *dmin__ * .25;
		if (std::abs(*dmin__ - *dn)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin__ + *dn)) {
		    gam = *dn;
		    a2 = 0.;
		    if (z__[nn - 5] > z__[nn - 7]) {
			return;
		    }
		    b2 = z__[nn - 5] / z__[nn - 7];
		    np = nn - 9;
		} else {
		    np = nn - (*pp << 1);
		    gam = *dn1;
		    if (z__[np - 4] > z__[np - 2]) {
			return;
		    }
		    a2 = z__[np - 4] / z__[np - 2];
		    if (z__[nn - 9] > z__[nn - 11]) {
			return;
		    }
		    b2 = z__[nn - 9] / z__[nn - 11];
		    np = nn - 13;
		}


		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = np; i4 >= i__1; i4 += -4) {
		    if (std::abs(b2)<PLUMED_GMX_FLOAT_MIN) {
			goto L20;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L20;
		    }
		}
L20:
		a2 *= 1.05;


		if (a2 < .563) {
		    s = gam * (1. -  std::sqrt(a2)) / (a2 + 1.);
		}
	    }
	} else if (std::abs(*dmin__ - *dn2)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin__ + *dn2)) {

	    *ttype = -5;
	    s = *dmin__ * .25;

	    np = nn - (*pp << 1);
	    b1 = z__[np - 2];
	    b2 = z__[np - 6];
	    gam = *dn2;
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
		return;
	    }
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);


	    if (*n0 - *i0 > 2) {
		b2 = z__[nn - 13] / z__[nn - 15];
		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
		    if (std::abs(b2)<PLUMED_GMX_FLOAT_MIN) {
			goto L40;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L40;
		    }
		}
L40:
		a2 *= 1.05;
	    }

	    if (a2 < .563) {
		s = gam * (1. -  std::sqrt(a2)) / (a2 + 1.);
	    }
	} else {

	    if (*ttype == -6) {
		g += (1. - g) * .333;
	    } else if (*ttype == -18) {
		g = .083250000000000005;
	    } else {
		g = .25;
	    }
	    s = g * *dmin__;
	    *ttype = -6;
	}

    } else if (*n0in == *n0 + 1) {

        if ( std::abs(*dmin1 - *dn1)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin1 + *dn1) &&
             std::abs(*dmin2 - *dn2)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin2 + *dn2)) {

	    *ttype = -7;
	    s = *dmin1 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (std::abs(b2)<PLUMED_GMX_FLOAT_MIN) {
		goto L60;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		a2 = b1;
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (((a2>b1) ? a2 : b1) * 100. < b2) {
		    goto L60;
		}
	    }
L60:
	    b2 =  std::sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
	    gap2 = *dmin2 * .5 - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
		*ttype = -8;
	    }
	} else {

	    s = *dmin1 * .25;
	    if (std::abs(*dmin1 - *dn1)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin1 + *dn1)) {
		s = *dmin1 * .5;
	    }
	    *ttype = -9;
	}

    } else if (*n0in == *n0 + 2) {

	if (std::abs(*dmin2 - *dn2)<PLUMED_GMX_FLOAT_EPS*std::abs(*dmin2 + *dn2) &&
        z__[nn - 5] * 2. < z__[nn - 7]) {
	    *ttype = -10;
	    s = *dmin2 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (std::abs(b2)<PLUMED_GMX_FLOAT_MIN) {
		goto L80;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (b1 * 100. < b2) {
		    goto L80;
		}
	    }
L80:
	    b2 =  std::sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
	    gap2 = z__[nn - 7] + z__[nn - 9] -  std::sqrt(z__[nn - 11]) * std::sqrt(z__[
		    nn - 9]) - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    s = *dmin2 * .25;
	    *ttype = -11;
	}
    } else if (*n0in > *n0 + 2) {

	s = 0.;
	*ttype = -12;
    }

    *tau = s;
    return;

}


}
}
#include <cmath>
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slasq5,SLASQ5)(int *i0, 
	int *n0,
	float *z__, 
	int *pp, 
	float *tau,
	float *dmin__, 
	float *dmin1, 
	float *dmin2, 
	float *dn,
	float *dnm1, 
	float *dnm2,
	int *ieee)
{
    int i__1;
    float d__1, d__2;

    float d__;
    int j4, j4p2;
    float emin, temp;

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4] - *tau;
    *dmin__ = d__;
    *dmin1 = -z__[j4];

    if (*ieee) {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		temp = z__[j4 + 1] / z__[j4 - 2];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4] = z__[j4 - 1] * temp;
		d__1 = z__[j4];
                if(d__1<emin)
                  emin = d__1;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		temp = z__[j4 + 2] / z__[j4 - 3];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4 - 1] = z__[j4] * temp;
		d__1 = z__[j4 - 1];
                if(d__1<emin)
                  emin = d__1;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
        if(*dn<*dmin__)
          *dmin__ = *dn;

    } else {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		    d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		    d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4 - 1];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	if (*dnm2 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	}
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	if (*dnm1 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	}
        if(*dn<*dmin__)
          *dmin__ = *dn;

    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;

}

}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasq6,SLASQ6)(int *i0, 
	int *n0, 
	float *z__, 
	int *pp, 
	float *dmin__, 
	float *dmin1, 
	float *dmin2,
	float *dn, 
	float *dnm1, 
	float *dnm2)
{
    int i__1;
    float d__1, d__2;

    /* Local variables */
    float d__;
    int j4, j4p2;
    float emin, temp;
    const float safemin = PLUMED_GMX_FLOAT_MIN*(1.0+PLUMED_GMX_FLOAT_EPS);

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;

    if (*pp == 0) {
	i__1 = 4*(*n0 - 3);
	for (j4 = *i0*4; j4 <= i__1; j4 += 4) {
	    z__[j4 - 2] = d__ + z__[j4 - 1];
	    if (std::abs(z__[j4 - 2])<PLUMED_GMX_FLOAT_MIN) {
		z__[j4] = 0.;
		d__ = z__[j4 + 1];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 1] < z__[j4 - 2] && safemin * z__[j4 
		    - 2] < z__[j4 + 1]) {
		temp = z__[j4 + 1] / z__[j4 - 2];
		z__[j4] = z__[j4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
	    }
	    if(d__<*dmin__)
	      *dmin__ = d__;

	    d__1 = emin, d__2 = z__[j4];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
    } else {
	i__1 = 4*(*n0 - 3);
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
	    z__[j4 - 3] = d__ + z__[j4];
	    if (std::abs(z__[j4 - 3])<PLUMED_GMX_FLOAT_MIN) {
		z__[j4 - 1] = 0.;
		d__ = z__[j4 + 2];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 2] < z__[j4 - 3] && safemin * z__[j4 
		    - 3] < z__[j4 + 2]) {
		temp = z__[j4 + 2] / z__[j4 - 3];
		z__[j4 - 1] = z__[j4] * temp;
		d__ *= temp;
	    } else {
		z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
	    }
	    if(d__<*dmin__)
	      *dmin__ = d__;
	    d__1 = emin, d__2 = z__[j4 - 1];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
    }

    *dnm2 = d__;
    *dmin2 = *dmin__;
    j4 = 4*(*n0 - 2) - *pp;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm2 + z__[j4p2];
    if (std::abs(z__[j4 - 2])<PLUMED_GMX_FLOAT_MIN) {
	z__[j4] = 0.;
	*dnm1 = z__[j4p2 + 2];
	*dmin__ = *dnm1;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dnm1 = *dnm2 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    if(*dnm1<*dmin__)
      *dmin__ = *dnm1;

    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (std::abs(z__[j4 - 2])<PLUMED_GMX_FLOAT_MIN) {
	z__[j4] = 0.;
	*dn = z__[j4p2 + 2];
	*dmin__ = *dn;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dn = *dnm1 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    if(*dn<*dmin__)
      *dmin__ = *dn;

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;


} 
}
}
#include <cmath>

#include "real.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasr,SLASR)(const char *side, 
       const char *pivot, 
       const char *direct, 
       int *m,
       int *n, 
       float *c__, 
       float *s, 
       float *a, 
       int *lda)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
    float temp;
    float ctemp, stemp;

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */

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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
		    if (std::abs(ctemp-1.0)>PLUMED_GMX_FLOAT_EPS || std::abs(stemp)>PLUMED_GMX_FLOAT_MIN) {
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
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)(const char *id, 
	int *n, 
	float *d__, 
	int *info)
{
    int i__1, i__2;

    int i__, j;
    float d1, d2, d3;
    int dir;
    float tmp;
    int endd;
    int stack[64];
    float dmnmx;
    int start;
    int stkpnt;

    --d__;

    *info = 0;
    dir = -1;
    if (*id=='D' || *id=='d') 
	dir = 0;
    else if (*id=='I' || *id=='i') 
	dir = 1;
   
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }
    if (*n <= 1) {
	return;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {


	if (dir == 0) {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L30;
		    }
		}
L30:
		;
	    }

	} else {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L50;
		    }
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return;

}
}
}
#include "lapack.h"
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;

void PLUMED_BLAS_F77_FUNC(slasrt2,SLASRT2)(const char *id, 
	      int *n, 
	      float *d__, 
	      int * key, 
	      int *info)
{
    int i__1, i__2;

    int i__, j;
    float d1, d2, d3;
    int dir;
    float tmp;
    int endd;
    int stack[64];
    float dmnmx;
    int start;
    int tmpkey, stkpnt;

    --key;
    --d__;

    *info = 0;
    dir = -1;
    if (*id=='D' || *id=='d')
	dir = 0;
    else if (*id=='I' || *id=='i')
	dir = 1;
    
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	return;
    }

    if (*n <= 1) {
	return;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start > 0) {

	if (dir == 0) {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			break;
		    }
		}
	    }

	} else {

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
			tmpkey = key[j];
			key[j] = key[j - 1];
			key[j - 1] = tmpkey;
		    } else {
			break;
		    }
		}
	    }

	}

    } else if (endd - start > 20) {

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		tmpkey = key[j];
		key[j] = key[i__];
		key[i__] = tmpkey;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }

    return;
}
}
}
#include <cmath>
#include "real.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slassq,SLASSQ)(int *n,
                        float *x,
                        int *incx,
                        float *scale,
                        float *sumsq)
{
  int ix;
  float absxi,t;

  if(*n>0) {
    for(ix=0;ix<=(*n-1)*(*incx);ix+=*incx) {
      if(std::abs(x[ix])>PLUMED_GMX_FLOAT_MIN) {
	absxi = std::abs(x[ix]);
	if(*scale<absxi) {
	  t = *scale/absxi;
	  t = t*t;
	  *sumsq = 1.0 + (*sumsq)*t;
	  *scale = absxi;
	} else {
	  t = absxi/(*scale);
	  *sumsq += t*t;
	}
      }
    }
  }
  return;
}
}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(slasv2,SLASV2)(float *f, 
                        float *g, 
                        float *h__, 
                        float *ssmin, 
                        float *ssmax, 
                        float *snr, 
                        float *csr, 
                        float *snl, 
                        float *csl)
{
    float d__1;

    float a, d__, l, m, r__, s, t, fa, ga, ha, ft, gt, ht, mm, tt,
	     clt, crt, slt, srt;
    int pmax;
    float temp;
    int swap;
    float tsign=1.0;
    int gasmal;

    ft = *f;
    fa = std::abs(ft);
    ht = *h__;
    ha = std::abs(*h__);

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

    }
    gt = *g;
    ga = std::abs(gt);
    if (std::abs(ga)<PLUMED_GMX_FLOAT_MIN) {

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = 1;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < PLUMED_GMX_FLOAT_EPS) {

		gasmal = 0;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

	    d__ = fa - ha;
	    if ( std::abs( fa - d__ )<PLUMED_GMX_FLOAT_EPS*std::abs( fa + d__ )) {
		l = 1.;
	    } else {
		l = d__ / fa;
	    }

	    m = gt / ft;
	    t = 2. - l;

	    mm = m * m;
	    tt = t * t;
	    s =  std::sqrt(tt + mm);

	    if ( std::abs(l)<PLUMED_GMX_FLOAT_MIN) {
		r__ = std::abs(m);
	    } else {
		r__ =  std::sqrt(l * l + mm);
	    }
	    a = (s + r__) * .5;

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if ( std::abs(mm)<PLUMED_GMX_FLOAT_MIN) {

		if (std::abs(l)<PLUMED_GMX_FLOAT_MIN) {
		    t = ( (ft>0) ? 2.0 : -2.0) * ( (gt>0) ? 1.0 : -1.0);
		} else {
		    t = gt / ( (ft>0) ? d__ : -d__) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
	    }
	    l =  std::sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

    if (pmax == 1) {
	tsign = ( (*csr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*f>0) ? 1.0 : -1.0);
    }
    if (pmax == 2) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*g>0) ? 1.0 : -1.0);
    }
    if (pmax == 3) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*snl>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    }
    if(tsign<0)
      *ssmax *= -1.0;
    d__1 = tsign * ( (*f>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    if(d__1<0)
      *ssmin *= -1.0;
    return;

}
}
}
#include "lapack.h"

/* LAPACK */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slaswp,SLASWP)(int *n,
	float *a,
	int *lda,
	int *k1,
	int *k2,
	int *ipiv,
	int *incx)
{
  int ix0,i1,i2,inc,n32;
  int ix,i,j,ip,k;
  float temp;

  if(*incx>0) {
    ix0 = *k1 - 1;
    i1 = *k1 - 1;
    i2 = *k2;
    inc = 1;
  } else if(*incx<0) {
    ix0 = *incx * (1- *k2);
    i1 = *k2 - 1;
    i2 = *k1;
    inc = -1;
  } else
    return;

  n32 = *n / 32;
  
  n32 *= 32;


  if(n32!=0) {
    for(j=0;j<n32;j+=32) {
      ix = ix0;
      for(i=i1;i<i2;i+=inc,ix+=*incx) {
	ip = ipiv[ix] - 1;
	if(ip != i) {
	  for(k=j;k<j+32;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	  }
	}
      }
    }
  }
  if(n32!=*n) {
    ix = ix0;
    for(i=i1;i<i2;i+=inc,ix+=*incx) {
      ip = ipiv[ix] - 1;
      if(ip != i) {
	for(k=n32;k<*n;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	}
      }
    }
  }
  return;
}
}
}
#include <cctype>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(slatrd,SLATRD)(const char *  uplo,
       int  *   n,
       int  *   nb,
       float * a,
       int *    lda,
       float * e,
       float * tau,
       float * w,
       int *    ldw)
{
  int i,iw;
  int ti1,ti2,ti3;
  float one,zero,minusone,alpha;
  const char ch=std::toupper(*uplo);

  one=1.0;
  minusone=-1.0;
  zero=0.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n;i>=(*n-*nb+1);i--) {
      iw = i -*n + *nb;
      
      if(i<*n) {
	ti1 = *n-i;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(a[ i*(*lda) + 0]),lda,&(w[iw*(*ldw)+(i-1)]),
	       ldw,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(w[ iw*(*ldw) + 0]),ldw,&(a[i*(*lda)+(i-1)]),
	       lda,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
      }

      if(i>1) {
	/*  Generate elementary reflector H(i) to annihilate
	 *              A(1:i-2,i) 
	 */
	ti1 = i-1;
	ti2 = 1;

	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i-2)]),&(a[(i-1)*(*lda)+0]),&ti2,&(tau[i-2]));
      
	e[i-2] = a[(i-1)*(*lda)+(i-2)];
	a[(i-1)*(*lda)+(i-2)] = 1.0;

	/* Compute W(1:i-1,i) */
	ti1 = i-1;
	ti2 = 1;

	/* BLAS */
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("U",&ti1,&one,a,lda,&(a[(i-1)*(*lda)+0]),&ti2,&zero,
	       &(w[(iw-1)*(*ldw)+0]),&ti2);
	if(i<*n) {
	  ti1 = i-1;
	  ti2 = *n-i;
	  ti3 = 1;
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[iw*(*ldw)+0]),ldw,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[i*(*lda)+0]),lda,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[i*(*lda)+0]),lda,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[iw*(*ldw)+0]),ldw,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	}
      
	ti1 = i-1;
	ti2 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-2]),&(w[(iw-1)*(*ldw)+0]),&ti2);
      
	alpha = -0.5*tau[i-2]*PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(w[(iw-1)*(*ldw)+0]),&ti2,
				    &(a[(i-1)*(*lda)+0]),&ti2);
      
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+0]),&ti2,&(w[(iw-1)*(*ldw)+0]),&ti2);

      }
    }
  } else {
    /* lower */
    for(i=1;i<=*nb;i++) {

      ti1 = *n-i+1;
      ti2 = i-1;
      ti3 = 1;
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(a[ i-1 ]),lda,&(w[ i-1 ]),
	       ldw,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);
      /* BLAS */
      PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(w[ i-1 ]),ldw,&(a[ i-1 ]),
	       lda,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);

      if(i<*n) {
	ti1 = *n - i;
	ti2 = (*n < i+2 ) ? *n : (i+2);
	ti3 = 1;
	/* LAPACK */
	PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+(ti2-1)]),&ti3,&(tau[i-1]));
	e[i-1] = a[(i-1)*(*lda)+(i)];
	a[(i-1)*(*lda)+(i)] = 1.0;
	
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("L",&ti1,&one,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),&ti2,
	       &zero,&(w[(i-1)*(*ldw)+i]),&ti2);
	ti1 = *n - i;
	ti2 = i-1;
	ti3 = 1;
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[ i ]),ldw,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[ i ]),lda,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[ i ]),lda,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	PLUMED_BLAS_F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[ i ]),ldw,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);

	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-1]),&(w[(i-1)*(*ldw)+i]),&ti3);
	alpha = -0.5*tau[i-1]*PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(w[(i-1)*(*ldw)+i]),&ti3,
				   &(a[(i-1)*(*lda)+i]),&ti3);
	
	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti3,&(w[(i-1)*(*ldw)+i]),&ti3);
      }
    }
  }
  return;
}
	


  
}
}
#include <cmath>

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sorg2r,SORG2R)(int *m, 
                        int *n,
                        int *k, 
                        float *a, 
                        int *lda,
                        float *tau,
                        float *work,
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    float r__1;
    int c__1 = 1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;

    if (*n <= 0) {
        return;
    }

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a[l + j * a_dim1] = 0.0;
	}
	a[j + j * a_dim1] = 1.0;
    }
    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    a[i__ + i__ * a_dim1] = 1.0;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__;
	    PLUMED_BLAS_F77_FUNC(slarf,SLARF)("L", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, 
                              &tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
	}
	if (i__ < *m) {
	    i__1 = *m - i__;
	    r__1 = -tau[i__];
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__1, &r__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[l + i__ * a_dim1] = 0.0;
	}
    }
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sorgbr,SORGBR)(const char *vect,
	int *m,
	int *n,
	int *k,
	float *a,
	int *lda,
	float *tau,
	float *work,
	int *lwork,
	int *info)
{
  int wantq,iinfo,j,i,i1,wrksz;
  int mn = (*m < *n) ? *m : *n;

  wantq = (*vect=='Q' || *vect=='q');

  *info = 0;
  wrksz = mn*DORGBR_BLOCKSIZE;
  if(*lwork==-1) {
    work[0] = wrksz;
    return;
  }
  
  if(*m==0 || *n==0)
    return;

  if(wantq) {
    if(*m>=*k)
      PLUMED_BLAS_F77_FUNC(sorgqr,SORGQR)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      for(j=*m;j>=2;j--) {
	a[(j-1)*(*lda)+0] = 0.0;
	for(i=j+1;i<=*m;i++)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-2)*(*lda)+(i-1)]; 
      }
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      if(*m>1) {
	i1 = *m-1;
	PLUMED_BLAS_F77_FUNC(sorgqr,SORGQR)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  } else {
    if(*k<*n)
      PLUMED_BLAS_F77_FUNC(sorglq,SORGLQ)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      for(j=2;j<=*n;j++) {
	for(i=j-1;i>=2;i--)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-1)*(*lda)+(i-2)]; 
	a[(j-1)*(*lda)+0] = 0.0;
      }
      if(*n>1) {
	i1 = *n-1;
	PLUMED_BLAS_F77_FUNC(sorglq,SORGLQ)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  }
  work[0] = wrksz;
  return;
}
 
}
}
#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sorgl2,SORGL2)(int *m,
                        int *n, 
                        int *k, 
                        float *a, 
                        int *lda, 
                        float *tau, 
                        float *work, 
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    float r__1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    i__ = (*m > 1) ? *m : 1;
    
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < i__) {
	*info = -5;
    }
    if (*info != 0) {
	return;
    }
    if (*m <= 0) {
	return;
    }

    if (*k < *m) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (l = *k + 1; l <= i__2; ++l) {
		a[l + j * a_dim1] = 0.0;
	    }
	    if (j > *k && j <= *m) {
		a[j + j * a_dim1] = 1.0;
	    }
	}
    }

    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.0;
		i__1 = *m - i__;
		i__2 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarf,SLARF)("R", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, 
               &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    }
	    i__1 = *n - i__;
	    r__1 = -tau[i__];
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__1, &r__1, &a[i__ + (i__ + 1) * a_dim1], lda);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[i__ + l * a_dim1] = 0.0;
	}
    }
    return;

}



}
}
#include "lapack.h"

#define SORGLQ_BLOCKSIZE    32
#define SORGLQ_MINBLOCKSIZE 2
#define SORGLQ_CROSSOVER    128


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sorglq,SORGLQ)(int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;

    int ldwork, lwkopt;
    int lquery;
    
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    ki = 0;
    nb = SORGLQ_BLOCKSIZE;
    lwkopt = (*m) * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*m) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

	nx = SORGLQ_CROSSOVER;
	if (nx < *k) {

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = SORGLQ_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = kk;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }
    if (kk < *m) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	PLUMED_BLAS_F77_FUNC(sorgl2,SORGL2)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *m) {

		i__2 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ - ib + 1;
		i__3 = *n - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)("Right", "Transpose", "Forward", "Rowwise", &i__2, &
			i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }

	    i__2 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(sorgl2,SORGL2)(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = i__ + ib - 1;
		for (l = i__; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (float) iws;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sorgqr,SORGQR)(int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;
    int lquery;
 
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    ki = 0;
    *info = 0;
    nb = DORGQR_BLOCKSIZE;
    lwkopt = (*n) * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*n <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

	nx = DORGQR_CROSSOVER;
	if (nx < *k) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DORGQR_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = *n;
	for (j = kk + 1; j <= i__1; ++j) {
	    i__2 = kk;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	PLUMED_BLAS_F77_FUNC(sorg2r,SORG2R)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *n) {

		i__2 = *m - i__ + 1;
		PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ + 1;
		i__3 = *n - i__ - ib + 1;
		PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
			1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &
			work[ib + 1], &ldwork);
	    }

	    i__2 = *m - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(sorg2r,SORG2R)(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ + ib - 1;
	    for (j = i__; j <= i__2; ++j) {
		i__3 = i__ - 1;
		for (l = 1; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (float) iws;
    return;

} 
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sorm2l,SORM2L)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	float *a,
	int *lda, 
	float *tau,
	float *c__,
	int *ldc, 
	float *work, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    int c__1 = 1;

    int i__, i1, i2, i3, mi, ni, nq;
    float aii;
    int left;
    int notran;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (*info != 0) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	return;
    }

    if ((left && notran) || (! left && ! notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
    } else {
	mi = *m;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - *k + i__;
	} else {

	    ni = *n - *k + i__;
	}

	aii = a[nq - *k + i__ + i__ * a_dim1];
	a[nq - *k + i__ + i__ * a_dim1] = 1.;
	PLUMED_BLAS_F77_FUNC(slarf,SLARF)(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[
		c_offset], ldc, &work[1]);
	a[nq - *k + i__ + i__ * a_dim1] = aii;
    }
    return;
}
}
}
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sorm2r,SORM2R)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc, 
	float *work, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    int i__, i1, i2, i3, ic, jc, mi, ni;
    float aii;
    int left;
    int notran;
    int c__1 = 1;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');

    ic = jc = 0;

    if (*m <= 0 || *n <= 0 || *k <= 0) {
	return;
    }

    if ((left && !notran) || (!left && notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

	    ni = *n - i__ + 1;
	    jc = i__;
	}


	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	PLUMED_BLAS_F77_FUNC(slarf,SLARF)(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1]);
	a[i__ + i__ * a_dim1] = aii;
    }
    return;

} 
}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sormbr,SORMBR)(const char *vect, 
	const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1;
 

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int notran;
    int applyq;
    char transt[1];
    int lwkopt;
    int lquery;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    applyq = (*vect=='Q' || *vect=='q');
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMQR_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (float) lwkopt;
    
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    work[1] = 1.;
    if (*m == 0 || *n == 0) {
	return;
    }

    if (applyq) {

	if (nq >= *k) {

	    PLUMED_BLAS_F77_FUNC(sormqr,SORMQR)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    PLUMED_BLAS_F77_FUNC(sormqr,SORMQR)(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
	}
    } else {

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}
	if (nq > *k) {

	    PLUMED_BLAS_F77_FUNC(sormlq,SORMLQ)(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    PLUMED_BLAS_F77_FUNC(sormlq,SORMLQ)(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo);
	}
    }
    work[1] = (float) lwkopt;
    return;


}


}
}
#include <cctype>
#include "real.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sorml2,SORML2)(const char *side,
                        const char *trans,
                        int *m,
                        int *n,
                        int *k,
                        float *a,
                        int *lda,
                        float *tau,
                        float *c,
                        int *ldc,
                        float *work,
                        int *info)
{
  const char xside=std::toupper(*side);
  const char xtrans=std::toupper(*trans);
  int i,i1,i2,i3,ni,mi,ic,jc;
  float aii;

  if(*m<=0 || *n<=0 || *k<=0)
    return;

  ic = jc = 0;

  if((xside=='L' && xtrans=='N') || (xside!='L' && xtrans!='N')) {
    i1 = 0;
    i2 = *k;
    i3 = 1;
  } else {
    i1 = *k-1;
    i2 = -1;
    i3 = -1;
  }
  
  if(xside=='L') {
    ni = *n;
    jc = 0;
  } else {
    mi = *m;
    ic = 0;
  }

  for(i=i1;i!=i2;i+=i3) {
    if(xside=='L') {
      mi = *m - i;
      ic = i;
    } else {
      ni = *n - i;
      jc = i;
    }
    aii = a[i*(*lda)+i];
    a[i*(*lda)+i] = 1.0;
    PLUMED_BLAS_F77_FUNC(slarf,SLARF)(side,&mi,&ni,&(a[i*(*lda)+i]),lda,tau+i,
	   &(c[jc*(*ldc)+ic]),ldc,work);
    a[i*(*lda)+i] = aii;
  }
  return;
}
	     
}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sormlq,SORMLQ)(const char *side, 
	const char *trans,
	int *m, 
	int *n, 
	int *k,
	float *a,
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, 
	    i__5;
  

    int i__;
    float t[4160]	/* was [65][64] */;
    int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork;
    char transt[1];
    int lwkopt;
    int lquery;
    int ldt = 65;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    ic = jc = 0;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMLQ_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (float) lwkopt;
    
    if (*info != 0) {
       	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMLQ_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {


	PLUMED_BLAS_F77_FUNC(sorml2,SORML2)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && notran) || (!left && !notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], 
		    lda, &tau[i__], t, &ldt);
	    if (left) {

		mi = *m - i__ + 1;
		ic = i__;
	    } else {

		ni = *n - i__ + 1;
		jc = i__;
	    }

	    PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ 
		    + i__ * a_dim1], lda, t, &ldt, &c__[ic + jc * c_dim1], 
		    ldc, &work[1], &ldwork);
	}
    }
    work[1] = (float) lwkopt;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sormql,SORMQL)(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *
	c__, int *ldc, float *work, int *lwork, int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;
    int c__65 = 65;

    int i__;
    float t[4160];
    int i1, i2, i3, ib, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork, lwkopt;
    int lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMQL_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (float) lwkopt;
    
    if (*info != 0) {
        return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMQL_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {

	PLUMED_BLAS_F77_FUNC(sorm2l,SORM2L)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && notran) || (! left && ! notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	} else {
	    mi = *m;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - *k + i__ + ib - 1;
	    PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], t, &c__65);
	    if (left) {

		mi = *m - *k + i__ + ib - 1;
	    } else {

		ni = *n - *k + i__ + ib - 1;
	    }

	    PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, t, &c__65, &c__[c_offset], ldc, &
		    work[1], &ldwork);
	}
    }
    work[1] = (float) lwkopt;
    return;

}


}
}
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void 
PLUMED_BLAS_F77_FUNC(sormqr,SORMQR)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc, 
	float *work, 
	int *lwork, 
	int *info)
{
   int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;

    int i__;
    float t[4160];
    int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork, lwkopt;
    int lquery;
    int ldt = 65;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

     ic = jc = 0;
     nb = DORMQR_BLOCKSIZE;
     lwkopt = nw * nb;
     work[1] = (float) lwkopt;

    if (*info != 0) {
	return;
    } else if (lquery) {
      return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMQR_MINBLOCKSIZE;
	}
    }

    if (nb < nbmin || nb >= *k) {

	PLUMED_BLAS_F77_FUNC(sorm2r,SORM2R)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && !notran) || (!left && notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(slarft,SLARFT)("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], t, &ldt);
	    if (left) {

		mi = *m - i__ + 1;
		ic = i__;
	    } else {
		ni = *n - i__ + 1;
		jc = i__;
	    }

	    PLUMED_BLAS_F77_FUNC(slarfb,SLARFB)(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, t, &ldt, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork);
	}
    }
    work[1] = (float) lwkopt;
    return;


}


}
}
#include "lapack.h"
#include "lapack_limits.h"


#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sormtr,SORMTR)(const char *side, 
	const char *uplo,
	const char *trans, 
	int *m, 
	int *n,
	float *a, 
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc,
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__2;

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int upper;
    int lwkopt;
    int lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    upper = (*uplo=='U' || *uplo=='u');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }


    nb = DORMQL_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (float) lwkopt;
    
    if (*info != 0) {
	i__2 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || nq == 1) {
	work[1] = 1.;
	return;
    }

    if (left) {
	mi = *m - 1;
	ni = *n;
    } else {
	mi = *m;
	ni = *n - 1;
    }

    if (upper) {
	i__2 = nq - 1;
	PLUMED_BLAS_F77_FUNC(sormql,SORMQL)(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &
		tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo);
    } else {
	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	i__2 = nq - 1;
	PLUMED_BLAS_F77_FUNC(sormqr,SORMQR)(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &
		c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
    }
    work[1] = (float) lwkopt;
    return;

}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sstebz,SSTEBZ)(const char *range, 
                        const char *order,
                        int *n,
                        float *vl, 
                        float *vu, 
                        int *il,
                        int *iu,
                        float *abstol, 
                        float *d__,
                        float *e,
                        int *m, 
                        int *nsplit, 
                        float *w,
                        int *iblock,
                        int *isplit,
                        float *work, 
                        int *iwork, 
                        int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2, d__3, d__4, d__5;
    int c__1 = 1;
    int c__3 = 3;
    int c__2 = 2;
    int c__0 = 0;

    int j, ib, jb, ie, je, nb;
    float gl;
    int im, in;
    float gu;
    int iw;
    float wl, wu;
    int nwl;
    float ulp, wlu, wul;
    int nwu;
    float tmp1, tmp2;
    int iend, ioff, iout, itmp1, jdisc;
    int iinfo;
    float atoli;
    int iwoff;
    float bnorm;
    int itmax;
    float wkill, rtoli, tnorm;
    int ibegin;
    int irange, idiscl;
    int idumma[1];
    int idiscu, iorder;
    int ncnvrg;
    float pivmin;
    int toofew;
    const float safemn = PLUMED_GMX_FLOAT_MIN*(1.0+PLUMED_GMX_FLOAT_EPS);

    --iwork;
    --work;
    --isplit;
    --iblock;
    --w;
    --e;
    --d__;

    *info = 0;

    if (*range=='A' || *range=='a') {
	irange = 1;
    } else if (*range=='V' || *range=='v') {
	irange = 2;
    } else if (*range=='I' || *range=='i') {
	irange = 3;
    } else {
	irange = 0;
    }

    if (*order=='B' || *order=='b') {
	iorder = 2;
    } else if (*order=='E' || *order=='e') {
	iorder = 1;
    } else {
	iorder = 0;
    }

    if (irange <= 0) {
	*info = -1;
    } else if (iorder <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (irange == 2) {
	if (*vl >= *vu) {
	    *info = -5;
	}
    } else if (irange == 3 && (*il < 1 || *il > (*n))) {
	*info = -6;
    } else if (irange == 3 && (*iu < ((*n<*il) ? *n : *il) || *iu > *n)) {
	*info = -7;
    }

    if (*info != 0) {
	return;
    }

    *info = 0;
    ncnvrg = 0;
    toofew = 0;

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (irange == 3 && *il == 1 && *iu == *n) {
	irange = 1;
    }

    ulp = 2*PLUMED_GMX_FLOAT_EPS;
    rtoli = ulp * 2.;
    nb = DSTEBZ_BLOCKSIZE;
    // cppcheck-suppress knownConditionTrueFalse
    if (nb <= 1) {
	nb = 0;
    }

    if (*n == 1) {
	*nsplit = 1;
	isplit[1] = 1;
	if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
	    *m = 0;
	} else {
	    w[1] = d__[1];
	    iblock[1] = 1;
	    *m = 1;
	}
	return;
    }

    *nsplit = 1;
    work[*n] = 0.;
    pivmin = 1.;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	d__1 = e[j - 1];
	tmp1 = d__1 * d__1;
	d__2 = ulp;
	if (std::abs(d__[j] * d__[j - 1]) * (d__2 * d__2) + safemn 
		> tmp1) {
	    isplit[*nsplit] = j - 1;
	    ++(*nsplit);
	    work[j - 1] = 0.;
	} else {
	    work[j - 1] = tmp1;
	    pivmin = (pivmin>tmp1) ? pivmin : tmp1;
	}
    }
    isplit[*nsplit] = *n;
    pivmin *= safemn;

    if (irange == 3) {

	gu = d__[1];
	gl = d__[1];
	tmp1 = 0.;

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    tmp2 =  std::sqrt(work[j]);
	    d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
	    gu = (d__1>d__2) ? d__1 : d__2;
	    d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
	    gl = (d__1<d__2) ? d__1 : d__2;
	    tmp1 = tmp2;
	}

	d__1 = gu, d__2 = d__[*n] + tmp1;
	gu = (d__1>d__2) ? d__1 : d__2;
	d__1 = gl, d__2 = d__[*n] - tmp1;
	gl = (d__1<d__2) ? d__1 : d__2;
	d__1 = std::abs(gl);
	d__2 = std::abs(gu);
	tnorm = (d__1>d__2) ? d__1 : d__2;
	gl = gl - tnorm * 2. * ulp * *n - pivmin * 4.;
	gu = gu + tnorm * 2. * ulp * *n + pivmin * 2.;

	itmax = (int) ((std::log(tnorm + pivmin) - std::log(pivmin)) / std::log(2.)) + 2;
	if (*abstol <= 0.) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = *abstol;
	}

	work[*n + 1] = gl;
	work[*n + 2] = gl;
	work[*n + 3] = gu;
	work[*n + 4] = gu;
	work[*n + 5] = gl;
	work[*n + 6] = gu;
	iwork[1] = -1;
	iwork[2] = -1;
	iwork[3] = *n + 1;
	iwork[4] = *n + 1;
	iwork[5] = *il - 1;
	iwork[6] = *iu;

	PLUMED_BLAS_F77_FUNC(slaebz,SLAEBZ)(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, 
		&d__[1], &e[1], &work[1], &iwork[5], &work[*n + 1], &work[*n 
		+ 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

	if (iwork[6] == *iu) {
	    wl = work[*n + 1];
	    wlu = work[*n + 3];
	    nwl = iwork[1];
	    wu = work[*n + 4];
	    wul = work[*n + 2];
	    nwu = iwork[4];
	} else {
	    wl = work[*n + 2];
	    wlu = work[*n + 4];
	    nwl = iwork[2];
	    wu = work[*n + 3];
	    wul = work[*n + 1];
	    nwu = iwork[3];
	}

	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
	    *info = 4;
	    return;
	}
    } else {


      /* avoid warnings for high gcc optimization */
      wlu = wul = 1.0;

	d__3 = std::abs(d__[1]) + std::abs(e[1]);
	d__4 = std::abs(d__[*n]) + std::abs(e[*n - 1]);
	tnorm = (d__3>d__4) ? d__3 : d__4;

	i__1 = *n - 1;
	for (j = 2; j <= i__1; ++j) {
	    d__4 = tnorm;
	    d__5 = std::abs(d__[j]) + std::abs(e[j - 1]) + std::abs(e[j]);
	    tnorm = (d__4>d__5) ? d__4 : d__5;
	}

	if (*abstol <= 0.) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = *abstol;
	}

	if (irange == 2) {
	    wl = *vl;
	    wu = *vu;
	} else {
	    wl = 0.;
	    wu = 0.;
	}
    }

    *m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;

    i__1 = *nsplit;
    for (jb = 1; jb <= i__1; ++jb) {
	ioff = iend;
	ibegin = ioff + 1;
	iend = isplit[jb];
	in = iend - ioff;

	if (in == 1) {

	    if (irange == 1 || wl >= d__[ibegin] - pivmin) {
		++nwl;
	    }
	    if (irange == 1 || wu >= d__[ibegin] - pivmin) {
		++nwu;
	    }
	    if (irange == 1 || ((wl < d__[ibegin] - pivmin) && (wu >= d__[ibegin] - pivmin))) {
		++(*m);
		w[*m] = d__[ibegin];
		iblock[*m] = jb;
	    }
	} else {

	    gu = d__[ibegin];
	    gl = d__[ibegin];
	    tmp1 = 0.;

	    i__2 = iend - 1;
	    for (j = ibegin; j <= i__2; ++j) {
		tmp2 = std::abs(e[j]);
		d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
		gu = (d__1>d__2) ? d__1 : d__2;
		d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
		gl = (d__1<d__2) ? d__1 : d__2;
		tmp1 = tmp2;
	    }

	    d__1 = gu, d__2 = d__[iend] + tmp1;
	    gu = (d__1>d__2) ? d__1 : d__2;
	    d__1 = gl, d__2 = d__[iend] - tmp1;
	    gl = (d__1<d__2) ? d__1 : d__2;
	    d__1 = std::abs(gl);
	    d__2 = std::abs(gu);
	    bnorm = (d__1>d__2) ? d__1 : d__2;
	    gl = gl - bnorm * 2. * ulp * in - pivmin * 2.;
	    gu = gu + bnorm * 2. * ulp * in + pivmin * 2.;

	    if (*abstol <= 0.) {
		d__1 = std::abs(gl);
		d__2 = std::abs(gu);
		atoli = ulp * ((d__1>d__2) ? d__1 : d__2);
	    } else {
		atoli = *abstol;
	    }

	    if (irange > 1) {
		if (gu < wl) {
		    nwl += in;
		    nwu += in;
		}
		gl = (gl>wl) ? gl : wl;
		gu = (gu<wu) ? gu : wu;
		if (gl >= gu) {
		}
		continue;
	    }

	    work[*n + 1] = gl;
	    work[*n + in + 1] = gu;
	    PLUMED_BLAS_F77_FUNC(slaebz,SLAEBZ)(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);

	    nwl += iwork[1];
	    nwu += iwork[in + 1];
	    iwoff = *m - iwork[1];

	    itmax = (int) ((log(gu - gl + pivmin) - log(pivmin)) / log(2.)
		    ) + 2;
	    PLUMED_BLAS_F77_FUNC(slaebz,SLAEBZ)(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);

	    i__2 = iout;
	    for (j = 1; j <= i__2; ++j) {
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

		if (j > iout - iinfo) {
		    ncnvrg = 1;
		    ib = -jb;
		} else {
		    ib = jb;
		}
		i__3 = iwork[j + in] + iwoff;
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
		    w[je] = tmp1;
		    iblock[je] = ib;
		}
	    }

	    *m += im;
	}
    }

    if (irange == 3) {
	im = 0;
	idiscl = *il - 1 - nwl;
	idiscu = nwu - *iu;

	if (idiscl > 0 || idiscu > 0) {
	    i__1 = *m;
	    for (je = 1; je <= i__1; ++je) {
		if (w[je] <= wlu && idiscl > 0) {
		    --idiscl;
		} else if (w[je] >= wul && idiscu > 0) {
		    --idiscu;
		} else {
		    ++im;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    *m = im;
	}
	if (idiscl > 0 || idiscu > 0) {

	    if (idiscl > 0) {
		wkill = wu;
		i__1 = idiscl;
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
		    iw = 0;
		    i__2 = *m;
		    for (je = 1; je <= i__2; ++je) {
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}
		    }
		    iblock[iw] = 0;
		}
	    }
	    if (idiscu > 0) {

		wkill = wl;
		i__1 = idiscu;
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
		    iw = 0;
		    i__2 = *m;
		    for (je = 1; je <= i__2; ++je) {
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}
		    }
		    iblock[iw] = 0;
		}
	    }
	    im = 0;
	    i__1 = *m;
	    for (je = 1; je <= i__1; ++je) {
		if (iblock[je] != 0) {
		    ++im;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    *m = im;
	}
	if (idiscl < 0 || idiscu < 0) {
	    toofew = 1;
	}
    }

    if (iorder == 1 && *nsplit > 1) {
	i__1 = *m - 1;
	for (je = 1; je <= i__1; ++je) {
	    ie = 0;
	    tmp1 = w[je];
	    i__2 = *m;
	    for (j = je + 1; j <= i__2; ++j) {
		if (w[j] < tmp1) {
		    ie = j;
		    tmp1 = w[j];
		}
	    }

	    if (ie != 0) {
		itmp1 = iblock[ie];
		w[ie] = w[je];
		iblock[ie] = iblock[je];
		w[je] = tmp1;
		iblock[je] = itmp1;
	    }
	}
    }

    *info = 0;
    if (ncnvrg) {
	++(*info);
    }
    if (toofew) {
	*info += 2;
    }
    return;

}


}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sstegr,SSTEGR)(const char *jobz, 
	const char *range, 
	int *n, 
	float *d__, 
	float *e, 
	float *vl, 
	float *vu, 
	int *il, 
	int *iu, 
	float *abstol, 
	int *m, 
	float *w, 
	float *z__, 
	int *ldz, 
	int *isuppz,
	float *work, 
	int *lwork, 
	int *iwork, 
	int *liwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2;
    float d__1, d__2;
    int c__1 = 1;

    int i__, j;
    int jj;
    float eps, tol, tmp, rmin, rmax;
    int itmp;
    float tnrm;
    float scale;
    int iinfo, iindw;
    int lwmin;
    int wantz;
    int iindbl;
    int valeig,alleig,indeig;
    float safmin,minval;
    float bignum;
    int iindwk, indgrs;
    float thresh;
    int iinspl, indwrk, liwmin, nsplit;
    float smlnum;
    int lquery;


    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    lquery = *lwork == -1 || *liwork == -1;
    lwmin = *n * 17;
    liwmin = *n * 10;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && (*il < 1 || *il > *n)) {
	*info = -8;
    } else if (indeig && (*iu < *il || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -14;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -17;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -19;
    }
    if (*info == 0) {
	work[1] = (float) lwmin;
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    minval = PLUMED_GMX_FLOAT_MIN;
    safmin = minval*(1.0+PLUMED_GMX_FLOAT_EPS);
    eps = PLUMED_GMX_FLOAT_EPS;
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin =  std::sqrt(smlnum);
    d__1 =  std::sqrt(bignum), d__2 = 1. / std::sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;
    scale = 1.;
    tnrm = PLUMED_BLAS_F77_FUNC(slanst,SLANST)("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin) {
	scale = rmin / tnrm;
    } else if (tnrm > rmax) {
	scale = rmax / tnrm;
    }
    if ( std::abs(scale-1.0)>PLUMED_GMX_FLOAT_EPS) {
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(n, &scale, &d__[1], &c__1);
	i__1 = *n - 1;
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__1, &scale, &e[1], &c__1);
	tnrm *= scale;
    }
    indgrs = 1;
    indwrk = (*n << 1) + 1;

    iinspl = 1;
    iindbl = *n + 1;
    iindw = (*n << 1) + 1;
    iindwk = *n * 3 + 1;

    thresh = eps * tnrm;
    PLUMED_BLAS_F77_FUNC(slarrex,SLARREX)(range, n, vl, vu, il, iu, &d__[1], &e[1], &thresh, &nsplit, &
	    iwork[iinspl], m, &w[1], &iwork[iindbl], &iwork[iindw], &work[
	    indgrs], &work[indwrk], &iwork[iindwk], &iinfo);
    
    if (iinfo != 0) {
	*info = 1;
	return;
    }

    if (wantz) {
	d__1 = *abstol, d__2 = (float) (*n) * eps;
	tol = (d__1>d__2) ? d__1 : d__2;
	PLUMED_BLAS_F77_FUNC(slarrvx,SLARRVX)(n, &d__[1], &e[1], &iwork[iinspl], m, &w[1], &iwork[iindbl], &
		iwork[iindw], &work[indgrs], &tol, &z__[z_offset], ldz, &
		isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
	if (iinfo != 0) {
	    *info = 2;
	    return;
	}
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	itmp = iwork[iindbl + j - 1];
	w[j] += e[iwork[iinspl + itmp - 1]];
    } 

    if (std::abs(scale-1.0)>PLUMED_GMX_FLOAT_EPS) {
	d__1 = 1. / scale;
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(m, &d__1, &w[1], &c__1);
    }
    if (nsplit > 1) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp) {
		    i__ = jj;
		    tmp = w[jj];
		}
	    }
	    if (i__ != 0) {
		w[i__] = w[j];
		w[j] = tmp;
		if (wantz) {
		    PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 
			    + 1], &c__1);
		    itmp = isuppz[(i__ << 1) - 1];
		    isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
		    isuppz[(j << 1) - 1] = itmp;
		    itmp = isuppz[i__ * 2];
		    isuppz[i__ * 2] = isuppz[j * 2];
		    isuppz[j * 2] = itmp;
		}
	    }
	}
    }

    work[1] = (float) lwmin;
    iwork[1] = liwmin;
    return;

} 
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

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

	onenrm = std::abs(d__[b1]) + std::abs(e[b1]);
	d__3 = onenrm;
	d__4 = std::abs(d__[bn]) + std::abs(e[bn - 1]);
	onenrm = (d__3>d__4) ? d__3 : d__4;
	i__2 = bn - 1;
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
	  d__4 = onenrm;
	  d__5 = std::abs(d__[i__]) + std::abs(e[i__ - 1]) + std::abs(e[i__]);
	    onenrm = (d__4>d__5) ? d__4 : d__5;
	}
	ortol = onenrm * .001;

	dtpcrt =  std::sqrt(.1 / blksiz);

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
		eps1 = std::abs(eps * xj);
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
	    d__3 = std::abs(work[indrv4 + blksiz]);
	    scl = blksiz * onenrm * ((d__2>d__3) ? d__2 : d__3) / PLUMED_BLAS_F77_FUNC(sasum,SASUM)(&blksiz, &work[
		    indrv1 + 1], &c__1);
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&blksiz, &scl, &work[indrv1 + 1], &c__1);

	    PLUMED_BLAS_F77_FUNC(slagts,SLAGTS)(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

	    if (jblk == 1) {
		goto L90;
	    }
	    if (std::abs(xj - xjm) > ortol) {
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
	    nrm = std::abs(work[indrv1 + jmax]);

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
#include <cmath>
#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssteqr,SSTEQR)(const char *    compz, 
                        int *     n, 
                        float *  d__, 
                        float *  e, 
                        float *  z__, 
                        int *     ldz, 
                        float *  work, 
                        int *     info)
{
    float c_b9 = 0.;
    float c_b10 = 1.;
    int c__0 = 0;
    int c__1 = 1;
    int c__2 = 2;
    int z_dim1, z_offset, i__1, i__2;
    float d__1, d__2;

    float b, c__, f, g;
    int i__, j, k, l, m;
    float p, r__, s;
    int l1, ii, mm, lm1, mm1, nm1;
    float rt1, rt2, eps;
    int lsv;
    float tst, eps2;
    int lend, jtot;
    float anorm;
    int lendm1, lendp1;
    int iscale;
    float safmin,minval;
    float safmax;
    int lendsv;
    float ssfmin;
    int nmaxit, icompz;
    float ssfmax;


    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    *info = 0;

    if (*compz=='N' || *compz=='n') {
	icompz = 0;
    } else if (*compz=='V' || *compz=='v') {
	icompz = 1;
    } else if (*compz=='I' || *compz=='i') {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < ((*n>1) ? *n : 1))) {
	*info = -6;
    }
    if (*info != 0) {
	return;
    }


    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    eps = PLUMED_GMX_FLOAT_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    minval = PLUMED_GMX_FLOAT_MIN;
    safmin = minval*(1.0+PLUMED_GMX_FLOAT_EPS);

    safmax = 1. / safmin;
    ssfmax =  std::sqrt(safmax) / 3.;
    ssfmin =  std::sqrt(safmin) / eps2;

    if (icompz == 2) {
	PLUMED_BLAS_F77_FUNC(slaset,SLASET)("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= i__1; ++m) {
	    tst = std::abs(e[m]);
	    if (std::abs(tst)<PLUMED_GMX_FLOAT_MIN) {
		goto L30;
	    }
	    if (tst <=  std::sqrt(std::abs(d__[m])) * std::sqrt(std::abs(d__[m + 1])) * eps) {
		e[m] = 0.;
		goto L30;
	    }
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    i__1 = lend - l + 1;
    anorm = PLUMED_BLAS_F77_FUNC(slanst,SLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (std::abs(anorm)<PLUMED_GMX_FLOAT_MIN) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    if (std::abs(d__[lend]) < std::abs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= i__1; ++m) {
  	        d__2 = std::abs(e[m]);
		tst = d__2 * d__2;
		if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m+ 1]) + safmin) {
		    goto L60;
		}
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L80;
	}

	if (m == l + 1) {
	    if (icompz > 0) {
		PLUMED_BLAS_F77_FUNC(slaev2,SLAEV2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
		work[l] = c__;
		work[*n - 1 + l] = s;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz);
	    } else {
		PLUMED_BLAS_F77_FUNC(slae2,SLAE2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
	    }
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l + 1] - p) / (e[l] * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l] / (g + ( (g>0) ? r__ : -r__ ) );

	s = 1.;
	c__ = 1.;
	p = 0.;

	mm1 = m - 1;
	i__1 = l;
	for (i__ = mm1; i__ >= i__1; --i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m - 1) {
		e[i__ + 1] = r__;
	    }
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = -s;
	    }
	}

	if (icompz > 0) {
	    mm = m - l + 1;
	    PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[l] = g;
	goto L40;

L80:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= i__1; --m) {
		d__2 = std::abs(e[m - 1]);
		tst = d__2 * d__2;
		if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m- 1]) + safmin) {
		    goto L110;
		}
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L130;
	}
	if (m == l - 1) {
	    if (icompz > 0) {
		PLUMED_BLAS_F77_FUNC(slaev2,SLAEV2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
		work[m] = c__;
		work[*n - 1 + m] = s;
		PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz);
	    } else {
		PLUMED_BLAS_F77_FUNC(slae2,SLAE2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
	    }
	    d__[l - 1] = rt1;
	    d__[l] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l - 1] / (g + ( (g>0) ? r__ : -r__ ));

	s = 1.;
	c__ = 1.;
	p = 0.;

	lm1 = l - 1;
	i__1 = lm1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    PLUMED_BLAS_F77_FUNC(slartg,SLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m) {
		e[i__ - 1] = r__;
	    }
	    g = d__[i__] - p;
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = s;
	    }
	}

	if (icompz > 0) {
	    mm = l - m + 1;
	    PLUMED_BLAS_F77_FUNC(slasr,SLASR)("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[lm1] = g;
	goto L90;

L130:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__])>PLUMED_GMX_FLOAT_MIN) {
	    ++(*info);
	}
    }
    goto L190;

L160:
    if (icompz == 0) {

	PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("I", n, &d__[1], info);

    } else {

	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
	    i__ = ii - 1;
	    k = i__;
	    p = d__[i__];
	    i__2 = *n;
	    for (j = ii; j <= i__2; ++j) {
		if (d__[j] < p) {
		    k = j;
		    p = d__[j];
		}
	    }
	    if (k != i__) {
		d__[k] = d__[i__];
		d__[i__] = p;
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

L190:
    return;
}


}
}
#include <cmath>
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssterf,SSTERF)(int *n, 
	float *d__, 
	float *e, 
	int *info)
{
    int i__1;
    float d__1;

    float c__;
    int i__, l, m;
    float p, r__, s;
    int l1;
    float bb, rt1, rt2, eps, rte;
    int lsv;
    float eps2, oldc;
    int lend, jtot;
    float gamma, alpha, sigma, anorm;
      int iscale;
    float oldgam;
    float safmax;
    int lendsv;
    float ssfmin;
    int nmaxit;
    float ssfmax;
    int c__0 = 0;
    int c__1 = 1;
    float c_b32 = 1.;
    const float safmin = PLUMED_GMX_FLOAT_MIN*(1.0+PLUMED_GMX_FLOAT_EPS);

    --e;
    --d__;

    *info = 0;

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	return;
    }
    if (*n <= 1) {
	return;
    }

    eps = PLUMED_GMX_FLOAT_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmax = 1. / safmin;
    ssfmax =  std::sqrt(safmax) / 3.;
    ssfmin =  std::sqrt(safmin) / eps2;

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

    l1 = 1;

L10:
    if (l1 > *n) {
      PLUMED_BLAS_F77_FUNC(slasrt,SLASRT)("I", n, &d__[1], info);
      return;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
	if (std::abs(e[m]) <=  std::sqrt(std::abs(d__[m])) * 
		 std::sqrt(std::abs(d__[m + 1])) * eps) {
	    e[m] = 0.;
	    goto L30;
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    i__1 = lend - l + 1;
    anorm = PLUMED_BLAS_F77_FUNC(slanst,SLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
	d__1 = e[i__];
	e[i__] = d__1 * d__1;
    }

    if (std::abs(d__[lend]) < std::abs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

L50:
	if (l != lend) {
	    i__1 = lend - 1;
	    for (m = l; m <= i__1; ++m) {
		if (std::abs(e[m]) <= eps2 * std::abs(d__[m] * d__[m + 1])) {
		    goto L70;
		}
	    }
	}
	m = lend;

L70:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L90;
	}
	if (m == l + 1) {
	    rte =  std::sqrt(e[l]);
	    PLUMED_BLAS_F77_FUNC(slae2,SLAE2)(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte =  std::sqrt(e[l]);
	sigma = (d__[l + 1] - p) / (rte * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l;
	for (i__ = m - 1; i__ >= i__1; --i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m - 1) {
		e[i__ + 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__ + 1] = oldgam + (alpha - gamma);
	    if (std::abs(c__)>PLUMED_GMX_FLOAT_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l] = s * p;
	d__[l] = sigma + gamma;
	goto L50;

L90:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

L100:
	i__1 = lend + 1;
	for (m = l; m >= i__1; --m) {
	    if (std::abs(e[m - 1]) <= eps2 * std::abs(d__[m] * d__[m - 1])) {
		goto L120;
	    }
	}
	m = lend;

L120:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L140;
	}

	if (m == l - 1) {
	    rte =  std::sqrt(e[l - 1]);
	    PLUMED_BLAS_F77_FUNC(slae2,SLAE2)(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l - 1] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte =  std::sqrt(e[l - 1]);
	sigma = (d__[l - 1] - p) / (rte * 2.);
	r__ = PLUMED_BLAS_F77_FUNC(slapy2,SLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l - 1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m) {
		e[i__ - 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__ + 1];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__] = oldgam + (alpha - gamma);
	    if (std::abs(c__)>PLUMED_GMX_FLOAT_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l - 1] = s * p;
	d__[l] = sigma + gamma;
	goto L100;

L140:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	PLUMED_BLAS_F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__])>PLUMED_GMX_FLOAT_MIN) {
	    ++(*info);
	}
    }
    return;
}


}
}
#include "lapack.h"


/* Normally, SSTEVR is the LAPACK wrapper which calls one
 * of the eigenvalue methods. However, our code includes a
 * version of SSTEGR which is never than LAPACK 3.0 and can
 * handle requests for a subset of eigenvalues/vectors too,
 * and it should not need to call SSTEIN.
 * Just in case somebody has a faster version in their lapack
 * library we still call the driver routine, but in our own
 * case this is just a wrapper to sstegr.
 */
#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(sstevr,SSTEVR)(const char *jobz, 
                        const char *range,
                        int *n,
                        float *d,
                        float *e,
                        float *vl, 
                        float *vu,
                        int *il, 
                        int *iu, 
                        float *abstol,
                        int *m,
                        float *w, 
                        float *z,
                        int *ldz,
                        int *isuppz, 
                        float *work, 
                        int *lwork, 
                        int *iwork,
                        int *liwork, 
                        int *info)
{
  PLUMED_BLAS_F77_FUNC(sstegr,SSTEGR)(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w,
	  z, ldz, isuppz, work, lwork, iwork, liwork, info);
  

    return;

}


}
}
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssyevr,SSYEVR)(const char *jobz, const char *range, const char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *
	il, int *iu, float *abstol, int *m, float *w, 
	float *z__, int *ldz, int *isuppz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    float d__1, d__2;

    /* Local variables */
    int c__1 = 1;
    int i__, j, nb, jj;
    float eps, tmp1;
    int indd, inde;
    float anrm;
    int imax;
    float rmin, rmax;
    int itmp1, inddd, indee;
    float sigma;
    int iinfo;
    int indwk;
    int lwmin;
    int lower, wantz;
    int alleig, indeig;
    int iscale, indibl, indifl;
    int valeig;
    float safmin,minval;
    float bignum;
    int indtau;
    int indwkn;
    int liwmin;
    int llwrkn, llwork;
    float smlnum;
    int lwkopt;
    int lquery;
    
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    lower = (*uplo=='L' || *uplo=='l');
    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    indibl = 0;
    lquery = *lwork == -1 || *liwork == -1;

    i__1 = 1;
    i__2 = *n * 26;

    if(*n>0) 
      lwmin = *n * 26;
    else
      lwmin = 1;

    if(*n>0) 
      liwmin = *n * 10;
    else
      liwmin = 1;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || (*uplo=='U' || *uplo=='u'))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < ((*n>1) ? *n : 1) ) {
	*info = -6;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -8;
	    }
	} else if (indeig) {
	  if (*il < 1 || *il > ((*n>1) ? *n : 1)) {
		*info = -9;
	    } else if (*iu < ((*n<*il) ? *n : *il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
      if (*ldz < 1 || (wantz && *ldz < *n)) {
	    *info = -15;
	} else if (*lwork < lwmin && ! lquery) {
	    *info = -18;
	} else if (*liwork < liwmin && ! lquery) {
	    *info = -20;
	}
    }

    if (*info == 0) {
      nb = 32;
      /* Computing MAX */
      i__1 = (nb + 1) * *n;
      lwkopt = (i__1>lwmin) ? i__1 : lwmin;
      work[1] = (float) lwkopt;
      iwork[1] = liwmin;
    } else 
      return;

    if (lquery) 
	return;
    
    *m = 0;
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    if (*n == 1) {
	work[1] = 7.;
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = a[a_dim1 + 1];
	} else {
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
		*m = 1;
		w[1] = a[a_dim1 + 1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }
    minval = PLUMED_GMX_FLOAT_MIN;
    safmin = minval*(1.0+PLUMED_GMX_FLOAT_EPS);
    eps = PLUMED_GMX_FLOAT_EPS;

    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin =  std::sqrt(smlnum);

    d__1 =  std::sqrt(bignum), d__2 = 1. / std::sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;

    iscale = 0;
    anrm = PLUMED_BLAS_F77_FUNC(slansy,SLANSY)("M", uplo, n, &a[a_offset], lda, &work[1]);
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm; 
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&j, &sigma, &a[j * a_dim1 + 1], &c__1);

	    }
	}
    }

    indtau = 1;
    inde = indtau + *n;
    indd = inde + *n;
    indee = indd + *n;
    inddd = indee + *n;
    indifl = inddd + *n;
    indwk = indifl + *n;
    llwork = *lwork - indwk + 1;
    PLUMED_BLAS_F77_FUNC(ssytrd,SSYTRD)(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwk], &llwork, &iinfo);

    i__1 = *n - 1;
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(&i__1, &work[inde], &c__1, &work[indee], &c__1);
    PLUMED_BLAS_F77_FUNC(scopy,SCOPY)(n, &work[indd], &c__1, &work[inddd], &c__1);

    PLUMED_BLAS_F77_FUNC(sstegr,SSTEGR)(jobz, range, n, &work[inddd], &work[indee], vl, vu, il, iu, 
	    abstol, m, &w[1], &z__[z_offset], ldz, &isuppz[1], 
	    &work[indwk], lwork, &iwork[1], liwork, info);
    if (wantz && *info == 0) {
      indwkn = inde;
      llwrkn = *lwork - indwkn + 1;
      PLUMED_BLAS_F77_FUNC(sormtr,SORMTR)("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
	      , &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
    }

    if (*info != 0) 
      return;

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&imax, &d__1, &w[1], &c__1);
    }

    if (wantz) {
	i__1 = *m - 1;
	
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		PLUMED_BLAS_F77_FUNC(sswap,SSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

    work[1] = (float) lwkopt;
    iwork[1] = liwmin;
    return;

}
}
}
#include <cctype>
#include <cmath>

#include "real.h"

#include "blas/blas.h"
#include "lapack.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssytd2,SSYTD2)(const char *    uplo,
	int *     n,
	float *  a,
	int *     lda,
	float *  d,
	float *  e,
	float *  tau,
    int *     info)
{
  float minusone,zero;
  float taui,alpha,tmp;
  int ti1,ti2,ti3,i;
  const char ch=std::toupper(*uplo);

  zero = 0.0;
  minusone = -1.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n-1;i>=1;i--) {

      ti1 = 1;
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&i,&(a[i*(*lda)+(i-1)]),&(a[i*(*lda)+0]),&ti1,&taui);
      e[i-1] = a[i*(*lda) + (i-1)];
      if(std::abs(taui)>PLUMED_GMX_FLOAT_MIN) {
	a[i*(*lda)+(i-1)] = 1.0;
      
	ti1 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)("U",&i,&taui,a,lda,&(a[i*(*lda)+0]),&ti1,&zero,tau,&ti1);

	tmp = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&i,tau,&ti1,&(a[i*(*lda)+0]),&ti1);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&i,&alpha,&(a[i*(*lda)+0]),&ti1,tau,&ti1);

	PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)("U",&i,&minusone,&(a[i*(*lda)+0]),&ti1,tau,&ti1,a,lda);

	a[i*(*lda)+(i-1)] = e[i-1]; 

      }
      d[i] = a[i*(*lda)+i];
      tau[i-1] = taui;
    }
    d[0] = a[0];
    
  } else {
    /* lower */

    for(i=1;i<*n;i++) {

      ti1 = *n - i;
      ti2 = ( *n < i+2) ? *n : i+2;
      ti3 = 1;
      PLUMED_BLAS_F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+ti2-1]),&ti3,&taui);

      e[i-1] = a[(i-1)*(*lda) + (i)];

      if(std::abs(taui)>PLUMED_GMX_FLOAT_MIN) {
	a[(i-1)*(*lda)+(i)] = 1.0;
      
	ti1 = *n - i;
	ti2 = 1;
	PLUMED_BLAS_F77_FUNC(ssymv,SSYMV)(uplo,&ti1,&taui,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),
	       &ti2,&zero,&(tau[i-1]),&ti2);
	
	tmp = PLUMED_BLAS_F77_FUNC(sdot,SDOT)(&ti1,&(tau[i-1]),&ti2,&(a[(i-1)*(*lda)+i]),&ti2);

	alpha = -0.5*taui*tmp;

	PLUMED_BLAS_F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2);

	PLUMED_BLAS_F77_FUNC(ssyr2,SSYR2)(uplo,&ti1,&minusone,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2,
	       &(a[(i)*(*lda)+i]),lda);

	a[(i-1)*(*lda)+(i)] = e[i-1]; 

      }
      d[i-1] = a[(i-1)*(*lda)+i-1];
      tau[i-1] = taui;
    }
    d[*n-1] = a[(*n-1)*(*lda)+(*n-1)];
 
  }
  return;
}
}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(ssytrd,SSYTRD)(const char *uplo, int *n, float *a, int *
	lda, float *d__, float *e, float *tau, float *
	work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, j, nb, kk, nx, iws;
    int nbmin, iinfo;
    int upper;
    int ldwork, lwkopt;
    int lquery;
    float c_b22 = -1.;
    float c_b23 = 1.;


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    upper = (*uplo=='U' || *uplo=='u');
    lquery = (*lwork == -1);

    if (! upper && ! (*uplo=='L' || *uplo=='l')) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ((1>*n) ? 1 : *n)) {
	*info = -4;
    } else if (*lwork < 1 && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {

      nb = DSYTRD_BLOCKSIZE;
      lwkopt = *n * nb;
      work[1] = (float) lwkopt;
    } else
      return;

    if (lquery) 
      return;
  
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    nx = *n;
    if (nb > 1 && nb < *n) {

	nx = DSYTRD_CROSSOVER;
	if (nx < *n) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		i__1 = *lwork / ldwork;
		nb = (i__1>1) ? i__1 : 1;
		nbmin = DSYTRD_MINBLOCKSIZE;
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

	    i__3 = i__ + nb - 1;
	    PLUMED_BLAS_F77_FUNC(slatrd,SLATRD)(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork);

	    i__3 = i__ - 1;
	    PLUMED_BLAS_F77_FUNC(ssyr2k,SSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 
		    + 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j - 1 + j * a_dim1] = e[j - 1];
		d__[j] = a[j + j * a_dim1];

	    }

	}

	PLUMED_BLAS_F77_FUNC(ssytd2,SSYTD2)(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
    } else {

	i__2 = *n - nx;
	i__1 = nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {


	    i__3 = *n - i__ + 1;
	    PLUMED_BLAS_F77_FUNC(slatrd,SLATRD)(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
		    tau[i__], &work[1], &ldwork);

	    i__3 = *n - i__ - nb + 1;
	    PLUMED_BLAS_F77_FUNC(ssyr2k,SSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + 
		    i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
		    i__ + nb + (i__ + nb) * a_dim1], lda);


	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j + 1 + j * a_dim1] = e[j];
		d__[j] = a[j + j * a_dim1];

	    }

	}


	i__1 = *n - i__ + 1;
	PLUMED_BLAS_F77_FUNC(ssytd2,SSYTD2)(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
		&tau[i__], &iinfo);
    }

    work[1] = (float) lwkopt;
    return;

}


}
}
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "blas/blas.h"
namespace PLMD{
namespace lapack{
using namespace blas;
void
PLUMED_BLAS_F77_FUNC(strti2,STRTI2)(const char *uplo,
	const char *diag, 
	int *n, 
	float *a,
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__2;

    int j;
    float ajj;
    int upper, nounit;
    int c__1 = 1;


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

    if (upper) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }

	    i__2 = j - 1;
	    PLUMED_BLAS_F77_FUNC(strmv,STRMV)("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1);
	    i__2 = j - 1;
	    PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
	}
    } else {

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(strmv,STRMV)("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1);
		i__1 = *n - j;
		PLUMED_BLAS_F77_FUNC(sscal,SSCAL)(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
	    }
	}
    }
    return;
}
}
}
#include <cmath>
#include "blas/blas.h"
#include "lapack.h"
#include "lapack_limits.h"

#include "real.h"

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
	    if (std::abs(a[*info + *info * a_dim1])<PLUMED_GMX_FLOAT_MIN) {
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
#endif

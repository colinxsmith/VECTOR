using System;
using System.Diagnostics;
namespace Blas
{
    public static class BlasLike
    {
        public static double lm_eps = Math.Abs((4.0 / 3 - 1) * 3 - 1);
        public static double lm_min = 2.2250738585072014e-308;
        public static double lm_rootmin = Math.Sqrt(lm_min);
        public static double lm_rooteps = Math.Sqrt(lm_eps);

        public static void daxpy(int n, double a, double[] x, int ix, double[] y, int iy)
        {
            if (a == 1)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] += x[iix];
            }
            else if (a == -1)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] -= x[iix];
            }
            else if (a != 0)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] += a * x[iix];
            }
        }


        public unsafe static void daxpy(int n, double a, double* x, int ix, double* y, int iy)
        {
            if (a == 1)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] += x[iix];
            }
            else if (a == -1)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] -= x[iix];
            }
            else if (a != 0)
            {
                for (int i = 0, iix = 0, iiy = 0; i < n; ++i, iiy += iy, iix += ix)
                    y[iiy] += a * x[iix];
            }
        }
        public static void daxpyvec(int n, double a, double[] x, double[] y)
        {
            daxpy(n, a, x, 1, y, 1);
        }
        public unsafe static void daxpyvec(int n, double a, double* x, double* y)
        {
            daxpy(n, a, x, 1, y, 1);
        }
        public static void dadd(int n, double[] x, int ix, double[] y, int iy, double[] z, int iz)
        {
            for (int i = 0, iix = 0, iiy = 0, iiz = 0; i < n; i++, iix += ix, iiy += iy, iiz += iz)
                z[iiz] = x[iix] + y[iiy];
        }
        public static void daddvec(int n, double[] x, double[] y, double[] z)
        {
            dadd(n, x, 1, y, 1, z, 1);
        }
        public static void dzero(int n, double[] x, int ix)
        {
            for (int i = 0, iix = 0; i < n; i++, iix += ix)
                x[iix] = 0;
        }
        public unsafe static void dzero(int n, double* x, int ix)
        {
            for (int i = 0, iix = 0; i < n; i++, iix += ix)
                x[iix] = 0;
        }
        public static void dzerovec(int n, double[] a)
        {
            for (int i = 0; i < n * sizeof(double); ++i)
            {
                Buffer.SetByte(a, i, 0);
            }
        }

        public unsafe static void dzerovec(int n, double* a)
        {
            while (n-- > 0)
            {
                *a++ = 0;
            }
        }
        public static void dcopy(int n, double[] x, int ix, double[] y, int iy)
        {
            if (ix == 1 && iy == 1)
            {
                dcopyvec(n, x, y);
            }
            else
            {
                for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0, iiy = iy < 0 ? -(n - 1) * iy : 0; i < n; i++, iix += ix, iiy += iy)
                    y[iiy] = x[iix];
            }
        }
        public static void dcopyvec(int n, double[] a, double[] b)
        {
            Buffer.BlockCopy(a, 0, b, 0, n * sizeof(double));
        }
        public static void dsub(int n, double[] x, int ix, double[] y, int iy, double[] z, int iz)
        {
            for (int i = 0, iix = 0, iiy = 0, iiz = 0; i < n; i++, iix += ix, iiy += iy, iiz += iz)
                z[iiz] = x[iix] - y[iiy];
        }
        public static void dsubvec(int n, double[] x, double[] y, double[] z)
        {
            dsub(n, x, 1, y, 1, z, 1);
        }
        public static double ddot(int n, double[] a, int ia, double[] b, int ib)
        {
            double back = 0;
            for (int i = 0, iia = 0, iib = 0; i < n; i++, iia += ia, iib += ib)
            {
                back += a[iia] * b[iib];
            }
            return back;
        }

        public unsafe static double ddot(int n, double* a, int ia, double* b, int ib)
        {
            double back = 0;
            for (int i = 0, iia = 0, iib = 0; i < n; i++, iia += ia, iib += ib)
            {
                back += a[iia] * b[iib];
            }
            return back;
        }
        public static double ddotvec(int n, double[] a, double[] b)
        {
            return ddot(n, a, 1, b, 1);
        }
        public unsafe static double ddotvec(int n, double* a, double* b)
        {
            return ddot(n, a, 1, b, 1);
        }
        public static void dset(int n, double a, double[] x, int ix)
        {
            for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
            {
                x[iix] = a;
            }
        }
        public static void dsetvec(int n, double a, double[] x)
        {
            dset(n, a, x, 1);
        }
        public static void dneg(int n, double[] x, int ix)
        {
            for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
                x[iix] = -x[iix];
        }
        public unsafe static void dneg(int n, double* x, int ix)
        {
            for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
                x[iix] = -x[iix];
        }
        public static void dnegvec(int n, double[] x)
        {
            for (int i = 0; i < n; ++i)
                x[i] = -x[i];
        }

        public static void dscal(int n, double a, double[] x, int ix)
        {
            if (a == 0)
            {
                dzero(n, x, ix);
            }
            else if (a == -1)
            {
                dneg(n, x, ix);
            }
            else
            {
                for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
                    x[iix] = a * x[iix];
            }
        }
        public unsafe static void dscal(int n, double da, double* dx,
    int incx)
        {
            /* System generated locals */
            int i__1, i__2;

            /* Local variables */
            int i__, m, mp1, nincx;


            /*  -- Reference BLAS level1 routine (version 3.8.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     November 2017 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */
            /* Parameter adjustments */
            --dx;
            /* Function Body */
            if (incx == 1)
            {
                /*        code for increment equal to 1 */
                /*        clean-up loop */
                m = n % 5;
                if (m != 0)
                {
                    i__1 = m;
                    if (da != 0)
                    {
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            dx[i__] = da * dx[i__];
                        }
                    }
                    else
                    {
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            dx[i__] = 0;
                        }
                    }
                    if (n < 5)
                    {
                        return;
                    }
                }
                mp1 = m + 1;
                i__1 = n;
                if (da != 0)
                {
                    if (da == 1)
                    {
                        for (i__ = mp1; i__ <= i__1; i__ += 5)
                        {
                            dx[i__] = dx[i__];
                            dx[i__ + 1] = dx[i__ + 1];
                            dx[i__ + 2] = dx[i__ + 2];
                            dx[i__ + 3] = dx[i__ + 3];
                            dx[i__ + 4] = dx[i__ + 4];
                        }
                    }
                    else
                    {
                        for (i__ = mp1; i__ <= i__1; i__ += 5)
                        {
                            dx[i__] = da * dx[i__];
                            dx[i__ + 1] = da * dx[i__ + 1];
                            dx[i__ + 2] = da * dx[i__ + 2];
                            dx[i__ + 3] = da * dx[i__ + 3];
                            dx[i__ + 4] = da * dx[i__ + 4];
                        }
                    }
                }
                else
                {
                    for (i__ = mp1; i__ <= i__1; i__ += 5)
                    {
                        dx[i__] = 0;
                        dx[i__ + 1] = 0;
                        dx[i__ + 2] = 0;
                        dx[i__ + 3] = 0;
                        dx[i__ + 4] = 0;
                    }
                }
            }
            else
            {
                /*        code for increment not equal to 1 */
                nincx = n * incx;
                i__1 = nincx;
                i__2 = incx;
                if (da == 0)
                {
                    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
                    {
                        dx[i__] = 0;
                    }
                }
                else
                {
                    if (da == 1)
                    {
                        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
                        {
                            dx[i__] = dx[i__];
                        }
                    }
                    else
                    {
                        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
                        {
                            dx[i__] = da * dx[i__];
                        }
                    }
                }
            }
        }
        public static void dscalvec(int n, double a, double[] x)
        {
            dscal(n, a, x, 1);
        }
        public unsafe static void dscalvec(int n, double a, double* x)
        {
            dscal(n, a, x, 1);
        }
        public static void dsccopy(int n, double a, double[] x, int ix, double[] y, int iy)
        {
            if (n > 0)
            {
                if (a == 0) dzero(n, y, iy);
                else if (a == 1) dcopy(n, x, ix, y, iy);
                else
                {
                    if (a == -1)
                    {
                        for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0, iiy = iy < 0 ? -(n - 1) * iy : 0; i < n; i++, iix += ix, iiy += iy)
                        {
                            y[iiy] = -x[iix];
                        }
                    }
                    else
                    {
                        for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0, iiy = iy < 0 ? -(n - 1) * iy : 0; i < n; i++, iix += ix, iiy += iy)
                        {
                            y[iiy] = a * x[iix];
                        }
                    }
                }
            }
        }
        public static void dsccopyvec(int n, double a, double[] x, double[] y)
        {
            dsccopy(n, a, x, 1, y, 1);
        }
        public unsafe static void dsssq(int n, double* x, int ix, double[] pscale, double[] psumsq, int px = 0)
        {
            /*
                  dsssqvec returns the values scl and smsq such that
                  ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
                  to be at least unity and the value of smsq will then satisfy
                  1.0 .le. smsq .le. ( sumsq + n ) .
                  scale is assumed to be non-negative and scl returns the value
                  scl = max( scale, abs( x( i ) ) ) .
                  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
                  scl and smsq are overwritten on SCALE and SUMSQ respectively.
                  The routine makes only one pass through the vector X.
          */
            if (n > 0)
            {
                double absxi, d, sumsq = psumsq[0], scale = pscale[0];
                Debug.Assert(scale >= 0);
                for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; ++i, iix += ix)
                {
                    absxi = x[iix + px];
                    if (absxi == 0) continue;
                    if (absxi < 0) absxi = -absxi;
                    if (scale < absxi)
                    {
                        d = scale / absxi;
                        sumsq = sumsq * (d * d) + 1;
                        scale = absxi;
                    }
                    else
                    {
                        d = absxi / scale;
                        sumsq += d * d;
                    }
                }
                pscale[0] = scale;
                psumsq[0] = sumsq;
            }
        }

        public static void dsssq(int n, double[] x, int ix, double[] pscale, double[] psumsq, int px = 0)
        {
            /*
                  dsssqvec returns the values scl and smsq such that
                  ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
                  to be at least unity and the value of smsq will then satisfy
                  1.0 .le. smsq .le. ( sumsq + n ) .
                  scale is assumed to be non-negative and scl returns the value
                  scl = max( scale, abs( x( i ) ) ) .
                  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
                  scl and smsq are overwritten on SCALE and SUMSQ respectively.
                  The routine makes only one pass through the vector X.
          */
            if (n > 0)
            {
                double absxi, d, sumsq = psumsq[0], scale = pscale[0];
                Debug.Assert(scale >= 0);
                for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; ++i, iix += ix)
                {
                    absxi = x[iix + px];
                    if (absxi == 0) continue;
                    if (absxi < 0) absxi = -absxi;
                    if (scale < absxi)
                    {
                        d = scale / absxi;
                        sumsq = sumsq * (d * d) + 1;
                        scale = absxi;
                    }
                    else
                    {
                        d = absxi / scale;
                        sumsq += d * d;
                    }
                }
                pscale[0] = scale;
                psumsq[0] = sumsq;
            }
        }

        public static void dsssqvec(int n, double[] x, double[] pscale, double[] psumsq, int px = 0)
        {
            if (n > 0)
            {
                double absxi, d, sumsq = psumsq[0], scale = pscale[0];
                Debug.Assert(scale >= 0);
                for (int i = 0; i < n; ++i)
                {
                    absxi = x[i + px];
                    if (absxi == 0) continue;
                    if (absxi < 0) absxi = -absxi;
                    if (scale < absxi)
                    {
                        d = scale / absxi;
                        sumsq = sumsq * (d * d) + 1;
                        scale = absxi;
                    }
                    else
                    {
                        d = absxi / scale;
                        sumsq += d * d;
                    }
                }
                pscale[0] = scale;
                psumsq[0] = sumsq;
            }
        }

        public static double dsum(int n, double[] x, int ix = 1, int px = 0)
        {
            double back = 0;
            for (int i = 0, iix = 0; i < n; i++, iix += ix)
            {
                back += x[iix + px];
            }
            return back;
        }
        public static double dsumvec(int n, double[] x, int px = 0)
        {
            return dsum(n, x, 1, px);
        }
        public unsafe static void dxminmax(int n, double* x, int ix, double* xmax, double* xmin)
        {
            if (n < 1) *xmin = *xmax = 0;
            else
            {
                double ax = Math.Abs(*x);
                double xm = ax, xn = ax;
                while (--n > 0)
                {
                    x += ix;
                    ax = *x;
                    if (ax < 0) ax = -ax;
                    if (ax > xm) xm = ax;
                    if (ax < xn) xn = ax;
                }
                *xmin = xn;
                *xmax = xm;
            }
        }
        public static void dxminmax(int n, double[] x, int ix, double[] xmax, double[] xmin, int px = 0)
        {
            if (n < 1) xmin[0] = xmax[0] = 0;
            else
            {
                double ax = Math.Abs(x[px]);
                double xm = ax, xn = ax;
                for (int i = 0, iix = 0; i < n; i++, iix += ix)
                {
                    ax = x[iix + px];
                    if (ax < 0) ax = -ax;
                    if (ax > xm) xm = ax;
                    if (ax < xn) xn = ax;
                }
                xmin[0] = xn;
                xmax[0] = xm;
            }
        }
        public unsafe static void detagen(int n, double* alpha, double* x, int ix, long* iswap, int* itrans)
        {
            long imax = 1000000000;
            int nzero;
            double xmax, absalf, tol, axi;
            double* v, vlim;
            *iswap = 0;
            *itrans = 0;
            if (n < 1) return;
            absalf = Math.Abs(*alpha);
            xmax = 0;
            for (v = x, vlim = x + n * ix; v != vlim; v += ix)
            {
                if (xmax < (axi = Math.Abs(*v)))
                {
                    xmax = axi;
                    imax = v - x;
                }
            }
            /* exit if  x  is very small */
            if (xmax <= lm_rootmin) return;
            /* see if an interchange is needed for stability */
            if (absalf < xmax)
            {
                *iswap = imax + 1;
                xmax = x[imax];
                x[imax] = *alpha;
                *alpha = xmax;
            }
            /*
                 form the multipliers in  x.  they will be no greater than one
                 in magnitude.  change negligible multipliers to zero
            */
            tol = Math.Abs(*alpha) * lm_eps;
            nzero = 0;
            for (v = x; v != vlim; v += ix)
            {
                if (Math.Abs(*v) > tol) *v = -*v / *alpha;
                else
                {
                    *v = 0;
                    ++nzero;
                }
            }
            /*z is zero only if nzero=n*/
            if (nzero < n) *itrans = 1;
        }
        public unsafe static void delm(int orthog, int n, double* x, int ix, double* y, int iy, double cs, double sn)
        {
            /*
                   If  orthog  is true, delm  applies a plane rotation.  otherwise,
                   elm computes the transformation (x y)*e  and returns the result
                   in  (x y),  where the 2 by 2 matrix  e  is defined by  cs  and  sn

                   as follows...
                   e  =    ( 1  sn )       if  cs>0 else   e  =    (     1 )
                           (     1 )                               ( 1  sn )
           */
            if (orthog == 0)
            {
                if (cs <= 0) dswap(n, x, ix, y, iy);
                if (sn != 0) daxpy(n, sn, x, ix, y, iy);
            }
            else dsymplanerotate(n, x, ix, y, iy, cs, sn);
        }
        public unsafe static void dsymplanerotate(int n, double* x, int ix, double* y, int iy, double c, double s)
        {
            int i__1, i__2;

            double temp1;
            int i, iix, iiy;
            --y;
            --x;

            if (n > 0 && s != 0)
            {
                if (c == 0 && s == 1)
                {
                    if (ix == iy && ix > 0)
                    {
                        i__1 = (n - 1) * ix + 1;
                        i__2 = ix;
                        for (iix = 1; (i__2 < 0 ? iix >= i__1 : iix <= i__1); iix += i__2)
                        {
                            temp1 = x[iix];
                            x[iix] = y[iix];
                            y[iix] = temp1;
                        }
                    }
                    else
                    {
                        if (iy >= 0)
                        {
                            iiy = 1;
                        }
                        else
                        {
                            iiy = 1 - (n - 1) * iy;
                        }
                        if (ix > 0)
                        {
                            i__2 = (n - 1) * ix + 1;
                            i__1 = ix;
                            for (iix = 1; (i__1 < 0 ? iix >= i__2 : iix <= i__2); iix += i__1)
                            {
                                temp1 = x[iix];
                                x[iix] = y[iiy];
                                y[iiy] = temp1;
                                iiy += iy;
                            }
                        }
                        else
                        {
                            iix = 1 - (n - 1) * ix;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[iix];
                                x[iix] = y[iiy];
                                y[iiy] = temp1;
                                iix += ix;
                                iiy += iy;
                            }
                        }
                    }
                }
                else if (c == 0 && s == -1)
                {
                    if (ix == iy && ix > 0)
                    {
                        i__1 = (n - 1) * ix + 1;
                        i__2 = ix;
                        for (iix = 1; (i__2 < 0 ? iix >= i__1 : iix <= i__1); iix += i__2)
                        {
                            temp1 = -x[iix];
                            x[iix] = -y[iix];
                            y[iix] = temp1;
                        }
                    }
                    else
                    {
                        if (iy >= 0)
                        {
                            iiy = 1;
                        }
                        else
                        {
                            iiy = 1 - (n - 1) * iy;
                        }
                        if (ix > 0)
                        {
                            i__2 = (n - 1) * ix + 1;
                            i__1 = ix;
                            for (iix = 1; (i__1 < 0 ? iix >= i__2 : iix <= i__2); iix += i__1)
                            {
                                temp1 = -x[iix];
                                x[iix] = -y[iiy];
                                y[iiy] = temp1;
                                iiy += iy;
                            }
                        }
                        else
                        {
                            iix = 1 - (n - 1) * ix;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = -x[iix];
                                x[iix] = -y[iiy];
                                y[iiy] = temp1;
                                iix += ix;
                                iiy += iy;
                            }
                        }
                    }
                }
                else
                {
                    if (ix == iy && ix > 0)
                    {
                        i__1 = (n - 1) * ix + 1;
                        i__2 = ix;
                        for (iix = 1; (i__2 < 0 ? iix >= i__1 : iix <= i__1); iix += i__2)
                        {
                            temp1 = x[iix];
                            x[iix] = c * temp1 + s * y[iix];
                            y[iix] = s * temp1 - c * y[iix];
                        }
                    }
                    else
                    {
                        if (iy >= 0)
                        {
                            iiy = 1;
                        }
                        else
                        {
                            iiy = 1 - (n - 1) * iy;
                        }
                        if (ix > 0)
                        {
                            i__2 = (n - 1) * ix + 1;
                            i__1 = ix;
                            for (iix = 1; (i__1 < 0 ? iix >= i__2 : iix <= i__2); iix += i__1)
                            {
                                temp1 = x[iix];
                                x[iix] = c * temp1 + s * y[iiy];
                                y[iiy] = s * temp1 - c * y[iiy];
                                iiy += iy;
                            }
                        }
                        else
                        {
                            iix = 1 - (n - 1) * ix;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[iix];
                                x[iix] = c * temp1 + s * y[iiy];
                                y[iiy] = s * temp1 - c * y[iiy];
                                iix += ix;
                                iiy += iy;
                            }
                        }
                    }
                }
            }
        }

        public unsafe static void dswap(int n, double* a, int ia, double* b, int ib)
        {
            for (; (n--) > 0; a += ia, b += ib)
            {
                double temp = *a;
                *a = *b;
                *b = temp;
            }
        }
        public unsafe static
        void dswapvec(int n, double* a, double* b)
        {
            while ((n--) > 0)
            {
                double temp = *a;
                *a++ = *b;
                *b++ = temp;
            }
        }
    }
}

using System;
using System.Diagnostics;
namespace Blas
{
    public static class BlasLike
    {
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
        public static void daxpyvec(int n, double a, double[] x, double[] y)
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
        public static void dzerovec(int n, double[] a)
        {
            for (int i = 0; i < n * sizeof(double); ++i)
            {
                Buffer.SetByte(a, i, 0);
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
        public static double ddotvec(int n, double[] a, double[] b)
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
        public static void dscalvec(int n, double a, double[] x)
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
    }
}

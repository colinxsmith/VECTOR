using System;
using System.Diagnostics;
namespace Blas
{
    public static class BlasLike
    {
        public static int baseref = 0;//Set to an array address for debug output
        public static double lm_eps = Math.Abs((((double)4) / 3 - 1) * 3 - 1);
        public static double lm_min = 2.2250738585072014e-308;
        public static double lm_rootmin = Math.Sqrt(lm_min);
        public static double lm_rooteps = Math.Sqrt(lm_eps);

        public static void daxpy1(int n, double a, double[] x, int ix, double[] y, int iy)
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

        public unsafe static void daxpy(int n, double da, double* dx,
            int incx, double* dy, int incy)
        {
            int i__, m, ix, iy, mp1;
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
            --dy;
            --dx;

            /* Function Body */
            if (n <= 0)
            {
                return;
            }
            if (da == 0)
            {
                return;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 4;
                if (m != 0)
                {
                    if (da == 1)
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dy[i__] += dx[i__];
                        }
                    }
                    else if (da == -1)
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dy[i__] -= dx[i__];
                        }
                    }
                    else
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dy[i__] += da * dx[i__];
                        }
                    }
                }
                if (n < 4)
                {
                    return;
                }
                mp1 = m + 1;
                if (da == 1)
                {
                    for (i__ = mp1; i__ <= n; i__ += 4)
                    {
                        dy[i__] += dx[i__];
                        dy[i__ + 1] += dx[i__ + 1];
                        dy[i__ + 2] += dx[i__ + 2];
                        dy[i__ + 3] += dx[i__ + 3];
                    }
                }
                else if (da == -1)
                {
                    for (i__ = mp1; i__ <= n; i__ += 4)
                    {
                        dy[i__] -= dx[i__];
                        dy[i__ + 1] -= dx[i__ + 1];
                        dy[i__ + 2] -= dx[i__ + 2];
                        dy[i__ + 3] -= dx[i__ + 3];
                    }
                }
                else
                {
                    for (i__ = mp1; i__ <= n; i__ += 4)
                    {
                        dy[i__] += da * dx[i__];
                        dy[i__ + 1] += da * dx[i__ + 1];
                        dy[i__ + 2] += da * dx[i__ + 2];
                        dy[i__ + 3] += da * dx[i__ + 3];
                    }
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                if (da == 1)
                {
                    for (i__ = 1; i__ <= n; ++i__)
                    {
                        dy[iy] += dx[ix];
                        ix += incx;
                        iy += incy;
                    }
                }
                else if (da == -1)
                {
                    for (i__ = 1; i__ <= n; ++i__)
                    {
                        dy[iy] -= dx[ix];
                        ix += incx;
                        iy += incy;
                    }
                }
                else
                {
                    for (i__ = 1; i__ <= n; ++i__)
                    {
                        dy[iy] += da * dx[ix];
                        ix += incx;
                        iy += incy;
                    }
                }
            }
        }


        public static void daxpy(int n, double da, double[] dx,
            int incx, double[] dy, int incy, int xstart = 0, int ystart = 0)
        {
            int i__, m, ix, iy;
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

            /* Function Body */
            if (n <= 0)
            {
                return;
            }
            if (da == 0)
            {
                return;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 4;
                if (m != 0)
                {
                    if (da == 1)
                    {
                        for (i__ = 0; i__ < m; ++i__)
                        {
                            dy[i__ + ystart] += dx[i__ + xstart];
                        }
                    }
                    else if (da == -1)
                    {
                        for (i__ = 0; i__ < m; ++i__)
                        {
                            dy[i__ + ystart] -= dx[i__ + xstart];
                        }
                    }
                    else
                    {
                        for (i__ = 0; i__ < m; ++i__)
                        {
                            dy[i__ + ystart] += da * dx[i__ + xstart];
                        }
                    }
                }
                if (n < 4)
                {
                    return;
                }
                if (da == 1)
                {
                    for (i__ = m; i__ < n; i__ += 4)
                    {
                        dy[i__ + ystart] += dx[i__ + xstart];
                        dy[i__ + 1 + ystart] += dx[i__ + 1 + xstart];
                        dy[i__ + 2 + ystart] += dx[i__ + 2 + xstart];
                        dy[i__ + 3 + ystart] += dx[i__ + 3 + xstart];
                    }
                }
                else if (da == -1)
                {
                    for (i__ = m; i__ < n; i__ += 4)
                    {
                        dy[i__ + ystart] -= dx[i__ + xstart];
                        dy[i__ + 1 + ystart] -= dx[i__ + 1 + xstart];
                        dy[i__ + 2 + ystart] -= dx[i__ + 2 + xstart];
                        dy[i__ + 3 + ystart] -= dx[i__ + 3 + xstart];
                    }
                }
                else
                {
                    for (i__ = m; i__ < n; i__ += 4)
                    {
                        dy[i__ + ystart] += da * dx[i__ + xstart];
                        dy[i__ + 1 + ystart] += da * dx[i__ + 1 + xstart];
                        dy[i__ + 2 + ystart] += da * dx[i__ + 2 + xstart];
                        dy[i__ + 3 + ystart] += da * dx[i__ + 3 + xstart];
                    }
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                if (da == 1)
                {
                    for (i__ = 0; i__ < n; ++i__)
                    {
                        dy[iy + ystart] += dx[ix + xstart];
                        ix += incx;
                        iy += incy;
                    }
                }
                else if (da == -1)
                {
                    for (i__ = 0; i__ < n; ++i__)
                    {
                        dy[iy + ystart] -= dx[ix + xstart];
                        ix += incx;
                        iy += incy;
                    }
                }
                else
                {
                    for (i__ = 0; i__ < n; ++i__)
                    {
                        dy[iy + ystart] += da * dx[ix + xstart];
                        ix += incx;
                        iy += incy;
                    }
                }
            }
        }


        public unsafe static void daxpy1(int n, double a, double* x, int ix, double* y, int iy)
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
        public static void daxpyvec(int n, double a, double[] x, double[] y, int xstart = 0, int ystart = 0)
        {
            daxpy(n, a, x, 1, y, 1, xstart, ystart);
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
        public static void dzero(int n, double[] x, int ix, int xstart = 0)
        {
            for (int i = 0, iix = 0; i < n; i++, iix += ix)
                x[iix + xstart] = 0;
        }
        public unsafe static void dzero(int n, double* x, int ix)
        {
            for (int i = 0, iix = 0; i < n; i++, iix += ix)
                x[iix] = 0;
        }
        public static void dzerovec(int n, double[] a, int astart = 0)
        {
            for (int i = 0; i < n; ++i)
            {
                a[i + astart] = 0;
            }
        }

        public unsafe static void dzerovec(int n, double* a)
        {
            while (n-- > 0)
            {
                *a++ = 0;
            }
        }
        public static void dcopy(int n, double[] dx, int incx, double[] dy, int incy)
        {
            int i__, m, ix, iy;
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


            /* Function Body */
            if (n <= 0)
            {
                return;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 7;
                if (m != 0)
                {
                    for (i__ = 0; i__ < m; ++i__)
                    {
                        dy[i__] = dx[i__];
                    }
                    if (n < 7)
                    {
                        return;
                    }
                }
                for (i__ = m; i__ < n; i__ += 7)
                {
                    dy[i__] = dx[i__];
                    dy[i__ + 1] = dx[i__ + 1];
                    dy[i__ + 2] = dx[i__ + 2];
                    dy[i__ + 3] = dx[i__ + 3];
                    dy[i__ + 4] = dx[i__ + 4];
                    dy[i__ + 5] = dx[i__ + 5];
                    dy[i__ + 6] = dx[i__ + 6];
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                for (i__ = 0; i__ < n; ++i__)
                {
                    dy[iy - 1] = dx[ix - 1];
                    ix += incx;
                    iy += incy;
                }
            }
        }

        public unsafe static void dcopy(int n, double* dx, int incx, double* dy, int incy)
        {
            int i__, m, ix, iy, mp1;
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
            --dy;
            --dx;

            /* Function Body */
            if (n <= 0)
            {
                return;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 7;
                if (m != 0)
                {
                    for (i__ = 1; i__ <= m; ++i__)
                    {
                        dy[i__] = dx[i__];
                    }
                    if (n < 7)
                    {
                        return;
                    }
                }
                mp1 = m + 1;
                for (i__ = mp1; i__ <= n; i__ += 7)
                {
                    dy[i__] = dx[i__];
                    dy[i__ + 1] = dx[i__ + 1];
                    dy[i__ + 2] = dx[i__ + 2];
                    dy[i__ + 3] = dx[i__ + 3];
                    dy[i__ + 4] = dx[i__ + 4];
                    dy[i__ + 5] = dx[i__ + 5];
                    dy[i__ + 6] = dx[i__ + 6];
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                for (i__ = 1; i__ <= n; ++i__)
                {
                    dy[iy] = dx[ix];
                    ix += incx;
                    iy += incy;
                }
            }
        }
        public static void dcopy1(int n, double[] x, int ix, double[] y, int iy)
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
        public static double ddot1(int n, double[] a, int ia, double[] b, int ib)

        {
            double back = 0;
            for (int i = 0, iia = 0, iib = 0; i < n; i++, iia += ia, iib += ib)
            {
                back += a[iia] * b[iib];
            }
            return back;
        }
        public unsafe static double ddot(int n, double* dx, int incx, double* dy,
            int incy)
        {   /* System generated locals */
            int i__1;
            double ret_val;

            /* Local variables */
            int i__, m, ix, iy, mp1;
            double dtemp;


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
            --dy;
            --dx;

            /* Function Body */
            ret_val = 0;
            dtemp = 0;
            if (n <= 0)
            {
                return ret_val;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 5;
                if (m != 0)
                {
                    i__1 = m;
                    for (i__ = 1; i__ <= i__1; ++i__)
                    {
                        dtemp += dx[i__] * dy[i__];
                    }
                    if (n < 5)
                    {
                        ret_val = dtemp;
                        return ret_val;
                    }
                }
                mp1 = m + 1;
                i__1 = n;
                for (i__ = mp1; i__ <= i__1; i__ += 5)
                {
                    dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] +
                        dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] +
                        dx[i__ + 4] * dy[i__ + 4];
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    dtemp += dx[ix] * dy[iy];
                    ix += incx;
                    iy += incy;
                }
            }
            ret_val = dtemp;
            return ret_val;

        } /* ddot_ */

        public static double ddot(int n, double[] dx, int incx, double[] dy,
            int incy, int dxstart = 0, int dystart = 0)
        {   /* System generated locals */
            int i__1;
            double ret_val;

            /* Local variables */
            int i__, m, ix, iy, mp1;
            double dtemp;


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

            /* Function Body */
            ret_val = 0;
            dtemp = 0;
            if (n <= 0)
            {
                return ret_val;
            }
            if (incx == 1 && incy == 1)
            {

                /*        code for both increments equal to 1 */


                /*        clean-up loop */

                m = n % 5;
                if (m != 0)
                {
                    i__1 = m;
                    for (i__ = 0; i__ < i__1; ++i__)
                    {
                        dtemp += dx[i__ + dxstart] * dy[i__ + dystart];
                    }
                    if (n < 5)
                    {
                        ret_val = dtemp;
                        return ret_val;
                    }
                }
                mp1 = m + 1;
                i__1 = n;
                for (i__ = mp1 - 1; i__ < i__1; i__ += 5)
                {
                    dtemp = dtemp + dx[i__ + dxstart] * dy[i__ + dystart] + dx[i__ + 1 + dxstart] * dy[i__ + 1 + dystart] +
                        dx[i__ + 2 + dxstart] * dy[i__ + 2 + dystart] + dx[i__ + 3 + dxstart] * dy[i__ + 3 + dystart] +
                        dx[i__ + 4 + dxstart] * dy[i__ + 4 + dystart];
                }
            }
            else
            {

                /*        code for unequal increments or equal increments */
                /*          not equal to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    dtemp += dx[ix - 1 + dxstart] * dy[iy - 1 + dystart];
                    ix += incx;
                    iy += incy;
                }
            }
            ret_val = dtemp;
            return ret_val;

        } /* ddot_ */


        public unsafe static double ddot1(int n, double* a, int ia, double* b, int ib)
        {
            double back = 0;
            for (int i = 0, iia = 0, iib = 0; i < n; i++, iia += ia, iib += ib)
            {
                back += a[iia] * b[iib];
            }
            return back;
        }
        public unsafe static double ddot2(int n, double* x, int incx, double* y, int incy)
        {
            double sum;
            sum = 0;
            while (n-- > 0)
            {
                sum += *x * *y;
                x += incx;
                y += incy;
            }
            return sum;
        }
        public static double ddotvec(int n, double[] a, double[] b, int astart = 0, int bstart = 0)
        {
            return ddot(n, a, 1, b, 1, astart, bstart);
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
        public static void dneg(int n, double[] x, int ix, int xstart = 0)
        {
            for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
                x[iix + xstart] = -x[iix + xstart];
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

        public unsafe static void dscal(int n, double a, double[] x, int ix, int xstart = 0)
        {
            if (baseref != 0)
            {
                fixed (double* xx = x)
                    dscal(n, a, xx + xstart, ix);
                return;
            }
            if (a == 0)
            {
                dzero(n, x, ix, xstart);
            }
            else if (a == -1)
            {
                dneg(n, x, ix, xstart);
            }
            else
            {
                for (int i = 0, iix = ix < 0 ? -(n - 1) * ix : 0; i < n; i++, iix += ix)
                    x[iix + xstart] = a * x[iix + xstart];
            }
        }
        public static void dger1(int m, int n, double alpha, double[] x, int incx, double[] y, int incy, double[] a, int lda, int xstart = 0, int ystart = 0, int astart = 0)
        {
            int a_dim1;
            //    size_t i, 
            int jy, kx, info;
            int j;
            //    double temp;
            //	vector px,pA,py;

            /*     .. Scalar Arguments .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  Purpose */
            /*  ======= */

            /*  DGER   performs the rank 1 operation */

            /*     A := alpha*x*y' + A, */

            /*  where alpha is a scalar, x is an m element vector, y is an n element */
            /*  vector and A is an m by n matrix. */

            /*  Parameters */
            /*  ========== */

            /*  M      - size_t. */
            /*           On entry, M specifies the number of rows of the matrix A. */
            /*           M must be at least zero. */
            /*           Unchanged on exit. */

            /*  N      - size_t. */
            /*           On entry, N specifies the number of columns of the matrix A. */
            /*           N must be at least zero. */
            /*           Unchanged on exit. */

            /*  ALPHA  - DOUBLE PRECISION. */
            /*           On entry, ALPHA specifies the scalar alpha. */
            /*           Unchanged on exit. */

            /*  X      - DOUBLE PRECISION array of dimension at least */
            /*           ( 1 + ( m - 1 )*abs( INCX ) ). */
            /*           Before entry, the incremented array X must contain the m */
            /*           element vector x. */
            /*           Unchanged on exit. */

            /*  INCX   - size_t. */
            /*           On entry, INCX specifies the increment for the elements of */
            /*           X. INCX must not be zero. */
            /*           Unchanged on exit. */

            /*  Y      - DOUBLE PRECISION array of dimension at least */
            /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
            /*           Before entry, the incremented array Y must contain the n */
            /*           element vector y. */
            /*           Unchanged on exit. */

            /*  INCY   - size_t. */
            /*           On entry, INCY specifies the increment for the elements of */
            /*           Y. INCY must not be zero. */
            /*           Unchanged on exit. */

            /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
            /*           Before entry, the leading m by n part of the array A must */
            /*           contain the matrix of coefficients. On exit, A is */
            /*           overwritten by the updated matrix. */

            /*  LDA    - size_t. */
            /*           On entry, LDA specifies the first dimension of A as declared */
            /*           in the calling (sub) program. LDA must be at least */
            /*           Math.Max( 1, m ). */
            /*           Unchanged on exit. */


            /*  Level 2 Blas routine. */

            /*  -- Written on 22-October-1986. */
            /*     Jack Dongarra, Argonne National Lab. */
            /*     Jeremy Du Croz, Nag Central Office. */
            /*     Sven Hammarling, Nag Central Office. */
            /*     Richard Hanson, Sandia National Labs. */


            /*     .. Parameters .. */
            /*     .. Local Scalars .. */
            /*     .. External Subroutines .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */
            /*     .. Executable Statements .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            a_dim1 = lda;

            /* Function Body */
            info = 0;
            if (m < 0)
            {
                info = 1;
            }
            else if (n < 0)
            {
                info = 2;
            }
            else if (incx == 0)
            {
                info = 5;
            }
            else if (incy == 0)
            {
                info = 7;
            }
            else if (lda < Math.Max(1, m))
            {
                info = 9;
            }
            if (info != 0)
            {
                Console.WriteLine($"dger1 error code {info}");
                return;
            }

            /*     Quick return if possible. */

            if (m == 0 || n == 0 || alpha == 0.0)
            {
                return;
            }

            /*     Start the operations. In this version the elements of A are */
            /*     accessed sequentially with one pass through A. */
            if (incy > 0)
            {
                jy = 0;
            }
            else
            {
                jy = 0 - (n - 1) * incy;
            }
            if (incx == 1)
            {
                if (m > 500)
                {
                    //#pragma omp parallel for private(j) schedule(dynamic)
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j + ystart] != 0.0)
                        {
                            daxpyvec(m, alpha * y[jy + incy * j + ystart], x, a, xstart, astart + j * a_dim1);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j + ystart] != 0.0)
                        {
                            daxpyvec(m, alpha * y[jy + incy * j + ystart], x, a, xstart, astart + j * a_dim1);
                        }
                    }
                }
            }
            else
            {
                if (incx > 0)
                {
                    kx = 0;
                }
                else
                {
                    kx = 0 - (m - 1) * incx;
                }
                if (m > 500)
                {
                    //#pragma omp parallel for private(j) schedule(dynamic)
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j + ystart] != 0.0)
                        {
                            daxpy(m, alpha * y[jy + incy * j + ystart], x, incx, a, 1, xstart + kx, astart + j * a_dim1);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j + ystart] != 0.0)
                        {
                            daxpy(m, alpha * y[jy + incy * j + ystart], x, incx, a, 1, xstart + kx, astart + j * a_dim1);
                        }
                    }
                }
            }
        }

        public unsafe static void dger(int m, int n, double alpha,
            double* x, int incx, double* y, int incy,
            double* a, int lda)
        {
            if (incx == 1 && incy == 1 && m > 1)
            {
                double oned = 1;
                int mm = (int)m, nn = (int)n, one = 1, ldaa = (int)lda;
                Console.WriteLine("Forward to dgemm from dger");
                char transa = 'N';
                char transb = 'T';
                dgemm(&transa, &transb, &mm, &nn, &one, &alpha, x, &mm, y, &nn, &oned, a, &ldaa);
            }
            else dger1(m, n, alpha, x, incx, y, incy, a, lda);
        }


        public unsafe static void dger(int m, int n, double alpha,
            double[] x, int incx, double[] y, int incy,
            double[] a, int lda, int xstart = 0, int ystart = 0, int astart = 0)
        {
            if (incx == 1 && incy == 1 && m > 1)
            {
                double oned = 1;
                int mm = (int)m, nn = (int)n, one = 1, ldaa = (int)lda;
                Console.WriteLine("Forward to dgemm from dger");
                char[] transa = { 'N' };
                char[] transb = { 'T' };
                dgemm(transa, transb, mm, nn, one, alpha, x, mm, y, nn, oned, a, ldaa, xstart, ystart, astart);
            }
            else
            {
                if (baseref != 0)
                {
                    fixed (double* px = x)
                    fixed (double* py = y)
                    fixed (double* pa = a)
                        dger1(m, n, alpha, px + xstart, incx, py + ystart, incy, pa + astart, lda);
                    return;
                }
                dger1(m, n, alpha, x, incx, y, incy, a, lda, xstart, ystart, astart);
            }
        }

        public unsafe static int dgemm(char* transa, char* transb, int* m, int*
n, int* k, double* alpha, double* a, int* lda,
double* b, int* ldb, double* beta, double* c__,
int* ldc)
        {
            /* System generated locals */
            int a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset;// i__1, i__2, 
                                                                     //    i__3;

            /* Local variables */
            int j, info, i__, l; double temp;
            short nota, notb;
            int ncola;
            int nrowa, nrowb;

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  Purpose */
            /*  ======= */

            /*  DGEMM  performs one of the matrix-matrix operations */

            /*     C := alpha*op( A )*op( B ) + beta*C, */

            /*  where  op( X ) is one of */

            /*     op( X ) = X   or   op( X ) = X', */

            /*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
            /*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

            /*  Arguments */
            /*  ========== */

            /*  TRANSA - CHARACTER*1. */
            /*           On entry, TRANSA specifies the form of op( A ) to be used in */
            /*           the matrix multiplication as follows: */

            /*              TRANSA = 'N' or 'n',  op( A ) = A. */

            /*              TRANSA = 'T' or 't',  op( A ) = A'. */

            /*              TRANSA = 'C' or 'c',  op( A ) = A'. */

            /*           Unchanged on exit. */

            /*  TRANSB - CHARACTER*1. */
            /*           On entry, TRANSB specifies the form of op( B ) to be used in */
            /*           the matrix multiplication as follows: */

            /*              TRANSB = 'N' or 'n',  op( B ) = B. */

            /*              TRANSB = 'T' or 't',  op( B ) = B'. */

            /*              TRANSB = 'C' or 'c',  op( B ) = B'. */

            /*           Unchanged on exit. */

            /*  M      - INTEGER. */
            /*           On entry,  M  specifies  the number  of rows  of the  matrix */
            /*           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
            /*           Unchanged on exit. */

            /*  N      - INTEGER. */
            /*           On entry,  N  specifies the number  of columns of the matrix */
            /*           op( B ) and the number of columns of the matrix C. N must be */
            /*           at least zero. */
            /*           Unchanged on exit. */

            /*  K      - INTEGER. */
            /*           On entry,  K  specifies  the number of columns of the matrix */
            /*           op( A ) and the number of rows of the matrix op( B ). K must */
            /*           be at least  zero. */
            /*           Unchanged on exit. */

            /*  ALPHA  - DOUBLE PRECISION. */
            /*           On entry, ALPHA specifies the scalar alpha. */
            /*           Unchanged on exit. */

            /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is */
            /*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
            /*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
            /*           part of the array  A  must contain the matrix  A,  otherwise */
            /*           the leading  k by m  part of the array  A  must contain  the */
            /*           matrix A. */
            /*           Unchanged on exit. */

            /*  LDA    - INTEGER. */
            /*           On entry, LDA specifies the first dimension of A as declared */
            /*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
            /*           LDA must be at least  Math.Max( 1, m ), otherwise  LDA must be at */
            /*           least  Math.Max( 1, k ). */
            /*           Unchanged on exit. */

            /*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is */
            /*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
            /*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
            /*           part of the array  B  must contain the matrix  B,  otherwise */
            /*           the leading  n by k  part of the array  B  must contain  the */
            /*           matrix B. */
            /*           Unchanged on exit. */

            /*  LDB    - INTEGER. */
            /*           On entry, LDB specifies the first dimension of B as declared */
            /*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
            /*           LDB must be at least  Math.Max( 1, k ), otherwise  LDB must be at */
            /*           least  Math.Max( 1, n ). */
            /*           Unchanged on exit. */

            /*  BETA   - DOUBLE PRECISION. */
            /*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
            /*           supplied as zero then C need not be set on input. */
            /*           Unchanged on exit. */

            /*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
            /*           Before entry, the leading  m by n  part of the array  C must */
            /*           contain the matrix  C,  except when  beta  is zero, in which */
            /*           case C need not be set on entry. */
            /*           On exit, the array  C  is overwritten by the  m by n  matrix */
            /*           ( alpha*op( A )*op( B ) + beta*C ). */

            /*  LDC    - INTEGER. */
            /*           On entry, LDC specifies the first dimension of C as declared */
            /*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
            /*           Math.Max( 1, m ). */
            /*           Unchanged on exit. */


            /*  Level 3 Blas routine. */

            /*  -- Written on 8-February-1989. */
            /*     Jack Dongarra, Argonne National Laboratory. */
            /*     Iain Duff, AERE Harwell. */
            /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
            /*     Sven Hammarling, Numerical Algorithms Group Ltd. */


            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. Parameters .. */
            /*     .. */

            /*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
            /*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
            /*     and  columns of  A  and the  number of  rows  of  B  respectively. */

            /* Parameter adjustments */
            a_dim1 = *lda;
            a_offset = 1 + a_dim1;
            a -= a_offset;
            b_dim1 = *ldb;
            b_offset = 1 + b_dim1;
            b -= b_offset;
            c_dim1 = *ldc;
            c_offset = 1 + c_dim1;
            c__ -= c_offset;

            /* Function Body */
            nota = (*transa == 'N') ? 1 : 0;//lsame_BITA(transa, "N", (ftnlen)1, (ftnlen)1);
            notb = (*transb == 'N') ? 1 : 0;//lsame_BITA(transb, "N", (ftnlen)1, (ftnlen)1);
            if (nota != 0)
            {
                nrowa = *m;
                ncola = *k;
            }
            else
            {
                nrowa = *k;
                ncola = *m;
            }
            if (notb != 0)
            {
                nrowb = *k;
            }
            else
            {
                nrowb = *n;
            }

            /*     Test the input parameters. */

            info = 0;
            if (nota == 0 && *transa != 'C' && *transa != 'T')
            {
                info = 1;
            }
            else if (notb == 0 && *transb != 'C' && *transb != 'T')
            {
                info = 2;
            }
            else if (*m < 0)
            {
                info = 3;
            }
            else if (*n < 0)
            {
                info = 4;
            }
            else if (*k < 0)
            {
                info = 5;
            }
            else if (*lda < Math.Max(1, nrowa))
            {
                info = 8;
            }
            else if (*ldb < Math.Max(1, nrowb))
            {
                info = 10;
            }
            else if (*ldc < Math.Max(1, *m))
            {
                info = 13;
            }
            if (info != 0)
            {
                Console.WriteLine($"DGEMM info {info}");
                //	xerbla_BITA("DGEMM ", &info, (ftnlen)6);
                return 0;
            }

            /*     Quick return if possible. */

            if ((*m == 0 || *n == 0 || *alpha == 0.0 || *k == 0) && *beta == 1.0)
            {
                return 0;
            }

            /*     And if  alpha.eq.zero. */

            if (*alpha == 0.0)
            {
                if (*beta == 0.0)
                {
                    //	    i__1 = *n;
                    //	    for (j = 1; j <= *n; ++j) {
                    //		i__2 = *m;
                    //   memset((double*)(c__ + c_dim1 + 1), 0, *n * *m * sizeof(double));
                    dzerovec(*n * *m, (c__ + c_dim1 + 1));
                    //		for (i__ = 1; i__ <= *m; ++i__) {
                    //		    c__[i__ + j * c_dim1] = 0.;
                    /* L10: */
                    //		}
                    /* L20: */
                    //	    }
                }
                else
                {
                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__) schedule(dynamic)
                    for (j = 1; j <= *n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1; i__ <= *m; ++i__)
                        {
                            c__[i__ + j * c_dim1] *= *beta;
                            /* L30: */
                        }
                        /* L40: */
                    }
                }
                return 0;
            }

            /*     Start the operations. */

            if (notb != 0)
            {
                if (nota != 0)
                {

                    /*           Form  C := alpha*A*B + beta*C. */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,l,temp) schedule(dynamic)
                    for (j = 1; j <= *n; ++j)
                    {
                        if (*beta == 0.0)
                        {
                            //		    i__2 = *m;
                            //		    for (i__ = 1; i__ <= *m; ++i__) {
                            //			c__[i__ + j * c_dim1] = 0.;
                            /* L50: */
                            //       memset((double*)(c__ + j * c_dim1 + 1), 0, *m * sizeof(double));
                            dzerovec(*m, (double*)(c__ + j * c_dim1 + 1));
                            //		    }
                        }
                        else if (*beta != 1.0)
                        {
                            //		    i__2 = *m;
                            for (i__ = 1; i__ <= *m; ++i__)
                            {
                                c__[i__ + j * c_dim1] *= *beta;
                                /* L60: */
                            }
                        }
                        //		i__2 = *k;
                        for (l = 1; l <= *k; ++l)
                        {
                            if (b[l + j * b_dim1] != 0.0)
                            {
                                temp = *alpha * b[l + j * b_dim1];
                                //			i__3 = *m;
                                if (temp > 0)
                                {
                                    for (i__ = 1; i__ <= *m; ++i__)
                                    {
                                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                                    }
                                    /* L70: */
                                }
                            }
                            /* L80: */
                        }
                        /* L90: */
                    }
                }
                else
                {

                    /*           Form  C := alpha*A'*B + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,temp,l)  schedule(dynamic)
                    for (j = 1; j <= *n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1; i__ <= *m; ++i__)
                        {
                            temp = 0.0;
                            //		    i__3 = *k;
                            for (l = 1; l <= *k; ++l)
                            {
                                temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                                /* L100: */
                            }
                            if (*beta == 0.0)
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1] = *alpha * temp;
                            }
                            else
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
                                else if (*beta != 1)
                                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                            }
                            /* L110: */
                        }
                        /* L120: */
                    }
                }
            }
            else
            {
                if (nota != 0)
                {

                    /*           Form  C := alpha*A*B' + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,l,temp) schedule(dynamic)
                    for (j = 1; j <= *n; ++j)
                    {
                        if (*beta == 0.0)
                        {
                            //  memset((double*)(c__ + j * c_dim1 + 1), 0, *m * sizeof(double));
                            dzerovec(*m, (double*)(c__ + j * c_dim1 + 1));
                        }
                        else
                        {
                            dscalvec(*m, *beta, c__ + 1 + j * c_dim1);
                        }
                        for (l = 1; l <= *k; ++l)
                        {
                            if (b[j + l * b_dim1] != 0)
                            {
                                temp = *alpha * b[j + l * b_dim1];
                                daxpyvec(*m, temp, a + 1 + l * a_dim1, c__ + 1 + j * c_dim1);
                            }
                        }
                    }
                }
                else
                {

                    /*           Form  C := alpha*A'*B' + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,l,temp) schedule(dynamic)
                    for (j = 1; j <= *n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1; i__ <= *m; ++i__)
                        {
                            temp = 0.0;
                            //		    i__3 = *k;
                            for (l = 1; l <= *k; ++l)
                            {
                                temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
                                /* L180: */
                            }
                            if (*beta == 0.0)
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1] = *alpha * temp;
                            }
                            else
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
                                else if (*beta != 1)
                                    c__[i__ + j * c_dim1] *= *beta;
                            }
                            /* L190: */
                        }
                        /* L200: */
                    }
                }
            }

            return 0;

            /*     End of DGEMM . */

        } /* dgemm_ */

        public static int dgemm(char[] transa, char[] transb, int m, int
n, int k, double alpha, double[] a, int lda,
double[] b, int ldb, double beta, double[] c__,
int ldc, int astart = 0, int bstart = 0, int cstart = 0)
        {
            /* System generated locals */
            int a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset;// i__1, i__2, 
                                                                     //    i__3;

            /* Local variables */
            int j, info, i__, l; double temp;
            short nota, notb;
            int ncola;
            int nrowa, nrowb;

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  Purpose */
            /*  ======= */

            /*  DGEMM  performs one of the matrix-matrix operations */

            /*     C := alpha*op( A )*op( B ) + beta*C, */

            /*  where  op( X ) is one of */

            /*     op( X ) = X   or   op( X ) = X', */

            /*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
            /*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

            /*  Arguments */
            /*  ========== */

            /*  TRANSA - CHARACTER*1. */
            /*           On entry, TRANSA specifies the form of op( A ) to be used in */
            /*           the matrix multiplication as follows: */

            /*              TRANSA = 'N' or 'n',  op( A ) = A. */

            /*              TRANSA = 'T' or 't',  op( A ) = A'. */

            /*              TRANSA = 'C' or 'c',  op( A ) = A'. */

            /*           Unchanged on exit. */

            /*  TRANSB - CHARACTER*1. */
            /*           On entry, TRANSB specifies the form of op( B ) to be used in */
            /*           the matrix multiplication as follows: */

            /*              TRANSB = 'N' or 'n',  op( B ) = B. */

            /*              TRANSB = 'T' or 't',  op( B ) = B'. */

            /*              TRANSB = 'C' or 'c',  op( B ) = B'. */

            /*           Unchanged on exit. */

            /*  M      - INTEGER. */
            /*           On entry,  M  specifies  the number  of rows  of the  matrix */
            /*           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
            /*           Unchanged on exit. */

            /*  N      - INTEGER. */
            /*           On entry,  N  specifies the number  of columns of the matrix */
            /*           op( B ) and the number of columns of the matrix C. N must be */
            /*           at least zero. */
            /*           Unchanged on exit. */

            /*  K      - INTEGER. */
            /*           On entry,  K  specifies  the number of columns of the matrix */
            /*           op( A ) and the number of rows of the matrix op( B ). K must */
            /*           be at least  zero. */
            /*           Unchanged on exit. */

            /*  ALPHA  - DOUBLE PRECISION. */
            /*           On entry, ALPHA specifies the scalar alpha. */
            /*           Unchanged on exit. */

            /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is */
            /*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
            /*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
            /*           part of the array  A  must contain the matrix  A,  otherwise */
            /*           the leading  k by m  part of the array  A  must contain  the */
            /*           matrix A. */
            /*           Unchanged on exit. */

            /*  LDA    - INTEGER. */
            /*           On entry, LDA specifies the first dimension of A as declared */
            /*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
            /*           LDA must be at least  Math.Max( 1, m ), otherwise  LDA must be at */
            /*           least  Math.Max( 1, k ). */
            /*           Unchanged on exit. */

            /*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is */
            /*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
            /*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
            /*           part of the array  B  must contain the matrix  B,  otherwise */
            /*           the leading  n by k  part of the array  B  must contain  the */
            /*           matrix B. */
            /*           Unchanged on exit. */

            /*  LDB    - INTEGER. */
            /*           On entry, LDB specifies the first dimension of B as declared */
            /*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
            /*           LDB must be at least  Math.Max( 1, k ), otherwise  LDB must be at */
            /*           least  Math.Max( 1, n ). */
            /*           Unchanged on exit. */

            /*  BETA   - DOUBLE PRECISION. */
            /*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
            /*           supplied as zero then C need not be set on input. */
            /*           Unchanged on exit. */

            /*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
            /*           Before entry, the leading  m by n  part of the array  C must */
            /*           contain the matrix  C,  except when  beta  is zero, in which */
            /*           case C need not be set on entry. */
            /*           On exit, the array  C  is overwritten by the  m by n  matrix */
            /*           ( alpha*op( A )*op( B ) + beta*C ). */

            /*  LDC    - INTEGER. */
            /*           On entry, LDC specifies the first dimension of C as declared */
            /*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
            /*           Math.Max( 1, m ). */
            /*           Unchanged on exit. */


            /*  Level 3 Blas routine. */

            /*  -- Written on 8-February-1989. */
            /*     Jack Dongarra, Argonne National Laboratory. */
            /*     Iain Duff, AERE Harwell. */
            /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
            /*     Sven Hammarling, Numerical Algorithms Group Ltd. */


            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. Parameters .. */
            /*     .. */

            /*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
            /*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
            /*     and  columns of  A  and the  number of  rows  of  B  respectively. */

            /* Parameter adjustments */
            a_dim1 = lda;
            a_offset = 1 + a_dim1;
            /*a -= a_offset; */
            astart -= a_offset;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            /*b -= b_offset;*/
            bstart -= b_offset;
            c_dim1 = ldc;
            c_offset = 1 + c_dim1;
            /*c__ -= c_offset;*/
            cstart -= c_offset;

            /* Function Body */
            nota = (transa[0] == 'N') ? 1 : 0;//lsame_BITA(transa, "N", (ftnlen)1, (ftnlen)1);
            notb = (transb[0] == 'N') ? 1 : 0;//lsame_BITA(transb, "N", (ftnlen)1, (ftnlen)1);
            if (nota != 0)
            {
                nrowa = m;
                ncola = k;
            }
            else
            {
                nrowa = k;
                ncola = m;
            }
            if (notb != 0)
            {
                nrowb = k;
            }
            else
            {
                nrowb = n;
            }

            /*     Test the input parameters. */

            info = 0;
            if (nota == 0 && transa[0] != 'C' && transa[0] != 'T')
            {
                info = 1;
            }
            else if (notb == 0 && transb[0] != 'C' && transb[0] != 'T')
            {
                info = 2;
            }
            else if (m < 0)
            {
                info = 3;
            }
            else if (n < 0)
            {
                info = 4;
            }
            else if (k < 0)
            {
                info = 5;
            }
            else if (lda < Math.Max(1, nrowa))
            {
                info = 8;
            }
            else if (ldb < Math.Max(1, nrowb))
            {
                info = 10;
            }
            else if (ldc < Math.Max(1, m))
            {
                info = 13;
            }
            if (info != 0)
            {
                Console.WriteLine($"DGEMM info {info}");
                //	xerbla_BITA("DGEMM ", &info, (ftnlen)6);
                return 0;
            }

            /*     Quick return if possible. */

            if ((m == 0 || n == 0 || alpha == 0.0 || k == 0) && beta == 1.0)
            {
                return 0;
            }

            /*     And if  alpha.eq.zero. */

            if (alpha == 0.0)
            {
                if (beta == 0.0)
                {
                    //	    i__1 = *n;
                    //	    for (j = 1; j <= *n; ++j) {
                    //		i__2 = *m;
                    //   memset((double*)(c__ + c_dim1 + 1), 0, *n * *m * sizeof(double));
                    dzerovec(n * m, (c__ /*+ c_dim1 + 1*/), c_dim1 + 1 + cstart);
                    //		for (i__ = 1; i__ <= *m; ++i__) {
                    //		    c__[i__ + j * c_dim1] = 0.;
                    /* L10: */
                    //		}
                    /* L20: */
                    //	    }
                }
                else
                {
                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__) schedule(dynamic)
                    for (j = 1; j <= n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1 + cstart; i__ <= m + cstart; ++i__)
                        {
                            c__[i__ + j * c_dim1] *= beta;
                            /* L30: */
                        }
                        /* L40: */
                    }
                }
                return 0;
            }

            /*     Start the operations. */

            if (notb != 0)
            {
                if (nota != 0)
                {

                    /*           Form  C := alpha*A*B + beta*C. */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,l,temp) schedule(dynamic)
                    for (j = 1; j <= n; ++j)
                    {
                        if (beta == 0.0)
                        {
                            //		    i__2 = *m;
                            //		    for (i__ = 1; i__ <= *m; ++i__) {
                            //			c__[i__ + j * c_dim1] = 0.;
                            /* L50: */
                            //       memset((double*)(c__ + j * c_dim1 + 1), 0, *m * sizeof(double));
                            dzerovec(m, (c__ /*+ j * c_dim1 + 1*/), j * c_dim1 + 1 + cstart);
                            //		    }
                        }
                        else if (beta != 1.0)
                        {
                            //		    i__2 = *m;
                            for (i__ = 1 + cstart; i__ <= m + cstart; ++i__)
                            {
                                c__[i__ + j * c_dim1] *= beta;
                                /* L60: */
                            }
                        }
                        //		i__2 = *k;
                        for (l = 1; l <= k; ++l)
                        {
                            if (b[l + j * b_dim1 + bstart] != 0.0)
                            {
                                temp = alpha * b[l + j * b_dim1 + bstart];
                                //			i__3 = *m;
                                if (temp > 0)
                                {
                                    for (i__ = 1; i__ <= m; ++i__)
                                    {
                                        c__[i__ + j * c_dim1 + cstart] += temp * a[i__ + l * a_dim1 + astart];
                                    }
                                    /* L70: */
                                }
                            }
                            /* L80: */
                        }
                        /* L90: */
                    }
                }
                else
                {

                    /*           Form  C := alpha*A'*B + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,temp,l)  schedule(dynamic)
                    for (j = 1; j <= n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            temp = 0.0;
                            //		    i__3 = *k;
                            for (l = 1; l <= k; ++l)
                            {
                                temp += a[l + i__ * a_dim1 + astart] * b[l + j * b_dim1 + bstart];
                                /* L100: */
                            }
                            if (beta == 0.0)
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1 + cstart] = alpha * temp;
                            }
                            else
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1 + cstart] = alpha * temp + beta * c__[i__ + j * c_dim1 + cstart];
                                else if (beta != 1)
                                    c__[i__ + j * c_dim1 + cstart] = beta * c__[i__ + j * c_dim1 + cstart];
                            }
                            /* L110: */
                        }
                        /* L120: */
                    }
                }
            }
            else
            {
                if (nota != 0)
                {

                    /*           Form  C := alpha*A*B' + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,l,temp) schedule(dynamic)
                    for (j = 1; j <= n; ++j)
                    {
                        if (beta == 0.0)
                        {
                            //  memset((double*)(c__ + j * c_dim1 + 1), 0, *m * sizeof(double));
                            dzerovec(m, (c__ /*+ j * c_dim1 + 1*/), j * c_dim1 + 1 + cstart);
                        }
                        else
                        {
                            dscalvec(m, beta, c__/* + 1 + j * c_dim1*/, 1 + j * c_dim1 + cstart);
                        }
                        for (l = 1; l <= k; ++l)
                        {
                            if (b[j + l * b_dim1] != 0)
                            {
                                temp = alpha * b[j + l * b_dim1 + bstart];
                                daxpyvec(m, temp, a /*+ 1 + l * a_dim1*/, c__ /*+ 1 + j * c_dim1*/, 1 + l * a_dim1 + astart, 1 + j * c_dim1 + cstart);
                            }
                        }
                    }
                }
                else
                {

                    /*           Form  C := alpha*A'*B' + beta*C */

                    //	    i__1 = *n;
                    //#pragma omp parallel for private(j,i__,l,temp) schedule(dynamic)
                    for (j = 1; j <= n; ++j)
                    {
                        //		i__2 = *m;
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            temp = 0.0;
                            //		    i__3 = *k;
                            for (l = 1; l <= k; ++l)
                            {
                                temp += a[l + i__ * a_dim1 + astart] * b[j + l * b_dim1 + bstart];
                                /* L180: */
                            }
                            if (beta == 0.0)
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1 + cstart] = alpha * temp;
                            }
                            else
                            {
                                if (temp != 0)
                                    c__[i__ + j * c_dim1 + cstart] = alpha * temp + beta * c__[i__ + j * c_dim1 + cstart];
                                else if (beta != 1)
                                    c__[i__ + j * c_dim1 + cstart] *= beta;
                            }
                            /* L190: */
                        }
                        /* L200: */
                    }
                }
            }

            return 0;

            /*     End of DGEMM . */

        } /* dgemm_ */

        public unsafe static int dgemv(char* trans, int m, int n, double alpha,
double* a, int lda, double* x, int incx,
double beta, double* y, int incy)
        {
            if (baseref != 0) Console.WriteLine($"DGEMV alpha {alpha} a {(int)a - baseref} lda {lda} x {(int)x - baseref} incy {incy}");
            /* System generated locals */
            int a_dim1, a_offset, i__1, i__2;

            /* Local variables */
            int i__, j, ix, iy, jx, jy, kx, ky, info;
            double temp;
            int lenx, leny;


            /*  -- Reference BLAS level2 routine (version 3.7.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Parameters .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            a_dim1 = lda;
            a_offset = 1 + a_dim1;
            a -= a_offset;
            --x;
            --y;

            /* Function Body */
            info = 0;
            if (*trans != 'N' && *trans != 'T' && *trans != 'C')
            {
                info = 1;
            }
            else if (m < 0)
            {
                info = 2;
            }
            else if (n < 0)
            {
                info = 3;
            }
            else if (lda < Math.Max(1, m))
            {
                info = 6;
            }
            else if (incx == 0)
            {
                info = 8;
            }
            else if (incy == 0)
            {
                info = 11;
            }
            if (info != 0)
            {
                Console.WriteLine($"In dgemv info is {info}");
                return info;
            }

            /*     Quick return if possible. */

            if (m == 0 || n == 0 || alpha == 0.0 && beta == 1.0)
            {
                return 0;
            }

            /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
            /*     up the start points in  X  and  Y. */

            if (*trans == 'N')
            {
                lenx = n;
                leny = m;
            }
            else
            {
                lenx = m;
                leny = n;
            }
            if (incx > 0)
            {
                kx = 1;
            }
            else
            {
                kx = 1 - (lenx - 1) * incx;
            }
            if (incy > 0)
            {
                ky = 1;
            }
            else
            {
                ky = 1 - (leny - 1) * incy;
            }

            /*     Start the operations. In this version the elements of A are */
            /*     accessed sequentially with one pass through A. */

            /*     First form  y := beta*y. */

            if (beta != 1.0)
            {
                if (incy == 1)
                {
                    if (beta == 0.0)
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[i__] = 0.0;
                            /* L10: */
                        }
                    }
                    else
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[i__] = beta * y[i__];
                            /* L20: */
                        }
                    }
                }
                else
                {
                    iy = ky;
                    if (beta == 0.0)
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[iy] = 0.0;
                            iy += incy;
                            /* L30: */
                        }
                    }
                    else
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[iy] = beta * y[iy];
                            iy += incy;
                            /* L40: */
                        }
                    }
                }
            }
            if (alpha == 0.0)
            {
                return 0;
            }
            if (*trans == 'N')
            {

                /*        Form  y := alpha*A*x + y. */

                jx = kx;
                if (incy == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = (alpha * x[jx]);
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            y[i__] += temp * a[i__ + j * a_dim1];
                            /* L50: */
                        }
                        jx += incx;
                        /* L60: */
                    }
                }
                else
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = (alpha * x[jx]);
                        iy = ky;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            y[iy] += temp * a[i__ + j * a_dim1];
                            iy += incy;
                            /* L70: */
                        }
                        jx += incx;
                        /* L80: */
                    }
                }
            }
            else
            {

                /*        Form  y := alpha*A**T*x + y. */

                jy = ky;
                if (incx == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = 0.0;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            temp += a[i__ + j * a_dim1] * x[i__];
                            /* L90: */
                        }
                        y[jy] += alpha * temp;
                        jy += incy;
                        /* L100: */
                    }
                }
                else
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = 0.0;
                        ix = kx;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            temp += a[i__ + j * a_dim1] * x[ix];
                            ix += incx;
                            /* L110: */
                        }
                        y[jy] += alpha * temp;
                        jy += incy;
                        /* L120: */
                    }
                }
            }

            return 0;

            /*     End of DGEMV . */

        }

        public unsafe static int dgemv(char[] trans, int m, int n, double alpha,
double[] a, int lda, double[] x, int incx,
double beta, double[] y, int incy, int astart = 0, int xstart = 0, int ystart = 0)
        {
            if (baseref != 0)
            {
                fixed (char* tt = trans)
                fixed (double* aa = a)
                fixed (double* xx = x)
                fixed (double* yy = y)
                    return dgemv(tt, m, n, alpha, aa + astart, lda, xx + xstart, incx, beta, yy + ystart, incy);
            }/* System generated locals */
            int a_dim1, a_offset, i__1, i__2;

            /* Local variables */
            int i__, j, ix, iy, jx, jy, kx, ky, info;
            double temp;
            int lenx, leny;


            /*  -- Reference BLAS level2 routine (version 3.7.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Parameters .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            a_dim1 = lda;
            a_offset = 1 + a_dim1;
            /*a -= a_offset;*/
            astart -= a_offset;
            /*--x;*/
            xstart -= 1;
            /*--y;*/
            ystart -= 1;

            /* Function Body */
            info = 0;
            if (trans[0] != 'N' && trans[0] != 'T' && trans[0] != 'C')
            {
                info = 1;
            }
            else if (m < 0)
            {
                info = 2;
            }
            else if (n < 0)
            {
                info = 3;
            }
            else if (lda < Math.Max(1, m))
            {
                info = 6;
            }
            else if (incx == 0)
            {
                info = 8;
            }
            else if (incy == 0)
            {
                info = 11;
            }
            if (info != 0)
            {
                Console.WriteLine($"In dgemv info is {info}");
                return info;
            }

            /*     Quick return if possible. */

            if (m == 0 || n == 0 || alpha == 0.0 && beta == 1.0)
            {
                return 0;
            }

            /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
            /*     up the start points in  X  and  Y. */

            if (trans[0] == 'N')
            {
                lenx = n;
                leny = m;
            }
            else
            {
                lenx = m;
                leny = n;
            }
            if (incx > 0)
            {
                kx = 1;
            }
            else
            {
                kx = 1 - (lenx - 1) * incx;
            }
            if (incy > 0)
            {
                ky = 1;
            }
            else
            {
                ky = 1 - (leny - 1) * incy;
            }

            /*     Start the operations. In this version the elements of A are */
            /*     accessed sequentially with one pass through A. */

            /*     First form  y := beta*y. */

            if (beta != 1.0)
            {
                if (incy == 1)
                {
                    if (beta == 0.0)
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[i__ + ystart] = 0.0;
                            /* L10: */
                        }
                    }
                    else
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[i__ + ystart] = beta * y[i__ + ystart];
                            /* L20: */
                        }
                    }
                }
                else
                {
                    iy = ky;
                    if (beta == 0.0)
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[iy + ystart] = 0.0;
                            iy += incy;
                            /* L30: */
                        }
                    }
                    else
                    {
                        i__1 = leny;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            y[iy + ystart] = beta * y[iy + ystart];
                            iy += incy;
                            /* L40: */
                        }
                    }
                }
            }
            if (alpha == 0.0)
            {
                return 0;
            }
            if (trans[0] == 'N')
            {

                /*        Form  y := alpha*A*x + y. */

                jx = kx;
                if (incy == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = (alpha * x[jx + xstart]);
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            y[i__ + ystart] += temp * a[i__ + j * a_dim1 + astart];
                            /* L50: */
                        }
                        jx += incx;
                        /* L60: */
                    }
                }
                else
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = (alpha * x[jx + xstart]);
                        iy = ky;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            y[iy + ystart] += temp * a[i__ + j * a_dim1 + astart];
                            iy += incy;
                            /* L70: */
                        }
                        jx += incx;
                        /* L80: */
                    }
                }
            }
            else
            {

                /*        Form  y := alpha*A**T*x + y. */

                jy = ky;
                if (incx == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = 0.0;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            temp += a[i__ + j * a_dim1 + astart] * x[i__ + xstart];
                            /* L90: */
                        }
                        y[jy + ystart] += alpha * temp;
                        jy += incy;
                        /* L100: */
                    }
                }
                else
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        temp = 0.0;
                        ix = kx;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            temp += a[i__ + j * a_dim1 + astart] * x[ix + xstart];
                            ix += incx;
                            /* L110: */
                        }
                        y[jy + ystart] += alpha * temp;
                        jy += incy;
                        /* L120: */
                    }
                }
            }

            return 0;

            /*     End of DGEMV . */

        }
        public unsafe static void dspr(char* uplo, int n, double alpha,
double* x, int incx, double* ap)
        {
            int i__, j, k, kk, ix, jx, kx = 0, info;
            double temp;
            /*  -- Reference BLAS level2 routine (version 3.7.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Parameters .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            //       --ap;
            //       --x;

            /* Function Body */
            info = 0;
            if (*uplo != 'U' && *uplo != 'L')
            {
                info = 1;
            }
            else if (n < 0)
            {
                info = 2;
            }
            else if (incx == 0)
            {
                info = 5;
            }
            if (info != 0)
            {
                Console.WriteLine($"dspr: Error {info}");
                return;
            }

            /*     Quick return if possible. */

            if (n == 0 || alpha == 0.0)
            {
                return;
            }

            /*     Set the start point in X if the increment is not unity. */

            if (incx <= 0)
            {
                kx = 1 - (n - 1) * incx;
            }
            else if (incx != 1)
            {
                kx = 1;
            }

            /*     Start the operations. In this version the elements of the array AP */
            /*     are accessed sequentially with one pass through AP. */

            kk = 1;
            if (*uplo == 'U')
            {

                /*        Form  A  when upper triangle is stored in AP. */

                if (incx == 1)
                {
                    for (j = 0; j < n; ++j)
                    {
                        if (x[j] != 0.0)
                        {
                            temp = alpha * x[j];
                            k = kk - 1;
                            for (i__ = 0; i__ <= j; ++i__)
                            {
                                ap[k] += x[i__] * temp;
                                ++k;
                                /* L10: */
                            }
                        }
                        kk += j + 1;
                        /* L20: */
                    }
                }
                else
                {
                    jx = kx - 1;
                    for (j = 0; j < n; ++j)
                    {
                        if (x[jx] != 0.0)
                        {
                            temp = alpha * x[jx];
                            ix = kx - 1;
                            for (k = kk - 1; k < kk + j; ++k)
                            {
                                ap[k] += x[ix] * temp;
                                ix += incx;
                                /* L30: */
                            }
                        }
                        jx += incx;
                        kk += j - 1;
                        /* L40: */
                    }
                }
            }
            else
            {

                /*        Form  A  when lower triangle is stored in AP. */

                if (incx == 1)
                {
                    for (j = 1; j <= n; ++j)
                    {
                        if (x[j - 1] != 0.0)
                        {
                            temp = alpha * x[j - 1];
                            k = kk - 1;
                            for (i__ = j - 1; i__ < n; ++i__)
                            {
                                ap[k] += x[i__] * temp;
                                ++k;
                                /* L50: */
                            }
                        }
                        kk = kk + n - j + 1;
                        /* L60: */
                    }
                }
                else
                {
                    jx = kx - 1;
                    for (j = 1; j <= n; ++j)
                    {
                        if (x[jx] != 0.0)
                        {
                            temp = alpha * x[jx];
                            ix = jx;
                            for (k = kk - 1; k < kk + n - j; ++k)
                            {
                                ap[k] += x[ix] * temp;
                                ix += incx;
                                /* L70: */
                            }
                        }
                        jx += incx;
                        kk = kk + n - j + 1;
                        /* L80: */
                    }
                }
            }
        }
        public static void dspr(char[] uplo, int n, double alpha,
        double[] x, int incx, double[] ap, int xstart = 0, int astart = 0)
        {
            int i__, j, k, kk, ix, jx, kx = 0, info;
            double temp;
            /*  -- Reference BLAS level2 routine (version 3.7.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Parameters .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. External Functions .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            //       --ap;
            //       --x;

            /* Function Body */
            info = 0;
            if (uplo[0] != 'U' && uplo[0] != 'L')
            {
                info = 1;
            }
            else if (n < 0)
            {
                info = 2;
            }
            else if (incx == 0)
            {
                info = 5;
            }
            if (info != 0)
            {
                Console.WriteLine($"dspr: Error {info}");
                return;
            }

            /*     Quick return if possible. */

            if (n == 0 || alpha == 0.0)
            {
                return;
            }

            /*     Set the start point in X if the increment is not unity. */

            if (incx <= 0)
            {
                kx = 1 - (n - 1) * incx;
            }
            else if (incx != 1)
            {
                kx = 1;
            }

            /*     Start the operations. In this version the elements of the array AP */
            /*     are accessed sequentially with one pass through AP. */

            kk = 1;
            if (uplo[0] == 'U')
            {

                /*        Form  A  when upper triangle is stored in AP. */

                if (incx == 1)
                {
                    for (j = 0; j < n; ++j)
                    {
                        if (x[j + xstart] != 0.0)
                        {
                            temp = alpha * x[j + xstart];
                            k = kk - 1;
                            for (i__ = 0; i__ <= j; ++i__)
                            {
                                ap[k + astart] += x[i__ + xstart] * temp;
                                ++k;
                                /* L10: */
                            }
                        }
                        kk += j + 1;
                        /* L20: */
                    }
                }
                else
                {
                    jx = kx - 1;
                    for (j = 0; j < n; ++j)
                    {
                        if (x[jx + xstart] != 0.0)
                        {
                            temp = alpha * x[jx + xstart];
                            ix = kx - 1;
                            for (k = kk - 1; k < kk + j; ++k)
                            {
                                ap[k + astart] += x[ix + xstart] * temp;
                                ix += incx;
                                /* L30: */
                            }
                        }
                        jx += incx;
                        kk += j - 1;
                        /* L40: */
                    }
                }
            }
            else
            {

                /*        Form  A  when lower triangle is stored in AP. */

                if (incx == 1)
                {
                    for (j = 1; j <= n; ++j)
                    {
                        if (x[j - 1 + xstart] != 0.0)
                        {
                            temp = alpha * x[j - 1 + xstart];
                            k = kk - 1;
                            for (i__ = j - 1; i__ < n; ++i__)
                            {
                                ap[k + astart] += x[i__ + xstart] * temp;
                                ++k;
                                /* L50: */
                            }
                        }
                        kk = kk + n - j + 1;
                        /* L60: */
                    }
                }
                else
                {
                    jx = kx - 1;
                    for (j = 1; j <= n; ++j)
                    {
                        if (x[jx + xstart] != 0.0)
                        {
                            temp = alpha * x[jx + xstart];
                            ix = jx;
                            for (k = kk - 1; k < kk + n - j; ++k)
                            {
                                ap[k + astart] += x[ix + xstart] * temp;
                                ix += incx;
                                /* L70: */
                            }
                        }
                        jx += incx;
                        kk = kk + n - j + 1;
                        /* L80: */
                    }
                }
            }
        }
        public unsafe static int idamaxvec1(int n, double* x)
        {
            double mxi = -1.0;
            int m = 1000000000, nn = n;
            double* pxi = x;

            while (n-- > 0)
            {
                double t = Math.Abs(*pxi++);
                if (t > mxi)
                {
                    mxi = t;
                    m = (int)(pxi - x);
                }
            }
            int back = m - 1;
            /*           for (int i = 0; i < nn; ++i)
                       {
                           Console.WriteLine($"iiii {back}   {x[i]}");
                       }*/
            return back;
        }
        public unsafe static int idamaxvec(int n, double* x)
        {
            int back = idamax(n, x, 1);
            /*           for (int i = 0; i < n; ++i)
                       {
                           Console.WriteLine($"iiiiiii {back}   {x[i]}");
                       }*/
            return back;
        }
        public static int idamaxvec(int n, double[] x, int xstart = 0)
        {
            int back = idamax(n, x, 1, xstart);
            /*           for (int i = 0; i < n; ++i)
                       {
                           Console.WriteLine($"iiiiiii {back}   {x[i]}");
                       }*/
            return back;
        }

        public static int idamax(int n, double[] x, int incx, int xstart = 0)
        {
            /* System generated locals */
            int ret_val;

            /* Local variables */
            int i__, ix;
            double dmax__;


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
            //--dx;

            /* Function Body */
            ret_val = 0;
            if (n < 1 || incx <= 0)
                return ret_val;
            ret_val = 0;
            if (n == 1)
                return ret_val;
            if (incx == 1)
            {
                /*        code for increment equal to 1 */
                dmax__ = Math.Abs(x[xstart]);
                for (i__ = 1; i__ < n; ++i__)
                {
                    if ((Math.Abs(x[i__ + xstart])) > dmax__)
                    {
                        ret_val = i__;
                        dmax__ = (Math.Abs(x[i__ + xstart]));
                    }
                }
            }
            else
            {
                /*        code for increment not equal to 1 */
                ix = 1;
                dmax__ = Math.Abs(x[xstart]);
                ix += incx;
                for (i__ = 1; i__ < n; ++i__)
                {
                    if ((Math.Abs(x[ix + xstart])) > dmax__)
                    {
                        ret_val = i__;
                        dmax__ = (Math.Abs(x[ix + xstart]));
                    }
                    ix += incx;
                }
            }
            return ret_val;
        }

        public unsafe static int idamax(int n, double* dx, int incx)
        {
            /* System generated locals */
            int ret_val;

            /* Local variables */
            int i__, ix;
            double dmax__;


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
            ret_val = 0;
            if (n < 1 || incx <= 0)
                return ret_val;
            ret_val = 1;
            if (n == 1)
                return ret_val - 1;// -1 because of FORTRAN style
            if (incx == 1)
            {
                /*        code for increment equal to 1 */
                dmax__ = Math.Abs(dx[1]);
                for (i__ = 2; i__ <= n; ++i__)
                {
                    if ((Math.Abs(dx[i__])) > dmax__)
                    {
                        ret_val = i__;
                        dmax__ = (Math.Abs(dx[i__]));
                    }
                }
            }
            else
            {
                /*        code for increment not equal to 1 */
                ix = 1;
                dmax__ = Math.Abs(dx[1]);
                ix += incx;
                for (i__ = 2; i__ <= n; ++i__)
                {
                    if ((Math.Abs(dx[ix])) > dmax__)
                    {
                        ret_val = i__;
                        dmax__ = (Math.Abs(dx[ix]));
                    }
                    ix += incx;
                }
            }
            return ret_val - 1;// -1 because of FORTRAN style
        }
        public unsafe static void dscal(int n, double da, double* dx, int incx)
        {
            if (baseref != 0) Console.WriteLine($"SCAL da {da} dx {(int)dx - baseref} incx {incx}");
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
            if (n <= 0 || incx <= 0)
            {
                return;
            }
            if (incx == 1)
            {
                m = n % 5;
                if (m != 0)
                {
                    if (da == 0)
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dx[i__] = 0;
                        }
                    }
                    else if (da == -1)
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dx[i__] = -dx[i__];
                        }
                    }
                    else
                    {
                        for (i__ = 1; i__ <= m; ++i__)
                        {
                            dx[i__] = da * dx[i__];
                        }
                    }
                    if (n < 5)
                    {
                        return;
                    }
                }
                mp1 = m + 1;
                if (da == 0)
                {
                    for (i__ = mp1; i__ <= n; i__ += 5)
                    {
                        dx[i__] = 0;
                        dx[i__ + 1] = 0;
                        dx[i__ + 2] = 0;
                        dx[i__ + 3] = 0;
                        dx[i__ + 4] = 0;
                    }
                }
                else if (da == -1)
                {
                    for (i__ = mp1; i__ <= n; i__ += 5)
                    {
                        dx[i__] = -dx[i__];
                        dx[i__ + 1] = -dx[i__ + 1];
                        dx[i__ + 2] = -dx[i__ + 2];
                        dx[i__ + 3] = -dx[i__ + 3];
                        dx[i__ + 4] = -dx[i__ + 4];
                    }
                }
                else
                {
                    for (i__ = mp1; i__ <= n; i__ += 5)
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
                /*        code for increment not equal to 1 */
                nincx = n * incx;
                if (da == 0)
                {
                    for (i__ = 1; incx < 0 ? i__ >= nincx : i__ <= nincx; i__ += incx)
                    {
                        dx[i__] = 0;
                    }
                }
                else if (da == -1)
                {
                    for (i__ = 1; incx < 0 ? i__ >= nincx : i__ <= nincx; i__ += incx)
                    {
                        dx[i__] = -dx[i__];
                    }
                }
                else
                {
                    for (i__ = 1; incx < 0 ? i__ >= nincx : i__ <= nincx; i__ += incx)
                    {
                        dx[i__] = da * dx[i__];
                    }
                }
            }
        }

        public static void dscalvec(int n, double a, double[] x, int xstart = 0)
        {
            dscal(n, a, x, 1, xstart);
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
        public unsafe static double dsum(int n, double* x, int ix = 1, int px = 0)
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

        public unsafe static void dswap(int n, double* dx, int incx,
    double* dy, int incy)
        {
            if (baseref != 0) Console.WriteLine($"SWAP dx {(int)dx - baseref} incx {incx} dy {(int)dy - baseref} incy {incy}");
            /* System generated locals */
            int i__1;

            /* Local variables */
            int i__, m, ix, iy, mp1;
            double dtemp;


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
            --dy;
            --dx;

            /* Function Body */
            if (n <= 0)
            {
                return;
            }
            if (incx == 1 && incy == 1)
            {

                /*       code for both increments equal to 1 */


                /*       clean-up loop */

                m = n % 3;
                if (m != 0)
                {
                    i__1 = m;
                    for (i__ = 1; i__ <= i__1; ++i__)
                    {
                        dtemp = dx[i__];
                        dx[i__] = dy[i__];
                        dy[i__] = dtemp;
                    }
                    if (n < 3)
                    {
                        return;
                    }
                }
                mp1 = m + 1;
                i__1 = n;
                for (i__ = mp1; i__ <= i__1; i__ += 3)
                {
                    dtemp = dx[i__];
                    dx[i__] = dy[i__];
                    dy[i__] = dtemp;
                    dtemp = dx[i__ + 1];
                    dx[i__ + 1] = dy[i__ + 1];
                    dy[i__ + 1] = dtemp;
                    dtemp = dx[i__ + 2];
                    dx[i__ + 2] = dy[i__ + 2];
                    dy[i__ + 2] = dtemp;
                }
            }
            else
            {

                /*       code for unequal increments or equal increments not equal */
                /*         to 1 */

                ix = 1;
                iy = 1;
                if (incx < 0)
                {
                    ix = (-(n) + 1) * incx + 1;
                }
                if (incy < 0)
                {
                    iy = (-(n) + 1) * incy + 1;
                }
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    dtemp = dx[ix];
                    dx[ix] = dy[iy];
                    dy[iy] = dtemp;
                    ix += incx;
                    iy += incy;
                }
            }
            return;
        }
        public unsafe static void dswap(int n, double[] a, int ia, double[] b, int ib, int astart = 0, int bstart = 0)
        {
            if (baseref != 0)
            {
                fixed (double* aa = a)
                fixed (double* bb = b)
                    dswap(n, aa + astart, ia, bb + bstart, ib); return;
            }
            for (int i = 0, iia = ia > 0 ? 0 : (1 - n) * ia, iib = ib > 0 ? 0 : (1 - n) * ib; i < n; i++, iia += ia, iib += ib)
            {
                double temp = a[iia + astart];
                a[iia + astart] = b[iib + bstart];
                b[iib + bstart] = temp;
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
        public unsafe static double didot(int n, double* x, int iix, double* y, int iy    /*increment for y*/    )
        {
            double sum = 0;
            if (n > 0)
            {
                if (iix != 0)
                {
                    if (iix > 0)
                    {
                        while (n-- > 0)
                        {
                            sum += (*x) * (*y);
                            y += iy;
                            x += iix++;
                        }
                    }
                    else
                    {
                        iix = -iix;
                        while (n-- > 0)
                        {
                            sum += (*x) * (*y);
                            y += iy;
                            x += iix--;
                        }
                    }
                }
                else sum = *x * dsum(n, y, iy);
            }
            return sum;
        }

        public static double didot(int n, double[] x, int iix, double[] y, int iy, int xstart = 0, int ystart = 0)
        {
            double sum = 0;
            if (n > 0)
            {
                if (iix != 0)
                {
                    if (iix > 0)
                    {
                        for (int i = 0, yi = 0, xi = 0; i < n; ++i, xi += iix++, yi += iy)
                        {
                            sum += x[xi + xstart] * y[yi + ystart];
                        }
                    }
                    else
                    {
                        iix = -iix;
                        for (int i = 0, yi = 0, xi = 0; i < n; ++i, xi += iix--, yi += iy)
                        {
                            sum += x[xi + xstart] * y[yi + ystart];
                        }
                    }
                }
                else sum = x[0] * dsum(n, y, iy, ystart);
            }
            return sum;
        }
        public unsafe static int dger1(int m, int n, double alpha,
    double* x, int incx, double* y, int incy,
    double* a, int lda)
        {

            if (baseref != 0) Console.WriteLine($"DGER1 alpha {alpha} x {(int)x - baseref} incx {incx} y {(int)y - baseref} incy {incy}");
            /* System generated locals */
            int a_dim1, a_offset, i__1, i__2;

            /* Local variables */
            int i__, j, ix, jy, kx, info;
            double temp;


            /*  -- Reference BLAS level2 routine (version 3.7.0) -- */
            /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  ===================================================================== */

            /*     .. Parameters .. */
            /*     .. */
            /*     .. Local Scalars .. */
            /*     .. */
            /*     .. External Subroutines .. */
            /*     .. */
            /*     .. Intrinsic Functions .. */
            /*     .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            --x;
            --y;
            a_dim1 = lda;
            a_offset = 1 + a_dim1;
            a -= a_offset;

            /* Function Body */
            info = 0;
            if (m < 0)
            {
                info = 1;
            }
            else if (n < 0)
            {
                info = 2;
            }
            else if (incx == 0)
            {
                info = 5;
            }
            else if (incy == 0)
            {
                info = 7;
            }
            else if (lda < Math.Max(1, m))
            {
                info = 9;
            }
            if (info != 0)
            {
                Console.WriteLine($"dger1: info {info}");
                return 0;
            }

            /*     Quick return if possible. */

            if (m == 0 || n == 0 || alpha == 0)
            {
                return 0;
            }

            /*     Start the operations. In this version the elements of A are */
            /*     accessed sequentially with one pass through A. */

            if (incy > 0)
            {
                jy = 1;
            }
            else
            {
                jy = 1 - (n - 1) * incy;
            }
            if (incx == 1)
            {
                i__1 = n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (y[jy] != 0)
                    {
                        temp = alpha * y[jy];
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            a[i__ + j * a_dim1] += x[i__] * temp;
                            /* L10: */
                        }
                    }
                    jy += incy;
                    /* L20: */
                }
            }
            else
            {
                if (incx > 0)
                {
                    kx = 1;
                }
                else
                {
                    kx = 1 - (m - 1) * incx;
                }
                i__1 = n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (y[jy] != 0)
                    {
                        temp = alpha * y[jy];
                        ix = kx;
                        i__2 = m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            a[i__ + j * a_dim1] += x[ix] * temp;
                            ix += incx;
                            /* L30: */
                        }
                    }
                    jy += incy;
                    /* L40: */
                }
            }

            return 0;

            /*     End of DGER  . */

        }
    }
}

using System;
using System.Diagnostics;
namespace Blas
{
    public static class Factorise
    {
        public static double lm_eps8 = 8 * Math.Abs((((double)4) / 3 - 1) * 3 - 1);
        public unsafe static void Writevec(int n, double* b, int ib)
        {
            if (BlasLike.baseref == 0) return;
            for (int i = 0, iib = 0; i < n; i++, iib += ib)
            {
                Console.WriteLine($"{i} -- {b[iib]}");
            }
        }
        public unsafe static void Writevec(int n, double[] b, int ib, int bstart = 0)
        {
            fixed (double* bb = b)
                Writevec(n, bb + bstart, ib);
        }
        public static int Solve1(char[] uplo, int n, int nrhs, double[] ap, int[] ipiv, double[] b, int ldb, int astart = 0, int pstart = 0, int bstart = 0, int root = 0)
        {
            int b_dim1, b_offset;
            int j, k;
            double ak, bk;
            int kc, kp;
            double akm1, bkm1;
            double akm1k;
            double denom;


            /*  -- LAPACK computational routine (version 3.7.0) -- */
            /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */
            /*  Handling for matrix square root Colin Smith February 2021 */
            //--ap;
            astart--;
            pstart--;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            bstart -= b_offset;

            int info = 0;
            if (uplo[0] != 'U' && uplo[0] != 'L')
            {
                info = -1;
            }
            else if (n < 0)
            {
                info = -2;
            }
            else if (nrhs < 0)
            {
                info = -3;
            }
            else if (ldb < Math.Max(1, n))
            {
                info = -7;
            }
            if (info != 0)
            {
                Console.WriteLine($"Solve: error code  {-info}");
                return -info;
            }

            /*     Quick return if possible */

            if (n == 0 || nrhs == 0)
            {
                return 0;
            }

            if (uplo[0] == 'U')
            {
                /*        Solve A*X = B, where A = U*D*U**T. */

                /*        First solve U*D*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */
                k = n;
                kc = n * (n + 1) / 2 + 1;
                while (k > 0)
                {
                    kc -= k;
                    if (ipiv[k + pstart] > 0)
                    {
                        /*           1 x 1 diagonal block */

                        /*           Interchange rows K and IPIV(K). */
                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            BlasLike.dger(k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by (U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 1, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);

                        }
                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 0)
                        {
                            var bot = ap[kc + k - 1 + astart] > 0 ? Math.Max(ap[kc + k - 1 + astart], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(ap[kc + k - 1 + astart], -BlasLike.lm_eps * BlasLike.lm_eps);
                            BlasLike.dscal(nrhs, 1.0 / bot, b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            BlasLike.dscal(nrhs, ap[kc + k - 1 + astart], b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 1)
                        {
                            var bot = Math.Max(ap[kc + k - 1 + astart], 0);
                            if (ap[kc + k - 1 + astart] < 0) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(bot), b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == -1)
                        {
                            var bot = Math.Max(ap[kc + k - 1 + astart], BlasLike.lm_eps * BlasLike.lm_eps);
                            if (ap[kc + k - 1 + astart] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(1.0 / bot), b, ldb, bstart + k + b_dim1);
                        }
                        --k;
                    }
                    else
                    {

                        /*           2 x 2 diagonal block */

                        /*           Interchange rows K-1 and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k - 1)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k - 1 + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in columns K-1 and K of A. */
                            BlasLike.dger(k - 2, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                            BlasLike.dger(k - 2, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc - (k - 1), bstart + k - 1 + b_dim1, bstart + b_dim1 + 1);
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in columns K-1 and K of A. */
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 2, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                            BlasLike.dgemv(TT, k - 2, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc - (k - 1), bstart + k - 1 + b_dim1);
                        }
                        if (root == 30)
                        {
                            /*           Multiply by the inverse of the diagonal block. */

                            akm1k = ap[kc + k - 2 + astart];
                            akm1 = ap[kc - 1 + astart] / akm1k;
                            ak = ap[kc + k - 1 + astart] / akm1k;
                            denom = akm1 * ak - 1.0;
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] / akm1k;
                                bk = b[k + j * b_dim1 + bstart] / akm1k;
                                b[k - 1 + j * b_dim1 + bstart] = (ak * bkm1 - bk) / denom;
                                b[k + j * b_dim1 + bstart] = (akm1 * bk - bkm1) / denom;
                            }
                        }
                        else if (root == 0)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            lambda[0] = lambda[0] > 0 ? Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[0], -BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = lambda[1] > 0 ? Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[1], -BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 /= lambda[0];
                                bk /= lambda[1];
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }

                        else if (root == 2)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 *= lambda[0];
                                bk *= lambda[1];
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 1)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < 0) return -10;
                            if (lambda[1] < 0) return -10;
                            lambda[0] = Math.Max(lambda[0], 0);
                            lambda[1] = Math.Max(lambda[1], 0);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 *= Math.Sqrt(lambda[0]);
                                bk *= Math.Sqrt(lambda[1]);
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == -1)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            if (lambda[1] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            lambda[0] = Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 /= Math.Sqrt(lambda[0]);
                                bk /= Math.Sqrt(lambda[1]);
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }

                        kc = kc - k + 1;
                        k += -2;
                    }
                }
                /*        Next solve U**T*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */
                k = 1;
                kc = 1;

                while (k <= n)
                {
                    if (ipiv[k + pstart] > 0)
                    {
                        /*           1 x 1 diagonal block */

                        /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                        /*           stored in column K of A. */
                        if (root == 0)   //*******
                        {
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                        }
                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        kc += k;
                        ++k;
                    }
                    else
                    {
                        if (root == 0)     //*******
                        {
                            /*           2 x 2 diagonal block */

                            /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                            /*           stored in columns K and K+1 of A. */
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap
                                , 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc + k, bstart + k + 1 + b_dim1);

                        }
                        else if (root == 2)
                        {
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + k, bstart + k + 1 + b_dim1, bstart + b_dim1 + 1);
                        }
                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        kc = kc + (k << 1) + 1;
                        k += 2;
                    }
                }
            }
            else
            {

                /*        Solve A*X = B, where A = L*D*L**T. */

                /*        First solve L*D*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;

                while (k <= n)
                {
                    if (ipiv[k + pstart] > 0)
                    {

                        /*           1 x 1 diagonal block */

                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }

                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by (L(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            }
                        }


                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 0)
                        {
                            var bot = ap[kc + astart] > 0 ? Math.Max(ap[kc + astart], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(ap[kc + astart], -BlasLike.lm_eps * BlasLike.lm_eps);
                            BlasLike.dscal(nrhs, 1.0 / bot, b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            BlasLike.dscal(nrhs, ap[kc + astart], b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 1)
                        {
                            var bot = Math.Max(ap[kc + astart], 0);
                            if (ap[kc + astart] < 0) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(bot), b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == -1)
                        {
                            var bot = Math.Max(ap[kc + astart], BlasLike.lm_eps * BlasLike.lm_eps);
                            if (ap[kc + astart] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(1.0 / bot), b, ldb, bstart + k + b_dim1);
                        }
                        kc = kc + n - k + 1;
                        ++k;
                    }
                    else
                    {

                        /*           2 x 2 diagonal block */

                        /*           Interchange rows K+1 and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k + 1)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + 1 + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in columns K and K+1 of A. */

                            if (k < n - 1)
                            {
                                BlasLike.dger(n - k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + 2, bstart + k + b_dim1, bstart + k + 2 + b_dim1);
                                BlasLike.dger(n - k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + n - k + 2, bstart + k + 1 + b_dim1, bstart + k + 2 + b_dim1);
                            }
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in columns K and K+1 of A. */

                            if (k < n - 1)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k - 1, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 2 + b_dim1, astart + kc + 2, bstart + k + b_dim1);
                                BlasLike.dgemv(TT, n - k - 1, nrhs, 1, b,
                                    ldb, ap, 1, 1, b, ldb, bstart + k + 2 + b_dim1, astart + kc + n - k + 2, bstart + k + 1 + b_dim1);
                            }
                        }

                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 30)
                        {
                            akm1k = ap[kc + 1 + astart];
                            akm1 = ap[kc + astart] / akm1k;
                            ak = ap[kc + n - k + 1 + astart] / akm1k;
                            denom = akm1 * ak - 1.0;
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] / akm1k;
                                bk = b[k + 1 + j * b_dim1 + bstart] / akm1k;
                                b[k + j * b_dim1 + bstart] = (ak * bkm1 - bk) / denom;
                                b[k + 1 + j * b_dim1 + bstart] = (akm1 * bk - bkm1) / denom;
                            }
                        }
                        else if (root == 0)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            lambda[0] = lambda[0] > 0 ? Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[0], -BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = lambda[1] > 0 ? Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[1], -BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 /= lambda[0];
                                bk /= lambda[1];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 1)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < 0) return -10;
                            if (lambda[1] < 0) return -10;
                            lambda[0] = Math.Max(lambda[0], 0);
                            lambda[1] = Math.Max(lambda[1], 0);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 *= Math.Sqrt(lambda[0]);
                                bk *= Math.Sqrt(lambda[1]);
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 2)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 *= lambda[0];
                                bk *= lambda[1];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == -1)
                        {
                            {
                                double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                                var lambda = new double[2];
                                var t = new double[4];
                                Factorise.Eigen2(S, lambda, t);
                                if (lambda[0] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                                if (lambda[1] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                                lambda[0] = Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps);
                                lambda[1] = Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps);
                                for (j = 1; j <= nrhs; ++j)
                                {
                                    bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                    bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                    bkm1 /= Math.Sqrt(lambda[0]);
                                    bk /= Math.Sqrt(lambda[1]);
                                    b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                    b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                                }
                            }
                        }
                        kc = kc + (n - k << 1) + 1;
                        k += 2;
                    }


                }


                /*        Next solve L**T*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;

                while (k >= 1)
                {
                    kc -= n - k + 1;
                    if (ipiv[k + pstart] > 0)
                    {
                        if (root == 0) //******
                        {
                            /*           1 x 1 diagonal block */

                            /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            }
                        }

                        else if (root == 2)
                        {
                            /*           1 x 1 diagonal block */

                            /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }
                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        --k;
                    }
                    else
                    {
                        if (root == 0) //*****
                        {
                            /*           2 x 2 diagonal block */

                            /*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
                            /*           stored in columns K-1 and K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b,
                                    ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc - (n - k), bstart + k - 1 + b_dim1);
                            }
                        }

                        else if (root == 2)
                        {
                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc - (n - k), bstart + k - 1 + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }
                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        kc -= n - k + 2;
                        k += -2;
                    }
                }
            }
            return 0;
        }


        public static int Solve(char[] uplo, int n, int nrhs, double[] ap, int[] ipiv, double[] b, int ldb, int astart = 0, int pstart = 0, int bstart = 0, int root = 0)
        {
            int b_dim1, b_offset;
            int j, k;
            double ak, bk;
            int kc, kp;
            double akm1, bkm1;
            double akm1k;
            double denom;


            /*  -- LAPACK computational routine (version 3.7.0) -- */
            /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
            /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
            /*     December 2016 */
            /*  Handling for matrix square root Colin Smith February 2021 */
            //--ap;
            astart--;
            pstart--;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            bstart -= b_offset;

            int info = 0;
            if (uplo[0] != 'U' && uplo[0] != 'L')
            {
                info = -1;
            }
            else if (n < 0)
            {
                info = -2;
            }
            else if (nrhs < 0)
            {
                info = -3;
            }
            else if (ldb < Math.Max(1, n))
            {
                info = -7;
            }
            if (info != 0)
            {
                Console.WriteLine($"Solve: error code  {-info}");
                return -info;
            }

            /*     Quick return if possible */

            if (n == 0 || nrhs == 0)
            {
                return 0;
            }

            if (uplo[0] == 'U')
            {
                /*        Solve A*X = B, where A = U*D*U**T. */

                /*        First solve U*D*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */
                k = n;
                kc = n * (n + 1) / 2 + 1;
                while (k > 0)
                {
                    kc -= k;
                    if (ipiv[k + pstart] > 0)
                    {
                        /*           1 x 1 diagonal block */

                        /*           Interchange rows K and IPIV(K). */
                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            //                    BlasLike.dger(k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                            //      x     y       a
                            for (int jj = 0, jy = 0; jj < nrhs; ++jj)
                            {
                                if (b[jy + ldb * jj + bstart + k + b_dim1] != 0.0)
                                {
                                    BlasLike.daxpyvec(k - 1, -b[jy + ldb * jj + bstart + k + b_dim1], ap, b, astart + kc, bstart + b_dim1 + 1 + jj * ldb);
                                }
                            }
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by (U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            char[] TT = { 'T' };
                            //BlasLike.dgemv(TT, k - 1, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                            for (int i2 = 0, iy = 0; i2 < ldb; ++i2, iy += ldb)
                            {
                                var sum = 0.0;
                                for (var i1 = 0; i1 < k - 1; ++i1)
                                {
                                    sum += b[bstart + b_offset + i1 + i2 * ldb] * ap[astart + kc + i1];
                                }
                                b[bstart + k + b_dim1 + iy] += sum;
                            }
                        }
                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 0)
                        {
                            var bot = ap[kc + k - 1 + astart] > 0 ? Math.Max(ap[kc + k - 1 + astart], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(ap[kc + k - 1 + astart], -BlasLike.lm_eps * BlasLike.lm_eps);
                            BlasLike.dscal(nrhs, 1.0 / bot, b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            BlasLike.dscal(nrhs, ap[kc + k - 1 + astart], b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 1)
                        {
                            var bot = Math.Max(ap[kc + k - 1 + astart], 0);
                            if (ap[kc + k - 1 + astart] < 0) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(bot), b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == -1)
                        {
                            var bot = Math.Max(ap[kc + k - 1 + astart], BlasLike.lm_eps * BlasLike.lm_eps);
                            if (ap[kc + k - 1 + astart] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(1.0 / bot), b, ldb, bstart + k + b_dim1);
                        }
                        --k;
                    }
                    else
                    {

                        /*           2 x 2 diagonal block */

                        /*           Interchange rows K-1 and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k - 1)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k - 1 + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in columns K-1 and K of A. */
                            BlasLike.dger(k - 2, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                            BlasLike.dger(k - 2, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc - (k - 1), bstart + k - 1 + b_dim1, bstart + b_dim1 + 1);
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in columns K-1 and K of A. */
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 2, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                            BlasLike.dgemv(TT, k - 2, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc - (k - 1), bstart + k - 1 + b_dim1);
                        }
                        if (root == 30)
                        {
                            /*           Multiply by the inverse of the diagonal block. */

                            akm1k = ap[kc + k - 2 + astart];
                            akm1 = ap[kc - 1 + astart] / akm1k;
                            ak = ap[kc + k - 1 + astart] / akm1k;
                            denom = akm1 * ak - 1.0;
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] / akm1k;
                                bk = b[k + j * b_dim1 + bstart] / akm1k;
                                b[k - 1 + j * b_dim1 + bstart] = (ak * bkm1 - bk) / denom;
                                b[k + j * b_dim1 + bstart] = (akm1 * bk - bkm1) / denom;
                            }
                        }
                        else if (root == 0)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            lambda[0] = lambda[0] > 0 ? Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[0], -BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = lambda[1] > 0 ? Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[1], -BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 /= lambda[0];
                                bk /= lambda[1];
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }

                        else if (root == 2)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 *= lambda[0];
                                bk *= lambda[1];
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 1)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < 0) return -10;
                            if (lambda[1] < 0) return -10;
                            lambda[0] = Math.Max(lambda[0], 0);
                            lambda[1] = Math.Max(lambda[1], 0);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 *= Math.Sqrt(lambda[0]);
                                bk *= Math.Sqrt(lambda[1]);
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == -1)
                        {
                            double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            if (lambda[1] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            lambda[0] = Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[1];
                                bk = b[k - 1 + j * b_dim1 + bstart] * t[2] + b[k + j * b_dim1 + bstart] * t[3];
                                bkm1 /= Math.Sqrt(lambda[0]);
                                bk /= Math.Sqrt(lambda[1]);
                                b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }

                        kc = kc - k + 1;
                        k += -2;
                    }
                }
                /*        Next solve U**T*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */
                k = 1;
                kc = 1;

                while (k <= n)
                {
                    if (ipiv[k + pstart] > 0)
                    {
                        /*           1 x 1 diagonal block */

                        /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                        /*           stored in column K of A. */
                        if (root == 0)   //*******
                        {
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            /*           Multiply by inv(U(K)), where U(K) is the transformation */
                            /*           stored in column K of A. */
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                        }
                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        kc += k;
                        ++k;
                    }
                    else
                    {
                        if (root == 0)     //*******
                        {
                            /*           2 x 2 diagonal block */

                            /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                            /*           stored in columns K and K+1 of A. */
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap
                                , 1, 1, b, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                            BlasLike.dgemv(TT, k - 1, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + b_offset, astart + kc + k, bstart + k + 1 + b_dim1);

                        }
                        else if (root == 2)
                        {
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                            BlasLike.dger(k - 1, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + k, bstart + k + 1 + b_dim1, bstart + b_dim1 + 1);
                        }
                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        kc = kc + (k << 1) + 1;
                        k += 2;
                    }
                }
            }
            else
            {

                /*        Solve A*X = B, where A = L*D*L**T. */

                /*        First solve L*D*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;

                while (k <= n)
                {
                    if (ipiv[k + pstart] > 0)
                    {

                        /*           1 x 1 diagonal block */

                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }

                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by (L(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            }
                        }


                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 0)
                        {
                            var bot = ap[kc + astart] > 0 ? Math.Max(ap[kc + astart], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(ap[kc + astart], -BlasLike.lm_eps * BlasLike.lm_eps);
                            BlasLike.dscal(nrhs, 1.0 / bot, b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 2)
                        {
                            BlasLike.dscal(nrhs, ap[kc + astart], b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == 1)
                        {
                            var bot = Math.Max(ap[kc + astart], 0);
                            if (ap[kc + astart] < 0) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(bot), b, ldb, bstart + k + b_dim1);
                        }
                        else if (root == -1)
                        {
                            var bot = Math.Max(ap[kc + astart], BlasLike.lm_eps * BlasLike.lm_eps);
                            if (ap[kc + astart] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                            BlasLike.dscal(nrhs, Math.Sqrt(1.0 / bot), b, ldb, bstart + k + b_dim1);
                        }
                        kc = kc + n - k + 1;
                        ++k;
                    }
                    else
                    {

                        /*           2 x 2 diagonal block */

                        /*           Interchange rows K+1 and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k + 1)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + 1 + b_dim1, bstart + kp + b_dim1);
                        }

                        if (root == 0 || root == -1)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in columns K and K+1 of A. */

                            if (k < n - 1)
                            {
                                BlasLike.dger(n - k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + 2, bstart + k + b_dim1, bstart + k + 2 + b_dim1);
                                BlasLike.dger(n - k - 1, nrhs, -1, ap, 1, b, ldb, b, ldb, astart + kc + n - k + 2, bstart + k + 1 + b_dim1, bstart + k + 2 + b_dim1);
                            }
                        }
                        else if (root == 1 || root == 2)
                        {
                            /*           Multiply by inv(L(K)), where L(K) is the transformation */
                            /*           stored in columns K and K+1 of A. */

                            if (k < n - 1)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k - 1, nrhs, 1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 2 + b_dim1, astart + kc + 2, bstart + k + b_dim1);
                                BlasLike.dgemv(TT, n - k - 1, nrhs, 1, b,
                                    ldb, ap, 1, 1, b, ldb, bstart + k + 2 + b_dim1, astart + kc + n - k + 2, bstart + k + 1 + b_dim1);
                            }
                        }

                        /*           Multiply by the inverse of the diagonal block. */
                        if (root == 30)
                        {
                            akm1k = ap[kc + 1 + astart];
                            akm1 = ap[kc + astart] / akm1k;
                            ak = ap[kc + n - k + 1 + astart] / akm1k;
                            denom = akm1 * ak - 1.0;
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] / akm1k;
                                bk = b[k + 1 + j * b_dim1 + bstart] / akm1k;
                                b[k + j * b_dim1 + bstart] = (ak * bkm1 - bk) / denom;
                                b[k + 1 + j * b_dim1 + bstart] = (akm1 * bk - bkm1) / denom;
                            }
                        }
                        else if (root == 0)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            lambda[0] = lambda[0] > 0 ? Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[0], -BlasLike.lm_eps * BlasLike.lm_eps);
                            lambda[1] = lambda[1] > 0 ? Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps) : Math.Min(lambda[1], -BlasLike.lm_eps * BlasLike.lm_eps);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 /= lambda[0];
                                bk /= lambda[1];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 1)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            if (lambda[0] < 0) return -10;
                            if (lambda[1] < 0) return -10;
                            lambda[0] = Math.Max(lambda[0], 0);
                            lambda[1] = Math.Max(lambda[1], 0);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 *= Math.Sqrt(lambda[0]);
                                bk *= Math.Sqrt(lambda[1]);
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == 2)
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 *= lambda[0];
                                bk *= lambda[1];
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                            }
                        }
                        else if (root == -1)
                        {
                            {
                                double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                                var lambda = new double[2];
                                var t = new double[4];
                                Factorise.Eigen2(S, lambda, t);
                                if (lambda[0] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                                if (lambda[1] < -BlasLike.lm_eps * BlasLike.lm_eps) return -10;
                                lambda[0] = Math.Max(lambda[0], BlasLike.lm_eps * BlasLike.lm_eps);
                                lambda[1] = Math.Max(lambda[1], BlasLike.lm_eps * BlasLike.lm_eps);
                                for (j = 1; j <= nrhs; ++j)
                                {
                                    bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[1];
                                    bk = b[k + j * b_dim1 + bstart] * t[2] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                    bkm1 /= Math.Sqrt(lambda[0]);
                                    bk /= Math.Sqrt(lambda[1]);
                                    b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[2];
                                    b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[1] + bk * t[3];
                                }
                            }
                        }
                        kc = kc + (n - k << 1) + 1;
                        k += 2;
                    }


                }


                /*        Next solve L**T*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;

                while (k >= 1)
                {
                    kc -= n - k + 1;
                    if (ipiv[k + pstart] > 0)
                    {
                        if (root == 0) //******
                        {
                            /*           1 x 1 diagonal block */

                            /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            }
                        }

                        else if (root == 2)
                        {
                            /*           1 x 1 diagonal block */

                            /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                            /*           stored in column K of A. */

                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }
                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        --k;
                    }
                    else
                    {
                        if (root == 0) //*****
                        {
                            /*           2 x 2 diagonal block */

                            /*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
                            /*           stored in columns K-1 and K of A. */

                            if (k < n)
                            {
                                char[] TT = { 'T' };
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b, ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                                BlasLike.dgemv(TT, n - k, nrhs, -1, b,
                                    ldb, ap, 1, 1, b, ldb, bstart + k + 1 + b_dim1, astart + kc - (n - k), bstart + k - 1 + b_dim1);
                            }
                        }

                        else if (root == 2)
                        {
                            if (k < n)
                            {
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                                BlasLike.dger(n - k, nrhs, 1, ap, 1, b, ldb, b, ldb, astart + kc - (n - k), bstart + k - 1 + b_dim1, bstart + k + 1 + b_dim1);
                            }
                        }
                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b, ldb, b, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        }
                        kc -= n - k + 2;
                        k += -2;
                    }
                }
            }
            return 0;
        }

        public static int Factor(char[] uplo, int n, double[] ap, int[] ipiv, int astart = 0, int pstart = 0)
        {
            int info;
            int i__, j, k;
            double t, r1, d11, d12, d21, d22;
            int kc, kk, kp;
            double wk;
            int kx, knc, kpc = -20, npp;// I put -20 where compiler said uninitialised
            double wkm1, wkp1;
            int imax = -20, jmax = -20;
            int kstep;
            double absakk;
            double colmax, rowmax, alpha;

            /*  -- LAPACK computational routine (version 3.7.0) -- */
            /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
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
            /*     .. Executable Statements .. */

            /*     Test the input parameters. */

            /* Parameter adjustments */
            /*--ipiv;
            --ap;*/
            astart--; pstart--;

            /* Function Body */
            info = 0;
            if (uplo[0] != 'U' && uplo[0] != 'L')
            {
                info = -1;
            }
            else if (n < 0)
            {
                info = -2;
            }
            if (info != 0)
            {
                Console.WriteLine($"Factor: Error {info}");
                return info;
            }

            /*     Initialize ALPHA for use in choosing pivot block size. */

            alpha = (Math.Sqrt(17) + 1) / 8;

            if (uplo[0] == 'U')
            {
                /*        Factorize A as U*D*U**T using the upper triangle of A */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2 */
                k = n;
                kc = (n - 1) * n / 2 + 1;
                while (k >= 1)
                {
                    knc = kc;
                    kstep = 1;

                    /*        Determine rows and columns to be interchanged and whether */
                    /*        a 1-by-1 or 2-by-2 pivot block will be used */

                    absakk = (Math.Abs(ap[kc + k - 1 + astart]));

                    /*        IMAX is the row-index of the largest off-diagonal element in */
                    /*        column K, and COLMAX is its absolute value */

                    if (k > 1)
                    {
                        imax = BlasLike.idamax(k - 1, ap, 1, kc + astart) + 1;
                        colmax = (Math.Abs(ap[kc + imax - 1 + astart]));
                    }
                    else
                    {
                        colmax = 0;
                    }
                    if (Math.Max(absakk, colmax) == 0)
                    {
                        /*           Column K is zero: set INFO and continue */
                        if (info == 0)
                        {
                            info = k;
                        }
                        kp = k;
                    }
                    else
                    {
                        if (absakk >= alpha * colmax)
                        {
                            /*              no interchange, use 1-by-1 pivot block */
                            kp = k;
                        }
                        else
                        {
                            rowmax = 0;
                            jmax = imax;
                            kx = imax * (imax + 1) / 2 + imax;
                            for (j = imax + 1; j <= k; ++j)
                            {
                                if ((Math.Abs(ap[kx + astart])) > rowmax)
                                {
                                    rowmax = (Math.Abs(ap[kx + astart]));
                                    jmax = j;
                                }
                                kx += j;
                            }
                            kpc = (imax - 1) * imax / 2 + 1;
                            if (imax > 1)
                            {
                                jmax = BlasLike.idamax(imax - 1, ap, 1, astart + kpc) + 1;
                                rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - 1 + astart]));
                            }

                            if (absakk >= alpha * colmax * (colmax / rowmax))
                            {
                                /*                 no interchange, use 1-by-1 pivot block */
                                kp = k;
                            }
                            else if ((Math.Abs(ap[kpc + imax - 1 + astart])) >= alpha * rowmax)
                            {
                                /*                 interchange rows and columns K and IMAX, use 1-by-1 */
                                /*                 pivot block */
                                kp = imax;
                            }
                            else
                            {
                                /*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
                                /*                 pivot block */
                                kp = imax;
                                kstep = 2;
                            }
                        }
                        kk = k - kstep + 1;
                        if (kstep == 2)
                        {
                            knc = knc - k + 1;
                        }
                        if (kp != kk)
                        {
                            /*              Interchange rows and columns KK and KP in the leading */
                            /*              submatrix A(1:k,1:k) */
                            BlasLike.dswap(kp - 1, ap, 1, ap, 1, astart + knc, astart + kpc);
                            kx = kpc + kp - 1;
                            for (j = kp + 1; j <= (kk - 1); ++j)
                            {
                                kx = kx + j - 1;
                                t = ap[knc + j - 1 + astart];
                                ap[knc + j - 1 + astart] = ap[kx + astart];
                                ap[kx + astart] = t;
                            }
                            t = ap[knc + kk - 1 + astart];
                            ap[knc + kk - 1 + astart] = ap[kpc + kp - 1 + astart];
                            ap[kpc + kp - 1 + astart] = t;
                            if (kstep == 2)
                            {
                                t = ap[kc + k - 2 + astart];
                                ap[kc + k - 2 + astart] = ap[kc + kp - 1 + astart];
                                ap[kc + kp - 1 + astart] = t;
                            }
                        }
                        /*           Update the leading submatrix */
                        if (kstep == 1)
                        {
                            /*              1-by-1 pivot block D(k): column k now holds */
                            /*              W(k) = U(k)*D(k) */
                            /*              where U(k) is the k-th column of U */
                            /*              Perform a rank-1 update of A(1:k-1,1:k-1) as */
                            /*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */
                            r1 = 1 / ap[kc + k - 1 + astart];
                            BlasLike.dspr(uplo, k - 1, -r1, ap, 1, ap, astart + kc, astart + 1);
                            /*              Store U(k) in column k */
                            BlasLike.dscal(k - 1, r1, ap, 1, astart + kc);
                        }
                        else
                        {
                            /*              2-by-2 pivot block D(k): columns k and k-1 now hold */
                            /*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                            /*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                            /*              of U */
                            /*              Perform a rank-2 update of A(1:k-2,1:k-2) as */
                            /*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
                            /*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */
                            if (k > 2)
                            {
                                d12 = ap[k - 1 + (k - 1) * k / 2 + astart];
                                d22 = ap[k - 1 + (k - 2) * (k - 1) / 2 + astart] / d12;
                                d11 = ap[k + (k - 1) * k / 2 + astart] / d12;
                                t = 1 / (d11 * d22 - 1);
                                d12 = t / d12;

                                for (j = k - 2; j >= 1; --j)
                                {
                                    wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2 + astart] -
                                                  ap[j + (k - 1) * k / 2 + astart]);
                                    wk = d12 * (d22 * ap[j + (k - 1) * k / 2 + astart] - ap[j + (k - 2) * (k - 1) / 2 + astart]);
                                    for (i__ = j; i__ >= 1; --i__)
                                    {
                                        ap[i__ + (j - 1) * j / 2 + astart] = ap[i__ + (j - 1) * j /
                                                                                 2 + astart] -
                                                                    ap[i__ + (k - 1) * k / 2 + astart] * wk - ap[i__ + (k - 2) * (k - 1) / 2 + astart] * wkm1;
                                    }
                                    ap[j + (k - 1) * k / 2 + astart] = wk;
                                    ap[j + (k - 2) * (k - 1) / 2 + astart] = wkm1;
                                }
                            }
                        }
                    }
                    /*        Store details of the interchanges in IPIV */
                    if (kstep == 1)
                    {
                        ipiv[k + pstart] = kp;
                    }
                    else
                    {
                        ipiv[k + pstart] = -kp;
                        ipiv[k - 1 + pstart] = -kp;
                    }
                    /*        Decrease K and return to the start of the main loop */
                    k -= kstep;
                    kc = knc - k;
                }
            }
            else
            {
                /*        Factorize A as L*D*L**T using the lower triangle of A */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2 */
                k = 1;
                kc = 1;
                npp = n * (n + 1) / 2;
                while (k <= n)
                {
                    knc = kc;
                    kstep = 1;
                    /*        Determine rows and columns to be interchanged and whether */
                    /*        a 1-by-1 or 2-by-2 pivot block will be used */
                    absakk = (Math.Abs(ap[kc + astart]));
                    /*        IMAX is the row-index of the largest off-diagonal element in */
                    /*        column K, and COLMAX is its absolute value */
                    if (k < n)
                    {
                        imax = k + BlasLike.idamax(n - k, ap, 1, astart + kc + 1) + 1;
                        colmax = (Math.Abs(ap[kc + imax - k + astart]));
                    }
                    else
                    {
                        colmax = 0;
                    }
                    if (Math.Max(absakk, colmax) == 0)
                    {
                        /*           Column K is zero: set INFO and continue */
                        if (info == 0)
                        {
                            info = k;
                        }
                        kp = k;
                    }
                    else
                    {
                        if (absakk >= alpha * colmax)
                        {
                            /*              no interchange, use 1-by-1 pivot block */
                            kp = k;
                        }
                        else
                        {
                            /*              JMAX is the column-index of the largest off-diagonal */
                            /*              element in row IMAX, and ROWMAX is its absolute value */
                            rowmax = 0;
                            kx = kc + imax - k;
                            for (j = k; j <= imax - 1; ++j)
                            {
                                if ((Math.Abs(ap[kx + astart])) > rowmax)
                                {
                                    rowmax = (Math.Abs(ap[kx + astart]));
                                    jmax = j;
                                }
                                kx = kx + n - j;
                            }
                            kpc = npp - (n - imax + 1) * (n - imax + 2) / 2 + 1;
                            if (imax < n)
                            {
                                jmax = imax + BlasLike.idamax(n - imax, ap, 1, astart + kpc + 1) + 1;
                                rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - imax + astart]));
                            }
                            if (absakk >= alpha * colmax * (colmax / rowmax))
                            {
                                /*                 no interchange, use 1-by-1 pivot block */
                                kp = k;
                            }
                            else if ((Math.Abs(ap[kpc + astart])) >= alpha * rowmax)
                            {
                                /*                 interchange rows and columns K and IMAX, use 1-by-1 */
                                /*                 pivot block */
                                kp = imax;
                            }
                            else
                            {
                                /*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
                                /*                 pivot block */
                                kp = imax;
                                kstep = 2;
                            }
                        }
                        kk = k + kstep - 1;
                        if (kstep == 2)
                        {
                            knc = knc + n - k + 1;
                        }
                        if (kp != kk)
                        {
                            /*              Interchange rows and columns KK and KP in the trailing */
                            /*              submatrix A(k:n,k:n) */
                            if (kp < n)
                            {
                                BlasLike.dswap(n - kp, ap, 1, ap, 1, astart + knc + kp - kk + 1, astart + kpc + 1);
                            }
                            kx = knc + kp - kk;
                            for (j = kk + 1; j <= kp - 1; ++j)
                            {
                                kx = kx + n - j + 1;
                                t = ap[knc + j - kk + astart];
                                ap[knc + j - kk + astart] = ap[kx + astart];
                                ap[kx + astart] = t;
                            }
                            t = ap[knc + astart];
                            ap[knc + astart] = ap[kpc + astart];
                            ap[kpc + astart] = t;
                            if (kstep == 2)
                            {
                                t = ap[kc + 1 + astart];
                                ap[kc + 1 + astart] = ap[kc + kp - k + astart];
                                ap[kc + kp - k + astart] = t;
                            }
                        }
                        /*           Update the trailing submatrix */
                        if (kstep == 1)
                        {
                            /*              1-by-1 pivot block D(k): column k now holds */
                            /*              W(k) = L(k)*D(k) */
                            /*              where L(k) is the k-th column of L */
                            if (k < n)
                            {
                                /*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */
                                /*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
                                r1 = 1 / ap[kc + astart];
                                BlasLike.dspr(uplo, n - k, -r1, ap, 1, ap, astart + kc + 1, astart + kc + n - k + 1);
                                /*                 Store L(k) in column K */
                                BlasLike.dscal(n - k, r1, ap, 1, astart + kc + 1);
                            }
                        }
                        else
                        {
                            /*              2-by-2 pivot block D(k): columns K and K+1 now hold */
                            /*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
                            /*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
                            /*              of L */
                            if (k < n - 1)
                            {
                                /*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */
                                /*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
                                /*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */
                                /*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
                                /*                 columns of L */
                                d21 = ap[k + 1 + (k - 1) * ((n << 1) - k) / 2 + astart];
                                d11 = ap[k + 1 + k * ((n << 1) - k - 1) / 2 + astart] / d21;
                                d22 = ap[k + (k - 1) * ((n << 1) - k) / 2 + astart] / d21;
                                t = 1 / (d11 * d22 - 1);
                                d21 = t / d21;
                                for (j = k + 2; j <= n; ++j)
                                {
                                    wk = d21 * (d11 * ap[j + (k - 1) * ((n << 1) - k) /
                                                                 2 + astart] -
                                                ap[j + k * ((n << 1) - k - 1) / 2 + astart]);
                                    wkp1 = d21 * (d22 * ap[j + k * ((n << 1) - k - 1) /
                                                                   2 + astart] -
                                                  ap[j + (k - 1) * ((n << 1) - k) / 2 + astart]);
                                    for (i__ = j; i__ <= n; ++i__)
                                    {
                                        ap[i__ + (j - 1) * ((n << 1) - j) / 2 + astart] = ap[i__ + (j - 1) * ((n << 1) - j) / 2 + astart] - ap[i__ + (k - 1) * ((n << 1) - k) / 2 + astart] * wk -
                                                                                  ap[i__ + k * ((n << 1) - k - 1) / 2 + astart] *
                                                                                      wkp1;
                                    }

                                    ap[j + (k - 1) * ((n << 1) - k) / 2 + astart] = wk;
                                    ap[j + k * ((n << 1) - k - 1) / 2 + astart] = wkp1;
                                }
                            }
                        }
                    }
                    /*        Store details of the interchanges in IPIV */
                    if (kstep == 1)
                    {
                        ipiv[k + pstart] = kp;
                    }
                    else
                    {
                        ipiv[k + pstart] = -kp;
                        ipiv[k + 1 + pstart] = -kp;
                    }
                    /*        Increase K and return to the start of the main loop */
                    k += kstep;
                    kc = knc + n - k + 2;
                }
            }

            return info;
        }
        public static void dsmxmulv(int n, double[] S, double[] x, double[] y, int ystart = 0)
        {
            int i, iS, ix;//This needed change to be compatable with BLAS ddot
            for (i = 1, iS = 0, ix = 0; i <= n; i++, ix++, iS += i)
                y[i - 1 + ystart] = BlasLike.ddot(i, S, -1, x, -1, iS + 1 - i) + BlasLike.didot(n - i, S, i + 1, x, 1, i + iS, 1 + ix);
        }
        public static void dsmxmulvT(int n, double[] S, double[] x, double[] y, int ystart = 0)
        {
            int i, iS, ix;
            for (i = 1, iS = 0, ix = 0; i <= n; i++, ix++, iS += n - i + 2)
                y[i - 1 + ystart] = BlasLike.ddot(n - i + 1, S, 1, x, 1, iS, ix) + BlasLike.didot(i - 1, S, -(n - 1), x, 1, i - 1);
        }
        public unsafe static void dsmxmulvT(int n, double* S, double* x, double* y)
        {
            int i, iS, ix;
            for (i = 1, iS = 0, ix = 0; i <= n; i++, ix++, iS += n - i + 2)
                y[i - 1] = BlasLike.ddot(n - i + 1, S + iS, 1, x + ix, 1) + BlasLike.didot(i - 1, S + i - 1, -(n - 1), x, 1);
        }
        public static void Eigen2(double[] S, double[] lambda, double[] t)
        {//Transformation for diagonalising a symmetric 2x2 matrix
            if (S[1] == 0)
            {
                lambda[0] = S[0];
                lambda[1] = S[2];
                t[0] = t[3] = 1;
                t[1] = t[2] = 0;
                return;
            }
            var d = Math.Sqrt((S[0] - S[2]) * (S[0] - S[2]) + 4 * S[1] * S[1]) / 2;
            var ab = (S[0] + S[2]) / 2;
            if (S[0] > S[2])
            {
                lambda[0] = ab + d;
                lambda[1] = ab - d;
            }
            else
            {
                lambda[0] = ab - d;
                lambda[1] = ab + d;
            }
            t[0] = ((lambda[0]) - S[2]) / S[1];
            t[2] = ((lambda[1]) - S[2]) / S[1];
            t[1] = t[3] = 1;
            var bot = Math.Sqrt(t[0] * t[0] + 1);
            t[0] /= bot;
            t[1] /= bot;
            bot = Math.Sqrt(t[2] * t[2] + 1);
            t[2] /= bot;
            t[3] /= bot;
        }
        public static void dmx_transpose(int n, int m, double[] a, double[] b, int astart = 0, int bstart = 0)
        {
            int size = m * n;
            if (b != a || bstart != astart)
            {
                int pa = astart, pb = bstart, pbmn = pb + size, panm = pa + size, paij;
                while (pb < pbmn) for (paij = pa++; paij < panm; paij += n)
                        b[pb++] = a[paij];
            }
            else if (size > 3)
            {//since a and b are identical bstart is ignored
                int row, column, current;
                size -= 2;
                for (int i = 1; i < size; i++)
                {
                    current = i;
                    do
                    {
                        column = current / m;
                        row = current % m;
                        current = n * row + column;
                    } while (current < i);
                    if (current > i)
                    {
                        double temp = a[i + astart];
                        a[i + astart] = a[current + astart];
                        a[current + astart] = temp;
                    }
                }
            }
        }
        public static void dmxmulv(int n, int m, double[] A, double[] x, double[] y, int astart = 0, int xstart = 0, int ystart = 0, bool atran = false)
        {
            int way = 2;//Which way is fastest?
            if (way == 1)
            {
                BlasLike.dzerovec(n, y, ystart);
                if (atran) for (int j = 0; j < m; ++j) BlasLike.daxpy(n, x[j + xstart], A, m, y, 1, astart + j, ystart);
                else for (int j = 0; j < m; ++j) BlasLike.daxpyvec(n, x[j + xstart], A, y, astart + j * n, ystart);
            }
            else if (way == 2)
            {
                BlasLike.dzerovec(n, y, ystart);
                var N = new char[1];
                N[0] = atran ? 'T' : 'N';
                if (atran) for (int j = 0; j < m; ++j) BlasLike.dgemv(N, m, n, 1, A, m, x, 1, 0, y, 1, astart, xstart, ystart);
                else for (int j = 0; j < m; ++j) BlasLike.dgemv(N, n, m, 1, A, n, x, 1, 0, y, 1, astart, xstart, ystart);
            }
            else
            {
                if (atran) for (int i = 0; i < n; i++) y[i + ystart] = BlasLike.ddot(m, x, 1, A, 1, xstart, i * m + astart);
                else for (int i = 0; i < n; i++) y[i + ystart] = BlasLike.ddot(m, x, 1, A, n, xstart, i + astart);
            }
        }
    }
}

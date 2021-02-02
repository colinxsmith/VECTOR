using System;
using System.Diagnostics;
namespace Blas
{
    public static class Factorise
    {
        public static double lm_eps8 = 8 * Math.Abs((((double)4) / 3 - 1) * 3 - 1);
        public unsafe static int dsptrfold(char* uplo, int n, double* ap, int* ipiv)
        {
            int info;


            int i, j, k;
            double t;
            double r1, d11, d12, d21, d22;
            int kc, kk, kp;
            double wk;
            int kx, knc, kpc = 0, npp;//kpc assigned to make this compile
            double wkm1, wkp1;
            int imax = 0, jmax;//imax assigned to make this compile
            double alpha;
            int kstep;
            double absakk;
            double colmax, rowmax;


            /*  -- LAPACK routine (version 3.0) -- */
            /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
            /*     Courant Institute, Argonne National Lab, and Rice University */
            /*     June 30, 1999 */

            /*     .. Scalar Arguments .. */
            /*     .. */
            /*     .. Array Arguments .. */
            /*     .. */

            /*  Purpose */
            /*  ======= */

            /*  DSPTRF computes the factorization of a real symmetric matrix A stored */
            /*  in packed format using the Bunch-Kaufman diagonal pivoting method: */

            /*     A = U*D*U**T  or  A = L*D*L**T */

            /*  where U (or L) is a product of permutation and unit upper (lower) */
            /*  triangular matrices, and D is symmetric and block diagonal with */
            /*  1-by-1 and 2-by-2 diagonal blocks. */

            /*  Arguments */
            /*  ========= */

            /*  UPLO    (input) CHARACTER*1 */
            /*          = 'U':  Upper triangle of A is stored; */
            /*          = 'L':  Lower triangle of A is stored. */

            /*  N       (input) size_t */
            /*          The order of the matrix A.  N >= 0. */

            /*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2) */
            /*          On entry, the upper or lower triangle of the symmetric matrix */
            /*          A, packed columnwise in a linear array.  The j-th column of A */
            /*          is stored in the array AP as follows: */
            /*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
            /*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */

            /*          On exit, the block diagonal matrix D and the multipliers used */
            /*          to obtain the factor U or L, stored as a packed triangular */
            /*          matrix overwriting A (see below for further details). */

            /*  IPIV    (output) size_t array, dimension (N) */
            /*          Details of the interchanges and the block structure of D. */
            /*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
            /*          interchanged and D(k,k) is a 1-by-1 diagonal block. */
            /*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
            /*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
            /*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) = */
            /*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
            /*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */

            /*  INFO    (output) size_t */
            /*          = 0: successful exit */
            /*          < 0: if INFO = -i, the i-th argument had an illegal value */
            /*          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization */
            /*               has been completed, but the block diagonal matrix D is */
            /*               exactly singular, and division by zero will occur if it */
            /*               is used to solve a system of equations. */

            /*  Further Details */
            /*  =============== */

            /*  5-96 - Based on modifications by J. Lewis, Boeing Computer Services */
            /*         Company */

            /*  If UPLO = 'U', then A = U*D*U', where */
            /*     U = P(n)*U(n)* ... *P(k)U(k)* ..., */
            /*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to */
            /*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
            /*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as */
            /*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such */
            /*  that if the diagonal block D(k) is of order s (s = 1 or 2), then */

            /*             (   I    v    0   )   k-s */
            /*     U(k) =  (   0    I    0   )   s */
            /*             (   0    0    I   )   n-k */
            /*                k-s   s   n-k */

            /*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k). */
            /*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), */
            /*  and A(k,k), and v overwrites A(1:k-2,k-1:k). */

            /*  If UPLO = 'L', then A = L*D*L', where */
            /*     L = P(1)*L(1)* ... *P(k)*L(k)* ..., */
            /*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to */
            /*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
            /*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as */
            /*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such */
            /*  that if the diagonal block D(k) is of order s (s = 1 or 2), then */

            /*             (   I    0     0   )  k-1 */
            /*     L(k) =  (   0    I     0   )  s */
            /*             (   0    v     I   )  n-k-s+1 */
            /*                k-1   s  n-k-s+1 */

            /*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k). */
            /*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k), */
            /*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1). */

            /*  ===================================================================== */

            /* Parameter adjustments */
            //       --ipiv;
            //       --ap;

            /* Function Body */
            info = 0;
            if (*uplo != 'U' && *uplo != 'L')
            {
                info = -1;
            }
            else if (n < 0)
            {
                info = -2;
            }
            if (info != 0)
            {
                Console.WriteLine($"dsptrf: error code {-info}");
                return info;
            }

            /*     Initialize ALPHA for use in choosing pivot block size. */

            alpha = (Math.Sqrt(17.0) + 1.0) / 8.0;

            if (*uplo == 'U')
            {

                /*        Factorize A as U*D*U' using the upper triangle of A */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2 */

                k = n;
                kc = (n - 1) * n / 2 + 1;
            L10:
                knc = kc;

                /*        If K < 1, exit from loop */

                if (k < 1)
                {
                    goto L110;
                }
                kstep = 1;

                /*        Determine rows and columns to be interchanged and whether */
                /*        a 1-by-1 or 2-by-2 pivot block will be used */

                absakk = Math.Abs(ap[kc + k - 2]);

                /*        imax is the row-index of the largest off-diagonal element in */
                /*        column K, and colmax is its absolute value */

                if (k > 1)
                {
                    imax = BlasLike.idamaxvec(k - 1, &ap[-1 + kc]) + 1;
                    colmax = Math.Abs(ap[kc + imax - 2]);
                }
                else
                {
                    colmax = 0.0;
                }

                if (Math.Max(absakk, colmax) == 0.0)
                {

                    /*           Column K is zero: set INFO and continue */

                    if (info == 0)
                    {
                        info = (int)k;
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

                        /*              jmax is the column-index of the largest off-diagonal */
                        /*              element in row imax, and rowmax is its absolute value */

                        rowmax = 0.0;
                        jmax = imax;
                        kx = (int)(imax * (imax + 1) / 2 + imax);
                        for (j = (int)(imax + 1); j <= k; ++j)
                        {
                            if (Math.Abs(ap[-1 + kx]) > rowmax)
                            {
                                rowmax = Math.Abs(ap[-1 + kx]);
                                jmax = j;
                            }
                            kx += j;
                            /* L20: */
                        }
                        kpc = (int)((imax - 1) * imax / 2 + 1);
                        if (imax > 1)
                        {
                            jmax = BlasLike.idamaxvec(imax - 1, &ap[-1 + kpc]) + 1;
                            /* Computing Math.Max */
                            rowmax = Math.Max(rowmax, Math.Abs(ap[-1 + kpc + jmax - 1]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if (Math.Abs(ap[-1 + kpc + imax - 1]) >= alpha *
                          rowmax)
                        {

                            /*                 interchange rows and columns K and imax, use 1-by-1 */
                            /*                 pivot block */

                            kp = (int)imax;
                        }
                        else
                        {

                            /*                 interchange rows and columns K-1 and imax, use 2-by-2 */
                            /*                 pivot block */

                            kp = (int)imax;
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

                        BlasLike.dswap((int)(kp - 1), &ap[-1 + knc], 1, &ap[-1 + kpc], 1);
                        kx = kpc + kp - 1;
                        for (j = kp + 1; j <= kk - 1; ++j)
                        {
                            kx = kx + j - 1;
                            t = ap[-1 + knc + j - 1];
                            ap[-1 + knc + j - 1] = ap[-1 + kx];
                            ap[-1 + kx] = t;
                            /* L30: */
                        }
                        t = ap[-1 + knc + kk - 1];
                        ap[-1 + knc + kk - 1] = ap[-1 + kpc + kp - 1];
                        ap[-1 + kpc + kp - 1] = t;
                        if (kstep == 2)
                        {
                            t = ap[-1 + kc + k - 2];
                            ap[-1 + kc + k - 2] = ap[-1 + kc + kp - 1];
                            ap[-1 + kc + kp - 1] = t;
                        }
                    }

                    /*           Update the leading submatrix */

                    if (kstep == 1)
                    {

                        /*              1-by-1 pivot block D(k): column k now holds */

                        /*              W(k) = U(k)*D(k) */

                        /*              where U(k) is the k-th column of U */

                        /*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

                        /*              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)' */

                        r1 = (Math.Abs(ap[-1 + kc + k - 1]) > lm_eps8 ? 1.0 / ap[-1 + kc + k - 1] : 0);
                        BlasLike.dspr(uplo, (int)(k - 1), -r1, &ap[-1 + kc], 1, &ap[-1 + 1]);

                        /*              Store U(k) in column k */

                        BlasLike.dscal((int)(k - 1), r1, &ap[-1 + kc], 1);
                    }
                    else
                    {

                        /*              2-by-2 pivot block D(k): columns k and k-1 now hold */

                        /*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

                        /*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                        /*              of U */

                        /*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

                        /*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )' */
                        /*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )' */

                        if (k > 2)
                        {
                            double desc;
                            d12 = ap[-1 + k - 1 + (k - 1) * k / 2];
                            d22 = ap[-1 + k - 1 + (k - 2) * (k - 1) / 2] / d12;
                            d11 = ap[-1 + k + (k - 1) * k / 2] / d12;
                            desc = (d11 * d22 - 1.0);
                            t = (Math.Abs(desc) > lm_eps8 ? 1.0 / desc : 0);
                            d12 = t / d12;

                            for (j = k - 2; j >= 1; --j)
                            {
                                wkm1 = d12 * (d11 * ap[-1 + j + (k - 2) * (k - 1) / 2] -
                                    ap[-1 + j + (k - 1) * k / 2]);
                                wk = d12 * (d22 * ap[-1 + j + (k - 1) * k / 2] - ap[-1 + j + (k
                                    - 2) * (k - 1) / 2]);
                                for (i = j; i >= 1; --i)
                                {
                                    ap[-1 + i + (j - 1) * j / 2] = ap[-1 + i + (j - 1) * j /
                                        2] - ap[-1 + i + (k - 1) * k / 2] * wk - ap[-1 +
                                        i + (k - 2) * (k - 1) / 2] * wkm1;
                                    /* L40: */
                                }
                                ap[-1 + j + (k - 1) * k / 2] = wk;
                                ap[-1 + j + (k - 2) * (k - 1) / 2] = wkm1;
                                /* L50: */
                            }

                        }

                    }
                }

                /*        Store details of the interchanges in IPIV */

                if (kstep == 1)
                {
                    ipiv[-1 + k] = kp;
                }
                else
                {
                    ipiv[-1 + k] = -kp;
                    ipiv[-1 + k - 1] = -kp;
                }

                /*        Decrease K and return to the start of the main loop */

                k -= kstep;
                kc = knc - k;
                goto L10;

            }
            else
            {

                /*        Factorize A as L*D*L' using the lower triangle of A */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2 */

                k = 1;
                kc = 1;
                npp = n * (n + 1) / 2;
            L60:
                knc = kc;

                /*        If K > N, exit from loop */

                if (k > n)
                {
                    goto L110;
                }
                kstep = 1;

                /*        Determine rows and columns to be interchanged and whether */
                /*        a 1-by-1 or 2-by-2 pivot block will be used */

                absakk = Math.Abs(ap[-1 + kc]);

                /*        imax is the row-index of the largest off-diagonal element in */
                /*        column K, and colmax is its absolute value */

                if (k < n)
                {
                    imax = k + BlasLike.idamaxvec(n - k, &ap[-1 + kc + 1]) + 1;
                    colmax = Math.Abs(ap[-1 + kc + imax - k]);
                }
                else
                {
                    colmax = 0.0;
                }

                if (Math.Max(absakk, colmax) == 0.0)
                {

                    /*           Column K is zero: set INFO and continue */

                    if (info == 0)
                    {
                        info = (int)k;
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

                        /*              jmax is the column-index of the largest off-diagonal */
                        /*              element in row imax, and rowmax is its absolute value */

                        rowmax = 0.0;
                        kx = kc + (int)imax - k;
                        for (j = k; j <= imax - 1; ++j)
                        {
                            if (Math.Abs(ap[-1 + kx]) > rowmax)
                            {
                                rowmax = Math.Abs(ap[-1 + kx]);
                                jmax = j;
                            }
                            kx = kx + n - j;
                            /* L70: */
                        }
                        kpc = npp - (n - (int)imax + 1) * (n - (int)imax + 2) / 2 + 1;
                        if (imax < n)
                        {
                            jmax = imax + BlasLike.idamaxvec(n - (int)imax, &ap[-1 + kpc + 1]) + 1;
                            /* Computing Math.Max */
                            rowmax = Math.Max(rowmax, Math.Abs(ap[-1 + kpc + jmax - imax]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if (Math.Abs(ap[-1 + kpc]) >= alpha * rowmax)
                        {

                            /*                 interchange rows and columns K and imax, use 1-by-1 */
                            /*                 pivot block */

                            kp = (int)imax;
                        }
                        else
                        {

                            /*                 interchange rows and columns K+1 and imax, use 2-by-2 */
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
                            BlasLike.dswap((int)(n - kp), &ap[-1 + knc + kp - kk + 1], 1, &ap[-1 + kpc + 1],
                                 1);
                        }
                        kx = knc + kp - kk;
                        for (j = kk + 1; j <= kp - 1; ++j)
                        {
                            kx = kx + n - j + 1;
                            t = ap[-1 + knc + j - kk];
                            ap[-1 + knc + j - kk] = ap[-1 + kx];
                            ap[-1 + kx] = t;
                            /* L80: */
                        }
                        t = ap[-1 + knc];
                        ap[-1 + knc] = ap[-1 + kpc];
                        ap[-1 + kpc] = t;
                        if (kstep == 2)
                        {
                            t = ap[-1 + kc + 1];
                            ap[-1 + kc + 1] = ap[-1 + kc + kp - k];
                            ap[-1 + kc + kp - k] = t;
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

                            /*                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)' */

                            r1 = (Math.Abs(ap[-1 + kc]) > lm_eps8 ? 1.0 / ap[-1 + kc] : 0);
                            BlasLike.dspr(uplo, (int)(n - k), -r1, &ap[-1 + kc + 1], 1, &ap[-1 + kc + n
                                - k + 1]);

                            /*                 Store L(k) in column K */

                            BlasLike.dscal((int)(n - k), r1, &ap[-1 + kc + 1], 1);
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

                            /*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )' */
                            /*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )' */

                            d21 = ap[-1 + k + 1 + (k - 1) * ((n << 1) - k) / 2];
                            d11 = ap[-1 + k + 1 + k * ((n << 1) - k - 1) / 2] / d21;
                            d22 = ap[-1 + k + (k - 1) * ((n << 1) - k) / 2] / d21;
                            t = d11 * d22 - 1;
                            t = (Math.Abs(t) > lm_eps8 ? 1.0 / t : 0);
                            d21 = t / d21;

                            for (j = k + 2; j <= n; ++j)
                            {
                                wk = d21 * (d11 * ap[-1 + j + (k - 1) * ((n << 1) - k) /
                                    2] - ap[-1 + j + k * ((n << 1) - k - 1) / 2]);
                                wkp1 = d21 * (d22 * ap[-1 + j + k * ((n << 1) - k - 1) /
                                    2] - ap[-1 + j + (k - 1) * ((n << 1) - k) / 2]);

                                for (i = j; i <= n; ++i)
                                {
                                    ap[-1 + i + (j - 1) * ((n << 1) - j) / 2] = ap[-1 + i
                                        + (j - 1) * ((n << 1) - j) / 2] - ap[-1 + i
                                        + (k - 1) * ((n << 1) - k) / 2] * wk -
                                        ap[-1 + i + k * ((n << 1) - k - 1) / 2] *
                                        wkp1;
                                    /* L90: */
                                }

                                ap[-1 + j + (k - 1) * ((n << 1) - k) / 2] = wk;
                                ap[-1 + j + k * ((n << 1) - k - 1) / 2] = wkp1;

                                /* L100: */
                            }
                        }
                    }
                }

                /*        Store details of the interchanges in IPIV */

                if (kstep == 1)
                {
                    ipiv[-1 + k] = kp;
                }
                else
                {
                    ipiv[-1 + k] = -kp;
                    ipiv[-1 + k + 1] = -kp;
                }

                /*        Increase K and return to the start of the main loop */

                k += kstep;
                kc = knc + n - k + 2;
                goto L60;

            }

        L110:
            return info;
        }
        public unsafe static int dsptrs(char* uplo, int n, int nrhs, double* ap, int* ipiv, double* b, int ldb)
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

            /* Parameter adjustments */
            --ap;
            --ipiv;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            b -= b_offset;

            /* Function Body */
            int info = 0;
            if (*uplo != 'U' && *uplo != 'L')
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
                Console.WriteLine($"dsptrs: error code  {-info}");
                return -info;
            }

            /*     Quick return if possible */

            if (n == 0 || nrhs == 0)
            {
                return 0;
            }

            if (*uplo == 'U')
            {

                /*        Solve A*X = B, where A = U*D*U**T. */

                /*        First solve U*D*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;
            L10:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L30;
                }

                kc -= k;
                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    BlasLike.dger(k - 1, nrhs, -1, &ap[kc], 1, &b[k + b_dim1], ldb, &b[
                        b_dim1 + 1], ldb);
                    Writevec(n, b + b_offset, nrhs);

                    /*           Multiply by the inverse of the diagonal block. */

                    BlasLike.dscal(nrhs, 1.0 / ap[kc + k - 1], &b[k + b_dim1], ldb);
                    --k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Interchange rows K-1 and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k - 1)
                    {
                        BlasLike.dswap(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in columns K-1 and K of A. */
                    BlasLike.dger(k - 2, nrhs, -1, &ap[kc], 1, &b[k + b_dim1], ldb, &b[
                        b_dim1 + 1], ldb);
                    BlasLike.dger(k - 2, nrhs, -1, &ap[kc - (k - 1)], 1, &b[k - 1 +
                        b_dim1], ldb, &b[b_dim1 + 1], ldb);
                    Writevec(n, b + b_offset, nrhs);

                    /*           Multiply by the inverse of the diagonal block. */

                    akm1k = ap[kc + k - 2];
                    akm1 = ap[kc - 1] / akm1k;
                    ak = ap[kc + k - 1] / akm1k;
                    denom = akm1 * ak - 1.0;
                    for (j = 1; j <= nrhs; ++j)
                    {
                        bkm1 = b[k - 1 + j * b_dim1] / akm1k;
                        bk = b[k + j * b_dim1] / akm1k;
                        b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                        /* L20: */
                    }
                    Writevec(n, b + b_offset, nrhs);
                    kc = kc - k + 1;
                    k += -2;
                }

                goto L10;
            L30:

                /*        Next solve U**T*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L40:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L50;
                }

                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    char TT = 'T';
                    BlasLike.dgemv(&TT, k - 1, nrhs, -1, &b[b_offset], ldb, &ap[kc]
                        , 1, 1, &b[k + b_dim1], ldb);
                    Writevec(n, b + b_offset, nrhs);

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }
                    kc += k;
                    ++k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                    /*           stored in columns K and K+1 of A. */
                    char TT = 'T';
                    BlasLike.dgemv(&TT, k - 1, nrhs, -1, &b[b_offset], ldb, &ap[kc]
                        , 1, 1, &b[k + b_dim1], ldb);
                    BlasLike.dgemv(&TT, k - 1, nrhs, -1, &b[b_offset], ldb, &ap[kc
                        + k], 1, 1, &b[k + 1 + b_dim1], ldb);

                    Writevec(n, b + b_offset, nrhs);
                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }
                    kc = kc + (k << 1) + 1;
                    k += 2;
                }

                goto L40;
            L50:

                ;
            }
            else
            {

                /*        Solve A*X = B, where A = L*D*L**T. */

                /*        First solve L*D*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L60:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L80;
                }

                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        BlasLike.dger(n - k, nrhs, -1, &ap[kc + 1], 1, &b[k + b_dim1],
                            ldb, &b[k + 1 + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

                    BlasLike.dscal(nrhs, 1.0 / ap[kc], &b[k + b_dim1], ldb);
                    Writevec(n, b + b_offset, nrhs);
                    kc = kc + n - k + 1;
                    ++k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Interchange rows K+1 and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k + 1)
                    {
                        BlasLike.dswap(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in columns K and K+1 of A. */

                    if (k < n - 1)
                    {
                        BlasLike.dger(n - k - 1, nrhs, -1, &ap[kc + 2], 1, &b[k + b_dim1],
                            ldb, &b[k + 2 + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                        BlasLike.dger(n - k - 1, nrhs, -1, &ap[kc + n - k + 2], 1, &b[k +
                            1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

                    akm1k = ap[kc + 1];
                    akm1 = ap[kc] / akm1k;
                    ak = ap[kc + n - k + 1] / akm1k;
                    denom = akm1 * ak - 1.0;
                    for (j = 1; j <= nrhs; ++j)
                    {
                        bkm1 = b[k + j * b_dim1] / akm1k;
                        bk = b[k + 1 + j * b_dim1] / akm1k;
                        b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                        /* L70: */
                    }
                    Writevec(n, b + b_offset, nrhs);
                    kc = kc + (n - k << 1) + 1;
                    k += 2;
                }

                goto L60;
            L80:

                /*        Next solve L**T*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;
            L90:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L100;
                }

                kc -= n - k + 1;
                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        char TT = 'T';
                        BlasLike.dgemv(&TT, n - k, nrhs, -1, &b[k + 1 + b_dim1],
                            ldb, &ap[kc + 1], 1, 1, &b[k + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }
                    --k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
                    /*           stored in columns K-1 and K of A. */

                    if (k < n)
                    {
                        char TT = 'T';
                        BlasLike.dgemv(&TT, n - k, nrhs, -1, &b[k + 1 + b_dim1],
                            ldb, &ap[kc + 1], 1, 1, &b[k + b_dim1], ldb);
                        BlasLike.dgemv(&TT, n - k, nrhs, -1, &b[k + 1 + b_dim1],
                            ldb, &ap[kc - (n - k)], 1, 1, &b[k - 1 +
                            b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }

                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                        Writevec(n, b + b_offset, nrhs);
                    }
                    kc -= n - k + 2;
                    k += -2;
                }

                goto L90;
            L100:
                ;
            }

            return 0;

            /*     End of DSPTRS */

        }

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
        public static int dsptrs(char[] uplo, int n, int nrhs, double[] ap, int[] ipiv, double[] b, int ldb, int astart = 0, int pstart = 0, int bstart = 0, int root = 0)
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

            /* Parameter adjustments */
            //--ap;
            astart--;
            //--ipiv;
            pstart--;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            //b -= b_offset;
            bstart -= b_offset;

            /* Function Body */
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
                Console.WriteLine($"dsptrs: error code  {-info}");
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
            L10:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L30;
                }

                kc -= k;
                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    BlasLike.dger(k - 1, nrhs, -1, ap/*[kc]*/, 1, b/*[k + b_dim1]*/, ldb, b/*[
                        b_dim1 + 1]*/, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);

                    Writevec(n, b, nrhs);
                    /*           Multiply by the inverse of the diagonal block. */
                    if (root == 0)
                    {
                        BlasLike.dscal(nrhs, 1.0 / ap[kc + k - 1 + astart], b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    else if (root == 1)
                    {
                        BlasLike.dscal(nrhs, Math.Sqrt(ap[kc + k - 1 + astart]), b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    else if (root == -1)
                    {
                        BlasLike.dscal(nrhs, Math.Sqrt(1.0 / ap[kc + k - 1 + astart]), b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    Writevec(n, b, nrhs);
                    --k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Interchange rows K-1 and -IPIV(K). */

                    kp = -ipiv[k + pstart];
                    if (kp != k - 1)
                    {
                        BlasLike.dswap(nrhs, b/*[k - 1 + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k - 1 + b_dim1, bstart + kp + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in columns K-1 and K of A. */
                    BlasLike.dger(k - 2, nrhs, -1, ap/*[kc]*/, 1, b/*[k + b_dim1]*/, ldb, b/*[
                        b_dim1 + 1]*/, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                    BlasLike.dger(k - 2, nrhs, -1, ap/*[kc - (k - 1)]*/, 1, b/*[k - 1 +
                        b_dim1]*/, ldb, b/*[b_dim1 + 1]*/, ldb, astart + kc - (k - 1), bstart + k - 1 + b_dim1, bstart + b_dim1 + 1);
                    Writevec(n, b, nrhs);
                    if (root == 0)
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
                            /* L20: */
                        }
                    }
                    else if (root == 1)
                    {
                        double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                        var lambda = new double[2];
                        var t = new double[4];
                        Factorise.Eigen2(S, lambda, t);
                        for (j = 1; j <= nrhs; ++j)
                        {
                            bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[2];
                            bk = b[k - 1 + j * b_dim1 + bstart] * t[1] + b[k + j * b_dim1 + bstart] * t[3];
                            bkm1 *= Math.Sqrt(lambda[0]);
                            bk *= Math.Sqrt(lambda[1]);
                            b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[1];
                            b[k + j * b_dim1 + bstart] = bkm1 * t[2] + bk * t[3];
                            /* L20: */
                        }
                    }

                    else if (root == -1)
                    {
                        double[] S = { ap[kc - 1 + astart], ap[kc + k - 2 + astart], ap[kc + k - 1 + astart] };
                        var lambda = new double[2];
                        var t = new double[4];
                        Factorise.Eigen2(S, lambda, t);
                        for (j = 1; j <= nrhs; ++j)
                        {
                            bkm1 = b[k - 1 + j * b_dim1 + bstart] * t[0] + b[k + j * b_dim1 + bstart] * t[2];
                            bk = b[k - 1 + j * b_dim1 + bstart] * t[1] + b[k + j * b_dim1 + bstart] * t[3];
                            bkm1 /= Math.Sqrt(lambda[0]);
                            bk /= Math.Sqrt(lambda[1]);
                            b[k - 1 + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[1];
                            b[k + j * b_dim1 + bstart] = bkm1 * t[2] + bk * t[3];
                            /* L20: */
                        }
                    }

                    Writevec(n, b, nrhs);
                    kc = kc - k + 1;
                    k += -2;
                }

                goto L10;
            L30:

                /*        Next solve U**T*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L40:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L50;
                }
                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    if (root == 0)
                    {
                        char[] TT = { 'T' };
                        BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc]*/
                            , 1, 1, b/*[k + b_dim1]*/, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);

                        Writevec(n, b, nrhs);
                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                            Writevec(n, b, nrhs);
                        }
                    }
                    kc += k;
                    ++k;
                }
                else
                {
                    if (root == 0)
                    {
                        /*           2 x 2 diagonal block */

                        /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                        /*           stored in columns K and K+1 of A. */
                        char[] TT = { 'T' };
                        BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc]*/
                            , 1, 1, b/*[k + b_dim1]*/, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                        BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc
                        + k]*/, 1, 1, b/*[k + 1 + b_dim1]*/, ldb, bstart + b_offset, astart + kc + k, bstart + k + 1 + b_dim1);

                        Writevec(n, b, nrhs);
                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                            Writevec(n, b, nrhs);
                        }
                    }
                    kc = kc + (k << 1) + 1;
                    k += 2;
                }

                goto L40;
            L50:

                ;
            }
            else
            {

                /*        Solve A*X = B, where A = L*D*L**T. */

                /*        First solve L*D*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L60:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L80;
                }

                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        BlasLike.dger(n - k, nrhs, -1, ap/*[kc + 1]*/, 1, b/*[k + b_dim1]*/,
                            ldb, b/*[k + 1 + b_dim1]*/, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by the inverse of the diagonal block. */
                    if (root == 0)
                    {
                        BlasLike.dscal(nrhs, 1.0 / ap[kc + astart], b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    else if (root == 1)
                    {
                        BlasLike.dscal(nrhs, Math.Sqrt(ap[kc + astart]), b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    else if (root == -1)
                    {
                        BlasLike.dscal(nrhs, Math.Sqrt(1.0 / ap[kc + astart]), b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    }
                    Writevec(n, b, nrhs);
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
                        BlasLike.dswap(nrhs, b/*[k + 1 + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, bstart + kp + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in columns K and K+1 of A. */

                    if (k < n - 1)
                    {
                        BlasLike.dger(n - k - 1, nrhs, -1, ap/*[kc + 2]*/, 1, b/*[k + b_dim1]*/,
                            ldb, b/*[k + 2 + b_dim1]*/, ldb, astart + kc + 2, bstart + k + b_dim1, bstart + k + 2 + b_dim1);
                        Writevec(n, b, nrhs);
                        BlasLike.dger(n - k - 1, nrhs, -1, ap/*[kc + n - k + 2]*/, 1, b/*[k +
                            1 + b_dim1]*/, ldb, b/*[k + 2 + b_dim1]*/, ldb, astart + kc + n - k + 2, bstart + k + 1 + b_dim1, bstart + k + 2 + b_dim1);
                        Writevec(n, b, nrhs);
                    }

                    /*           Multiply by the inverse of the diagonal block. */
                    if (root == 0)
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
                            /* L70: */
                        }
                    }
                    else if (root == 1)
                    {
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[2];
                                bk = b[k + j * b_dim1 + bstart] * t[1] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 *= Math.Sqrt(lambda[0]);
                                bk *= Math.Sqrt(lambda[1]);
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[1];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[2] + bk * t[3];
                                /* L20: */
                            }
                        }
                    }
                    else if (root == -1)
                    {
                        {
                            double[] S = { ap[kc + astart], ap[kc + 1 + astart], ap[kc + n - k + 1 + astart] };
                            var lambda = new double[2];
                            var t = new double[4];
                            Factorise.Eigen2(S, lambda, t);
                            for (j = 1; j <= nrhs; ++j)
                            {
                                bkm1 = b[k + j * b_dim1 + bstart] * t[0] + b[k + 1 + j * b_dim1 + bstart] * t[2];
                                bk = b[k + j * b_dim1 + bstart] * t[1] + b[k + 1 + j * b_dim1 + bstart] * t[3];
                                bkm1 /= Math.Sqrt(lambda[0]);
                                bk /= Math.Sqrt(lambda[1]);
                                b[k + j * b_dim1 + bstart] = bkm1 * t[0] + bk * t[1];
                                b[k + 1 + j * b_dim1 + bstart] = bkm1 * t[2] + bk * t[3];
                                /* L20: */
                            }
                        }
                   }
                    Writevec(n, b, nrhs);
                    kc = kc + (n - k << 1) + 1;
                    k += 2;
                }

                goto L60;
            L80:

                /*        Next solve L**T*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;
            L90:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L100;
                }

                kc -= n - k + 1;

                if (ipiv[k + pstart] > 0)
                {
                    if (root == 0)
                    {
                        /*           1 x 1 diagonal block */

                        /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                        /*           stored in column K of A. */

                        if (k < n)
                        {
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                                ldb, ap/*[kc + 1]*/, 1, 1, b/*[k + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            Writevec(n, b, nrhs);
                        }

                        /*           Interchange rows K and IPIV(K). */

                        kp = ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                            Writevec(n, b, nrhs);
                        }
                    }
                    --k;
                }
                else
                {
                    if (root == 0)
                    {
                        /*           2 x 2 diagonal block */

                        /*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
                        /*           stored in columns K-1 and K of A. */

                        if (k < n)
                        {
                            char[] TT = { 'T' };
                            BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                                ldb, ap/*[kc + 1]*/, 1, 1, b/*[k + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                            BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                                ldb, ap/*[kc - (n - k)]*/, 1, 1, b/*[k - 1 +
                            b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc - (n - k), bstart + k - 1 + b_dim1);
                            Writevec(n, b, nrhs);
                        }

                        /*           Interchange rows K and -IPIV(K). */

                        kp = -ipiv[k + pstart];
                        if (kp != k)
                        {
                            BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                            Writevec(n, b, nrhs);
                        }
                    }
                    kc -= n - k + 2;
                    k += -2;
                }

                goto L90;
            L100:
                ;
            }

            return 0;

            /*     End of DSPTRS */

        }

        public static int dsptrs1(char[] uplo, int n, int nrhs, double[] ap, int[] ipiv, double[] b, int ldb, int astart = 0, int pstart = 0, int bstart = 0)
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

            /* Parameter adjustments */
            /*--ap;*/
            astart--;
            /*--ipiv;*/
            pstart--;
            b_dim1 = ldb;
            b_offset = 1 + b_dim1;
            /*b -= b_offset;*/
            bstart -= b_offset;

            /* Function Body */
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
                Console.WriteLine($"dsptrs: error code  {-info}");
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
            L10:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L30;
                }

                kc -= k;
                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    BlasLike.dger(k - 1, nrhs, -1, ap/*[kc]*/, 1, b/*[k + b_dim1]*/, ldb, b/*[
                        b_dim1 + 1]*/, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);

                    /*           Multiply by the inverse of the diagonal block. */

                    BlasLike.dscal(nrhs, 1.0 / ap[kc + k - 1 + astart], b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
                    --k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Interchange rows K-1 and -IPIV(K). */

                    kp = -ipiv[k + pstart];
                    if (kp != k - 1)
                    {
                        BlasLike.dswap(nrhs, b/*[k - 1 + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k - 1 + b_dim1, bstart + kp + b_dim1);
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in columns K-1 and K of A. */
                    BlasLike.dger(k - 2, nrhs, -1, ap/*[kc]*/, 1, b/*[k + b_dim1]*/, ldb, b/*[
                        b_dim1 + 1]*/, ldb, astart + kc, bstart + k + b_dim1, bstart + b_dim1 + 1);
                    BlasLike.dger(k - 2, nrhs, -1, ap/*[kc - (k - 1)]*/, 1, b/*[k - 1 +
                        b_dim1]*/, ldb, b/*[b_dim1 + 1]*/, ldb, astart + kc - (k - 1), bstart + k - 1 + b_dim1, bstart + b_dim1 + 1);

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
                        /* L20: */
                    }
                    kc = kc - k + 1;
                    k += -2;
                }

                goto L10;
            L30:

                /*        Next solve U**T*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L40:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L50;
                }

                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */
                    char[] TT = { 'T' };
                    BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc]*/
                        , 1, 1, b/*[k + b_dim1]*/, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }
                    kc += k;
                    ++k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                    /*           stored in columns K and K+1 of A. */
                    char[] TT = { 'T' };
                    BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc]*/
                        , 1, 1, b/*[k + b_dim1]*/, ldb, bstart + b_offset, astart + kc, bstart + k + b_dim1);
                    BlasLike.dgemv(TT, k - 1, nrhs, -1, b/*[b_offset]*/, ldb, ap/*[kc
                        + k]*/, 1, 1, b/*[k + 1 + b_dim1]*/, ldb, bstart + b_offset, astart + kc + k, bstart + k + 1 + b_dim1);

                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }
                    kc = kc + (k << 1) + 1;
                    k += 2;
                }

                goto L40;
            L50:

                ;
            }
            else
            {

                /*        Solve A*X = B, where A = L*D*L**T. */

                /*        First solve L*D*X = B, overwriting B with X. */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = 1;
                kc = 1;
            L60:

                /*        If K > N, exit from loop. */

                if (k > n)
                {
                    goto L80;
                }

                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        BlasLike.dger(n - k, nrhs, -1, ap/*[kc + 1]*/, 1, b/*[k + b_dim1]*/,
                            ldb, b/*[k + 1 + b_dim1]*/, ldb, astart + kc + 1, bstart + k + b_dim1, bstart + k + 1 + b_dim1);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

                    BlasLike.dscal(nrhs, 1.0 / ap[kc + astart], b/*[k + b_dim1]*/, ldb, bstart + k + b_dim1);
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
                        BlasLike.dswap(nrhs, b/*[k + 1 + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, bstart + kp + b_dim1);
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in columns K and K+1 of A. */

                    if (k < n - 1)
                    {
                        BlasLike.dger(n - k - 1, nrhs, -1, ap/*[kc + 2]*/, 1, b/*[k + b_dim1]*/,
                            ldb, b/*[k + 2 + b_dim1]*/, ldb, astart + kc + 2, bstart + k + b_dim1, k + 2 + b_dim1);
                        BlasLike.dger(n - k - 1, nrhs, -1, ap/*[kc + n - k + 2]*/, 1, b/*[k +
                            1 + b_dim1]*/, ldb, b/*[k + 2 + b_dim1]*/, ldb, astart + kc + n - k + 2, bstart + k + 1 + b_dim1, bstart + k + 2 + b_dim1);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

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
                        /* L70: */
                    }
                    kc = kc + (n - k << 1) + 1;
                    k += 2;
                }

                goto L60;
            L80:

                /*        Next solve L**T*X = B, overwriting B with X. */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2, depending on the size of the diagonal blocks. */

                k = n;
                kc = n * (n + 1) / 2 + 1;
            L90:

                /*        If K < 1, exit from loop. */

                if (k < 1)
                {
                    goto L100;
                }

                kc -= n - k + 1;
                if (ipiv[k + pstart] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        char[] TT = { 'T' };
                        BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                            ldb, ap/*[kc + 1]*/, 1, 1, b/*[k + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + b_dim1);
                    }

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }
                    --k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
                    /*           stored in columns K-1 and K of A. */

                    if (k < n)
                    {
                        char[] TT = { 'T' };
                        BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                            ldb, ap/*[kc + 1]*/, 1, 1, b/*[k + b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc + 1, bstart + k + 1 + b_dim1);
                        BlasLike.dgemv(TT, n - k, nrhs, -1, b/*[k + 1 + b_dim1]*/,
                            ldb, ap/*[kc - (n - k)]*/, 1, 1, b/*[k - 1 +
                            b_dim1]*/, ldb, bstart + k + 1 + b_dim1, astart + kc - (n - k), bstart + k - 1 + b_dim1);
                    }

                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k + pstart];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, b/*[k + b_dim1]*/, ldb, b/*[kp + b_dim1]*/, ldb, bstart + k + b_dim1, bstart + kp + b_dim1);
                    }
                    kc -= n - k + 2;
                    k += -2;
                }

                goto L90;
            L100:
                ;
            }

            return 0;

            /*     End of DSPTRS */

        }
        public unsafe static int dsptrf(char* uplo, int n, double* ap, int* ipiv)
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
            --ipiv;
            --ap;

            /* Function Body */
            info = 0;
            if (*uplo != 'U' && *uplo != 'L')
            {
                info = -1;
            }
            else if (n < 0)
            {
                info = -2;
            }
            if (info != 0)
            {
                Console.WriteLine($"dsptrf: Error {info}");
                return info;
            }

            /*     Initialize ALPHA for use in choosing pivot block size. */

            alpha = (Math.Sqrt(17) + 1) / 8;

            if (*uplo == 'U')
            {

                /*        Factorize A as U*D*U**T using the upper triangle of A */

                /*        K is the main loop index, decreasing from N to 1 in steps of */
                /*        1 or 2 */

                k = n;
                kc = (n - 1) * n / 2 + 1;
            L10:
                knc = kc;

                /*        If K < 1, exit from loop */

                if (k < 1)
                {
                    goto L110;
                }
                kstep = 1;

                /*        Determine rows and columns to be interchanged and whether */
                /*        a 1-by-1 or 2-by-2 pivot block will be used */

                absakk = (Math.Abs(ap[kc + k - 1]));

                /*        IMAX is the row-index of the largest off-diagonal element in */
                /*        column K, and COLMAX is its absolute value */

                if (k > 1)
                {
                    imax = BlasLike.idamax(k - 1, &ap[kc], 1) + 1;
                    colmax = (Math.Abs(ap[kc + imax - 1]));
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
                            if ((Math.Abs(ap[kx])) > rowmax)
                            {
                                rowmax = (Math.Abs(ap[kx]));
                                jmax = j;
                            }
                            kx += j;
                            /* L20: */
                        }
                        kpc = (imax - 1) * imax / 2 + 1;
                        if (imax > 1)
                        {
                            jmax = BlasLike.idamax(imax - 1, &ap[kpc], 1) + 1;
                            rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - 1]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if ((Math.Abs(ap[kpc + imax - 1])) >= alpha * rowmax)
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
                        BlasLike.dswap(kp - 1, &ap[knc], 1, &ap[kpc], 1);
                        kx = kpc + kp - 1;
                        for (j = kp + 1; j <= (kk - 1); ++j)
                        {
                            kx = kx + j - 1;
                            t = ap[knc + j - 1];
                            ap[knc + j - 1] = ap[kx];
                            ap[kx] = t;
                            /* L30: */
                        }
                        t = ap[knc + kk - 1];
                        ap[knc + kk - 1] = ap[kpc + kp - 1];
                        ap[kpc + kp - 1] = t;
                        if (kstep == 2)
                        {
                            t = ap[kc + k - 2];
                            ap[kc + k - 2] = ap[kc + kp - 1];
                            ap[kc + kp - 1] = t;
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

                        r1 = 1 / ap[kc + k - 1];
                        BlasLike.dspr(uplo, k - 1, -r1, &ap[kc], 1, &ap[1]);

                        /*              Store U(k) in column k */
                        BlasLike.dscal(k - 1, r1, &ap[kc], 1);
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

                            d12 = ap[k - 1 + (k - 1) * k / 2];
                            d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
                            d11 = ap[k + (k - 1) * k / 2] / d12;
                            t = 1 / (d11 * d22 - 1);
                            d12 = t / d12;

                            for (j = k - 2; j >= 1; --j)
                            {
                                wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] -
                                              ap[j + (k - 1) * k / 2]);
                                wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k - 2) * (k - 1) / 2]);
                                for (i__ = j; i__ >= 1; --i__)
                                {
                                    ap[i__ + (j - 1) * j / 2] = ap[i__ + (j - 1) * j /
                                                                             2] -
                                                                ap[i__ + (k - 1) * k / 2] * wk - ap[i__ + (k - 2) * (k - 1) / 2] * wkm1;
                                    /* L40: */
                                }
                                ap[j + (k - 1) * k / 2] = wk;
                                ap[j + (k - 2) * (k - 1) / 2] = wkm1;
                                /* L50: */
                            }
                        }
                    }
                }

                /*        Store details of the interchanges in IPIV */

                if (kstep == 1)
                {
                    ipiv[k] = kp;
                }
                else
                {
                    ipiv[k] = -kp;
                    ipiv[k - 1] = -kp;
                }

                /*        Decrease K and return to the start of the main loop */

                k -= kstep;
                kc = knc - k;
                goto L10;
            }
            else
            {

                /*        Factorize A as L*D*L**T using the lower triangle of A */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2 */

                k = 1;
                kc = 1;
                npp = n * (n + 1) / 2;
            L60:
                knc = kc;

                /*        If K > N, exit from loop */

                if (k > n)
                {
                    goto L110;
                }
                kstep = 1;

                /*        Determine rows and columns to be interchanged and whether */
                /*        a 1-by-1 or 2-by-2 pivot block will be used */

                absakk = (Math.Abs(ap[kc]));

                /*        IMAX is the row-index of the largest off-diagonal element in */
                /*        column K, and COLMAX is its absolute value */

                if (k < n)
                {
                    imax = k + BlasLike.idamax(n - k, &ap[kc + 1], 1) + 1;
                    colmax = (Math.Abs(ap[kc + imax - k]));
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
                            if ((Math.Abs(ap[kx])) > rowmax)
                            {
                                rowmax = (Math.Abs(ap[kx]));
                                jmax = j;
                            }
                            kx = kx + n - j;
                            /* L70: */
                        }
                        kpc = npp - (n - imax + 1) * (n - imax + 2) / 2 + 1;
                        if (imax < n)
                        {
                            jmax = imax + BlasLike.idamax(n - imax, &ap[kpc + 1], 1) + 1;
                            rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - imax]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if ((Math.Abs(ap[kpc])) >= alpha * rowmax)
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
                            BlasLike.dswap(n - kp, &ap[knc + kp - kk + 1], 1, &ap[kpc + 1], 1);
                        }
                        kx = knc + kp - kk;
                        for (j = kk + 1; j <= kp - 1; ++j)
                        {
                            kx = kx + n - j + 1;
                            t = ap[knc + j - kk];
                            ap[knc + j - kk] = ap[kx];
                            ap[kx] = t;
                            /* L80: */
                        }
                        t = ap[knc];
                        ap[knc] = ap[kpc];
                        ap[kpc] = t;
                        if (kstep == 2)
                        {
                            t = ap[kc + 1];
                            ap[kc + 1] = ap[kc + kp - k];
                            ap[kc + kp - k] = t;
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

                            r1 = 1 / ap[kc];
                            BlasLike.dspr(uplo, n - k, -r1, &ap[kc + 1], 1, &ap[kc + n - k + 1]);

                            /*                 Store L(k) in column K */
                            BlasLike.dscal(n - k, r1, &ap[kc + 1], 1);
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

                            d21 = ap[k + 1 + (k - 1) * ((n << 1) - k) / 2];
                            d11 = ap[k + 1 + k * ((n << 1) - k - 1) / 2] / d21;
                            d22 = ap[k + (k - 1) * ((n << 1) - k) / 2] / d21;
                            t = 1 / (d11 * d22 - 1);
                            d21 = t / d21;
                            for (j = k + 2; j <= n; ++j)
                            {
                                wk = d21 * (d11 * ap[j + (k - 1) * ((n << 1) - k) /
                                                             2] -
                                            ap[j + k * ((n << 1) - k - 1) / 2]);
                                wkp1 = d21 * (d22 * ap[j + k * ((n << 1) - k - 1) /
                                                               2] -
                                              ap[j + (k - 1) * ((n << 1) - k) / 2]);

                                for (i__ = j; i__ <= n; ++i__)
                                {
                                    ap[i__ + (j - 1) * ((n << 1) - j) / 2] = ap[i__ + (j - 1) * ((n << 1) - j) / 2] - ap[i__ + (k - 1) * ((n << 1) - k) / 2] * wk -
                                                                              ap[i__ + k * ((n << 1) - k - 1) / 2] *
                                                                                  wkp1;
                                    /* L90: */
                                }

                                ap[j + (k - 1) * ((n << 1) - k) / 2] = wk;
                                ap[j + k * ((n << 1) - k - 1) / 2] = wkp1;

                                /* L100: */
                            }
                        }
                    }
                }

                /*        Store details of the interchanges in IPIV */

                if (kstep == 1)
                {
                    ipiv[k] = kp;
                }
                else
                {
                    ipiv[k] = -kp;
                    ipiv[k + 1] = -kp;
                }

                /*        Increase K and return to the start of the main loop */

                k += kstep;
                kc = knc + n - k + 2;
                goto L60;
            }

        L110:
            return info;
        }

        public static int dsptrf(char[] uplo, int n, double[] ap, int[] ipiv, int astart = 0, int pstart = 0)
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
                Console.WriteLine($"dsptrf: Error {info}");
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
            L10:
                knc = kc;

                /*        If K < 1, exit from loop */

                if (k < 1)
                {
                    goto L110;
                }
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
                            /* L20: */
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
                            /* L30: */
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
                                    /* L40: */
                                }
                                ap[j + (k - 1) * k / 2 + astart] = wk;
                                ap[j + (k - 2) * (k - 1) / 2 + astart] = wkm1;
                                /* L50: */
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
                goto L10;
            }
            else
            {

                /*        Factorize A as L*D*L**T using the lower triangle of A */

                /*        K is the main loop index, increasing from 1 to N in steps of */
                /*        1 or 2 */

                k = 1;
                kc = 1;
                npp = n * (n + 1) / 2;
            L60:
                knc = kc;

                /*        If K > N, exit from loop */

                if (k > n)
                {
                    goto L110;
                }
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
                            /* L70: */
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
                            /* L80: */
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
                                    /* L90: */
                                }

                                ap[j + (k - 1) * ((n << 1) - k) / 2 + astart] = wk;
                                ap[j + k * ((n << 1) - k - 1) / 2 + astart] = wkp1;

                                /* L100: */
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
                goto L60;
            }

        L110:
            return info;
        }
        public unsafe static void dsmxmulv(int n, double* S, double* x, double* y)
        {
            int i;//This needed change to be compatable with BLAS ddot
            for (i = 1; i <= n; i++, x++, S += i)
                *y++ = BlasLike.ddot(i, S + 1 - i, -1, x + 1 - i, -1) + BlasLike.didot(n - i, S + i, i + 1, x + 1, 1);
        }
        public static void dsmxmulv(int n, double[] S, double[] x, double[] y)
        {
            int i, iS, ix;//This needed change to be compatable with BLAS ddot
            for (i = 1, iS = 0, ix = 0; i <= n; i++, ix++, iS += i)
                y[i - 1] = BlasLike.ddot(i, S, -1, x, -1, iS + 1 - i) + BlasLike.didot(n - i, S, i + 1, x, 1, i + iS, 1 + ix);
        }
        public static void dsmxmulvT(int n, double[] S, double[] x, double[] y)
        {
            int i, iS, ix;
            for (i = 1, iS = 0, ix = 0; i <= n; i++, ix++, iS += n - i + 2)
                y[i - 1] = BlasLike.ddot(n - i + 1, S, 1, x, 1, iS, ix) + BlasLike.didot(i - 1, S, -(n - 1), x, 1, i - 1);
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
    }
}

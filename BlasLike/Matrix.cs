using System;
using System.Diagnostics;
namespace Blas
{
    public static class Factorise
    {
        public static double lm_eps8 = 8 * Math.Abs((4.0 / 3 - 1) * 3 - 1);
        public unsafe static int dsptrf(char* uplo, int n, double* ap, int* ipiv)
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

                absakk = Math.Abs(ap[kc + k - 1]);

                /*        imax is the row-index of the largest off-diagonal element in */
                /*        column K, and colmax is its absolute value */

                if (k > 1)
                {
                    imax = idamaxvec(k - 1, &ap[kc]) + 1;
                    colmax = Math.Abs(ap[kc + imax - 1]);
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
                            if (Math.Abs(ap[kx]) > rowmax)
                            {
                                rowmax = Math.Abs(ap[kx]);
                                jmax = j;
                            }
                            kx += j;
                            /* L20: */
                        }
                        kpc = (int)((imax - 1) * imax / 2 + 1);
                        if (imax > 1)
                        {
                            jmax = idamaxvec(imax - 1, &ap[kpc]) + 1;
                            /* Computing Math.Max */
                            rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - 1]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if (Math.Abs(ap[kpc + imax - 1]) >= alpha *
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

                        BlasLike.dswap((int)(kp - 1), &ap[knc], 1, &ap[kpc], 1);
                        kx = kpc + kp - 1;
                        for (j = kp + 1; j <= kk - 1; ++j)
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

                        /*              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)' */

                        r1 = (Math.Abs(ap[kc + k - 1]) > lm_eps8 ? 1.0 / ap[kc + k - 1] : 0);
                        dspr(uplo, (int)(k - 1), -r1, &ap[kc], 1, &ap[1]);

                        /*              Store U(k) in column k */

                        BlasLike.dscal((int)(k - 1), r1, &ap[kc], 1);
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
                            d12 = ap[k - 1 + (k - 1) * k / 2];
                            d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
                            d11 = ap[k + (k - 1) * k / 2] / d12;
                            desc = (d11 * d22 - 1.0);
                            t = (Math.Abs(desc) > lm_eps8 ? 1.0 / desc : 0);
                            d12 = t / d12;

                            for (j = k - 2; j >= 1; --j)
                            {
                                wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] -
                                    ap[j + (k - 1) * k / 2]);
                                wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k
                                    - 2) * (k - 1) / 2]);
                                for (i = j; i >= 1; --i)
                                {
                                    ap[i + (j - 1) * j / 2] = ap[i + (j - 1) * j /
                                        2] - ap[i + (k - 1) * k / 2] * wk - ap[
                                        i + (k - 2) * (k - 1) / 2] * wkm1;
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

                absakk = Math.Abs(ap[kc]);

                /*        imax is the row-index of the largest off-diagonal element in */
                /*        column K, and colmax is its absolute value */

                if (k < n)
                {
                    imax = k + idamaxvec(n - k, &ap[kc + 1]) + 1;
                    colmax = Math.Abs(ap[kc + imax - k]);
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
                            if (Math.Abs(ap[kx]) > rowmax)
                            {
                                rowmax = Math.Abs(ap[kx]);
                                jmax = j;
                            }
                            kx = kx + n - j;
                            /* L70: */
                        }
                        kpc = npp - (n - (int)imax + 1) * (n - (int)imax + 2) / 2 + 1;
                        if (imax < n)
                        {
                            jmax = imax + idamaxvec(n - (int)imax, &ap[kpc + 1]) + 1;
                            /* Computing Math.Max */
                            rowmax = Math.Max(rowmax, Math.Abs(ap[kpc + jmax - imax]));
                        }

                        if (absakk >= alpha * colmax * (colmax / rowmax))
                        {

                            /*                 no interchange, use 1-by-1 pivot block */

                            kp = k;
                        }
                        else if (Math.Abs(ap[kpc]) >= alpha * rowmax)
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
                            BlasLike.dswap((int)(n - kp), &ap[knc + kp - kk + 1], 1, &ap[kpc + 1],
                                 1);
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

                            /*                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)' */

                            r1 = (Math.Abs(ap[kc]) > lm_eps8 ? 1.0 / ap[kc] : 0);
                            dspr(uplo, (int)(n - k), -r1, &ap[kc + 1], 1, &ap[kc + n
                                - k + 1]);

                            /*                 Store L(k) in column K */

                            BlasLike.dscal((int)(n - k), r1, &ap[kc + 1], 1);
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

                            d21 = ap[k + 1 + (k - 1) * ((n << 1) - k) / 2];
                            d11 = ap[k + 1 + k * ((n << 1) - k - 1) / 2] / d21;
                            d22 = ap[k + (k - 1) * ((n << 1) - k) / 2] / d21;
                            t = d11 * d22 - 1;
                            t = (Math.Abs(t) > lm_eps8 ? 1.0 / t : 0);
                            d21 = t / d21;

                            for (j = k + 2; j <= n; ++j)
                            {
                                wk = d21 * (d11 * ap[j + (k - 1) * ((n << 1) - k) /
                                    2] - ap[j + k * ((n << 1) - k - 1) / 2]);
                                wkp1 = d21 * (d22 * ap[j + k * ((n << 1) - k - 1) /
                                    2] - ap[j + (k - 1) * ((n << 1) - k) / 2]);

                                for (i = j; i <= n; ++i)
                                {
                                    ap[i + (j - 1) * ((n << 1) - j) / 2] = ap[i
                                        + (j - 1) * ((n << 1) - j) / 2] - ap[i
                                        + (k - 1) * ((n << 1) - k) / 2] * wk -
                                        ap[i + k * ((n << 1) - k - 1) / 2] *
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
        public unsafe static void dspr(char* uplo, int n, double alpha,
    double* x, int incx, double* ap)
        {
            /* System generated locals */
            int i__1, i__2;

            /* Local variables */
            int i__, j, k, kk, ix, jx, kx=0, info;
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
            --ap;
            --x;

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
                return ;
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
            if (*uplo=='U')
            {

                /*        Form  A  when upper triangle is stored in AP. */

                if (incx == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (x[j] != 0.0)
                        {
                            temp = alpha * x[j];
                            k = kk;
                            i__2 = j;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                ap[k] += x[i__] * temp;
                                ++k;
                                /* L10: */
                            }
                        }
                        kk += j;
                        /* L20: */
                    }
                }
                else
                {
                    jx = kx;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (x[jx] != 0.0)
                        {
                            temp = alpha * x[jx];
                            ix = kx;
                            i__2 = kk + j - 1;
                            for (k = kk; k <= i__2; ++k)
                            {
                                ap[k] += x[ix] * temp;
                                ix += incx;
                                /* L30: */
                            }
                        }
                        jx += incx;
                        kk += j;
                        /* L40: */
                    }
                }
            }
            else
            {

                /*        Form  A  when lower triangle is stored in AP. */

                if (incx == 1)
                {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (x[j] != 0.0)
                        {
                            temp = alpha * x[j];
                            k = kk;
                            i__2 = n;
                            for (i__ = j; i__ <= i__2; ++i__)
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
                    jx = kx;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (x[jx] != 0.0)
                        {
                            temp = alpha * x[jx];
                            ix = jx;
                            i__2 = kk + n - j;
                            for (k = kk; k <= i__2; ++k)
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

            return ;

            /*     End of DSPR  . */

        }
        public unsafe static int idamaxvec(int n, double* x)
        {
            double mxi = -1.0;
            int m = 1000000000;
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
            return m - 1;
        }
        public unsafe static int dsptrs(char* uplo, int n, int nrhs,
    double* ap, int* ipiv, double* b, int ldb)
        {
            double c_b7 = -1.0;
            int c__1 = 1;
            double c_b19 = 1.0;

            /* System generated locals */
            int b_dim1, b_offset, i__1;
            double d__1;

            /* Local variables */
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
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */

                    i__1 = k - 1;
                    dger(i__1, nrhs, c_b7, &ap[kc], c__1, &b[k + b_dim1], ldb, &b[
                        b_dim1 + 1], ldb);

                    /*           Multiply by the inverse of the diagonal block. */

                    d__1 = 1.0 / ap[kc + k - 1];
                    BlasLike.dscal(nrhs, d__1, &b[k + b_dim1], ldb);
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
                    }

                    /*           Multiply by inv(U(K)), where U(K) is the transformation */
                    /*           stored in columns K-1 and K of A. */

                    i__1 = k - 2;
                    dger(i__1, nrhs, c_b7, &ap[kc], c__1, &b[k + b_dim1], ldb, &b[
                        b_dim1 + 1], ldb);
                    i__1 = k - 2;
                    dger(i__1, nrhs, c_b7, &ap[kc - (k - 1)], c__1, &b[k - 1 +
                        b_dim1], ldb, &b[b_dim1 + 1], ldb);

                    /*           Multiply by the inverse of the diagonal block. */

                    akm1k = ap[kc + k - 2];
                    akm1 = ap[kc - 1] / akm1k;
                    ak = ap[kc + k - 1] / akm1k;
                    denom = akm1 * ak - 1.0;
                    i__1 = nrhs;
                    for (j = 1; j <= i__1; ++j)
                    {
                        bkm1 = b[k - 1 + j * b_dim1] / akm1k;
                        bk = b[k + j * b_dim1] / akm1k;
                        b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
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

                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(U**T(K)), where U(K) is the transformation */
                    /*           stored in column K of A. */

                    i__1 = k - 1;
                    char TT = 'T';
                    dgemv(&TT, i__1, nrhs, c_b7, &b[b_offset], ldb, &ap[kc]
                        , c__1, c_b19, &b[k + b_dim1], ldb);

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                    }
                    kc += k;
                    ++k;
                }
                else
                {

                    /*           2 x 2 diagonal block */

                    /*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
                    /*           stored in columns K and K+1 of A. */

                    i__1 = k - 1;
                    char TT = 'T';
                    dgemv(&TT, i__1, nrhs, c_b7, &b[b_offset], ldb, &ap[kc]
                        , c__1, c_b19, &b[k + b_dim1], ldb);
                    i__1 = k - 1;
                    dgemv(&TT, i__1, nrhs, c_b7, &b[b_offset], ldb, &ap[kc
                        + k], c__1, c_b19, &b[k + 1 + b_dim1], ldb);

                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        i__1 = n - k;
                        dger(i__1, nrhs, c_b7, &ap[kc + 1], c__1, &b[k + b_dim1],
                            ldb, &b[k + 1 + b_dim1], ldb);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

                    d__1 = 1.0 / ap[kc];
                    BlasLike.dscal(nrhs, d__1, &b[k + b_dim1], ldb);
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
                    }

                    /*           Multiply by inv(L(K)), where L(K) is the transformation */
                    /*           stored in columns K and K+1 of A. */

                    if (k < n - 1)
                    {
                        i__1 = n - k - 1;
                        dger(i__1, nrhs, c_b7, &ap[kc + 2], c__1, &b[k + b_dim1],
                            ldb, &b[k + 2 + b_dim1], ldb);
                        i__1 = n - k - 1;
                        dger(i__1, nrhs, c_b7, &ap[kc + n - k + 2], c__1, &b[k +
                            1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
                    }

                    /*           Multiply by the inverse of the diagonal block. */

                    akm1k = ap[kc + 1];
                    akm1 = ap[kc] / akm1k;
                    ak = ap[kc + n - k + 1] / akm1k;
                    denom = akm1 * ak - 1.0;
                    i__1 = nrhs;
                    for (j = 1; j <= i__1; ++j)
                    {
                        bkm1 = b[k + j * b_dim1] / akm1k;
                        bk = b[k + 1 + j * b_dim1] / akm1k;
                        b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
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
                if (ipiv[k] > 0)
                {

                    /*           1 x 1 diagonal block */

                    /*           Multiply by inv(L**T(K)), where L(K) is the transformation */
                    /*           stored in column K of A. */

                    if (k < n)
                    {
                        i__1 = n - k;
                        char TT = 'T';
                        dgemv(&TT, i__1, nrhs, c_b7, &b[k + 1 + b_dim1],
                            ldb, &ap[kc + 1], c__1, c_b19, &b[k + b_dim1], ldb);
                    }

                    /*           Interchange rows K and IPIV(K). */

                    kp = ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                        i__1 = n - k;
                        char TT = 'T';
                        dgemv(&TT, i__1, nrhs, c_b7, &b[k + 1 + b_dim1],
                            ldb, &ap[kc + 1], c__1, c_b19, &b[k + b_dim1], ldb);
                        i__1 = n - k;
                        dgemv(&TT, i__1, nrhs, c_b7, &b[k + 1 + b_dim1],
                            ldb, &ap[kc - (n - k)], c__1, c_b19, &b[k - 1 +
                            b_dim1], ldb);
                    }

                    /*           Interchange rows K and -IPIV(K). */

                    kp = -ipiv[k];
                    if (kp != k)
                    {
                        BlasLike.dswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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

        } /* dsptrs_ */
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
        public unsafe static void dger1(int m, int n, double alpha,
    double* x, int incx, double* y, int incy,
    double* a, int lda)
        {
            int a_dim1;
            //    size_t i, 
            int jy, kx, info;
            long j;
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
                Console.WriteLine($"dger error code {info}");
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
                        if (y[jy + incy * j] != 0.0)
                        {
                            BlasLike.daxpyvec(m, alpha * y[jy + incy * j], x, a + j * a_dim1);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j] != 0.0)
                        {
                            BlasLike.daxpyvec(m, alpha * y[jy + incy * j], x, a + j * a_dim1);
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
                        if (y[jy + incy * j] != 0.0)
                        {
                            BlasLike.daxpy(m, alpha * y[jy + incy * j], x + kx, incx, a + j * a_dim1, 1);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < n; ++j/*,py+=incy*/)
                    {
                        if (y[jy + incy * j] != 0.0)
                        {
                            BlasLike.daxpy(m, alpha * y[jy + incy * j], x + kx, incx, a + j * a_dim1, 1);
                        }
                    }
                }
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
            if (nota == 0 && (*transa == 'C' ? 1 : 0) == 0/*lsame_BITA(transa, "C", (ftnlen)1, (ftnlen)1) */&& *transa != 'T'/*lsame_BITA(
	    transa, "T", (ftnlen)1, (ftnlen)1)*/)
            {
                info = 1;
            }
            else if (notb == 0 && (*transb == 'C' ? 1 : 0) == 0 /*lsame_BITA(transb, "C", (ftnlen)1, (ftnlen)1)*/ && *transb != 'T'
              /*lsame_BITA(transb, "T", (ftnlen)1, (ftnlen)1)*/)
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
                    BlasLike.dzerovec(*n * *m, (c__ + c_dim1 + 1));
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
                            BlasLike.dzerovec(*m, (double*)(c__ + j * c_dim1 + 1));
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
                            BlasLike.dzerovec(*m, (double*)(c__ + j * c_dim1 + 1));
                        }
                        else
                        {
                            BlasLike.dscalvec(*m, *beta, c__ + 1 + j * c_dim1);
                        }
                        for (l = 1; l <= *k; ++l)
                        {
                            if (b[j + l * b_dim1] != 0)
                            {
                                temp = *alpha * b[j + l * b_dim1];
                                BlasLike.daxpyvec(*m, temp, a + 1 + l * a_dim1, c__ + 1 + j * c_dim1);
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
        public unsafe static int dgemv(char* trans, int m, int n, double alpha,
        double* a, int lda, double* x, int incx,
            double beta, double* y, int incy)
        {
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

        } /* dgemv_ */


    }
}

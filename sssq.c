
/* 
                    Verify conversion from c to c#
Make on linux with
gcc -O -I ~/safeqp sssq.c -L ~/safeqp -lsafeqp -o sssq

Make on cygwin in safeqp with
gcc -O -I. ../VECTOR/sssq.c -DMSDOSS -L. -lsafeqp -o sssq 

Make 64 bit version on Windows in safeqp64 with
cl -D__SYSNT__ -I . safeqp.lib ..\VECTOR\sssq.c

Result
scale=  9.0000000000000000, sumsq=  2.0370370370370372
scale= 10.0000000000000000, sumsq=  2.2000000000000002
scale= 10.0000000000000000, sumsq=  3.8499999999999996
scale= 10.0000000000000000, sumsq=  3.8500000000000001
  4.0571891388307382  -3.0885621722338525
  5.4686269665968865  -3.1771243444677046
  6.8800647943630340  -3.2656865167015567
  8.2915026221291814  -3.3542486889354093
piv[0] -13584
a[0][0]   1.0000000000000000
piv[1] 0
a[1][0]   2.0000000000000000  a[1][1]   3.0000000000000000
piv[2] -2145620944
a[2][0]   4.0000000000000000  a[2][1]   5.0000000000000000  a[2][2]   6.0000000000000000
piv[0] 1
a[0][0]  -0.1428571428571426
piv[1] 2
a[1][0]   1.1428571428571432  a[1][1]  -1.1666666666666661
piv[2] 3
a[2][0]   0.6666666666666666  a[2][1]   0.8333333333333333  a[2][2]   6.0000000000000000
  3.0000000000000049  -3.0000000000000067   1.0000000000000022
  1.0000000000000004   2.0000000000000009   2.9999999999999982
piv[0] 1
a[0][0]   1.0000000000000000  a[0][1]   2.0000000000000000  a[0][2]   4.0000000000000000
piv[1] 2
a[1][1]   3.0000000000000000  a[1][2]   5.0000000000000000
piv[2] 3
a[2][2]   6.0000000000000000
piv[0] 3
a[0][0]   6.0000000000000000  a[0][1]   0.8333333333333333  a[0][2]   0.6666666666666666
piv[1] 2
a[1][1]  -1.1666666666666661  a[1][2]   1.1428571428571432
piv[2] 3
a[2][2]  -0.1428571428571426
  3.0000000000000049  -3.0000000000000067   1.0000000000000022
  1.0000000000000004   2.0000000000000009   2.9999999999999991

*/
#include <ldefns.h>

int main()
{
    {
        double a[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int n = 10;
        double scale = 1;
        double sumsq = 0;
        dsssq(n / 2, a, -2, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
        scale = 1;
        sumsq = 0;
        dsssq(n / 2, a + 1, -2, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
        scale = 1;
        sumsq = 0;
        dsssq(n, a, -1, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
        scale = 1;
        sumsq = 0;
        dsssqvec(n, a, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
    }
    {
        int orthog = 1, n = 4, i;
        double x[] = {1, 2, 3, 4};
        double y[] = {5, 6, 7, 8};
        double cs = 0.75, sn = sqrt(1 - cs * cs);

        delm(orthog, n, x, 1, y, 1, cs, sn);
        for (i = 0; i < n; ++i)
            printf("%20.16f %20.16f\n", x[i], y[i]);
    }
    {
        int n = 3, i, j, ij;
        double a[] = {1, 2, 3, 4, 5, 6};
        double acopy[6];
        dcopyvec(6, a, acopy);
        short_scl piv[3];
        char *U = "U";
        int back = 10;
        for (i = 0, ij = 0; i < n; ++i)
        {
            printf("piv[%d] %ld\n", i, piv[i]);
            for (j = 0; j <= i; ++j, ++ij)
            {
                printf("a[%d][%d] %20.16f  ", i, j, a[ij]);
            }
            printf("\n");
        }
        back = dsptrf(U, n, a, piv);
        for (i = 0, ij = 0; i < n; ++i)
        {
            printf("piv[%d] %ld\n", i, piv[i]);
            for (j = 0; j <= i; ++j, ++ij)
            {
                printf("a[%d][%d] %20.16f  ", i, j, a[ij]);
            }
            printf("\n");
        }
        double b[] = {1, 2, 3};
        dsptrs(U, n, 1, a, piv, b, n);
        double lm_eps = fabs((((double)4) / 3 - 1) * 3 - 1);
        double c[3];
        dsmxmulv(n, acopy, b, c);
        printf("%20.16f %20.16f %20.16f\n", b[0], b[1], b[2]);
        printf("%20.16f %20.16f %20.16f\n", c[0], c[1], c[2]);
        double at[] = {1, 2, 4, 3, 5, 6};
        dcopyvec(6, at, a);
        U = "L";
        double bb[] = {1, 2, 3};
        for (i = 0, ij = 0; i < n; ++i)
        {
            printf("piv[%d] %ld\n", i, piv[i]);
            for (j = i; j < n ; ++j, ++ij)
            {
                printf("a[%d][%d] %20.16f  ", i, j, a[ij]);
            }
            printf("\n");
        }
        dcopyvec(n, bb, b);
        back = dsptrf(U, n, a, piv);
        for (i = 0, ij = 0; i < n; ++i)
        {
            printf("piv[%d] %ld\n", i, piv[i]);
            for (j = i; j < n ; ++j, ++ij)
            {
                printf("a[%d][%d] %20.16f  ", i, j, a[ij]);
            }
            printf("\n");
        }
        dsptrs(U, n, 1, a, piv, b, n);
        dsmxmulvT(n, at, b, c);
        printf("%20.16f %20.16f %20.16f\n", b[0], b[1], b[2]);
        printf("%20.16f %20.16f %20.16f\n", c[0], c[1], c[2]);
    }
}
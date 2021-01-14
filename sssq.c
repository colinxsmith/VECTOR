
/* 
                    Verify conversion from c to c#
Make on linux with
gcc -O -I ~/safeqp sssq.c -L ~/safeqp -lsafeqp -o sssq

Make on Windows in safeqp64 with
cl -D__SYSNT__ -I . safeqp.lib ..\VECTOR\sssq.c

Result
./sssq 
scale=   9.000000000000000, sumsq=   2.037037037037037
scale=  10.000000000000000, sumsq=   2.200000000000000
scale=  10.000000000000000, sumsq=   3.850000000000000
scale=  10.000000000000000, sumsq=   3.850000000000000
*/
#include <ldefns.h>

int main()
{
    {
        double a[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int n = 10;
        double scale = 1;
        double sumsq = 0;
        dsssq(n / 2, a, 2, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
        scale = 1;
        sumsq = 0;
        dsssq(n / 2, a + 1, 2, &scale, &sumsq);
        printf("scale=%20.16f, sumsq=%20.16f\n", scale, sumsq);
        scale = 1;
        sumsq = 0;
        dsssq(n, a, 1, &scale, &sumsq);
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
}

/* 
                    Verify conversion from c to c#
Make on linux with
gcc -O -I ~/safeqp sssq.c -L ~/safeqp -lsafeqp -o sssq

Make on Windows in safeqp64 with
cl -D__SYSNT__ -I . safeqp.lib ..\VECTOR\sssq.c

Result
./sssq 
scale=  9.0000000000000000, sumsq=  2.0370370370370372
scale= 10.0000000000000000, sumsq=  2.2000000000000002
scale= 10.0000000000000000, sumsq=  3.8499999999999996
scale= 10.0000000000000000, sumsq=  3.8500000000000001
  4.0571891388307382  -3.0885621722338525
  5.4686269665968865  -3.1771243444677046
  6.8800647943630340  -3.2656865167015567
  8.2915026221291814  -3.3542486889354093
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
}
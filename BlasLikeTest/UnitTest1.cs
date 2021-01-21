using Microsoft.VisualStudio.TestTools.UnitTesting;
using Blas;
using System;
namespace BlasLikeTest
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void Test_daxpyvec()
        {
            int n = 2;
            double a = 0.5;
            double[] x = { 1, 2 };
            double[] y = { 2, 3 };
            BlasLike.daxpyvec(n, a, x, y);
            Assert.IsTrue(y[0] == 2.5 && y[1] == 4);
        }
        [TestMethod]
        public void Test_daddvec()
        {
            int n = 2;
            double[] x = { 1, 2 };
            double[] y = { 2, 3 };
            double[] z = { 2, 3 };
            BlasLike.daddvec(n, x, y, z);
            Assert.IsTrue(z[0] == 3 && z[1] == 5);
        }
        [TestMethod]
        public void Test_dsubvec()
        {
            int n = 2;
            double[] z = { 2, 3 };
            double[] x = { 1, 2 };
            double[] y = { 2, 3 };
            BlasLike.dsubvec(n, x, y, z);
            Assert.IsTrue(z[0] == -1 && z[1] == -1);
        }
        [TestMethod]
        public void Test_ddotvec()
        {
            int n = 2;
            double[] x = { 1, 2 };
            double[] y = { 2, 3 };
            Assert.IsTrue(BlasLike.ddotvec(n, x, y) == 8);
        }
        [TestMethod]
        public void Test_dcopyvec()
        {
            int n = 2;
            double[] x = { 1, 2 };
            double[] y = { 2, 3 };
            BlasLike.dcopyvec(n, x, y);
            Assert.IsTrue(y[0] == 1 && y[1] == 2);
        }
        [TestMethod]
        public void Test_dcopy1()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 2, 3 };
            BlasLike.dcopy(n, x, 3, y, 1);
            Assert.IsTrue(y[0] == 1 && y[1] == 4);
        }
        [TestMethod]
        public void Test_dcopy2()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 2, 3 };
            BlasLike.dcopy(n, x, -3, y, 1);
            Assert.IsTrue(y[0] == 4 && y[1] == 1);
        }
        [TestMethod]
        public void Test_dzerovec()
        {
            int n = 2;
            double[] x = { 1, 2 };
            BlasLike.dzerovec(n, x);
            Assert.IsTrue(x[0] == 0 && x[1] == 0);
        }
        [TestMethod]
        public void Test_dsetvec1()
        {
            int n = 2;
            double[] x = { 1, 2 };
            BlasLike.dsetvec(n, 10, x);
            Assert.IsTrue(x[0] == 10 && x[1] == 10);
        }
        [TestMethod]
        public void Test_dset1()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dset(n, 10, x, 3);
            Assert.IsTrue(x[0] == 10 && x[1] == 2 && x[2] == 3 && x[3] == 10 && x[4] == 5 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dset2()
        {
            int n = 3;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dset(n, 10, x, -2);
            Assert.IsTrue(x[0] == 10 && x[1] == 2 && x[2] == 10 && x[3] == 4 && x[4] == 10 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dsetvec2()
        {
            int n = 6;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dsetvec(n, 10, x);
            Assert.IsTrue(x[0] == 10 && x[1] == 10 && x[2] == 10 && x[3] == 10 && x[4] == 10 && x[5] == 10);
        }
        [TestMethod]
        public void Test_dneg1()
        {
            int n = 3;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dneg(n, x, 2);
            Assert.IsTrue(x[0] == -1 && x[1] == 2 && x[2] == -3 && x[3] == 4 && x[4] == -5 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dneg2()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dneg(n, x, -3);
            Assert.IsTrue(x[0] == -1 && x[1] == 2 && x[2] == 3 && x[3] == -4 && x[4] == 5 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dnegvec()
        {
            int n = 6;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dnegvec(n, x);
            Assert.IsTrue(x[0] == -1 && x[1] == -2 && x[2] == -3 && x[3] == -4 && x[4] == -5 && x[5] == -6);
        }
        [TestMethod]
        public void Test_dscal1()
        {
            int n = 3;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dscal(n, 3, x, 2);
            Assert.IsTrue(x[0] == 3 && x[1] == 2 && x[2] == 9 && x[3] == 4 && x[4] == 15 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dscal2()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dscal(n, 2, x, 3);
            Assert.IsTrue(x[0] == 2 && x[1] == 2 && x[2] == 3 && x[3] == 8 && x[4] == 5 && x[5] == 6);
        }
        [TestMethod]
        public void Test_dscalvec()
        {
            int n = 6;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            BlasLike.dscalvec(n, 3, x);
            Assert.IsTrue(x[0] == 3 && x[1] == 6 && x[2] == 9 && x[3] == 12 && x[4] == 15 && x[5] == 18);
        }
        [TestMethod]
        public void Test_dsccopy1()
        {
            int n = 3;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 1, 1, 1 };
            BlasLike.dsccopy(n, 2, x, 2, y, 1);
            Assert.IsTrue(y[0] == 2 && y[1] == 6 && y[2] == 10);
        }
        [TestMethod]
        public void Test_dsccopy2()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 1, 1 };
            BlasLike.dsccopy(n, -1, x, -3, y, 1);
            Assert.IsTrue(y[0] == -4 && y[1] == -1);
        }
        [TestMethod]
        public void Test_dsccopy3()
        {
            int n = 2;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 1, 1 };
            BlasLike.dsccopy(n, -2, x, -3, y, 1);
            Assert.IsTrue(y[0] == -8 && y[1] == -2);
        }
        [TestMethod]
        public void Test_dsccopyvec()
        {
            int n = 6;
            double[] x = { 1, 2, 3, 4, 5, 6 };
            double[] y = { 1, 1, 1, 1, 1, 1 };
            BlasLike.dsccopyvec(n, -2, x, y);
            Assert.IsTrue(y[0] == -2 && y[1] == -4 && y[2] == -6 && y[3] == -8 && y[4] == -10 && y[5] == -12);
        }
        [TestMethod]
        public void Test_dsssq1()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            BlasLike.dsssq(a.Length / 2, a, 2, scale, sumsq);
            Assert.IsTrue(scale[0] == 9 && sumsq[0] == 2.037037037037037, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }
        [TestMethod]
        public unsafe void Test_dsssq2()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            fixed (double* aa = a)
                BlasLike.dsssq(a.Length / 2, aa + 1, -2, scale, sumsq);
            Assert.IsTrue(scale[0] == 10 && sumsq[0] == 2.2, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }
        [TestMethod]
        public unsafe void Test_dsssq2a()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            fixed (double* aa = &(a[1]))
                BlasLike.dsssq(a.Length / 2, aa, -2, scale, sumsq);
            Assert.IsTrue(scale[0] == 10 && sumsq[0] == 2.2, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }

        [TestMethod]
        public void Test_dsssq3()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            BlasLike.dsssq(a.Length / 2, a, -2, scale, sumsq, 1);
            Assert.IsTrue(scale[0] == 10 && sumsq[0] == 2.2, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }
        [TestMethod]
        public void Test_dsssqvec()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            BlasLike.dsssqvec(a.Length, a, scale, sumsq);
            Assert.IsTrue(scale[0] == 10 && sumsq[0] == 3.85, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }
        [TestMethod]
        public void Test_dsum1()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double sum = BlasLike.dsum(a.Length / 2, a, 2);
            Assert.IsTrue(sum == 25, $"sum is {sum}");
        }
        [TestMethod]
        public void Test_dsum2()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double sum = BlasLike.dsum(a.Length / 2, a, 2, 1);
            Assert.IsTrue(sum == 30, $"sum is {sum}");
        }
        [TestMethod]
        public void Test_dsumvec()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double sum = BlasLike.dsumvec(a.Length, a);
            Assert.IsTrue(sum == 55, $"sum is {sum}");
        }
        [TestMethod]
        public unsafe void Test_dxminmax1()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double min, max;
            fixed (double* aa = a)
                BlasLike.dxminmax(a.Length, aa, 1, &max, &min);
            Assert.IsTrue(max == 10 && min == 1, $"max is {max}, min is {min}");
        }
        [TestMethod]
        public void Test_dxminmax1a()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double[] min = { 1 }, max = { 1 };
            BlasLike.dxminmax(a.Length, a, 1, max, min);
            Assert.IsTrue(max[0] == 10 && min[0] == 1, $"max is {max}, min is {min}");
        }
        [TestMethod]
        public unsafe void Test_detagen()
        {
            double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double alpha = 0;
            long iswap = 0;
            int itrans = 0;
            fixed (double* xx = x)
                BlasLike.detagen(x.Length / 2, &alpha, xx + 1, 2, &iswap, &itrans);
            Assert.IsTrue(alpha == 10 && iswap == 9 && itrans == 1, $"alpha={alpha}, iswap={iswap}, itrans={itrans}, lm_rooteps={BlasLike.lm_rooteps}");
        }
        [TestMethod]
        public unsafe void Test_dswap()
        {
            double[] a = { 1, 2, 3, 4 };
            double[] b = { 11, 12, 13, 14 };
            fixed (double* pa = a)
            fixed (double* pb = b)
                BlasLike.dswap(a.Length / 2, pa, 2, pb, 1);
            Assert.IsTrue(a[0] == 11 && a[1] == 2 && a[2] == 12 && a[3] == 4 && b[0] == 1 && b[1] == 3 && b[2] == 13 && b[3] == 14, $"{a[0]},{a[1]},{a[2]},{a[3]}  {b[0]},{b[1]},{b[2]},{b[3]}");
        }
        [TestMethod]
        public unsafe void Test_dswapvec()
        {
            double[] a = { 1, 2, 3, 4 };
            double[] b = { 11, 12, 13, 14 };
            fixed (double* pa = a)
            fixed (double* pb = b)
                BlasLike.dswapvec(a.Length, pa, pb);
            Assert.IsTrue(a[0] == 11 && a[1] == 12 && a[2] == 13 && a[3] == 14 && b[0] == 1 && b[1] == 2 && b[2] == 3 && b[3] == 4, $"{a[0]},{a[1]},{a[2]},{a[3]}  {b[0]},{b[1]},{b[2]},{b[3]}");
        }
        [TestMethod]
        public unsafe void Test_delm1()
        {
            int orthog = 0;
            double[] a = { 1, 2, 3, 4 };
            double[] b = { 5, 6, 7, 8 };
            double cs = 0.75, sn = Math.Sqrt(1 - cs * cs);
            fixed (double* px = a)
            fixed (double* py = b)
                BlasLike.delm(orthog, a.Length, px, 1, py, 1, cs, sn);
            Assert.IsTrue(a[0] == 1 && a[1] == 2 && a[2] == 3 && a[3] == 4 && b[0] == 5.6614378277661475 && b[1] == 7.322875655532295 && b[2] == 8.984313483298443 && b[3] == 10.64575131106459, $"{a[0]},{a[1]},{a[2]},{a[3]}  {b[0]},{b[1]},{b[2]},{b[3]}");
        }
        [TestMethod]
        public unsafe void Test_delm2()
        {
            int orthog = 1;
            double[] a = { 1, 2, 3, 4 };
            double[] b = { 5, 6, 7, 8 };
            double cs = 0.75, sn = Math.Sqrt(1 - cs * cs);
            fixed (double* px = a)
            fixed (double* py = b)
                BlasLike.delm(orthog, a.Length, px, 1, py, 1, cs, sn);
            Assert.IsTrue(a[0] == 4.0571891388307382 && a[1] == 5.4686269665968865 && a[2] == 6.880064794363034 && a[3] == 8.2915026221291814 && b[0] == -3.0885621722338525 && b[1] == -3.1771243444677046 && b[2] == -3.2656865167015567 && b[3] == -3.3542486889354093, $"{a[0]},{a[1]},{a[2]},{a[3]}  {b[0]},{b[1]},{b[2]},{b[3]}");
        }
        [TestMethod]
        public unsafe void Test_MatrixFactoriseAndSolveU()
        {
            int n = 3;
            double[] a = {  1,
                            2,3,
                            4,5,6};
            int[] piv = { 1, 2, 3 };
            char[] U = { 'U' };
            int back = 10;
            fixed (double* ap = a)
            fixed (int* ipiv = piv)
            fixed (char* UP = U)
                back = Factorise.dsptrf(UP, n, ap, ipiv);
            Assert.IsTrue(back == 0, $"{back} {piv[0]} {piv[1]} {piv[2]}  {a[0]} {a[1]} {a[2]} {a[3]} {a[4]} {a[5]} ");
            double[] b = { 1, 0, 1 };
            double[] bcopy = (double[])b.Clone();
            fixed (double* ap = a)
            fixed (double* bp = b)
            fixed (int* ipiv = piv)
            fixed (char* UP = U)
                Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
            double[] a1 = { 1, 2, 4 };
            double[] a2 = { 2, 3, 5 };
            double[] a3 = { 4, 5, 6 };
            double c1 = BlasLike.ddotvec(n, b, a1);
            double c2 = BlasLike.ddotvec(n, b, a2);
            double c3 = BlasLike.ddotvec(n, b, a3);
            double okerror = BlasLike.lm_eps * 64;
            Assert.IsTrue(back == 0 && Math.Abs(c1 - bcopy[0]) < okerror && Math.Abs(c2 - bcopy[1]) < okerror && Math.Abs(c3 - bcopy[2]) < okerror, $"back={back} {c1},{c2},{c3} ");
        }
        [TestMethod]
        public unsafe void Test_MatrixFactoriseAndSolveL()
        {
            int n = 3;
            double[] a = { 1,2,4,
                             3,5,
                               6};
            int[] piv = { 1, 2, 3 };
            char[] U = { 'L' };
            int back = 10;
            fixed (double* ap = a)
            fixed (int* ipiv = piv)
            fixed (char* UP = U)
                back = Factorise.dsptrf(UP, n, ap, ipiv);
            Assert.IsTrue(back == 0, $"{back} {piv[0]} {piv[1]} {piv[2]}  {a[0]} {a[1]} {a[2]} {a[3]} {a[4]} {a[5]} ");
            double[] b = { 1, 0, 1 };
            double[] bcopy = (double[])b.Clone();
            fixed (double* ap = a)
            fixed (double* bp = b)
            fixed (int* ipiv = piv)
            fixed (char* UP = U)
                Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
            double[] a1 = { 1, 2, 4 };
            double[] a2 = { 2, 3, 5 };
            double[] a3 = { 4, 5, 6 };
            double c1 = BlasLike.ddotvec(n, b, a1);
            double c2 = BlasLike.ddotvec(n, b, a2);
            double c3 = BlasLike.ddotvec(n, b, a3);
            double okerror = BlasLike.lm_eps * 64;
            Assert.IsTrue(back == 0 && Math.Abs(c1 - bcopy[0]) < okerror && Math.Abs(c2 - bcopy[1]) < okerror && Math.Abs(c3 - bcopy[2]) < okerror, $"back={back} {c1},{c2},{c3} ");
        }
    }
}

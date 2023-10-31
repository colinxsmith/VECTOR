using Microsoft.VisualStudio.TestTools.UnitTesting;
using Blas;
using Solver;
using System;
using ActiveSet;
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
            double scale = s1;
            double sumsq = 0;
            BlasLike.dsssq(a.Length / 2, a, 2, ref scale, ref sumsq);
            Assert.IsTrue(scale == 9 && sumsq == 2.037037037037037, $"scale is {scale}, sumsq is {sumsq}");
        }
        [TestMethod]
        public unsafe void Test_dsssq2()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double s1 = 1;
            double[] scale = { s1 };
            double[] sumsq = { 0 };
            fixed (double* aa = a)
            fixed (double* sscale = scale)
            fixed (double* ssumsq = sumsq)
                BlasLike.dsssq(a.Length / 2, aa + 1, -2, sscale, ssumsq);
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
            fixed (double* sscale = scale)
            fixed (double* ssumsq = sumsq)
                BlasLike.dsssq(a.Length / 2, aa, -2, sscale, ssumsq);
            Assert.IsTrue(scale[0] == 10 && sumsq[0] == 2.2, $"scale is {scale[0]}, sumsq is {sumsq[0]}");
        }

        [TestMethod]
        public void Test_dsssq3()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double scale = 1;
            double sumsq = 0;
            BlasLike.dsssq(a.Length / 2, a, -2, ref scale, ref sumsq, 1);
            Assert.IsTrue(scale == 10 && sumsq == 2.2, $"scale is {scale}, sumsq is {sumsq}");
        }
        [TestMethod]
        public void Test_dsssqvec()
        {
            double[] a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

            {
                double scale = 1;
                double sumsq = 0;
                BlasLike.dsssqvec(a.Length, a, ref scale, ref sumsq);
                Assert.IsTrue(scale == 10 && sumsq == 3.85, $"scale is {scale}, sumsq is {sumsq}");
            }
            {
                double scale = 1;
                double sumsq = 0;
                BlasLike.dsssqvec(a.Length, a, ref scale, ref sumsq);
                Assert.IsTrue(scale == 10 && sumsq == 3.85, $"scale is {scale}, sumsq is {sumsq}");
            }
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
            double min = 1, max = 1;
            BlasLike.dxminmax(a.Length, a, 1, ref max, ref min);
            Assert.IsTrue(max == 10 && min == 1, $"max is {max}, min is {min}");
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
            double[] y = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            alpha = 0;
            int iswap1 = 0;
            itrans = 0;
            BlasLike.detagen(x.Length / 2, ref alpha, y, 2, ref iswap1, ref itrans, 1);
            Assert.IsTrue(alpha == 10 && iswap == 9 && itrans == 1, $"alpha={alpha}, iswap1={iswap1}, itrans={itrans}, lm_rooteps={BlasLike.lm_rooteps}");
            Assert.IsTrue(x[0] == y[0] && x[1] == y[1] && x[2] == y[2] && x[3] == y[3] && x[4] == y[4] && x[5] == y[5] && x[6] == y[6] && x[7] == y[7] && x[8] == y[8] && x[9] == y[9],
            $"{x[0]} {x[1]} {x[2]} {x[3]} {x[4]} {x[5]} {x[6]} {x[7]} {x[8]}\n{x[0]} {x[1]} {x[2]} {x[3]} {x[4]} {x[5]} {x[6]} {x[7]} {x[8]} ");
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
            int n = 4;
            double[] a = {  1,
                            2,3,
                            4,5,6,
                            7,8,9,10};
            double[] acopy = (double[])a.Clone();
            int[] piv = { 1, 2, 3, 4 };
            char U = 'U';
            int back = 10;
            back = Factorise.Factor(U, n, a, piv);
            Assert.IsTrue(back == 0);
            double[] b = { 1, 2, 3, 4 };
            double[] bcopy = (double[])b.Clone();
            Factorise.Solve(U, n, 1, a, piv, b, n);
            double okerror = BlasLike.lm_eps * 16;
            double[] c = new double[n];
            Factorise.dsmxmulv(n, acopy, b, c);
            double[] errorvec = new double[n];
            BlasLike.dsubvec(n, bcopy, c, errorvec);
            double error = Math.Sqrt(BlasLike.ddotvec(n, errorvec, errorvec)) / n;
            Assert.IsTrue(back == 0 && error < okerror, $"back={back} error={error} {c[0]},{c[1]},{c[2]},{c[3]} ");
        }
        [TestMethod]
        public unsafe void Test_MatrixFactoriseAndSolveL()
        {
            int n = 4;
            double[] a ={1,2,4,7,
                            3,5,8,
                            6,9,
                            10};
            double[] acopy = (double[])a.Clone();
            int[] piv = { 1, 2, 3, 4 };
            char U = 'L';
            int back = 10;
            back = Factorise.Factor(U, n, a, piv);
            Assert.IsTrue(back == 0);
            double[] b = { 1, 2, 3, 4 };
            double[] bcopy = (double[])b.Clone();
            Factorise.Solve(U, n, 1, a, piv, b, n);
            double okerror = BlasLike.lm_eps * 16;
            double[] c = new double[n];
            Factorise.dsmxmulvT(n, acopy, b, c);
            double[] errorvec = new double[n];
            BlasLike.dsubvec(n, bcopy, c, errorvec);
            double error = Math.Sqrt(BlasLike.ddotvec(n, errorvec, errorvec)) / n;
            Assert.IsTrue(back == 0 && error < okerror, $"back={back} error={error} {c[0]},{c[1]},{c[2]},{c[3]} ");
        }
        [TestMethod]
        public void Test_MatrixFactorisationsAndUse()
        {
            /* 
            Generate 2 random symmetric matrices such that the lower packed version of
            one is equal to the upper packed version of the other.
            Test that upper and lower of the solver are working, but this shows that the
            working is not identical!
            */
            var n = 2000;
            var cov = new double[n * (n + 1) / 2];
            for (var i = 0; i < n; ++i)
            {
                for (var j = i; j < n; j++)
                {
                    var cc = new Random();
                    cov[j * (j + 1) / 2 + i] = cc.NextDouble();
                }
            }
            var S = new double[n * (n + 1) / 2];
            var ST = new double[n * (n + 1) / 2];
            for (var i = 0; i < n; ++i)
            {
                for (var j = i; j < n; j++)
                {
                    S[j * (j + 1) / 2 + i] = ST[i * n - i * (i - 1) / 2 + j - i] = cov[j * (j + 1) / 2 + i];
                }
            }
            var unit1 = new double[n];
            var unit1T = new double[n];
            for (var i = 0; i < n; ++i) unit1T[i] = unit1[i] = 1;
            char U = 'U';
            char L = 'L';
            var ipiv = new int[n];
            var Sbefore = (double[])S.Clone();
            var back = Factorise.Factor(U, n, S, ipiv);
            Factorise.Solve(U, n, 1, S, ipiv, unit1, n);
            var ipivT = new int[n];
            var STbefore = (double[])ST.Clone();
            var backT = Factorise.Factor(L, n, ST, ipivT);
            Factorise.Solve(L, n, 1, ST, ipivT, unit1T, n);
            var c = new double[n];
            Factorise.dsmxmulv(n, Sbefore, unit1, c);
            var cT = new double[n];
            Factorise.dsmxmulvT(n, STbefore, unit1T, cT);
            var diff = new double[n];
            int negpiv = 0, negpivT = 0;
            for (var i = 0; i < n; ++i)
            {
                if (ipiv[i] < 0) negpiv++;
                if (ipivT[i] < 0) negpivT++;
            }
            BlasLike.dsubvec(n, unit1T, unit1, diff);
            var error = Math.Sqrt(BlasLike.ddotvec(n, diff, diff) / n);
            Assert.IsTrue(error < BlasLike.lm_rooteps, $"{error} back={back} backT={backT} negpiv={negpiv} negpivT={negpivT}\n {unit1[0]},{unit1[1]},{unit1[2]},{unit1[3]} \n {unit1T[0]},{unit1T[1]},{unit1T[2]},{unit1T[3]} \n {c[0]},{c[1]},{c[2]},{c[3]} \n {cT[0]},{cT[1]},{cT[2]},{cT[3]}");
        }
        [TestMethod]
        public void Test_Eigen2()
        {
            double[] M = { 1, 2, 3 };
            var lambda = new double[2];
            var t = new double[4];
            Factorise.Eigen2(M, lambda, t);
            var lt = new double[4];

            lt[0] = (M[0] * t[0] + M[1] * t[1]) / lambda[0];
            lt[1] = (M[1] * t[0] + M[2] * t[1]) / lambda[0];
            lt[2] = (M[0] * t[2] + M[1] * t[3]) / lambda[1];
            lt[3] = (M[1] * t[2] + M[2] * t[3]) / lambda[1];
            Assert.IsTrue(Math.Abs(lt[0] - t[0]) < 4 * BlasLike.lm_eps, $"{lt[0]} {t[0]}\n{M[0]} {M[1]} {lambda[0]} {t[0]} {t[2]}\n{M[1]} {M[2]} {lambda[1]} {t[1]} {t[3]}");
        }
        [TestMethod]
        public void Test_transpose()
        {
            double[] am = { 11, 12, 13, 21, 22, 23 };
            var bm = new double[6];
            Factorise.dmx_transpose(3, 2, am, bm);
            Factorise.dmx_transpose(3, 2, am, am);
            Assert.IsTrue(am[0] == bm[0] && am[1] == bm[1] && am[2] == bm[2] && am[3] == bm[3] && am[4] == bm[4] && am[5] == bm[5]);
            Factorise.dmx_transpose(2, 3, am, bm);
            Factorise.dmx_transpose(2, 3, am, am);
            Assert.IsTrue(am[0] == bm[0] && am[1] == bm[1] && am[2] == bm[2] && am[3] == bm[3] && am[4] == bm[4] && am[5] == bm[5]);
        }
        [TestMethod]
        public void Test_MatrixTimesVector()
        {
            double[] am = { 11, 12,
                                21,22,
                                31,32 };
            var xx = new double[3];
            var yy = new double[2];
            xx[1] = 1;
            Factorise.dmxmulv(2, 3, am, xx, yy);
            Assert.IsTrue(yy[0] == am[2] && yy[1] == am[3]);
            Factorise.dmx_transpose(2, 3, am, am);
            Factorise.dmxmulv(2, 3, am, xx, yy, 0, 0, 0, true);
            Assert.IsTrue(yy[0] == am[1] && yy[1] == am[4]);
            var xxx = new double[2];
            var yyy = new double[3];
            xxx[0] = 1;
            Factorise.dmxmulv(3, 2, am, xxx, yyy);
            Assert.IsTrue(yyy[0] == am[0] && yyy[1] == am[1] && yyy[2] == am[2]);
        }
        [TestMethod]
        public void Test_SemiDefinateInverse()
        {
            var n = 40;
            var cov = new double[n * (n + 1) / 2];
            for (int i = 0; i < cov.Length; ++i) cov[i] = i + 1;
            cov[4] = -cov[4];
            var S = new double[n * (n + 1) / 2];
            var ST = new double[n * (n + 1) / 2];
            for (var i = 0; i < n; ++i)
            {
                for (var j = i; j < n; j++)
                {
                    ST[i * n - i * (i - 1) / 2 + j - i] = S[j * (j + 1) / 2 + i] = cov[j * (j + 1) / 2 + i];
                }
            }

            char way = 'L';
            var piv = new int[n];
            var back = way == 'L' ? Factorise.Factor(way, n, ST, piv) : Factorise.Factor(way, n, S, piv);
            var Sback = new double[n * n];
            for (int i = 0; i < n; ++i) Sback[i * n + i] = 1;
            var whichroot = 2;
            var info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
            whichroot = 0;
            info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
            for (int i = 0; i < n; ++i) Sback[n * i + i] += -1;
            var error = Math.Sqrt(BlasLike.ddotvec(n * n, Sback, Sback)) / n;
            Assert.IsTrue(error < BlasLike.lm_eps * 128, $"error is {error}");
            way = 'U';
            back = way == 'L' ? Factorise.Factor(way, n, ST, piv) : Factorise.Factor(way, n, S, piv);
            BlasLike.dzerovec(Sback.Length, Sback);
            for (int i = 0; i < n; ++i) Sback[i * n + i] = 1;
            whichroot = 2;
            info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
            whichroot = 0;
            info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
            for (int i = 0; i < n; ++i) Sback[n * i + i] += -1;
            error = Math.Sqrt(BlasLike.ddotvec(n * n, Sback, Sback)) / n;
            Assert.IsTrue(error < BlasLike.lm_eps * 128, $"error is {error}");
        }
        [TestMethod]
        public void Test_TriangleInverse()
        {
            var n = 3;
            double[] U = {1,
                             1,1,
                             2,3,1};
            double[] L = {1,1,2,
                             1,3,
                             1};
            var way = 'L';
            for (int itest = 0; itest < 2; ++itest)
            {
                way = itest == 0 ? 'L' : 'U';
                var unit = new double[n * n * 2];
                var piv = new int[n];
                for (int i = 0; i < n; i++)
                {
                    piv[i] = i + 1;
                    unit[i + n * i] = 1;
                    unit[i + n * i + n * n] = 1;
                }
                var symdef = Factorise.Solve(way, n, n, way == 'L' ? L : U, piv, unit, n, 0, 0, 0, 1);
                symdef = Factorise.Solve(way, n, n, way == 'L' ? L : U, piv, unit, n, 0, 0, n * n, -1);
                var start = 0;
                start = n * n;
                Factorise.dmx_transpose(n, n, unit, unit, n * n, n * n);
                if (way == 'L') Assert.IsTrue(unit[start + 6] == (-unit[6] + unit[3] * unit[7]), $"{unit[start + 6]} {-unit[6] + unit[3] * unit[7]}");
                else if (way == 'U') Assert.IsTrue(unit[start + 2] == (-unit[2] + unit[1] * unit[5]), $"{unit[start + 2]} {-unit[2] + unit[1] * unit[5]}");
            }
        }
        [TestMethod]
        public void Test_FixSingular()
        {
            var n = 3;
            char way = 'U';
            double[] M ={1,
                           1,1,
                           1,1,1};//singular
            var S = new double[n * (n + 1) / 2];
            BlasLike.dcopyvec(S.Length, M, S);
            var piv = new int[n];
            var info = Factorise.Factor(way, n, S, piv);
            var resolve = new double[n * n];
            BlasLike.dzerovec(n * n, resolve);
            for (int ii = 0; ii < n; ++ii) resolve[ii + n * ii] = 1;
            Factorise.Solve(way, n, n, S, piv, resolve, n, 0, 0, 0, 0, true);
            Assert.IsTrue(resolve[n * n - 1] == 1, $"\n{resolve[0]} {resolve[1]} {resolve[2]} \n{resolve[3]} {resolve[4]} {resolve[5]} \n{resolve[6]} {resolve[7]} {resolve[8]} ");
            info = Factorise.Factor(way, n, S, piv);

            way = 'L';
            BlasLike.dcopyvec(S.Length, M, S);
            info = Factorise.Factor(way, n, S, piv);
            BlasLike.dzerovec(n * n, resolve);
            for (int ii = 0; ii < n; ++ii) resolve[ii + n * ii] = 1;
            Factorise.Solve(way, n, n, S, piv, resolve, n, 0, 0, 0, 0, true);
            Assert.IsTrue(resolve[0] == 1, $"\n{resolve[0]} {resolve[1]} {resolve[2]} \n{resolve[3]} {resolve[4]} {resolve[5]} \n{resolve[6]} {resolve[7]} {resolve[8]} ");
        }
        [TestMethod]
        public void Test_ordering()
        {
            double[] S = { 1.2, 4, -2.1, 5, 45 };
            var order = new int[S.Length];
            var inverse = new int[S.Length];
            var dropbad = new byte[S.Length];
            Ordering.Order.getorder(S.Length, S, order);
            Ordering.Order.Reorder_gen(S.Length, order, S);
            Assert.IsTrue(S[0] == 45 && S[4] == -2.1, $"\n{order[0]} {order[1]} {order[2]} {order[3]} {order[4]} \n{S[0]} {S[1]} {S[2]} {S[3]} {S[4]} ");

            for (int i = 0; i < S.Length; ++i) inverse[order[i]] = i;
            Ordering.Order.Reorder_gen(S.Length, inverse, S);
            Assert.IsTrue(S[0] == 1.2 && S[4] == 45, $"\n{S[0]} {S[1]} {S[2]} {S[3]} {S[4]} ");

            dropbad[2] = 1;
            Ordering.Order.getorder(S.Length, S, order, dropbad);
            Ordering.Order.Reorder_gen(S.Length, order, S);
            Assert.IsTrue(S[0] == -2.1 && S[4] == 1.2, $"\n{order[0]} {order[1]} {order[2]} {order[3]} {order[4]} \n{S[0]} {S[1]} {S[2]} {S[3]} {S[4]} ");
            for (int i = 0; i < S.Length; ++i) inverse[order[i]] = i;
            Ordering.Order.Reorder_gen(S.Length, inverse, S);

            Ordering.Order.getorderabs(S.Length, S, order);
            Ordering.Order.Reorder_gen(S.Length, order, S);
            Assert.IsTrue(S[0] == 45 && S[4] == 1.2, $"\n{S[0]} {S[1]} {S[2]} {S[3]} {S[4]} ");
            string[] symmU ={"11",
                              "21","22",
                              "31","32","33",
                              "41","42","43","44",
                              "51","52","53","54","55"};
            Ordering.Order.ReorderSymm(5, order, symmU, 'U');
            Assert.IsTrue(symmU[0] == "55" && symmU[14] == "11", $"\n{order[0]} {order[1]} {order[2]} {order[3]} {order[4]}\n{symmU[0]},{symmU[14]}");
            string[] symmL ={"11","12","13","14","15",
                                  "22","23","24","25",
                                       "33","34","35",
                                            "44","45",
                                                 "55"};
            Assert.IsTrue(symmL[1] == "12" && symmL[8] == "25", $"{symmL[1]},{symmL[8]}");
            Ordering.Order.ReorderSymm(5, order, symmL, 'L');
            Assert.IsTrue(symmL[1] == "45" && symmL[8] == "14", $"\n{order[0]} {order[1]} {order[2]} {order[3]} {order[4]}\n{symmL[1]},{symmL[8]}");
        }
        [TestMethod]
        public unsafe void Test_LPand_QP()
        {
            var n = 10;
            var m = 2;
            var x = new double[n];
            double[] c = { 1, 2, 3, 4, 5, 6, 17, 8, 9, 10 };
            double[] A = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,
                               0, 0, 1, 1, 1, 0, 0, 0, 0, 0};
            double[] L = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.1 };
            double[] U = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5 };
            Factorise.dmx_transpose(n, m, A, A);
            double[] hess = new double[n * (n + 1) / 2];
            var tdata = 2 * n;
            var timeD = new double[n, tdata];

            for (int i = 0; i < n; ++i)
            {
                for (int time = 0; time < tdata; ++time)
                {
                    var dat = new Random();
                    timeD[i, time] = dat.NextDouble();
                }
            }
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    hess[i * (i + 1) / 2 + j] = 0;
                    var ti = 0.0;
                    var tj = 0.0;
                    for (int time = 0; time < tdata; ++time)
                    {
                        ti += timeD[i, time];
                        tj += timeD[j, time];
                        hess[i * (i + 1) / 2 + j] += timeD[i, time] * timeD[j, time];
                    }
                    hess[i * (i + 1) / 2 + j] = hess[i * (i + 1) / 2 + j] / tdata - ti / tdata * tj / tdata;
                }
            }
            BlasLike.dscalvec(hess.Length, 1e3, hess);
            var obj = -0.9;
            var iter = 0;
            short back;
            var implied = new double[n];
            for (int i = 0; i < 2; ++i)
            {
                if (i == 0)
                {
                    BlasLike.dsetvec(x.Length, 1.0 / n, x);
                    var budget = 1.0;
                    var constraintVal = new double[m];
                    back = ActiveSet.Optimise.LPopt(n, m, x, L, U, A, c, ref obj, ref iter);
                    Factorise.dmxmulv(m, n, A, x, constraintVal);
                    Assert.IsTrue(back == 0 && Math.Abs(constraintVal[0] - budget) < BlasLike.lm_eps * 16, $"LP back is {back} {BlasLike.ddotvec(n, x, c)} {obj} {iter} {constraintVal[1]}");
                }
                else
                {
                    BlasLike.dsetvec(x.Length, 1.0 / n, x);
                    var budget = 1.0;
                    var opt = new ActiveSet.Optimise();
                    back = opt.QPopt(n, m, x, L, U, A, c, hess, ref obj, ref iter);
                    Factorise.dsmxmulv(n, hess, x, implied);
                    var constraintVal = new double[m];
                    Factorise.dmxmulv(m, n, A, x, constraintVal);
                    Assert.IsTrue(back == 0 && Math.Abs(constraintVal[0] - budget) < BlasLike.lm_eps * 16, $"QP back is {back} {BlasLike.ddotvec(n, x, c) + 0.5 * BlasLike.ddotvec(n, implied, x)} {obj} {iter} {constraintVal[1]}");
                }
            }
        }
        [TestMethod]
        public void Test_SOCP()
        {
            //Simplest SOCP optimisation; minimise c.x so that sum(x*x)=1
            int n = 12;
            var m = 1;
            var x = new double[n];
            double[] b = { 1 };
            double[] c = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0 };
            double[] A = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
            BlasLike.dscalvec(c.Length, 1e-3, c);
            int[] cone = { n };
            int nvar = 0;
            foreach (int ic in cone) nvar += ic;
            int[] typecone = { (int)InteriorPoint.conetype.SOCP };

            Factorise.dmx_transpose(n, m, A, A);// Does nothing here
            var opt1 = new InteriorPoint.Optimise(n, m, x, A, b, c);
            var back = opt1.Opt("SOCP", cone, typecone, true);
            var util = BlasLike.ddotvec(n, c, x);
            Assert.IsTrue(back != 6 && Math.Abs(x[n - 1] - 1) < 1e-8, $"util = {util}; cone check{x[n - 1] - 1}; back={back}");
        }
        [TestMethod]
        public void Test_TwoDimensionalMatrix()
        {
            //This shows how 2 dimension matrices in BITA20 are converted
            //to 1 dimensional arrays for use in the optimiser
            var n = 3;
            var nf = 2;
            double[] FL ={ 1, 3, 4,
                          3, 4, 7 };
            double[] f1 = { 1, 3, 4 };
            double[] f2 = { 3, 4, 7 };
            double[] s1 = { 1, 3 };
            double[] s2 = { 3, 4 };
            double[] s3 = { 4, 7 };
            var FLas2D = Portfolio.Portfolio.oneD2twoD(n, nf, FL);
            Assert.IsTrue(s1[0] == FLas2D[0][0] && s1[1] == FLas2D[0][1]);
            Assert.IsTrue(s2[0] == FLas2D[1][0] && s2[1] == FLas2D[1][1]);
            Assert.IsTrue(s3[0] == FLas2D[2][0] && s3[1] == FLas2D[2][1]);
            Assert.IsTrue(f1[0] == FLas2D[0][0] && f1[1] == FLas2D[1][0] && f1[2] == FLas2D[2][0]);
            Assert.IsTrue(f2[0] == FLas2D[0][1] && f2[1] == FLas2D[1][1] && f2[2] == FLas2D[2][1]);
            var FLback = Portfolio.Portfolio.twoD2oneD(n, nf, FLas2D);
            BlasLike.dsubvec(n * nf, FL, FLback, FLback);
            var test = BlasLike.ddotvec(n * nf, FLback, FLback);
            Assert.IsTrue(test <= BlasLike.lm_eps, $"test = {test}");
            n = 5;
            var m = 2;
            double[] a1 = { 1, 1, 1, 1, 1 };
            double[] a2 = { 1, 2, 3, 4, 5 };
            double[][] Aas2D = { a1, a2 };
            var A = Portfolio.Portfolio.twoD2oneD(m, n, Aas2D);
            var Aas2Dback = Portfolio.Portfolio.oneD2twoD(m, n, A);
            BlasLike.dsubvec(n, Aas2D[0], Aas2Dback[0], Aas2Dback[0]);
            BlasLike.dsubvec(n, Aas2D[1], Aas2Dback[1], Aas2Dback[1]);
            test = BlasLike.ddotvec(n, Aas2Dback[0], Aas2Dback[0]) + BlasLike.ddotvec(n, Aas2Dback[1], Aas2Dback[1]);
            Assert.IsTrue(test <= 2 * BlasLike.lm_eps, $"test = {test}");
            //A bit of overkill to show that the transpose option works;
            //if we transpose in twoD2oneD then we must transpose back somewhere else
            A = Portfolio.Portfolio.twoD2oneD(m, n, Aas2D, transpose: true);
            Factorise.dmx_transpose(n, m, A, A);
            Aas2Dback = Portfolio.Portfolio.oneD2twoD(m, n, A);
            BlasLike.dsubvec(n, Aas2D[0], Aas2Dback[0], Aas2Dback[0]);
            BlasLike.dsubvec(n, Aas2D[1], Aas2Dback[1], Aas2Dback[1]);
            test = BlasLike.ddotvec(n, Aas2Dback[0], Aas2Dback[0]) + BlasLike.ddotvec(n, Aas2Dback[1], Aas2Dback[1]);
            Assert.IsTrue(test <= 2 * BlasLike.lm_eps, $"test = {test}");
        }
        [TestMethod]
        public void Test_FMP()
        {
            //FMP calculates data for compressed risk model
            var nfac = 2;
            var n = 3;
            double[] SV = { 1, 2, 3 };
            var Q = new double[(nfac + 1) * n];
            double[] FC = { 2,
                                1, 3 };
            double[] FL = { 1, 0, 1,
                            0, 1, 1 }; //Factors by assets
            var uplow = 'L';
            var back = Factorise.FMP(n, nfac, FC, SV, FL, Q, uplow);
            double[] correct = { 2, 1, 3, 1, 3, 4, 3, 4, 7 };
            var result = new double[n * n];
            for (var i = 0; i < n; ++i) //Multiply out the factor part of the compressed model
            {
                Factorise.dmxmulv(n, nfac, Q, Q, result, n, n + nfac * i, i * n, true);
            }
            var test = 0.0;
            for (var i = 0; i < result.Length; ++i)
            {
                test = Math.Max(Math.Abs(result[i] - correct[i]), test);
            }
            Assert.IsTrue(back == 0 && test <= BlasLike.lm_eps * 8, $"{test}:{result[0]},{result[1]},{result[2]},{result[3]},{result[4]},{result[5]},{result[6]},{result[7]},{result[8]},");
            var w = new double[n];
            w[0] = 1;
            var Qw = new double[n];
            Factorise.FacMul(n, nfac, Q, w, Qw);
            result[0] += Q[0];
            BlasLike.dsubvec(n, Qw, result, Qw);
            test = BlasLike.ddotvec(n, Qw, Qw);
            Assert.IsTrue(Math.Abs(test) < BlasLike.lm_eps2, $"test={test}");
            var Qinv = (double[])Q.Clone();
            //Invert the factorised risk model
            back = Factorise.FMPinverse(n, nfac, Qinv, uplow);
            var Qwback = new double[n];
            w[0] = 1;
            w[1] = 2;
            w[2] = 3;
            Factorise.FacMul(n, nfac, Q, w, Qw);
            Factorise.FacMul(n, nfac, Qinv, Qw, Qwback, inverse: true);
            BlasLike.dsubvec(n, w, Qwback, Qwback);
            test = BlasLike.ddotvec(n, Qwback, Qwback);
            Assert.IsTrue(back == 0 && test <= BlasLike.lm_eps, $"Inverse test check {test}");
        }
        [TestMethod]
        public void Test_digit()
        {
            var d1 = 2.0 + 1e-13;
            Assert.IsTrue(Portfolio.Portfolio.check_digit(d1) == 2.0, $"{Portfolio.Portfolio.check_digit(d1)}");
            d1 = 2.0 - 1e-13;
            Assert.IsTrue(Portfolio.Portfolio.check_digit(d1) == 2.0, $"{Portfolio.Portfolio.check_digit(d1)}");
        }
        [TestMethod]
        public void Test_round()
        {
            var d1 = 0.010099999999999998;
            Assert.IsTrue(Portfolio.Portfolio.digitRound(d1, 1e7) == 0.0101, $"{Portfolio.Portfolio.digitRound(d1, 1e7)}");
            d1 = 0.008900000000000002;
            Assert.IsTrue(Portfolio.Portfolio.digitRound(d1) == 0.0089, $"{Portfolio.Portfolio.digitRound(d1)}");
            d1 = 0.0001234556789;
            Assert.IsTrue(Portfolio.Portfolio.digitRound(d1) == 0.0001234556789, $"{Portfolio.Portfolio.digitRound(d1)}");
        }
        [TestMethod]
        public void Test_sigfig()
        {
            var d1 = 23123.45567;
            Assert.IsTrue(Portfolio.Portfolio.rounder(d1, 4) == 23120, $"{Portfolio.Portfolio.rounder(d1, 4)}");
            d1 = 0.001234567;
            Assert.IsTrue(Portfolio.Portfolio.rounder(d1, 2) == 0.0012, $"{Portfolio.Portfolio.rounder(d1, 2)}");
            d1 = 1.283456789;
            Assert.IsTrue(Portfolio.Portfolio.rounder(d1) == 1.28, $"{Portfolio.Portfolio.rounder(d1)}");
        }
        [TestMethod]
        public void Test_MCTR()
        {
            var opt = new Portfolio.Portfolio("");
            opt.n = 5;
            opt.ntrue = opt.n;
            double[] Q = { 1,
                            1e-2,2,
                            1e-2,1e-2,3,
                            1e-2,1e-2,1e-2,4,
                            1e-2,1e-2,1e-2,1e-2,5 };
            double[] w = { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double[] bench = { 1.0, 1.0, 1.0, 1.0, 1.0 };
            opt.Q = Q;
            BlasLike.dscalvec(opt.n, 1.0 / 15.0, w);//make portfolio and benchmark weights sum to 1
            BlasLike.dscalvec(opt.n, 1.0 / opt.n, bench);
            var breakdown = (double[])w.Clone();
            opt.hessmull(w.Length, opt.Q, w, breakdown);
            var risk = Math.Sqrt(BlasLike.ddotvec(w.Length, w, breakdown));//portfolio risk
            opt.RiskBreakdown(w, null, breakdown);//breakdown is MCTR
            var risktest = BlasLike.ddotvec(w.Length, w, breakdown);//check we get risk
            Assert.IsTrue(Math.Abs(risk - risktest) < BlasLike.lm_eps8, $"{risk - risktest}");

            var breakdownr = (double[])w.Clone();
            var beta = (double[])w.Clone();
            opt.hessmull(opt.n, opt.Q, bench, beta);
            var benchvar = BlasLike.ddotvec(opt.n, bench, beta);//benchmark variance
            var benchrisk = Math.Sqrt(benchvar);//benchmark risk

            BlasLike.dscalvec(opt.n, 1.0 / benchvar, beta);//asset beta
            var portbeta = BlasLike.ddotvec(opt.n, w, beta);//portfolio beta

            opt.RiskBreakdown(w, bench, breakdownr, beta);
            BlasLike.daxpyvec(w.Length, -portbeta, bench, w);//residual weights
            var resrisktest = BlasLike.ddotvec(w.Length, w, breakdownr);
            BlasLike.daxpyvec(w.Length, portbeta, bench, w);//restore weights
            // total risk squared = residual risk squared +systematic risk squared (i.e. beta*benchrisk squared)
            Assert.IsTrue(Math.Abs(risk * risk - portbeta * portbeta * benchrisk * benchrisk - resrisktest * resrisktest) < BlasLike.lm_eps8, $"{risk * risk - portbeta * portbeta * benchrisk * benchrisk - resrisktest * resrisktest}");
        }
        [TestMethod]
        public void Test_LOSS()
        {
            double[] s = { 1, 2, 5, 4, 3, 7, 6, 5, 8, 9, -1, 7, -3, 8, 9, 6 };
            double[] target = new double[s.Length];
            bool[] breakIndex = new bool[s.Length];
            BlasLike.dsetvec(s.Length, 4.3, target);
            var loss = Portfolio.Portfolio.LOSS(s, target, breakIndex);
            var numberOfLosses = 0;
            foreach (var k in breakIndex) if (k) numberOfLosses++;
            Assert.IsTrue(numberOfLosses == 6, $"number of losses {numberOfLosses}");
            Assert.IsTrue(loss == 19.8, $"Loss is {loss}");
        }
        [TestMethod]
        public void Test_ETL()
        {
            double[] s = { 1, 2, 15, 4, 3, 7, 16, 5, 18, 9, -1, 7, -3,
             8, 9, 6, 5, 3, 8, 9, 1, 4, 8, 9, 4, 0, 12 };
            double tail = 0.05;
            bool[] breakIndex = new bool[s.Length];
            var VAR = 0.0;
            var VARindex = -1;
            var etl = Portfolio.Portfolio.ETL(s, tail, ref VAR, ref VARindex, breakIndex);
            var numberInTail = 0;
            foreach (var k in breakIndex) if (k) numberInTail++;
            Assert.IsTrue(Math.Abs(etl - 17.481481481481506) < BlasLike.lm_eps, $"{etl - 17.481481481481506}");
            Assert.IsTrue(Math.Abs(VAR - 16.000000000000092) < BlasLike.lm_eps, $"{VAR - 16.000000000000092}");
            Assert.IsTrue(VARindex == 6, $"{VARindex}");
            Assert.IsTrue(numberInTail == 1, $"{numberInTail}");
        }
        [TestMethod]
        public void Test_Composite()
        {
            var n = 5;
            var ncomp = 2;
            var ntrue = n - ncomp;
            var nfac = 2;
            double[] compositeWeights = { 0.333, 0.333, 0.333,
                                            0.4, 0.6, 0 };
            double[] SV = { 1, 2, 3 };
            double[] FC = { 0.5,
                            0.1, 0.6 };
            double[] FL = { 2, 3, 4,
                            4, -1, 2 };
            var port = new Portfolio.FPortfolio("");
            //        port.n = n;
            port.ntrue = ntrue;
            port.nfac = nfac;
            port.SV = SV;
            port.FL = FL;
            port.FC = FC;
            port.nfixed = 0;
            port.ncomp = ncomp;
            port.compw = compositeWeights;
            port.makeQ();
            var testq = (double[])port.compQ.Clone();
            var piv = new int[ncomp];
            var back = Factorise.Factor('U', ncomp, testq, piv);
            Assert.IsFalse(back != 0, $"compQ should be positive definite! {back} is not zero");
            double[] w = { 0.4, 0.6, 0, 0, -1 };//A fully hedged portfolio
            var implied = (double[])w.Clone();
            port.hessmull(n, port.Q, w, implied);
            var variance = BlasLike.ddotvec(n, w, implied) * 0.5;
            Assert.IsTrue(Math.Abs(variance) <= BlasLike.lm_eps, $"variance should be zero, not {variance}");
            variance = port.Variance(w);
            //Now test for full covariance case. Generate covariances from risk model
            var C = new double[ntrue * (ntrue + 1) / 2];
            Factorise.Fac2Cov(ntrue, nfac, port.Q, C);
            var portC = new Portfolio.Portfolio("");
            //        portC.n = n;
            portC.ntrue = ntrue;
            portC.nfixed = 0;
            portC.Q = C;
            portC.ncomp = ncomp;
            portC.compw = compositeWeights;
            portC.makeQ();
            portC.hessmull(n, port.Q, w, implied);
            variance = BlasLike.ddotvec(n, w, implied) * 0.5;
            Assert.IsTrue(Math.Abs(variance) <= BlasLike.lm_eps, $"variance should be zero, not {variance}");
            int[] order = { 2, 0, 3, 1, 4 };//correct dropping order for basket or trade number
            var inverse = (int[])order.Clone();
            for (var i = 0; i < n; ++i) inverse[order[i]] = i;

            portC.mainorder = (int[])order.Clone();
            portC.mainorderInverse = (int[])inverse.Clone();
            portC.createMainOrderTrue(n);
            Ordering.Order.Reorder(n, portC.mainorder, w);
            Ordering.Order.Reorder(n, portC.mainorderInverse, w);
            portC.hessmull(n, portC.Q, w, implied);
            Ordering.Order.Reorder(n, portC.mainorder, w);
            Ordering.Order.Reorder(n, portC.mainordertrue, w);
            Ordering.Order.Reorder(n, portC.mainordertrueInverse, w);
            var fixedW = new double[n];
            var fixedSecondOrder = new double[n];
            var nfixed = 2;
            BlasLike.dcopyvec(n, w, fixedW);
            BlasLike.dzerovec(n - nfixed, fixedW);
            Ordering.Order.Reorder(n, portC.mainordertrue, w);
            Ordering.Order.Reorder(n, portC.mainordertrue, fixedW);
            Ordering.Order.Reorder(n, portC.mainordertrue, order);//combined order
            Ordering.Order.ReorderSymm(ntrue, order, portC.Q);
            var orderPortC = new int[ncomp];
            var orderPortCInverse = new int[ncomp];
            for (var i = 0; i < ncomp; ++i) orderPortC[i] = order[i + ntrue] - ntrue;
            for (var i = 0; i < ncomp; ++i) orderPortCInverse[orderPortC[i]] = i;
            Ordering.Order.ReorderSymm(ncomp, orderPortC, portC.compQ);
            //Need to reorder over assets and composites in compw and compImplied
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, order, portC.compw, i * ntrue);
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, order, portC.compImplied, i * ntrue);
            Ordering.Order.Reorder_gen(ncomp, orderPortC, portC.compw, ntrue, columns: true);
            Ordering.Order.Reorder_gen(ncomp, orderPortC, portC.compImplied, ntrue, columns: true);
            portC.nfixed = 0;
            portC.hessmull(n, portC.Q, w, implied);
            portC.hessmull(n, portC.Q, fixedW, fixedSecondOrder);
            variance = BlasLike.ddotvec(n, w, implied) * 0.5;
            Assert.IsTrue(variance < BlasLike.lm_eps, $"reordering without fixing gives non zero variance; ({variance})");
            //When using this in optimisation w and implied have order; mainorder i.e. the fixed weights are at the bottom
            //fixedW, fixedSecondOrder, compw and compImplied etc have order mainordertrue*mainorder since the entries must be compatible with the risk model
            portC.nfixed = nfixed;
            portC.fixedW = fixedW;
            portC.nfixedComp = 1;
            portC.nfixedTrue = 1;
            var fixedVariance = BlasLike.ddotvec(n, fixedSecondOrder, fixedW) * 0.5;
            portC.fixedSecondOrder = fixedSecondOrder;
            //BlasLike.dsubvec(n,w,fixedW,w);
            BlasLike.dzerovec(n, implied);
            Ordering.Order.Reorder(n, portC.mainordertrueInverse, w);
            Ordering.Order.Reorder(n, portC.mainordertrueInverse, implied);
            portC.hessmull(n - nfixed, portC.Q, w, implied);
            //Change to order with fixed assets at the end
            Ordering.Order.Reorder(n, portC.mainordertrueInverse, fixedSecondOrder);
            Ordering.Order.Reorder(n, portC.mainordertrueInverse, fixedW);
            variance = BlasLike.ddotvec(n - nfixed, w, implied) * 0.5 + BlasLike.ddotvec(n - nfixed, fixedSecondOrder, w) + BlasLike.ddotvec(n, implied, fixedW) + fixedVariance;
            Assert.IsTrue(Math.Abs(variance) <= BlasLike.lm_eps, $"reordering with fixed and non-fixed combined gives non zero variance; ({variance})");


            Ordering.Order.Reorder(n, portC.mainordertrueInverse, order);//Reset order
            Ordering.Order.Reorder(n, portC.mainorderInverse, w);

            //Drop order
            Ordering.Order.Reorder(n, portC.mainorder, w);
            port.mainorder = (int[])order.Clone();
            var orderInverse = (int[])order.Clone();
            for (var i = 0; i < n; i++) orderInverse[order[i]] = i;
            port.mainorderInverse = (int[])orderInverse.Clone();
            port.createMainOrderTrue(n);
            BlasLike.dcopyvec(n, w, fixedW);
            BlasLike.dzerovec(n - nfixed, fixedW);

            //Apply to previous factor risk model in correct order
            Ordering.Order.Reorder(n, portC.mainordertrue, w);
            Ordering.Order.Reorder(n, portC.mainordertrue, fixedW);
            Ordering.Order.Reorder(n, portC.mainordertrue, order);//combined order

            Ordering.Order.Reorder(ntrue, order, port.Q);
            Ordering.Order.Reorder_gen(ntrue, order, port.Q, nfac, 1, astart: ntrue, columns: true);
            orderPortC = new int[ncomp];
            orderPortCInverse = new int[ncomp];
            for (var i = 0; i < ncomp; ++i) orderPortC[i] = order[i + ntrue] - ntrue;
            for (var i = 0; i < ncomp; ++i) orderPortCInverse[orderPortC[i]] = i;
            Ordering.Order.ReorderSymm(ncomp, orderPortC, port.compQ);
            //Need to reorder over assets and composites in compw and compImplied
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, order, port.compw, i * ntrue);
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, order, port.compImplied, i * ntrue);
            Ordering.Order.Reorder_gen(ncomp, orderPortC, port.compw, ntrue, columns: true);
            Ordering.Order.Reorder_gen(ncomp, orderPortC, port.compImplied, ntrue, columns: true);
            port.nfixed = 0;
            port.hessmull(n, port.Q, fixedW, fixedSecondOrder);
            port.fixedW = fixedW;
            port.fixedSecondOrder = fixedSecondOrder;


            port.nfixed = nfixed;
            port.nfixedComp = 1;
            port.nfixedTrue = 1;


            //Change to order with fixed assets at the end
            BlasLike.dzerovec(n, implied);
            Ordering.Order.Reorder(n, port.mainordertrueInverse, w);
            Ordering.Order.Reorder(n, port.mainordertrueInverse, implied);
            port.hessmull(n - nfixed, port.Q, w, implied);
            Ordering.Order.Reorder(n, port.mainordertrueInverse, fixedSecondOrder);
            Ordering.Order.Reorder(n, port.mainordertrueInverse, fixedW);
            variance = BlasLike.ddotvec(n - nfixed, w, implied) * 0.5 + BlasLike.ddotvec(n - nfixed, fixedSecondOrder, w) + BlasLike.ddotvec(n, implied, fixedW) + fixedVariance;
            Assert.IsTrue(Math.Abs(variance) <= BlasLike.lm_eps, $"reordering with fixed and non-fixed combined gives non zero variance; ({variance})");
            Ordering.Order.Reorder(n, port.mainorderInverse, w);
            for (var i = 0; i < n; ++i) inverse[order[i]] = i;
            Ordering.Order.Reorder(ntrue, inverse, port.Q);
            Ordering.Order.Reorder_gen(ntrue, inverse, port.Q, nfac, 1, astart: ntrue, columns: true);


            Ordering.Order.ReorderSymm(ncomp, orderPortC, port.compQ);
            //Need to reorder over assets and composites in compw and compImplied
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, inverse, port.compw, i * ntrue);
            for (var i = 0; i < ncomp; ++i) Ordering.Order.Reorder(ntrue, inverse, port.compImplied, i * ntrue);
            Ordering.Order.Reorder_gen(ncomp, orderPortCInverse, port.compw, ntrue, columns: true);
            Ordering.Order.Reorder_gen(ncomp, orderPortCInverse, port.compImplied, ntrue, columns: true);
            port.nfixed = 0;
            variance = port.Variance(w);
            Assert.IsTrue(Math.Abs(variance) <= BlasLike.lm_eps, $"reordering back to start gives non zero variance; ({variance})");
        }
        [TestMethod]
        public void Test_readLicence()
        {
            ColourConsole.WriteLine($"1 int has length {(double)sizeof(UInt32) / (double)sizeof(byte)} bytes");
            var licence = new Licensing.Licence();
            Assert.IsTrue(licence.fromRegistry(false) > 0);
            var testlicence = licence.licenceByteValue;
            if (testlicence == null) return;
            var vv = new Licensing.validator_t();
            for (var i = 0; i < vv.b.Length; ++i)
            {
                vv.b[i] = testlicence[i];
            }
            Assert.IsTrue(licence.check_valid(ref vv));
            var aStart = licence.UnixTimeStampToDateTime((double)vv.start);
            ColourConsole.WriteLine($"Start time {aStart}");
            var aStop = licence.UnixTimeStampToDateTime((double)vv.stop);
            ColourConsole.WriteLine($"Start time {aStop}");
            var uStart = new DateTimeOffset(aStart).ToUnixTimeSeconds();
            var uStop = new DateTimeOffset(aStop).ToUnixTimeSeconds();
            ColourConsole.WriteLine($"{vv.start} {uStart}");
            ColourConsole.WriteLine($"{vv.stop} {uStop}");
            Assert.IsTrue(uStart == vv.start);
            Assert.IsTrue(uStop == vv.stop);
        }
        [TestMethod]
        public void Test_convert_licence()
        {
            var licence = new Licensing.Licence();
            Assert.IsTrue(licence.fromRegistry(false) > 0);
            var testlicence = licence.licenceByteValue;
            if (testlicence == null) return;
            Licensing.byteint curveKeys = new Licensing.byteint();
            curveKeys.byte1 = testlicence[16];
            curveKeys.byte2 = testlicence[17];
            curveKeys.byte3 = testlicence[18];
            curveKeys.byte4 = testlicence[19];
            ColourConsole.WriteLine($"Curve Keys integer {curveKeys.mainint}");
            int hid = 0, start = 0, stop = 0;
            licence.convert(testlicence, ref hid, ref start, ref stop);
            var newbytes = new byte[16];
            int newhid = hid;
            int newstart = start;
            int newstop = stop;
            licence.convert(newbytes, ref newhid, ref newstart, ref newstop);
            var check = true;
            for (var i = 0; i < 16; i++) check &= newbytes[i] == testlicence[i];
            Assert.IsTrue(check);
        }
        [TestMethod]
        public void Test_adapters()
        {
            var licence = new Licensing.Licence();
            var testhid = licence.VolId();
            ColourConsole.WriteInfo($"Volid {testhid}  {testhid:x}");
            Assert.IsTrue(licence.fromRegistry(false) > 0);
            var testlicence = licence.licenceByteValue;
            if (testlicence == null) return;
            int hid = 0, start = 0, stop = 0;
            licence.convert(testlicence, ref hid, ref start, ref stop);
            Licensing.byteint curveKeys = new Licensing.byteint();
            curveKeys.byte1 = testlicence[16];
            curveKeys.byte2 = testlicence[17];
            curveKeys.byte3 = testlicence[18];
            curveKeys.byte4 = testlicence[19];
            ColourConsole.WriteLine($"Curve Keys integer {curveKeys.mainint}");
            hid -= curveKeys.mainint;
            DateTimeOffset now = new DateTimeOffset(DateTime.Now);
            var timenow = now.ToUnixTimeSeconds();
            var pass = true;
            pass = pass && (timenow < stop);
            pass = pass && (timenow > start);
            pass = pass && (hid == testhid || hid == 0x13101955);
            Assert.IsTrue(pass);
            start = (int)timenow - 23;
            hid = (int)testhid;
            curveKeys.mainint = (int)Convert.ToInt32("111101110", 2);
            testlicence[16] = curveKeys.byte1;
            testlicence[17] = curveKeys.byte2;
            testlicence[18] = curveKeys.byte3;
            testlicence[19] = curveKeys.byte4;
            hid += curveKeys.mainint;
            licence.convert(testlicence, ref hid, ref start, ref stop);
            Assert.IsTrue(licence.toRegistry(false));
        }
        [TestMethod]
        public void Test_CheckLicence()
        {
            var licence = new Licensing.Licence();
            var vers = licence.CheckLicence(true);
            Assert.IsTrue(vers);
        }
        [TestMethod]
        public void Test_Version()
        {
            var lic = new Licensing.Licence();
            var version = lic.version;
            var revision = lic.revision;
            var vers = lic.CheckLicence(true);
            Assert.IsTrue(vers, $"Licence is invalid; version is {version}: revision {revision}");
        }
        [TestMethod]
        public void Test_DeleteKey()
        {
            var licence = new Licensing.Licence();
            var vers = licence.deleteKey();
            //    Assert.IsTrue(vers);
            //Now we put a licence back
            DateTimeOffset now = new DateTimeOffset(DateTime.Now);
            var start = (int)now.ToUnixTimeSeconds();
            var bmonth = 10;
            var extra = now.Month >= bmonth ? 1 : 0;
            var stopl = new DateTimeOffset(new DateTime(now.Year + extra, bmonth, 26, 12, 0, 0));
            int stop = (int)stopl.ToUnixTimeSeconds();
            int hid = 0x13101955;
            int keys = Convert.ToInt32("101010011", 2);//0x153;// binary 101010011
            hid += keys;
            Licensing.byteint curveKeys = new Licensing.byteint();
            curveKeys.mainint = keys;
            licence.licenceByteValue = new byte[20];
            licence.licenceByteValue[16] = curveKeys.byte1;
            licence.licenceByteValue[17] = curveKeys.byte2;
            licence.licenceByteValue[18] = curveKeys.byte3;
            licence.licenceByteValue[19] = curveKeys.byte4;
            licence.convert(licence.licenceByteValue, ref hid, ref start, ref stop);
            Assert.IsTrue(licence.toRegistry(false));
        }
    }
}

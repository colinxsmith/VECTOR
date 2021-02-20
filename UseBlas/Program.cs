using System;
using System.Runtime.InteropServices;
using Blas;
using Solver;
using Microsoft.Win32;

namespace UseBlas
{
    class Program
    {
        static unsafe void Main(string[] args)
        {
            var a = 4.0;
            double[] x = { 1, 2, 3 };
            double[] y = { 4, 5, 6 };
            BlasLike.daxpyvec(x.Length, a, x, y);
            var iy = 0;
            foreach (var yy in y)
            {
                Console.WriteLine($"y[{iy}]={y.GetValue(iy)} {y[iy++]}");
            }
            {
                double[] aref = new double[1];
                fixed (double* aq = aref)
                    BlasLike.baseref = 0;
                int n = 4;
                double[] aa ={1,2,4,7,
                            3,5,8,
                            6,9,
                            10};
                double[] acopy = (double[])aa.Clone();
                int[] piv = { 1, 2, 3, 4 };
                char[] U = { 'L' };
                int back = 10;
                back = Factorise.Factor(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    //            Factorise.Solve(UP, n, 1, ap, ipiv, bp, n);
                    Factorise.Solve(U, n, 1, aa, piv, b, n);
                double[] c = new double[n];
                fixed (double* acp = acopy)
                fixed (double* bp = b)
                fixed (double* cp = c)
                    //                   Factorise.dsmxmulvT(n, acp, bp, cp);
                    Factorise.dsmxmulvT(n, acopy, b, c);
                Console.WriteLine($"back={back} safe dsmxmulvT {c[0]},{c[1]},{c[2]},{c[3]} ");
                Factorise.dsmxmulvT(n, acopy, b, c);
                double[] errorvec = new double[n];
                BlasLike.dsubvec(n, bcopy, c, errorvec);
                double error = Math.Sqrt(BlasLike.ddotvec(n, errorvec, errorvec)) / n;
                Console.WriteLine($"b {b[0]},{b[1]},{b[2]},{b[3]} ");
                Console.WriteLine($"back={back} error={error.ToString("e1")} dsmxmulvT {c[0]},{c[1]},{c[2]},{c[3]} ");
            }
            {
                double[] aref = new double[1];
                fixed (double* aq = aref)
                    BlasLike.baseref = 0;
                int n = 4;
                double[] aa = { 1,
                                2,3,
                                4,5,6,
                                7,8,9,10};
                double[] acopy = (double[])aa.Clone();
                int[] piv = { 1, 2, 3, 4 };
                char[] U = { 'U' };
                int back = 10;
                back = Factorise.Factor(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    //               Factorise.Solve(UP, n, 1, ap, ipiv, bp, n);
                    Factorise.Solve(U, n, 1, aa, piv, b, n);
                double[] c = new double[n];
                fixed (double* acp = acopy)
                fixed (double* bp = b)
                fixed (double* cp = c)
                    //                   Factorise.dsmxmulv(n, acp, bp, cp);
                    Factorise.dsmxmulv(n, acopy, b, c);
                Console.WriteLine($"back={back} safe dsmxmulv {c[0]},{c[1]},{c[2]},{c[3]} ");
                Factorise.dsmxmulv(n, acopy, b, c);
                double[] errorvec = new double[n];
                BlasLike.dsubvec(n, bcopy, c, errorvec);
                double error = Math.Sqrt(BlasLike.ddotvec(n, errorvec, errorvec)) / n;
                Console.WriteLine($"b {b[0]},{b[1]},{b[2]},{b[3]} ");
                Console.WriteLine($"back={back} error={error.ToString("e1")} dsmxmulv {c[0]},{c[1]},{c[2]},{c[3]} ");
            }
            {
                double[] aref = new double[1];
                fixed (double* aq = aref)
                    BlasLike.baseref = 0;
                int n = 4;
                double[] S = new double[n * (n + 1) / 2];
                int[] ji = new int[n * (n + 1) / 2];
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = i; j < n; j++, ij++)
                    {
                        //S[ij] = i * n - i * (i - 1) / 2 + j - i;
                        //S[ij] = j * (j + 1) / 2 + i;
                        //S[ij]=ij;
                        ji[j * (j + 1) / 2 + i] = ij;
                        //S[j * (j + 1) / 2 + i] = i * n - i * (i - 1) / 2 + j - i;
                    }
                }
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = i; j < n; j++, ij++)
                    {
                        S[j * (j + 1) / 2 + i] = ji[ij];
                    }
                }

                double[] SS = (double[])S.Clone();
                double[] c = new double[n];
                for (int i = 1; i <= c.Length; ++i)
                {
                    c[i - 1] = i;
                }
                double[] cc = (double[])c.Clone();
                double[] b = new double[n];
                int[] piv = new int[n];
                char[] U = { 'U' };
                int back;
                fixed (double* SSS = S)
                fixed (char* UU = U)
                fixed (int* pv = piv)
                fixed (double* ccc = c)
                {
                    // back = Factorise.Factor(UU, n, SSS, pv);
                    // Factorise.Solve(UU, n, 1, SSS, pv, ccc, n);
                }
                back = Factorise.Factor(U, n, S, piv);
                Factorise.Solve(U, n, 1, S, piv, c, n);
                Factorise.dsmxmulv(n, SS, c, b);
                BlasLike.dsubvec(n, b, cc, b);
                double error = Math.Sqrt(BlasLike.ddotvec(n, b, b) / n);
                Console.WriteLine($"{U[0]}\t\t{back}\tError: {error}");
            }
            {
                double[] aref = new double[1];
                fixed (double* aq = aref)
                    BlasLike.baseref = 0;
                int n = 4;
                double[] S = new double[n * (n + 1) / 2];
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; j++, ij++)
                    {
                        //S[ij] = i * (i + 1) / 2 + j ;
                        //S[ij] = i * n - i * (i s- 1) / 2 + j - i;
                        //S[ij]=ij;
                        S[j * n - j * (j - 1) / 2 + i - j] = i * (i + 1) / 2 + j;
                        //S[i * (i + 1) / 2 + j] = j * n - j * (j - 1) / 2 + i - j;
                    }
                }
                double[] SS = (double[])S.Clone();
                double[] c = new double[n];
                for (int i = 1; i <= c.Length; ++i)
                {
                    c[i - 1] = i;
                }
                double[] cc = (double[])c.Clone();
                double[] b = new double[n];
                int[] piv = new int[n];
                char[] U = { 'L' };
                int back;
                fixed (double* SSS = S)
                fixed (char* UU = U)
                fixed (int* pv = piv)
                    //              back = Factorise.Factor(UU, n, SSS, pv);
                    back = Factorise.Factor(U, n, S, piv);
                fixed (double* SSS = S)
                fixed (int* pv = piv)
                fixed (char* UU = U)
                fixed (double* ccc = c)
                    //Factorise.Solve(UU, n, 1, SSS, pv, ccc, n);
                    Factorise.Solve(U, n, 1, S, piv, c, n);
                Factorise.dsmxmulvT(n, SS, c, b);
                BlasLike.dsubvec(n, b, cc, b);
                double error = Math.Sqrt(BlasLike.ddotvec(n, b, b) / n);
                Console.WriteLine($"{U[0]}\t\t{back}\tError: {error}");
            }

            {
                /* 
                Generate 2 random symmetric matrices such that the lower packed version of
                one is equal to the upper packed version of the other.
                Test that upper and lower of the solver are workijng, but this shows that the
                workijng is not identical!
                */
                var n = 100;
                var cov = new double[n * (n + 1) / 2];
                fixed (double* cv = cov)
                    BlasLike.baseref = 0 * (int)cv;
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
                        ST[i * n - i * (i - 1) / 2 + j - i] = S[j * (j + 1) / 2 + i] = cov[j * (j + 1) / 2 + i];
                    }
                }

                var unit1 = new double[n];
                var unit1T = new double[n];
                //           for (var i = 0; i < n; ++i) unit1T[i] =unit1[i] = 1;
                unit1T[0] = unit1[0] = 1;
                var diff = new double[n];
                Factorise.dsmxmulv(n, S, unit1, diff);
                Console.WriteLine($"{diff[0]},{diff[1]},{diff[2]},{diff[3]}");
                Factorise.dsmxmulvT(n, ST, unit1T, diff);
                Console.WriteLine($"{diff[0]},{diff[1]},{diff[2]},{diff[3]}");
                char[] U = { 'U' };
                char[] L = { 'L' };
                var ipiv = new int[n];
                var Sbefore = (double[])S.Clone();
                var back = Factorise.Factor(U, n, S, ipiv);
                Factorise.Solve(U, n, 1, S, ipiv, unit1, n);
                int[] ipivT = new int[n];
                var STbefore = (double[])ST.Clone();
                var c = new double[n];
                Factorise.dsmxmulv(n, Sbefore, unit1, c);
                var backT = Factorise.Factor(L, n, ST, ipivT);
                Factorise.Solve(L, n, 1, ST, ipivT, unit1T, n);

                var cT = new double[n];
                Factorise.dsmxmulvT(n, STbefore, unit1T, cT);
                int negpiv = 0, negpivT = 0;
                for (var i = 0; i < n; ++i)
                {
                    if (ipiv[i] < 0) negpiv++;
                    if (ipivT[i] < 0) negpivT++;
                }
                BlasLike.dsubvec(n, unit1T, unit1, diff);
                var error = Math.Sqrt(BlasLike.ddotvec(n, diff, diff) / n);
                Console.WriteLine($"{error} back={back} backT={backT} negpiv={negpiv} negpivT={negpivT}\n {unit1[0]},{unit1[1]},{unit1[2]},{unit1[3]} \n {unit1T[0]},{unit1T[1]},{unit1T[2]},{unit1T[3]} \n {c[0]},{c[1]},{c[2]},{c[3]} \n {cT[0]},{cT[1]},{cT[2]},{cT[3]}");
            }
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
                char[] way = { 'U' };
                var piv = new int[n];
                var back = way[0] == 'L' ? Factorise.Factor(way, n, ST, piv) : Factorise.Factor(way, n, S, piv);
                var Sback = new double[n * n];
                for (int i = 0; i < n; ++i) Sback[i * n + i] = 1;
                var whichroot = 2;
                var info = way[0] == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                whichroot = 0;
                info = way[0] == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                for (int i = 0; i < n; ++i) Sback[n * i + i] += -1;
                var error = Math.Sqrt(BlasLike.ddotvec(n * n, Sback, Sback)) / n;
                Console.WriteLine($"{error}");
                way[0] = 'L';
                back = way[0] == 'L' ? Factorise.Factor(way, n, ST, piv) : Factorise.Factor(way, n, S, piv);
                BlasLike.dzerovec(Sback.Length, Sback);
                for (int i = 0; i < n; ++i) Sback[i * n + i] = 1;
                whichroot = 2;
                info = way[0] == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                whichroot = 0;
                info = way[0] == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                for (int i = 0; i < n; ++i) Sback[n * i + i] += -1;
                error = Math.Sqrt(BlasLike.ddotvec(n * n, Sback, Sback)) / n;
                Console.WriteLine($"{error}");
            }
            {
                var n = 500;
                var tdata = 40;
                char[] way = { 'U' };
                var cov = new double[n * (n + 1) / 2];
                var M = new double[n * (n + 1) / 2];
                var MT = new double[n * (n + 1) / 2];
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
                        cov[i * (i + 1) / 2 + j] = 0;
                        var ti = 0.0;
                        var tj = 0.0;
                        for (int time = 0; time < tdata; ++time)
                        {
                            ti += timeD[i, time];
                            tj += timeD[j, time];
                            cov[i * (i + 1) / 2 + j] += timeD[i, time] * timeD[j, time];
                        }
                        cov[i * (i + 1) / 2 + j] = cov[i * (i + 1) / 2 + j] / tdata - ti / tdata * tj / tdata;
                    }
                }
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        M[i * (i + 1) / 2 + j] = cov[i * (i + 1) / 2 + j];
                        MT[j * n - j * (j - 1) / 2 + i - j] = cov[i * (i + 1) / 2 + j];
                    }
                }
                var piv = new int[n];
                var start = new double[n * n];
                var xxx = new double[n];
                for (int i = 0; i < n; ++i)
                {
                    xxx[i] = 1;
                    if (way[0] == 'U') Factorise.dsmxmulv(n, M, xxx, start, i * n);
                    else Factorise.dsmxmulvT(n, MT, xxx, start, i * n);
                    xxx[i] = 0;
                }
                var back = (way[0] == 'U') ? Factorise.Factor(way, n, M, piv) : Factorise.Factor(way, n, MT, piv);
                var r = new double[n * n];
                for (int i = 0; i < n; ++i) r[i * n + i] = 1;
                Console.WriteLine($"{r[0]} {r[1]} {r[2]}");
                Console.WriteLine($"{r[3]} {r[4]} {r[5]}");
                Console.WriteLine($"{r[6]} {r[7]} {r[8]}");
                var rI = new double[n * n];
                for (int i = 0; i < n; ++i) rI[i * n + i] = 1;
                var rBack = new double[n * n];
                for (int i = 0; i < n; ++i) rBack[i * n + i] = 1;
                var whichroot = 2;
                var symback = (way[0] == 'U') ? Factorise.Solve(way, n, n, M, piv, rBack, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, rBack, n, 0, 0, 0, whichroot);
                whichroot = 0;
                symback = (way[0] == 'U') ? Factorise.Solve(way, n, n, M, piv, rI, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, rI, n, 0, 0, 0, whichroot);
                whichroot = 1;
                symback = (way[0] == 'U') ? Factorise.Solve(way, n, n, M, piv, r, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, r, n, 0, 0, 0, whichroot);
                if (whichroot != 0 && way[0] == 'L') Factorise.dmx_transpose(n, n, r, r);
                Console.WriteLine($"{r[0]} {r[1]} {r[2]}");
                Console.WriteLine($"{r[3]} {r[4]} {r[5]}");
                Console.WriteLine($"{r[6]} {r[7]} {r[8]}");

                var rr2 = new double[n * (n + 1) / 2];
                double[] lower = new double[n * n];
                var negpiv = 0;
                for (int i = 0; i < n; ++i)
                {
                    if (piv[i] < 0) negpiv++;
                }
                Console.WriteLine($"{negpiv} negative pivots");
                Factorise.dmx_transpose(n, n, r, lower);
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; j++, ij++)
                    {
                        rr2[ij] = way[0] == 'L' ? BlasLike.ddotvec(n, lower, lower, i * n, j * n) : BlasLike.ddotvec(n, r, r, i * n, j * n);
                    }
                }
                var diff = new double[n * (n + 1) / 2];
                BlasLike.dsubvec(n * (n + 1) / 2, cov, rr2, diff);
                var error = BlasLike.ddotvec(n * (n + 1) / 2, diff, diff) / (n * (n + 1) / 2);
                Console.WriteLine($"error in new covariance {error}");
            }
            {
                var n = 3;
                var tdata = 30;
                char[] way = { 'U' };
                var cov = new double[n * (n + 1) / 2];
                var M = new double[n * (n + 1) / 2];
                var MT = new double[n * (n + 1) / 2];
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
                        cov[i * (i + 1) / 2 + j] = 0;
                        var ti = 0.0;
                        var tj = 0.0;
                        for (int time = 0; time < tdata; ++time)
                        {
                            ti += timeD[i, time];
                            tj += timeD[j, time];
                            cov[i * (i + 1) / 2 + j] += timeD[i, time] * timeD[j, time];
                        }
                        cov[i * (i + 1) / 2 + j] = cov[i * (i + 1) / 2 + j] / tdata - ti / tdata * tj / tdata;
                    }
                }
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        M[i * (i + 1) / 2 + j] = cov[i * (i + 1) / 2 + j];
                        MT[j * n - j * (j - 1) / 2 + i - j] = cov[i * (i + 1) / 2 + j];
                    }
                }
                var piv = new int[n];
                var back = Factorise.Factor(way, n, M, piv);
                var FM = new double[n * n];
                var rootFM = new double[n * n];
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        rootFM[i * n + j] = FM[i * n + j] = (i != j) ? M[i * (i + 1) / 2 + j] : 1;
                    }
                    for (int j = i + 1; j < n; ++j)
                    {
                        rootFM[i * n + j] = FM[i * n + j] = 0;
                    }
                }
                if (n == 3)
                {
                    Console.WriteLine($"{FM[0]} {FM[1]} {FM[2]} ");
                    Console.WriteLine($"{FM[3]} {FM[4]} {FM[5]} ");
                    Console.WriteLine($"{FM[6]} {FM[7]} {FM[8]} ");
                }
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        FM[i * n + j] *= M[i * (i + 3) / 2];
                        rootFM[i * n + j] *= Math.Sqrt(M[i * (i + 3) / 2]);
                    }
                }
                if (n == 3)
                {
                    Console.WriteLine($"\n{FM[0]} {FM[1]} {FM[2]} ");
                    Console.WriteLine($"{FM[3]} {FM[4]} {FM[5]} ");
                    Console.WriteLine($"{FM[6]} {FM[7]} {FM[8]} ");
                }
                for (int i = 0; i < n; ++i)
                {
                    for (int k = i + 1; k < n; ++k)
                    {
                        BlasLike.daxpy(n, M[k * (k + 1) / 2 + i], FM, n, FM, n, k, i);
                    }
                }

                if (n == 3)
                {
                    Console.WriteLine($"\n{FM[0]} {FM[1]} {FM[2]} ");
                    Console.WriteLine($"{FM[3]} {FM[4]} {FM[5]} ");
                    Console.WriteLine($"{FM[6]} {FM[7]} {FM[8]} ");

                    Console.WriteLine($"\n{rootFM[0]} {rootFM[1]} {rootFM[2]} ");
                    Console.WriteLine($"{rootFM[3]} {rootFM[4]} {rootFM[5]} ");
                    Console.WriteLine($"{rootFM[6]} {rootFM[7]} {rootFM[8]} ");
                }

                var checkFM = new double[n * (n + 1) / 2];
                Factorise.dmx_transpose(n, n, rootFM, rootFM);
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        checkFM[i * (i + 1) / 2 + j] = BlasLike.ddotvec(n, rootFM, rootFM, i * n, j * n);
                    }
                }
            }
            {
                double[] am = { 11, 12,
                                21,22,
                                31,32 };
                var xx = new double[3];
                var yy = new double[2];
                for (int i = 0; i < 3; ++i)
                {
                    xx[i] = 1;
                    Factorise.dmxmulv(2, 3, am, xx, yy);
                    Console.WriteLine($"{xx[0]},{xx[1]},{xx[2]}   {yy[0]},{yy[1]}");
                    xx[i] = 0;
                }
                Factorise.dmx_transpose(2, 3, am, am);
                var xxx = new double[2];
                var yyy = new double[3];
                for (int i = 0; i < 2; ++i)
                {
                    xxx[i] = 1;
                    Factorise.dmxmulv(3, 2, am, xxx, yyy);
                    Console.WriteLine($"{xxx[0]},{xxx[1]}   {yyy[0]},{yyy[1]},{yyy[2]}");
                    xxx[i] = 0;
                }
            }
            {
                var n = 3;
                var tdata = 4000;
                var cov = new double[n * (n + 1) / 2];
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
                        cov[i * (i + 1) / 2 + j] = 0;
                        var ti = 0.0;
                        var tj = 0.0;
                        for (int time = 0; time < tdata; ++time)
                        {
                            ti += timeD[i, time];
                            tj += timeD[j, time];
                            cov[i * (i + 1) / 2 + j] += timeD[i, time] * timeD[j, time];
                        }
                        cov[i * (i + 1) / 2 + j] = cov[i * (i + 1) / 2 + j] / tdata - ti / tdata * tj / tdata;
                    }
                }

                var piv = new int[n];
                var S = new double[n * (n + 1)];
                BlasLike.dcopyvec(cov.Length, cov, S);
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        S[j * n - j * (j - 1) / 2 + i - j + cov.Length] = cov[i * (i + 1) / 2 + j];
                    }
                }
                char[] way = { 'U' };
                var R = new double[n * n * 2];
                var R_Inverse = new double[n * n + n * (n + 1) / 2];
                var whichroot = 0;
                var back = Factorise.Factor(way, n, S, piv, way[0] == 'U' ? 0 : cov.Length);
                var npiv = 0;
                for (int i = 0; i < n; ++i) if (piv[i] < 0) npiv++;
                Console.WriteLine($"{npiv} negative pivots");
                BlasLike.dzerovec(R.Length, R);
                for (int i = 0; i < n; i++) R[i * n + i] = 1;
                for (int i = 0; i < n; i++) R_Inverse[i * n + i] = 1;
                for (int i = 0; i < n; i++) R[i * n + i + n * n] = 1;
                whichroot = 1;
                var symprob = Factorise.Solve(way, n, n, S, piv, R, n, way[0] == 'U' ? 0 : cov.Length, 0, 0, whichroot);
                whichroot = -1;
                symprob = Factorise.Solve(way, n, n, S, piv, R, n, way[0] == 'U' ? 0 : cov.Length, 0, n * n, whichroot);
                if (symprob == -10) Console.WriteLine("not positive definite!");
                whichroot = 0;
                symprob = Factorise.Solve(way, n, n, S, piv, R_Inverse, n, way[0] == 'U' ? 0 : cov.Length, 0, 0, whichroot);
                var unit = new double[n * n];
                for (int i = 0; i < n; ++i)
                {//Check inv(R)*inv(RT) is inverse of cov
                    for (int j = 0; j <= i; ++j)
                    {
                        R_Inverse[n * n + i * (i + 1) / 2 + j] = BlasLike.ddotvec(n, R, R, i * n + n * n, j * n + n * n);
                    }
                }
                Factorise.dmx_transpose(n, n, R, R, n * n, n * n);
                for (int i = 0; i < n; ++i)//Check R*inv(RT)=I
                    Factorise.dmxmulv(n, n, R, R, unit, 0, n * n + i * n, i * n);
                var testunit = BlasLike.ddotvec(unit.Length, unit, unit);
                Console.WriteLine($"test unit {testunit}");
            }
            {
                var n = 3;
                double[] U ={1,
                             1,1,
                             2,3,1};
                double[] L ={1,1,2,
                             1,3,
                             1};
                char[] way = { 'U' };
                var unit = new double[n * n * 2];
                var piv = new int[n];
                for (int i = 0; i < n; i++)
                {
                    piv[i] = i + 1;
                    unit[i + n * i] = 1;
                    unit[i + n * i + n * n] = 1;
                }
                var symdef = Factorise.Solve(way, n, n, way[0] == 'L' ? L : U, piv, unit, n, 0, 0, 0, 1);
                symdef = Factorise.Solve(way, n, n, way[0] == 'L' ? L : U, piv, unit, n, 0, 0, n * n, -1);
                var start = 0;
                Console.WriteLine($"{unit[start + 0]} {unit[start + 1]} {unit[start + 2]} ");
                Console.WriteLine($"{unit[start + 3]} {unit[start + 4]} {unit[start + 5]} ");
                Console.WriteLine($"{unit[start + 6]} {unit[start + 7]} {unit[start + 8]} ");
                start = n * n;
                Factorise.dmx_transpose(n, n, unit, unit, n * n, n * n);
                Console.WriteLine($"\n{unit[start + 0]} {unit[start + 1]} {unit[start + 2]} ");
                Console.WriteLine($"{unit[start + 3]} {unit[start + 4]} {unit[start + 5]} ");
                Console.WriteLine($"{unit[start + 6]} {unit[start + 7]} {unit[start + 8]} ");
                if (way[0] == 'L') Console.WriteLine($"{unit[start + 6]} {-unit[6] + unit[3] * unit[7]}");
                else if (way[0] == 'U') Console.WriteLine($"{unit[start + 2]} {-unit[2] + unit[1] * unit[5]}");
            }
            {
                //double[] S = { 3, 20, 2, 0.1, -0.1, 1 }; //not positive definite
                double[] S = { 10, 0.1, 2, 0.1, 0.1, 3 }; //positive definite
                var n = 3;
                var piv = new int[n];
                char[] way = { 'U' };
                var info = Factorise.Factor(way, n, S, piv);
                var root = new double[n * n];
                for (int i = 0; i < n; ++i) root[i + i * n] = 1;
                var rootinv = new double[n * n];
                var result = new double[n * n];
                for (int i = 0; i < n; ++i) rootinv[i + i * n] = 1;
                var symback = Factorise.Solve(way, n, n, S, piv, root, n, 0, 0, 0, 1, true);
                symback = Factorise.Solve(way, n, n, S, piv, rootinv, n, 0, 0, 0, -1, true);
                if (symback == 0)
                {
                    for (int i = 0; i < n; ++i)
                    {
                        Factorise.dmxmulv(n, n, root, rootinv, result, 0, i * n, i * n, true);
                    }
                    Console.WriteLine($"\n{result[0]},{result[1]},{result[2]}");
                    Console.WriteLine($"{result[3]},{result[4]},{result[5]}");
                    Console.WriteLine($"{result[6]},{result[7]},{result[8]}");
                    for (int i = 0; i < n; ++i)
                    {
                        Factorise.dmxmulv(n, n, rootinv, root, result, 0, i * n, i * n, true);
                    }
                    Console.WriteLine($"\n{result[0]},{result[1]},{result[2]}");
                    Console.WriteLine($"{result[3]},{result[4]},{result[5]}");
                    Console.WriteLine($"{result[6]},{result[7]},{result[8]}");
                }
            }
            {
                var n = 3;
                char[] way = { 'L' };
                double[] M ={1,
                           1,1,
                           1,1,1};//singular
                var S = new double[n * (n + 1) / 2];
                BlasLike.dcopyvec(S.Length, M, S);
                var piv = new int[n];
                var info = Factorise.Factor(way, n, S, piv);
                var resolve = new double[n * n];
                for (int i = 0; i < 2; ++i)
                {
                    BlasLike.dzerovec(n * n, resolve);
                    for (int ii = 0; ii < n; ++ii) resolve[ii + n * ii] = 1;
                    Factorise.Solve(way, n, n, S, piv, resolve, n, 0, 0, 0, 0, i == 0 ? false : true);
                    Console.WriteLine($"{resolve[0]} {resolve[1]} {resolve[2]} ");
                    Console.WriteLine($"{resolve[3]} {resolve[4]} {resolve[5]} ");
                    Console.WriteLine($"{resolve[6]} {resolve[7]} {resolve[8]} ");
                }

                way[0] = 'U';
                BlasLike.dcopyvec(S.Length, M, S);
                info = Factorise.Factor(way, n, S, piv);
                for (int i = 0; i < 2; ++i)
                {
                    BlasLike.dzerovec(n * n, resolve);
                    for (int ii = 0; ii < n; ++ii) resolve[ii + n * ii] = 1;
                    Factorise.Solve(way, n, n, S, piv, resolve, n, 0, 0, 0, 0, i == 0 ? false : true);
                    Console.WriteLine($"{resolve[0]} {resolve[1]} {resolve[2]} ");
                    Console.WriteLine($"{resolve[3]} {resolve[4]} {resolve[5]} ");
                    Console.WriteLine($"{resolve[6]} {resolve[7]} {resolve[8]} ");
                }
            }
            {
                double[] S = { 1.2, 4, -2.1, 5, 45 };
                var order = new int[S.Length];
                Ordering.Order.getorder(S.Length, S, order);
                Ordering.Order.Display(order, "New Order");
                Ordering.Order.Reorder_gen(S.Length, order, S);
                Ordering.Order.Display(S, "S big at the start");
                var inverse = new int[S.Length];
                for (int i = 0; i < S.Length; ++i) inverse[order[i]] = i;
                Ordering.Order.Display(inverse, "Inverse Order");
                Ordering.Order.Reorder_gen(S.Length, inverse, S);
                Ordering.Order.Display(S, "Original S");
                Ordering.Order.getorderabs(S.Length, S, order);
                Ordering.Order.Reorder_gen(S.Length, order, S);
                Ordering.Order.Display(S, "abs(S) big at the start");
                double[] MN = {11,12,13,14,15,
                            21,22,23,24,25};
                var m = 2;
                int[] ord = { 0,1,2,4,3 };
                Ordering.Order.Display(MN, "Two rows", m);
                Ordering.Order.Reorder_gen(MN.Length/m,ord,MN,m,MN.Length/m);
                Ordering.Order.Display(MN, "Two rows", m);
                Factorise.dmx_transpose(MN.Length/m,m,MN,MN);
                Ordering.Order.Display(MN, "Two Columns", MN.Length/m);
                Factorise.dmx_transpose(m,MN.Length/m,MN,MN);
                Ordering.Order.Reorder_gen(MN.Length/m,ord,MN,m,MN.Length/m);
                Ordering.Order.Display(MN, "Two Rows", m);
                Factorise.dmx_transpose(MN.Length/m,m,MN,MN);
                Ordering.Order.Reorder_gen(MN.Length/m,ord,MN,m);
                Ordering.Order.Display(MN, "Two Columns", MN.Length/m);
            }
            var isWindows = RuntimeInformation.IsOSPlatform(OSPlatform.Windows);
            if (isWindows) //Show how to read and write to Windows registry
            {
                string ourkey = "Software\\safeqp";
                if (args.Length == 1)
                {
                    ourkey = "Software\\" + args[0];
                }
                Byte[] lic = null;
                string licence = "";
                try
                {
                    RegistryKey safekey = Registry.CurrentUser, newkey;
                    if (safekey != null)
                    {
                        newkey = safekey.OpenSubKey(ourkey);
                        if (newkey == null)
                        {
                            safekey.Dispose();
                            throw new Exception("No key");
                        }
                        lic = (Byte[])newkey.GetValue(ourkey);
                        if (lic == null)
                        {
                            safekey.Dispose();
                            throw new Exception("No lic");
                        }
                        for (int i = 0; i < lic.Length; ++i)
                        {
                            licence += string.Format("{0:x2};", lic[i]);
                        }
                    }
                    Console.WriteLine($"Our Key = \t\t{ourkey}");
                    Console.WriteLine($"Current Licence: \t{licence}");
                    safekey.Dispose();
                }
                catch (Exception prob)
                {
                    Console.WriteLine("exception" + prob);
                }
                lic = new Byte[20];
                var fiddlelic = "18;ed;58;7a;e8;46;d1;6e;1d;5a;04;ae;0b;ad;66;83;ff;03;00;00";
                var il = 0;
                foreach (string ll in fiddlelic.Split(';'))
                {
                    lic.SetValue((Byte)(Convert.ToInt32(ll, 16)), il++);
                }
                licence = "";
                for (int i = 0; i < (lic != null ? lic.Length : 0); ++i)
                {
                    licence += string.Format("{0:x2};", lic[i]);
                }
                Console.WriteLine($"New Licence to Write: \t{licence}");
                try
                {
                    RegistryKey safekey = Registry.CurrentUser, newkey;
                    if (safekey != null)
                    {
                        newkey = safekey.CreateSubKey(ourkey);
                        if (newkey == null)
                        {
                            safekey.Dispose();
                            throw new Exception("No key");
                        }
                        if (lic != null && lic.Length > 0)
                            newkey.SetValue(ourkey, lic);
                        else if (newkey.GetValue(ourkey) != null)
                            newkey.DeleteValue(ourkey);
                    }
                    safekey.Dispose();
                    safekey = Registry.CurrentUser;
                    newkey = safekey.OpenSubKey(ourkey);
                    lic = (Byte[])newkey.GetValue(ourkey);
                    if (lic == null)
                    {
                        safekey.Dispose();
                        throw new Exception("No lic");
                    }
                    safekey.Dispose();
                    licence = "";
                    for (int i = 0; i < lic.Length; ++i)
                    {
                        licence += string.Format("{0:x2};", lic[i]);
                    }
                    Console.WriteLine($"Newly Written Licence: \t{licence}");
                }
                catch (Exception prob)
                {
                    Console.WriteLine("exception" + prob);
                }
                Console.WriteLine($"{ourkey} cleared");
            }
        }
    }
}


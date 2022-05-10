using System;
using System.Runtime.InteropServices;
using Blas;
using Solver;
using DataFile;
using Portfolio;
using Microsoft.Win32;

namespace UseBlas
{
    class ETLpass
    {
        public double[] returns;
        public int T;
        public double inc;
    }
    class Program
    {
        static unsafe void Main(string[] args)
        {
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
            }
            {
                var aref = new double[1];
                fixed (double* aq = aref)
                    BlasLike.baseref = 0;
                int n = 4;
                double[] aa ={1,2,4,7,
                            3,5,8,
                            6,9,
                            10};
                double[] acopy = (double[])aa.Clone();
                int[] piv = { 1, 2, 3, 4 };
                char U = 'L';
                int back = 10;
                back = Factorise.Factor(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
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
                char U = 'U';
                int back = 10;
                back = Factorise.Factor(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone(); Factorise.Solve(U, n, 1, aa, piv, b, n);
                double[] c = new double[n];
                //fixed (double* acp = acopy)
                //fixed (double* bp = b)
                //fixed (double* cp = c)
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
                char U = 'U';
                int back;
                back = Factorise.Factor(U, n, S, piv);
                Factorise.Solve(U, n, 1, S, piv, c, n);
                Factorise.dsmxmulv(n, SS, c, b);
                BlasLike.dsubvec(n, b, cc, b);
                double error = Math.Sqrt(BlasLike.ddotvec(n, b, b) / n);
                Console.WriteLine($"{U}\t\t{back}\tError: {error}");
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
                char U = 'L';
                int back;    //              back = Factorise.Factor(UU, n, SSS, pv);
                back = Factorise.Factor(U, n, S, piv); //Factorise.Solve(UU, n, 1, SSS, pv, ccc, n);
                Factorise.Solve(U, n, 1, S, piv, c, n);
                Factorise.dsmxmulvT(n, SS, c, b);
                BlasLike.dsubvec(n, b, cc, b);
                double error = Math.Sqrt(BlasLike.ddotvec(n, b, b) / n);
                Console.WriteLine($"{U}\t\t{back}\tError: {error}");
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
                char U = 'U';
                char L = 'L';
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
                char way = 'U';
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
                Console.WriteLine($"{error}");
                way = 'L';
                back = way == 'L' ? Factorise.Factor(way, n, ST, piv) : Factorise.Factor(way, n, S, piv);
                BlasLike.dzerovec(Sback.Length, Sback);
                for (int i = 0; i < n; ++i) Sback[i * n + i] = 1;
                whichroot = 2;
                info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                whichroot = 0;
                info = way == 'L' ? Factorise.Solve(way, n, n, ST, piv, Sback, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, S, piv, Sback, n, 0, 0, 0, whichroot);
                for (int i = 0; i < n; ++i) Sback[n * i + i] += -1;
                error = Math.Sqrt(BlasLike.ddotvec(n * n, Sback, Sback)) / n;
                Console.WriteLine($"{error}");
            }
            {
                var n = 500;
                var tdata = 40;
                char way = 'U';
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
                    if (way == 'U') Factorise.dsmxmulv(n, M, xxx, start, 0, i * n);
                    else Factorise.dsmxmulvT(n, MT, xxx, start, 0, i * n);
                    xxx[i] = 0;
                }
                var back = (way == 'U') ? Factorise.Factor(way, n, M, piv) : Factorise.Factor(way, n, MT, piv);
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
                var symback = (way == 'U') ? Factorise.Solve(way, n, n, M, piv, rBack, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, rBack, n, 0, 0, 0, whichroot);
                whichroot = 0;
                symback = (way == 'U') ? Factorise.Solve(way, n, n, M, piv, rI, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, rI, n, 0, 0, 0, whichroot);
                whichroot = 1;
                symback = (way == 'U') ? Factorise.Solve(way, n, n, M, piv, r, n, 0, 0, 0, whichroot) : Factorise.Solve(way, n, n, MT, piv, r, n, 0, 0, 0, whichroot);
                if (whichroot != 0 && way == 'L') Factorise.dmx_transpose(n, n, r, r);
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
                        rr2[ij] = way == 'L' ? BlasLike.ddotvec(n, lower, lower, i * n, j * n) : BlasLike.ddotvec(n, r, r, i * n, j * n);
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
                char way = 'U';
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
                char way = 'U';
                var R = new double[n * n * 2];
                var R_Inverse = new double[n * n + n * (n + 1) / 2];
                var whichroot = 0;
                var back = Factorise.Factor(way, n, S, piv, way == 'U' ? 0 : cov.Length);
                var npiv = 0;
                for (int i = 0; i < n; ++i) if (piv[i] < 0) npiv++;
                Console.WriteLine($"{npiv} negative pivots");
                BlasLike.dzerovec(R.Length, R);
                for (int i = 0; i < n; i++) R[i * n + i] = 1;
                for (int i = 0; i < n; i++) R_Inverse[i * n + i] = 1;
                for (int i = 0; i < n; i++) R[i * n + i + n * n] = 1;
                whichroot = 1;
                var symprob = Factorise.Solve(way, n, n, S, piv, R, n, way == 'U' ? 0 : cov.Length, 0, 0, whichroot);
                whichroot = -1;
                symprob = Factorise.Solve(way, n, n, S, piv, R, n, way == 'U' ? 0 : cov.Length, 0, n * n, whichroot);
                if (symprob == -10) Console.WriteLine("not positive definite!");
                whichroot = 0;
                symprob = Factorise.Solve(way, n, n, S, piv, R_Inverse, n, way == 'U' ? 0 : cov.Length, 0, 0, whichroot);
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
                char way = 'U';
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
                Console.WriteLine($"{unit[start + 0]} {unit[start + 1]} {unit[start + 2]} ");
                Console.WriteLine($"{unit[start + 3]} {unit[start + 4]} {unit[start + 5]} ");
                Console.WriteLine($"{unit[start + 6]} {unit[start + 7]} {unit[start + 8]} ");
                start = n * n;
                Factorise.dmx_transpose(n, n, unit, unit, n * n, n * n);
                Console.WriteLine($"\n{unit[start + 0]} {unit[start + 1]} {unit[start + 2]} ");
                Console.WriteLine($"{unit[start + 3]} {unit[start + 4]} {unit[start + 5]} ");
                Console.WriteLine($"{unit[start + 6]} {unit[start + 7]} {unit[start + 8]} ");
                if (way == 'L') Console.WriteLine($"{unit[start + 6]} {-unit[6] + unit[3] * unit[7]}");
                else if (way == 'U') Console.WriteLine($"{unit[start + 2]} {-unit[2] + unit[1] * unit[5]}");
            }
            {
                //double[] S = { 3, 20, 2, 0.1, -0.1, 1 }; //not positive definite
                double[] S = { 10, 0.1, 2, 0.1, 0.1, 3 }; //positive definite
                var n = 3;
                var piv = new int[n];
                char way = 'U';
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
                char way = 'L';
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

                way = 'U';
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
                var dropbad = new byte[S.Length];
                dropbad[2] = 1;
                Ordering.Order.getorder(S.Length, S, order, dropbad);
                Ordering.Order.Display(order, "New Order");
                Ordering.Order.Reorder_gen(S.Length, order, S);
                Ordering.Order.Display(S, "S big at the start");
                var inverse = new int[S.Length];
                for (int i = 0; i < S.Length; ++i) inverse[order[i]] = i;
                Ordering.Order.Display(inverse, "Inverse Order");
                Ordering.Order.Reorder_gen(S.Length, inverse, S);
                Ordering.Order.Display(S, "Original S");
                dropbad[2] = 0;
                dropbad[0] = 1;
                Ordering.Order.getorderabs(S.Length, S, order, dropbad);
                Ordering.Order.Reorder_gen(S.Length, order, S);
                Ordering.Order.Display(S, "abs(S) big at the start");
                string[] MN = { "11", "12", "13", "14", "15",
                                "21", "22", "23", "24", "25" ,
                                "31", "32", "33", "34", "35" };
                var m = 3;
                int[] ord = { 3, 4, 2, 0, 1 };
                Ordering.Order.Display(MN, "Two Rows", m);
                Ordering.Order.Reorder_gen(MN.Length / m, ord, MN, m, MN.Length / m);
                Ordering.Order.Display(MN, "Two Rows", m);
                Factorise.dmx_transpose(MN.Length / m, m, MN, MN);
                Ordering.Order.Display(MN, "Two Columns", MN.Length / m);
                Factorise.dmx_transpose(m, MN.Length / m, MN, MN);
                Ordering.Order.Reorder_gen(MN.Length / m, ord, MN, m, MN.Length / m);
                Ordering.Order.Display(MN, "Two Rows", m);
                Factorise.dmx_transpose(MN.Length / m, m, MN, MN);
                Ordering.Order.Display(MN, "Two Columns", MN.Length / m);
                Ordering.Order.Reorder_gen(MN.Length / m, ord, MN, m, 1, true);
                Ordering.Order.Display(MN, "Two Columns", MN.Length / m);
                Factorise.dmx_transpose(m, MN.Length / m, MN, MN);
                Ordering.Order.Display(MN, "Two Rows", m);
                Ordering.Order.Reorder_gen(MN.Length / m, ord, MN, m, MN.Length / m);
                Ordering.Order.Display(MN, "Two Rows", m);
                string[] symm ={"11",
                              "21","22",
                              "31","32","33",
                              "41","42","43","44",
                              "51","52","53","54","55"};
                Ordering.Order.Display(5, symm, "Symmetric", 'U');
                Ordering.Order.ReorderSymm(5, ord, symm);
                Ordering.Order.Display(5, symm, "Symmetric");
                string[] symmL ={"11","12","13","14","15",
                                   "22","23","24","25",
                                        "33","34","35",
                                             "44","45",
                                                  "55"};
                Ordering.Order.Display(5, symmL, "Symmetric", 'L');
                Ordering.Order.ReorderSymm(5, ord, symmL, 'L');
                Ordering.Order.Display(5, symmL, "Symmetric", 'L');
            }
            {
                double[] xx = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
                var n = 7;
                var nn = 4;
                var m = 2;
                Ordering.Order.Display(xx, "Before");
                Ordering.Order.bound_reorganise(1, n, nn, m, xx);
                Ordering.Order.Display(xx, "After");
                Ordering.Order.bound_reorganise(0, n, nn, m, xx);
                Ordering.Order.Display(xx, "Reset");
            }

            {
                var n = 5;
                var alpha = 2.3;
                double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
                var iswap = 0;
                var itrans = 0;
                BlasLike.detagen(n - 1, ref alpha, x, 2, ref iswap, ref itrans, 1);
                foreach (var p in x) Console.WriteLine(p);
                alpha = 2.3;
                iswap = 0;
                itrans = 0;
                double[] y = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
                BlasLike.detagen(n - 1, ref alpha, y, 2, ref iswap, ref itrans, 1);
                foreach (var p in x) Console.WriteLine(p);
            }
            {
                Console.WriteLine("--------------------ActiveSet---------------");
                var n = 10;
                var m = 2;
                var x = new double[n];
                double[] c = { 1, 2, 3, 4, 5, 6, 17, 8, 9, 10 };
                double[] A = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,
                                   0, 0, 1, 1, 1, 0, 0, 0, 0, 0};
                double[] L = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.3 };
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
                //   foreach (var a in hess) Console.WriteLine(a + ",");
                double[] hesst ={0.07622384475840693,
                                    -0.0016365991417207626,
                                    0.08604333317714313,
                                    -0.011316860823247121,
                                    -0.030434772808639987,
                                    0.11211083919582854,
                                    -0.005442735249764186,
                                    0.004464673502705685,
                                    0.02884916845203872,
                                    0.07030671358796253,
                                    0.011709823440838929,
                                    0.02086299454025381,
                                    -0.015170796975208789,
                                    -0.0005067292359139386,
                                    0.059822876920856805,
                                    0.030315371151201392,
                                    0.0034082127351833247,
                                    -0.022933127136383458,
                                    0.0016861298409444059,
                                    0.012600705278736524,
                                    0.08316673826220947,
                                    0.0115607066361883,
                                    -0.0034806816621613668,
                                    -0.0010519351552628065,
                                    -0.0028576346296797506,
                                    0.008573189337515441,
                                    -0.025115408356197216,
                                    0.08293672602298946,
                                    0.014699977532219966,
                                    -0.019735088940840084,
                                    0.05551477101220931,
                                    -0.019246450916202917,
                                    -0.01412419584221336,
                                    -0.0025011076826836065,
                                    -0.013239288762594226,
                                    0.0907984789770602,
                                    -0.013526556333058826,
                                    0.0098303768770443,
                                    -0.042483694444829995,
                                    -0.035685449490480525,
                                    0.015581944899869915,
                                    -0.0008483921768689395,
                                    -0.006279985059017501,
                                    -0.00375529364246191,
                                    0.10266907898364636,
                                    0.03999356589115988,
                                    0.045353250040677695,
                                    -0.0019274282309346968,
                                    0.010712719440722496,
                                    0.024816489058402724,
                                    0.03416835576827826,
                                    -0.001184721225989449,
                                    -0.011735959703553567,
                                    -0.04135384837511391,
                                    0.12585651441173307};
                hess = hesst;
                var obj = -9.6;
                var iter = 0;
                short back;
                BlasLike.dsetvec(x.Length, 1.0 / n, x);
                back = ActiveSet.Optimise.LPopt(n, m, x, L, U, A, c, ref obj, ref iter);
                Console.WriteLine($"back is {back} {BlasLike.ddotvec(n, x, c)} {obj} {iter} iterations");
                ActiveSet.Optimise.printV("x", x);
                var implied = new double[m];
                Factorise.dmxmulv(m, n, A, x, implied);
                foreach (var cc in implied) Console.WriteLine($"Constraint value {cc}");
                BlasLike.dsetvec(x.Length, 1.0 / n, x);
                BlasLike.dscalvec(hess.Length, 1e3, hess);
                var opt = new ActiveSet.Optimise();
                if (opt.h == null) opt.h = opt.qphess1;
                back = opt.QPopt(n, m, x, L, U, A, c, hess, ref obj, ref iter);
                implied = new double[n];
                Factorise.dsmxmulv(n, hess, x, implied);
                Console.WriteLine($"back is {back} {BlasLike.ddotvec(n, x, c) + 0.5 * BlasLike.ddotvec(n, implied, x)} {obj}  {iter} iterations");
                ActiveSet.Optimise.printV("x", x);
                implied = new double[m];
                Factorise.dmxmulv(m, n, A, x, implied);
                foreach (var cc in implied) Console.WriteLine($"Constraint value {cc}");
            }
            {
                Console.WriteLine("----------------------InteriorPoint---------------");
                int nh = 10;
                double[] H ={0.07622384475840693,
                                    -0.0016365991417207626,
                                    0.08604333317714313,
                                    -0.011316860823247121,
                                    -0.030434772808639987,
                                    0.11211083919582854,
                                    -0.005442735249764186,
                                    0.004464673502705685,
                                    0.02884916845203872,
                                    0.07030671358796253,
                                    0.011709823440838929,
                                    0.02086299454025381,
                                    -0.015170796975208789,
                                    -0.0005067292359139386,
                                    0.059822876920856805,
                                    0.030315371151201392,
                                    0.0034082127351833247,
                                    -0.022933127136383458,
                                    0.0016861298409444059,
                                    0.012600705278736524,
                                    0.08316673826220947,
                                    0.0115607066361883,
                                    -0.0034806816621613668,
                                    -0.0010519351552628065,
                                    -0.0028576346296797506,
                                    0.008573189337515441,
                                    -0.025115408356197216,
                                    0.08293672602298946,
                                    0.014699977532219966,
                                    -0.019735088940840084,
                                    0.05551477101220931,
                                    -0.019246450916202917,
                                    -0.01412419584221336,
                                    -0.0025011076826836065,
                                    -0.013239288762594226,
                                    0.0907984789770602,
                                    -0.013526556333058826,
                                    0.0098303768770443,
                                    -0.042483694444829995,
                                    -0.035685449490480525,
                                    0.015581944899869915,
                                    -0.0008483921768689395,
                                    -0.006279985059017501,
                                    -0.00375529364246191,
                                    0.10266907898364636,
                                    0.03999356589115988,
                                    0.045353250040677695,
                                    -0.0019274282309346968,
                                    0.010712719440722496,
                                    0.024816489058402724,
                                    0.03416835576827826,
                                    -0.001184721225989449,
                                    -0.011735959703553567,
                                    -0.04135384837511391,
                                    0.12585651441173307};
                int n = 12;
                BlasLike.dscalvec(H.Length, 1e3, H);
                var m = 3;
                var x = new double[n];
                var nslack = 2;
                double[] b = { 1, 0.5, 0.3 };
                double[] c = { 1, 2, 3, 4, 5, 6, 17, 8, 9, 10, 0, 0 };
                double[] A = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                               0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0,
                               0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1};
                Factorise.dmx_transpose(n, m, A, A);
                var opt1 = new InteriorPoint.Optimise(n, m, x, A, b, c, nh, H);
                var back = opt1.Opt("QP", null, null, false);
                Console.WriteLine($"{back}");
                var implied = new double[m];
                var truex = (double[])x.Clone();
                BlasLike.dzerovec(nslack, truex, n - nslack);
                Factorise.dmxmulv(m, n, A, truex, implied);
                foreach (var cc in implied) Console.WriteLine($"Constraint value {cc}");
                var cx = BlasLike.ddotvec(c.Length, c, truex);
                implied = new double[nh];
                Factorise.dsmxmulv(nh, H, x, implied);
                var xHx = BlasLike.ddotvec(nh, implied, truex);
                Console.WriteLine($"Linear {cx}\nQuadratic {cx + xHx / 2}");
            }

            {
                Console.WriteLine("SOCP");
                int n = 12;
                var m = 2;
                var x = new double[n];
                double[] b = { 1, 1 };
                double[] c = { 1, 2, 3, 4, 5, 6000, 7, 8, 9, 10, 11, 0 };
                double[] A = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0 ,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

                int[] cone = { n };
                int nvar = 0;
                foreach (int ic in cone) nvar += ic;
                int[] typecone = { (int)InteriorPoint.conetype.SOCP };

                Factorise.dmx_transpose(n, m, A, A);
                var opt1 = new InteriorPoint.Optimise(n, m, x, A, b, c);
                var back = opt1.Opt("SOCP", cone, typecone, true);
                Console.WriteLine($"{back}");
                var implied = new double[m];
                var truex = (double[])x.Clone();
                Factorise.dmxmulv(m, n, A, truex, implied);
                foreach (var cc in implied) Console.WriteLine($"Constraint value {cc}");
                var cx = BlasLike.ddotvec(c.Length, c, truex);
                Console.WriteLine($"Linear {cx}");
                Console.WriteLine($"SOCP x check {Math.Sqrt(BlasLike.ddotvec(x.Length - 1, x, x))} {x[x.Length - 1]}");
            }
            {
                Console.WriteLine("--------------FMP-------------");
                var nfac = 2;
                var n = 3;
                double[] SV = { 1, 2, 3 };
                var Q = new double[(nfac + 1) * n];
                double[] FC = { 2,
                                1, 3 };
                double[] FL = { 1, 0, 1,
                                 0, 1, 1 }; //Factors by assets
                bool needTransposeInFMP = true; //Set FL by factors here, don't transpose in FMP
                if (needTransposeInFMP) Factorise.dmx_transpose(n, nfac, FL, FL);
                var back = Factorise.FMP(n, nfac, FC, SV, FL, Q, 'U', !needTransposeInFMP);
                if (back > 0) Console.WriteLine($"unstable FC: n={back}");
                else if (back == -10) Console.WriteLine($"FC is not positive definite!!!!!!");
                else Console.WriteLine($"Condition {back} from solver");
                var result = new double[n * n];
                /*    for (var i = 0; i < n; ++i) // Same as below, shows how dmxmulv works!
                    {
                        for (var k = 0; k < nfac; ++k)
                        {
                            BlasLike.daxpy(n, Q[n + i * nfac + k], Q, nfac, result, 1, n + k, i * n);

                        }
                    }*/
                for (var i = 0; i < n; ++i) //Multiply out the factor part of the compressed risk model
                {
                    Factorise.dmxmulv(n, nfac, Q, Q, result, n, n + nfac * i, i * n, true);
                }
                ActiveSet.Optimise.printV("Compressed risk model", Q);
                ActiveSet.Optimise.printV("Check the factor part 2 1 3 1 3 4 3 4 7", result);
                var order = new int[n];
                Ordering.Order.getorder(n, Q, order, null, 0, 1, 0);
                Ordering.Order.Display(order, "New Order");
                Ordering.Order.Reorder(n, order, Q);
                bool trans = false; // Transpose to use columns
                if (!trans) Factorise.dmx_transpose(nfac, n, Q, Q, n, n);
                Ordering.Order.Reorder_gen(n, order, Q, nfac, n, trans, n);
                if (!trans) Factorise.dmx_transpose(n, nfac, Q, Q, n, n);
                ActiveSet.Optimise.printV("Compressed risk model", Q);
                //    BlasLike.dzerovec(n * n, result);
                for (var i = 0; i < n; ++i) //Multiply out the factor part of the compressed risk model
                {
                    Factorise.dmxmulv(n, nfac, Q, Q, result, n, n + nfac * i, i * n, true);
                }
                ActiveSet.Optimise.printV("Check the factor part 7 4 3 4 3 1 3 1 2", result);
            }
            {
                Console.WriteLine("--------------Data File-------------");
                using (var TestData = new InputSomeData()) // We can use using because DataFile.InputSomeData has Dispose() method
                {
                    TestData.stringFields = "names";
                    TestData.intFields = "n m";
                    try
                    {
                        TestData.Read("./UseBlas/testData");
                    }
                    catch
                    {
                        TestData.Read("testData");
                    }
                    TestData.PrintField("n");
                    TestData.PrintField("m");
                    TestData.PrintField("c");
                    TestData.PrintField("alpha");
                    TestData.PrintField("names");
                }
                using (var TestData = new InputSomeData())
                {
                    Console.WriteLine("--------------Test FMP with real data from file-------------");
                    TestData.doubleFields = "FC SV FL";
                    TestData.intFields = "n nfac";
                    TestData.stringFields = "names";
                    try
                    {
                        TestData.Read("./pylog.log");
                    }
                    catch
                    {
                        TestData.Read("../pylog.log");
                    }
                    var n = TestData.mapInt["n"][0];
                    var nfac = TestData.mapInt["nfac"][0];
                    var SV = TestData.mapDouble["SV"];
                    var FL = TestData.mapDouble["FL"];
                    var FC = TestData.mapDouble["FC"];
                    var FCl = (double[])FC.Clone();
                    for (int ij = 0, i = 0; i < nfac; ++i)
                    {
                        for (var j = i; j < nfac; ++j, ij++)
                        {
                            FCl[ij] = FC[j * (j + 1) / 2 + i];
                        }
                    }
                    var names = TestData.mapString == null ? null : TestData.mapString["names"];// All unwanted data in the file will go into names
                    if (names != null)
                    {
                        Console.WriteLine($"{names.Length} names");
                        Array.Resize(ref names, n);
                        Console.WriteLine($"{names.Length} names");
                        TestData.mapString["names"] = names;
                    }
                    TestData.Write();
                    var Q = new double[(nfac + 1) * n];
                    var result = (double[])new double[n * n];
                    var back = Factorise.FMP(n, nfac, FCl, SV, FL, Q, 'L');
                    if (back == -10) Console.WriteLine("Factor Covraince matrix is not positive definite");
                    var upto = 30;
                    for (var i = 0; i < n; ++i) //Multiply out the factor part of the compressed risk model
                    {
                        Factorise.dmxmulv(n, nfac, Q, Q, result, n, n + nfac * i, i * n, true);
                    }
                    for (var i = 0; i < n; ++i)
                    {
                        result[i * n + i] += Q[i];
                    }
                    ActiveSet.Optimise.printV("Asset COV via lower case", result, upto);
                    back = Factorise.FMP(n, nfac, FC, SV, FL, Q);
                    if (back == -10) Console.WriteLine("Factor Covraince matrix is not positive definite");
                    for (var i = 0; i < n; ++i) //Multiply out the factor part of the compressed risk model
                    {
                        Factorise.dmxmulv(n, nfac, Q, Q, result, n, n + nfac * i, i * n, true);
                    }
                    for (var i = 0; i < n; ++i)
                    {
                        result[i * n + i] += Q[i];
                    }
                    ActiveSet.Optimise.printV("Asset COV via upper case", result, upto);
                    var w = new double[n];
                    var Qw = new double[n];
                    w[0] = 1;
                    Factorise.FacMul(n, nfac, Q, w, Qw);
                    var variance = BlasLike.ddotvec(n, w, Qw);
                    Console.WriteLine($"variance = {variance:F8}");
                    ActiveSet.Optimise.printV("Test first row of assembled covariance matrix", Qw, upto);
                    var COV = new double[n * (n + 1) / 2];
                    Factorise.Fac2Cov(n, nfac, Q, COV, 0, 0, 'L');
                    Factorise.CovMul(n, COV, w, Qw, 0, 0, 0, 'L');
                    ActiveSet.Optimise.printV("Test first row of generated covariance matrix 'L'", Qw, upto);
                    Factorise.Fac2Cov(n, nfac, Q, COV);
                    Factorise.CovMul(n, COV, w, Qw);
                    ActiveSet.Optimise.printV("Test first row of generated covariance matrix 'U'", Qw, upto);
                }
            }
            {
                Console.WriteLine("-------------------------Portfolio----------------");
                FPortfolio port = new FPortfolio("");
                using (var OptData = new InputSomeData())
                {
                    OptData.intFields = "n m nfac";
                    OptData.doubleFields = "delta kappa gamma alpha Q L U A buy sell initial bench FC FL SV";
                    OptData.stringFields = "names";
                    string inFile = "pylog.log";
                    try { OptData.Read(inFile); }
                    catch { OptData.Read("../" + inFile); }
                    port.n = OptData.mapInt["n"][0];
                    port.m = OptData.mapInt["m"][0];
                    port.gamma = OptData.mapDouble["gamma"][0];
                    port.kappa = OptData.mapDouble["kappa"][0];
                    port.delta = OptData.mapDouble["delta"][0];
                    port.c = OptData.mapDouble["alpha"];
                    BlasLike.dnegvec(port.c.Length, port.c);
                    port.initial = OptData.mapDouble["initial"];
                    port.bench = OptData.mapDouble["bench"];
                    port.L = OptData.mapDouble["L"];
                    port.U = OptData.mapDouble["U"];
                    port.A = OptData.mapDouble["A"];
                    port.buy = OptData.mapDouble["buy"];
                    port.sell = OptData.mapDouble["sell"];
                    port.names = OptData.mapString["names"];
                    port.nfac = OptData.mapInt["nfac"][0];
                    port.FC = OptData.mapDouble["FC"];
                    port.FL = OptData.mapDouble["FL"];
                    port.SV = OptData.mapDouble["SV"];
                    port.names = OptData.mapString["names"];
                    if (port.names != null) Array.Resize(ref port.names, port.n);
                }
                port.OptimiseTest();
            }
            {
                ColourConsole.WriteLine("Conditional Value at Risk", ConsoleColor.DarkYellow);
                Portfolio.Portfolio opt = new Portfolio.Portfolio("");
                double[] DATA;
                string[] names;
                int nstocks, tlen;
                using (var CVarData = new InputSomeData())
                {
                    CVarData.doubleFields = "DATA";//The DATA is LOSS i.e. - returns
                    CVarData.intFields = "n tlen";
                    CVarData.stringFields = "names";
                    try
                    {
                        CVarData.Read("./GL");
                    }
                    catch
                    {
                        CVarData.Read("../GL");
                    }
                    DATA = CVarData.mapDouble["DATA"];
                    names = CVarData.mapString["names"];
                    nstocks = CVarData.mapInt["n"][0];
                    tlen = CVarData.mapInt["tlen"][0];
                }
                var testw = new double[nstocks];
                BlasLike.dnegvec(DATA.Length, DATA);
                //double[] ox = { 0.17942402, 0.01108059, 0.01129700, 0.05966976, 0.19277110, 0.01090517, 0.22144947, 0.30223059, 0.00843431, 0.00273799 };
                //testw = (double[])ox.Clone();
                BlasLike.dsetvec(nstocks, 1.0 / nstocks, testw);
                var s = new double[tlen];
                double tail = 0.05;//confidence is (1-tail) I think, so 95% confident here
                ETLpass info = new ETLpass();
                info.inc = tail;
                info.T = tlen;
                info.returns = s;
                double smax = 1, smin = -1;
                Factorise.dmxmulv(tlen, nstocks, DATA, testw, s);
                BlasLike.dxminmax(DATA.Length, DATA, 1, ref smax, ref smin);
                smax *= (double)nstocks;
                smin = -smax;
                ///<summary>Traditional way to get VAR and CVAR</summary>
                double varorder(ref double cvar, int include)
                {
                    var order = new int[tlen];
                    var iorder = new int[tlen];
                    Ordering.Order.getorder(tlen, s, order, null, 0, -1);
                    for (var i = 0; i < tlen; ++i) iorder[order[i]] = i;
                    Ordering.Order.Reorder(tlen, order, s);
                    cvar = BlasLike.dsumvec(include, s, tlen - include) / include;
                    Ordering.Order.Reorder(tlen, iorder, s);
                    var back = s[order[tlen - include - 1]];
                    return back;
                }
                double cvarord1 = 0, cvarord2 = 0;
                var VARord1 = varorder(ref cvarord1, (int)(tail * tlen));
                info.inc = (double)((int)(tail * tlen)) / tlen;
                var cavrordo1 = cvar(VARord1, info);
                var VARord2 = varorder(ref cvarord2, (int)(tail * tlen + 1));
                info.inc = (double)((int)(tail * tlen + 1)) / tlen;
                var cavrordo2 = cvar(VARord2, info);

                var VARinterord = (VARord2 - VARord1) * (tail * tlen - (double)(int)(tail * tlen)) + VARord1;
                var ETLinterord = (cvarord2 - cvarord1) * (tail * tlen - (double)(int)(tail * tlen)) + cvarord1;

                ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Traditional Method[/darkyellow]\n[green]Interpolated VAR {VARinterord,16:E8}[/green] [darkyellow]Interpolated ETL {ETLinterord,16:E8}[/darkyellow]");
                double cvar(double X, object kk)
                {
                    ETLpass info = (ETLpass)kk;
                    double r = 0;
                    for (var i = 0; i < info.T; ++i)
                    {
                        r += Math.Max(0, info.returns[i] - X);
                    }
                    var back = X + r / Portfolio.Portfolio.check_digit(info.inc * tlen);
                    return back;
                }
                double cvar1d(ref double var1, object kk)
                {
                    ETLpass info = (ETLpass)kk;
                    var tailkeep = tail;
                    var1 = ActiveSet.Optimise.PathMin(cvar, smin, smax, BlasLike.lm_eps8, 0, kk);
                    var back = cvar(var1, kk);
                    ColourConsole.WriteEmbeddedColourLine($"[yellow]Value at Risk {var1,16:E8}[/yellow], [green]Expected Tail Loss {back,16:E8}[/green] [magenta]at {info.inc}[/magenta]");
                    tail = tailkeep;
                    return back;
                }
                ColourConsole.WriteInfo("Don't need to reorder x any more, use optimisation instead");
                double VAR = -1;
                info.inc = tail;
                var ETL = cvar1d(ref VAR, info);
                var tail1 = (double)(long)((tail * tlen)) / tlen;
                double VAR1 = -1;
                info.inc = tail1;
                var ETL1 = cvar1d(ref VAR1, info);
                var tail2 = (double)(long)((tail * tlen + 1)) / tlen;
                double VAR2 = -1;
                info.inc = tail2;
                var ETL2 = cvar1d(ref VAR2, info);
                var VARinter = (VAR2 - VAR1) * (tail - tail1) * tlen + VAR1;
                var ETLinter = (ETL2 - ETL1) * (tail - tail1) * tlen + ETL1;

                ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Optimisation Method[/darkyellow]\n[green]Interpolated VAR {VARinter,16:E8}[/green] [darkyellow]Interpolated ETL {ETLinter,16:E8}[/darkyellow]");
                info.inc = tail;
                var ETLinter2 = cvar(VARinter, info);
                //We find ETLinter2 is the same as ETL, i.e. true optimised CVAR using inferred VAR at 0.05
                //is the same as CVAR calculated using 0.05 with the interpolated VAR from the actual 
                //time variables above and below 0.05.
                ColourConsole.WriteEmbeddedColourLine($"[green]CVAR({VARinter,16:E8})[/green] [darkyellow]gives ETL {ETLinter2,16:E8}[/darkyellow] ([red]{ETL - ETLinter2:e16}[/red])");

                //Now try LP
                var m = tlen;
                var n = tlen + 1;
                double[] x = new double[n];
                double[] LL = new double[n + m];
                double[] L = new double[n + m];
                double[] U = new double[n + m];
                double[] A = new double[n * m];
                double[] c = new double[n];
                BlasLike.dsetvec(tlen, 1.0 / (tail * tlen), c);
                BlasLike.dsetvec(1, 1.0, c, tlen);
                BlasLike.dsetvec(tlen, 0, L);
                BlasLike.dsetvec(tlen, smax, U);
                BlasLike.dsetvec(1, smin, L, tlen);
                BlasLike.dsetvec(1, smax, U, tlen);
                BlasLike.dsetvec(m, smax, U, n);
                for (var i = 0; i < m; ++i)
                {
                    L[n + i] = Math.Max(s[i], 0);
                }
                for (var i = 0; i < tlen; ++i)
                {
                    BlasLike.dset(1, 1, A, m, i + m * i);
                    BlasLike.dset(1, 1, A, m, i + m * tlen);
                }
                opt.n = n;
                opt.m = m;
                opt.L = L;
                opt.U = U;
                opt.A = A;
                opt.c = c;
                opt.names = new string[n];
                for (var i = 0; i < tlen; ++i)
                {
                    opt.names[i] = "time" + (i + 1);
                }
                opt.names[tlen] = "VAR";
                BlasLike.dsetvec(n, 1.0 / n, x);
                double back = -12;
                bool useIP = false;
                if (!useIP) back = opt.ActiveOpt(1, x, LL);
                else back = opt.InteriorOpt(1e-9, x, LL);
                var ccc = new double[m];
                Factorise.dmxmulv(m, n, A, x, ccc);
                var cvarLP = x[tlen] + BlasLike.dsumvec(m, x) / (tail * tlen);
                var cvarLPVAR = cvar(x[tlen], info);
                ColourConsole.WriteEmbeddedColourLine($"[green]LP objective[/green]\t\t\t\t\t[darkyellow]{BlasLike.ddotvec(n, c, x),16:E8}[/darkyellow]");
                ColourConsole.WriteEmbeddedColourLine($"[green]LP gives VAR={x[tlen],16:E8}[/green] and [darkyellow]CVAR =\t{cvarLP,16:E8}[/darkyellow]");
                ColourConsole.WriteEmbeddedColourLine($"[green]CVAR({x[tlen],16:E8})[/green] =\t\t\t[darkyellow]{cvarLP,16:E8}[/darkyellow]");
                ColourConsole.WriteInfo("Optimise to get best portfolio weights");
                opt.n = nstocks + tlen + 1;
                opt.m = tlen + 1;
                opt.A = new double[opt.n * opt.m];
                opt.L = new double[opt.n + opt.m];
                opt.U = new double[opt.n + opt.m];
                opt.c = new double[opt.n];
                opt.ntrue = nstocks;
                opt.mtrue = 1;
                BlasLike.dsetvec(tlen, 1.0 / Portfolio.Portfolio.check_digit(info.inc * tlen), opt.c, nstocks);
                BlasLike.dsetvec(1, 1.0, opt.c, nstocks + tlen);
                BlasLike.dsetvec(nstocks, 0, opt.L);
                BlasLike.dsetvec(nstocks, 1, opt.U);
                BlasLike.dsetvec(tlen, 0, opt.L, nstocks);
                BlasLike.dsetvec(tlen, smax, opt.U, nstocks);
                BlasLike.dsetvec(1, smin, opt.L, nstocks + tlen);
                BlasLike.dsetvec(1, smax, opt.U, nstocks + tlen);
                BlasLike.dsetvec(1, 1, opt.L, opt.n);
                BlasLike.dsetvec(1, 1, opt.U, opt.n);
                BlasLike.dsetvec(tlen, 0, opt.L, opt.n + 1);
                BlasLike.dsetvec(tlen, smax, opt.U, opt.n + 1);

                BlasLike.dset(nstocks, 1, opt.A, opt.m);
                for (var i = 0; i < tlen; ++i)
                {//These constraints ensure that the optimisation variables are constrained so that (LOSS above VAR)+VAR-portfolio LOSS>=0
                    BlasLike.dsccopy(nstocks, -1, DATA, tlen, opt.A, opt.m, i, i + 1);
                    BlasLike.dset(1, 1, opt.A, opt.m, 1 + i + opt.m * (i + nstocks));
                    BlasLike.dset(1, 1, opt.A, opt.m, 1 + i + opt.m * (tlen + nstocks));
                }
                x = new double[opt.n];
                LL = new double[opt.n + opt.m];
                opt.names = (string[])names.Clone();
                Array.Resize(ref opt.names, opt.n);
                for (var i = 0; i < tlen; ++i)
                {
                    opt.names[i + nstocks] = "time" + (i + 1);
                }
                opt.names[opt.n - 1] = "VAR";

                if (!useIP) back = opt.ActiveOpt(1, x, LL);
                else back = opt.InteriorOpt(1e-9, x, LL);
                var VARopt = x[x.Length - 1];
                var CVARopt = BlasLike.ddotvec(opt.n, x, opt.c);
                var xopt = (double[])x.Clone();
                Array.Resize(ref xopt, nstocks);

                Factorise.dmxmulv(tlen, nstocks, DATA, xopt, s);

                info.inc = tail;
                var ETLcheck = cvar(VARopt, info);
                ETL = cvar1d(ref VAR, info);
                ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Using optimised portfolio weights[/darkyellow]\n[cyan]VAR\t{VARopt,16:E8} [/cyan] [yellow]Check using cavr1d {VAR,16:E8}[/yellow]\n[cyan]CVAR\t{CVARopt,16:E8} [/cyan] [yellow]Check using cvar1d {ETL,16:E8}[/yellow] [green]Check using cvar({VARopt,16:E8}) {ETLcheck,16:E8}[/green]");
                return;
            }

            {
                Console.WriteLine("GAIN/LOSS");
                Portfolio.Portfolio opt = new Portfolio.Portfolio("");
                double[] DATA;
                double lambda, R;
                string[] names;
                int n, tlen;
                using (var GainLossData = new InputSomeData())
                {
                    GainLossData.doubleFields = "DATA lambda R";
                    GainLossData.intFields = "n tlen";
                    GainLossData.stringFields = "names";
                    try
                    {
                        GainLossData.Read("./GL");
                    }
                    catch
                    {
                        GainLossData.Read("../GL");
                    }
                    DATA = GainLossData.mapDouble["DATA"];
                    names = GainLossData.mapString["names"];
                    n = GainLossData.mapInt["n"][0];
                    lambda = GainLossData.mapDouble["lambda"][0];
                    R = GainLossData.mapDouble["R"][0];
                    tlen = GainLossData.mapInt["tlen"][0];
                }
                bool useIP = false;
                opt.GainLossSetUp(n, tlen, DATA, names, R, lambda, useIP);
            }
            {
                var filename = "costlog";
                Console.WriteLine("BUY/SELL and LONG/SHORT");
                double[] SV = null, FC = null, FL = null, L = null, U = null, alpha = null, initial = null, A = null;
                double[] buy = null, sell = null, bench = null, Q = null, A_abs = null, Abs_U = null, Abs_L = null;
                double gamma, delta, kappa, value, valuel, rmin, rmax, min_holding, min_trade, minRisk;
                double[] minhold = new double[0], mintrade = new double[0], min_lot = null, size_lot = null;
                int n, nfac, m, nabs, mabs, round, basket, tradenum;
                int[] I_A = null;
                string[] names;
                using (var buysell = new InputSomeData())
                {
                    buysell.doubleFields = "Q SV FC FL L U alpha initial A buy sell gamma delta kappa bench value valuel rmin rmax Abs_U Abs_L A_abs min_holding min_trade min_lot size_lot minRisk";
                    buysell.intFields = "n nfac m nabs mabs I_A round basket tradenum";
                    buysell.stringFields = "names";
                    try
                    {
                        buysell.Read("./" + filename);
                    }
                    catch
                    {
                        buysell.Read("../" + filename);
                    }
                    Q = buysell.mapDouble["Q"];
                    try { SV = buysell.mapDouble["SV"]; }
                    catch { SV = null; }
                    try { FC = buysell.mapDouble["FC"]; }
                    catch { FC = null; }
                    try { FL = buysell.mapDouble["FL"]; }
                    catch { FL = null; }
                    L = buysell.mapDouble["L"];
                    U = buysell.mapDouble["U"];
                    alpha = buysell.mapDouble["alpha"];
                    initial = buysell.mapDouble["initial"];
                    bench = buysell.mapDouble["bench"];
                    A = buysell.mapDouble["A"];
                    buy = buysell.mapDouble["buy"];
                    sell = buysell.mapDouble["sell"];
                    n = buysell.mapInt["n"][0];
                    m = buysell.mapInt["m"][0];
                    delta = buysell.mapDouble["delta"][0];
                    kappa = buysell.mapDouble["kappa"][0];
                    gamma = buysell.mapDouble["gamma"][0];
                    minRisk = buysell.mapDouble["minRisk"][0];
                    value = buysell.mapDouble["value"][0];
                    valuel = buysell.mapDouble["valuel"][0];
                    rmin = buysell.mapDouble["rmin"][0];
                    rmax = buysell.mapDouble["rmax"][0];
                    nfac = buysell.mapInt["nfac"][0];
                    names = buysell.mapString["names"];
                    nabs = buysell.mapInt["nabs"][0];
                    mabs = buysell.mapInt["mabs"][0];
                    round = buysell.mapInt["round"][0];
                    basket = buysell.mapInt["basket"][0];
                    tradenum = buysell.mapInt["tradenum"][0];
                    A_abs = buysell.mapDouble["A_abs"];
                    Abs_U = buysell.mapDouble["Abs_U"];
                    Abs_L = buysell.mapDouble["Abs_L"];
                    I_A = buysell.mapInt["I_A"];
                    min_holding = buysell.mapDouble["min_holding"][0];
                    min_trade = buysell.mapDouble["min_trade"][0];
                    min_lot = buysell.mapDouble["min_lot"];
                    size_lot = buysell.mapDouble["size_lot"];
                }
                /*   L[0] = -0.017583279913421082;
                      U[0] = -0.017583279913421082;
                      L[n - 10] = 0;//-1e-3;
                      U[n - 10] = 0;//1e-3;
                      L[n - 20] = 0;//-1e-3;
                      U[n - 20] = 0;//1e-3;*/
                bool useIp = true;
                double[] minlot = null;
                double[] sizelot = null;
                var roundw = new double[n];
                var shake = new int[n];
                for (var i = 0; i < n; ++i) shake[i] = -1;
                if (round == 1)
                {
                    if (min_lot != null && min_lot.Length == 1)
                    {
                        minlot = new double[n];
                        BlasLike.dsetvec(n, min_lot[0], minlot);
                    }
                    else if (min_lot != null && min_lot.Length == n)
                        minlot = min_lot;
                    if (size_lot != null && size_lot.Length == 1)
                    {
                        sizelot = new double[n];
                        BlasLike.dsetvec(n, size_lot[0], sizelot);
                    }
                    else if (size_lot != null && size_lot.Length == n)
                        sizelot = size_lot;
                }
                if (min_holding >= 0)
                {
                    minhold = new double[n];
                    BlasLike.dsetvec(n, min_holding, minhold);
                }
                else minhold = null;
                if (min_trade >= 0)
                {
                    mintrade = new double[n];
                    BlasLike.dsetvec(n, min_trade, mintrade);
                }
                else mintrade = null;
                var trades = tradenum;
                useIp = false;
                Portfolio.Portfolio.OptParamRound Op = new Portfolio.Portfolio.OptParamRound();
                Portfolio.Portfolio.INFO sendInput = new Portfolio.Portfolio.INFO(); Op.basket = basket;
                Op.trades = trades;
                Op.lower = L;
                Op.m = m;
                Op.n = n;
                Op.upper = U;
                Op.minholdlot = null;
                Op.mintradelot = null;
                sendInput.n = n;
                sendInput.m = m;
                sendInput.nfac = nfac;
                sendInput.A = A;
                sendInput.L = L;
                sendInput.U = U; sendInput.delta = delta;
                sendInput.value = value;
                sendInput.valuel = valuel;
                sendInput.rmin = rmin;
                sendInput.rmax = rmax;
                sendInput.alpha = alpha;
                sendInput.initial = initial;
                sendInput.buy = buy;
                sendInput.sell = sell;
                sendInput.names = names;
                sendInput.useIP = useIp;
                sendInput.nabs = nabs;
                sendInput.A_abs = A_abs;
                sendInput.L_abs = Abs_L;
                sendInput.U_abs = Abs_U;
                sendInput.mabs = mabs;
                sendInput.I_a = I_A;
                sendInput.bench = bench;
                if (nfac > -1)
                {
                    FPortfolio opt = new FPortfolio("");
                    opt.SV = SV;
                    opt.FL = FL;
                    opt.FC = FC;
                    opt.nfac = nfac;
                    opt.bench = bench;
                    var targetRisk = minRisk;
                    sendInput.target = targetRisk;
                    opt.gamma = opt.kappa = 5e-2;
                    opt.CalcRisk(opt.gamma, sendInput);
                    opt.BoundsSetToSign(n, sendInput.L, sendInput.U, initial, opt.wback);
                    //  opt.DropRisk(basket, trades, targetRisk, sendInput);
                    sendInput.useIP = false;
                    if (round == 1)
                    {
                        opt.Rounding(basket, trades, initial, minlot, sizelot, roundw, null, null, sendInput);
                        opt.roundcheck(n, roundw, initial, minlot, sizelot, shake);
                    }
                    else
                    {
                        Op.x = opt.wback; Op.MoreInfo = sendInput;
                        opt.Thresh(Op, mintrade == null ? null : initial, mintrade == null ? minhold : mintrade, roundw, mintrade == null ? null : minhold);
                        opt.thresh_check(n, roundw, mintrade == null ? null : initial, L, U, mintrade == null ? minhold : mintrade, mintrade == null ? null : minhold, BlasLike.lm_eps8, shake);
                    }
                    foreach (var i in shake)
                    {
                        if (i != -1) ColourConsole.WriteEmbeddedColourLine($"[green]{names[i]}[/green][red] was not rounded properly! {roundw[i],26:e16}[/red]");
                    }
                }
                else
                {
                    Portfolio.Portfolio opt = new Portfolio.Portfolio("");
                    opt.Q = Q;
                    opt.bench = bench;
                    var targetRisk = minRisk;
                    sendInput.target = targetRisk;
                    opt.gamma = opt.kappa = 0.5;
                    opt.CalcRisk(opt.gamma, sendInput);
                    opt.BoundsSetToSign(n, sendInput.L, sendInput.U, initial, opt.wback);
                    // opt.DropRisk(basket, trades, targetRisk, sendInput);
                    sendInput.useIP = false;
                    if (round == 1)
                    {
                        opt.Rounding(basket, trades, initial, minlot, sizelot, roundw, null, null, sendInput);
                        opt.roundcheck(n, roundw, initial, minlot, sizelot, shake);
                    }
                    else
                    {
                        Op.x = opt.wback; Op.MoreInfo = sendInput;
                        opt.Thresh(Op, mintrade == null ? null : initial, mintrade == null ? minhold : mintrade, roundw, mintrade == null ? null : minhold);
                        opt.thresh_check(n, roundw, mintrade == null ? null : initial, L, U, mintrade == null ? minhold : mintrade, mintrade == null ? null : minhold, BlasLike.lm_eps8, shake);
                    }
                    foreach (var i in shake)
                    {
                        if (i != -1) ColourConsole.WriteEmbeddedColourLine($"[green]{names[i]}[/green][red] was not rounded properly! {roundw[i],26:e16}[/red]");
                    }
                }
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
                var fiddlelic = "9e;d2;25;3f;29;c5;56;d6;ef;7e;21;82;13;31;f6;d6;03;02;00;00;";
                if (fiddlelic.EndsWith(';'))
                {
                    fiddlelic = fiddlelic.Remove(fiddlelic.Length - 1, 1);
                }
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
                        throw new Exception("No licence!");
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


using System;
using System.Runtime.InteropServices;
using Blas;
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
                back = Factorise.dsptrf(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    //            Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
                    Factorise.dsptrs(U, n, 1, aa, piv, b, n);
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
                back = Factorise.dsptrf(U, n, aa, piv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]} {piv[3]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} {aa[6]} {aa[7]} {aa[8]} {aa[9]} ");
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    //               Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
                    Factorise.dsptrs(U, n, 1, aa, piv, b, n);
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
                    // back = Factorise.dsptrf(UU, n, SSS, pv);
                    // Factorise.dsptrs(UU, n, 1, SSS, pv, ccc, n);
                }
                back = Factorise.dsptrf(U, n, S, piv);
                Factorise.dsptrs(U, n, 1, S, piv, c, n);
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
                    //              back = Factorise.dsptrf(UU, n, SSS, pv);
                    back = Factorise.dsptrf(U, n, S, piv);
                fixed (double* SSS = S)
                fixed (int* pv = piv)
                fixed (char* UU = U)
                fixed (double* ccc = c)
                    //Factorise.dsptrs(UU, n, 1, SSS, pv, ccc, n);
                    Factorise.dsptrs(U, n, 1, S, piv, c, n);
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
                var back = Factorise.dsptrf(U, n, S, ipiv);
                Factorise.dsptrs(U, n, 1, S, ipiv, unit1, n);
                int[] ipivT = new int[n];
                var STbefore = (double[])ST.Clone();
                var c = new double[n];
                Factorise.dsmxmulv(n, Sbefore, unit1, c);
                var backT = Factorise.dsptrf(L, n, ST, ipivT);
                Factorise.dsptrs(L, n, 1, ST, ipivT, unit1T, n);

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
                var n = 500;
                var tdata = 300;
                char[] way = { 'L' };
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
                var back = (way[0] == 'U') ? Factorise.dsptrf(way, n, M, piv) : Factorise.dsptrf(way, n, MT, piv);
                var r = new double[n * n];
                for (int i = 0; i < n; ++i) r[i * n + i] = 1;
                Console.WriteLine($"{r[0]} {r[1]} {r[2]}");
                Console.WriteLine($"{r[3]} {r[4]} {r[5]}");
                Console.WriteLine($"{r[6]} {r[7]} {r[8]}");
                var whichroot = 0;
                var rI = new double[n * n];
                for (int i = 0; i < n; ++i) rI[i * n + i] = 1;
                var symback = (way[0] == 'U') ? Factorise.dsptrs(way, n, n, M, piv, rI, n, 0, 0, 0, whichroot) : Factorise.dsptrs(way, n, n, MT, piv, rI, n, 0, 0, 0, whichroot);
                whichroot = 1;
                symback = (way[0] == 'U') ? Factorise.dsptrs(way, n, n, M, piv, r, n, 0, 0, 0, whichroot) : Factorise.dsptrs(way, n, n, MT, piv, r, n, 0, 0, 0, whichroot);
                if (symback != -10)
                {
                    if (whichroot != 0 && way[0] == 'L') Factorise.dmx_transpose(n, n, r, r);
                    Console.WriteLine($"{r[0]} {r[1]} {r[2]}");
                    Console.WriteLine($"{r[3]} {r[4]} {r[5]}");
                    Console.WriteLine($"{r[6]} {r[7]} {r[8]}");
                }
                var rr2 = new double[n * (n + 1) / 2];
                double[] lower = new double[n * n];

                Factorise.dmx_transpose(n, n, r, lower);
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; j++, ij++)
                    {
                        rr2[ij] = way[0] == 'L' ? BlasLike.ddotvec(n, lower, lower, i * n, j * n) : BlasLike.ddotvec(n, r, r, i * n, j * n);
                    }
                }
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
                var back = Factorise.dsptrf(way, n, M, piv);
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
                    for (int j = 0; j < n; ++j)
                    {
                        for (int k = i + 1; k < n; ++k)
                        {
                            FM[i * n + j] += M[k * (k + 1) / 2 + i] * FM[k * n + j];
                        }
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
                lic = new Byte[3];
                lic[0] = 110;
                lic[1] = 111;
                lic[2] = 112;
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


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
                Test that upper and lower of the solver are working, but this shows that the
                working is not identical!
                */
                var n = 4;
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
                        ST[i * n - i * (i - 1) / 2 + j - i]=S[j * (j + 1) / 2 + i] = cov[j * (j + 1) / 2 + i];
                    }
                }

                var unit1 = new double[n];
                //           for (var i = 0; i < n; ++i) unit1[i] = 1;
                unit1[0] = 1;
                var unit1T = new double[n];
                //           for (var i = 0; i < n; ++i) unit1T[i] = 1;
                unit1T[0] = 1;
                char[] U = { 'U' };
                char[] L = { 'L' };
                var ipiv = new int[n];
                var Sbefore = (double[])S.Clone();
                var back = Factorise.dsptrf(U, n, S, ipiv);
                Factorise.dsptrs(U, n, 1, S, ipiv, unit1, n);
                int[] ipivT = new int[n];
                var STbefore = (double[])ST.Clone();
                var backT = Factorise.dsptrf(L, n, ST, ipivT);
                Factorise.dsptrs(L, n, 1, ST, ipivT, unit1T, n);
                var c = new double[n];
                Factorise.dsmxmulv(n, Sbefore, unit1, c);

                var cT = new double[n];
                Factorise.dsmxmulvT(n, STbefore, unit1T, cT);
                var diff = new double[n];
                Factorise.dsmxmulv(n, Sbefore, c, diff);
                Factorise.dsmxmulvT(n, STbefore, cT, diff);
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


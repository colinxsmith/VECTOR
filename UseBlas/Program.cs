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
                int n = 4;
                double[] aa ={1,2,4,7,
                            3,5,8,
                            6,9,
                            10};
                double[] acopy = (double[])aa.Clone();
                int[] piv = { 1, 2, 3, 4 };
                char[] U = { 'L' };
                int back = 10;
                fixed (double* ap = aa)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    back = Factorise.dsptrf(UP, n, ap, ipiv);
                double[] b = { 1, 2, 3, 4 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
                double okerror = BlasLike.lm_eps * 16;
                double[] c = new double[n];
                Factorise.dsmxmulvT(n, acopy, b, c);
                double[] errorvec = new double[n];
                BlasLike.dsubvec(n, bcopy, c, errorvec);
                double error = Math.Sqrt(BlasLike.ddotvec(n, errorvec, errorvec)) / n;
                Console.WriteLine($"back={back} error={error} {c[0]},{c[1]},{c[2]},{c[3]} ");
            }


            {
                int n = 3;
                double[] aa = { 1,
                            2,3,
                            4,5,6};
                double[] acopy = (double[])aa.Clone();
                int[] piv = { 1, 2, 3 };
                char[] U = { 'U' };
                int back = 10;
                fixed (double* ap = aa)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    back = Factorise.dsptrf(UP, n, ap, ipiv);
                Console.WriteLine($"{back} {piv[0]} {piv[1]} {piv[2]}  {aa[0]} {aa[1]} {aa[2]} {aa[3]} {aa[4]} {aa[5]} ");
                double[] b = { 1, 2, 3 };
                double[] bcopy = (double[])b.Clone();
                fixed (double* ap = aa)
                fixed (double* bp = b)
                fixed (int* ipiv = piv)
                fixed (char* UP = U)
                    Factorise.dsptrs(UP, n, 1, ap, ipiv, bp, n);
                double[] a1 = { 1, 2, 4 };
                double[] a2 = { 2, 3, 5 };
                double[] a3 = { 4, 5, 6 };
                double[] c = new double[n];
                Factorise.dsmxmulv(n, acopy, b, c);
                Console.WriteLine($"back={back} dsmxmulv {c[0]},{c[1]},{c[2]} ");
                double c1 = BlasLike.ddotvec(n, b, a1);
                double c2 = BlasLike.ddotvec(n, b, a2);
                double c3 = BlasLike.ddotvec(n, b, a3);
                double okerror = BlasLike.lm_eps * 64;
                Console.WriteLine($"back={back} {c1},{c2},{c3} ");
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


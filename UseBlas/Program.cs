using System;
using System.Runtime.InteropServices;
using Blas;
using Microsoft.Win32;

namespace UseBlas
{
    class Program
    {
        static void Main(string[] args)
        {
            double a = 4;
            double[] x = { 1, 2, 3 };
            double[] y = { 4, 5, 6 };
            BlasLike.daxpyvec(x.Length, a, x, y);
            for (int i = 0; i < x.Length; ++i)
            {
                Console.WriteLine($"y[{i}]={y.GetValue(i)} {y[i]}");
            }
            bool isWindows = RuntimeInformation.IsOSPlatform(OSPlatform.Windows);
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


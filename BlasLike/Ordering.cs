using System;
using System.Collections.Generic;
namespace Ordering
{
    public class compareitem
    {
        public double x;
        public byte bad;
        public compareitem()
        {

        }
        public compareitem(double x, byte bad)
        {
            this.x = x;
            this.bad = bad;
        }
    }
    class compare : IComparer<compareitem>
    {// I used to think return 0 was correct but got some wrong answers so changed to return 1
        public int Compare(compareitem x, compareitem y)
        {
            if ((y.x - x.x) > 0 && (y.bad == x.bad)) return 1;//smallest have largest index
            else if ((y.bad == x.bad)) return -1;
            else if (y.bad > 0) return 1;
            else return -1;
        }
    }
    class compareAbs : IComparer<compareitem>
    {// I used to think return 0 was correct but got some wrong answers so changed to return 1
        public int Compare(compareitem x, compareitem y)
        {
            if ((Math.Abs(y.x) - Math.Abs(x.x)) > 0 && (y.bad == x.bad)) return 1;//smallest have largest index
            else if ((y.bad == x.bad)) return -1;
            else if (y.bad > 0) return 1;
            else return -1;
        }
    }
    public class Order
    {
        public static void swap<Q>(ref Q a, ref Q b)
        {
            var k2 = a;
            a = b;
            b = k2;
        }
        public static void getorder(int n, double[] x, int[] order = null, byte[] dropbad = null, double init = 0, short sign = 1)
        {
            var cmp = new compare();
            var xx = new compareitem[n];
            for (int i = 0; i < n; ++i)
            {
                xx[i] = new compareitem((x[i] - init) * sign, (dropbad == null) ? (byte)0 : dropbad[i]);
            }
            if (order == null) Array.Sort(xx, cmp);
            else
            {
                var xxx = (compareitem[])xx.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xxx, order, cmp);
            }
        }
        public static void getorderabs(int n, double[] x, int[] order = null, byte[] dropbad = null)
        {
            var cmp = new compareAbs();
            var xx = new compareitem[n];

            for (int i = 0; i < n; ++i)
            {
                xx[i] = new compareitem(x[i], dropbad == null ? (byte)0 : dropbad[i]);
            }
            if (order == null) Array.Sort(xx, cmp);
            else
            {
                var xxx = (compareitem[])xx.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xxx, order, cmp);
            }
        }
        public static void Reorder<T>(int n, int[] order, T[] array)
        {
            if (array == null || order == null) return;
            int i;
            var marked = new bool[n];//This will set marked[i] = false for all i
            int j, k;
            for (i = 0; i < n; ++i) marked[i] = false;
            for (i = 0; i < n; i++)
            {
                if (!marked[i])
                {
                    for (j = i, k = order[j]; k != i; k = order[j = k])
                    {
                        swap(ref array[k], ref array[j]);
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }


        public static void ReorderSymm<T>(int n, int[] order, T[] array, char way = 'U')
        {
            if (array == null || order == null) return;
            var marked = new bool[n];//This will set marked[i] = false for all i
            if (way == 'U')
            {
                for (var i = 0; i < n; i++)
                {
                    if (!marked[i])
                    {
                        for (int j = i, k = order[j]; k != i; k = order[j = k])
                        {
                            var jj = j * (j + 1) / 2;
                            var kk = k * (k + 1) / 2;
                            for (int l = 0; l < n; ++l)
                            {
                                var ll = l * (l + 1) / 2;
                                if (l == j)
                                {
                                    swap(ref array[jj + j], ref array[kk + k]);
                                }
                                else if (l < j && l < k)
                                {
                                    swap(ref array[jj + l], ref array[kk + l]);
                                }
                                else if (l > j && l > k)
                                {
                                    swap(ref array[ll + j], ref array[ll + k]);
                                }
                                else if (l > j && l < k)
                                {
                                    swap(ref array[ll + j], ref array[kk + l]);
                                }
                                else if (l < j && l > k)
                                {
                                    swap(ref array[jj + l], ref array[ll + k]);
                                }
                            }
                            marked[k] = true;
                        }
                        marked[i] = true;
                    }
                }
            }
            else if (way == 'L')
            {
                for (var i = 0; i < n; i++)
                {
                    if (!marked[i])
                    {
                        for (int j = i, k = order[j]; k != i; k = order[j = k])
                        {
                            var jj = j * n - j * (j - 1) / 2;
                            var kk = k * n - k * (k - 1) / 2;
                            for (int l = 0; l < n; ++l)
                            {
                                var ll = l * n - l * (l - 1) / 2;
                                if (l == j)
                                {
                                    swap(ref array[jj], ref array[kk]);
                                }
                                else if (l > j && l > k)
                                {
                                    swap(ref array[jj + l - j], ref array[kk + l - k]);
                                }
                                else if (l < j && l < k)
                                {
                                    swap(ref array[ll + j - l], ref array[ll + k - l]);
                                }
                                else if (l > j && l < k)
                                {
                                    swap(ref array[jj + l - j], ref array[ll + k - l]);
                                }
                                else if (l < j && l > k)
                                {
                                    swap(ref array[jj + l - j], ref array[kk + l - k]);
                                }
                            }
                            marked[k] = true;
                        }
                        marked[i] = true;
                    }
                }
            }
            else return;
        }
        public static void Reorder_gen<T>(int n, int[] order, T[] array, int m = 1, int im = 1)
        {
            if (m == 1 && im == 1) { Reorder(n, order, array); return; }
            if (array == null || order == null) return;
            var marked = new bool[n];//This will set marked[i] = false for all i
            for (var i = 0; i < n; i++)
            {
                if (!marked[i])
                {
                    for (int j = i, k = order[j]; k != i; k = order[j = k])
                    {
                        if (im == 1)
                        {
                            for (int l = 0; l < m; l++)
                            {
                                swap(ref array[k * m + l], ref array[j * m + l]);
                            }
                        }
                        else
                        {
                            for (int l = 0, ll = 0; l < m; l++, ll += im)
                            {
                                swap(ref array[k + ll], ref array[j + ll]);
                            }
                        }
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }
        public static void Display<T>(int n, T[] arr1, string lab = "", char way = 'U')
        {
            if (lab != "") Console.WriteLine(lab);
            if (way == 'U')
            {
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = 0; j <= i; ++j)
                    {
                        Console.Write(arr1[ij++] + " ");
                    }
                    Console.Write("\n");
                }
                Console.Write("\n");
            }
            else if (way == 'L')
            {
                for (int i = 0, ij = 0; i < n; ++i)
                {
                    for (int j = i; j < n; ++j)
                    {
                        Console.Write(arr1[ij++] + " ");
                    }
                    Console.Write("\n");
                }
                Console.Write("\n");
            }
            else return;
        }
        public static void Display<T>(T[] arr1, string lab = "", int m = 1)
        {
            if (lab != "") Console.WriteLine(lab);
            var columns = arr1.Length / m;
            if (columns == 1)
            {
                for (var i = 0; i < arr1.Length; i++)
                {
                    Console.Write(arr1[i] + " ");
                }
            }
            else
            {
                for (var i = 0; i < m; ++i)
                {
                    for (var j = 0; j < columns; ++j)
                    {
                        Console.Write(arr1[i * columns + j] + " ");
                    }
                    Console.Write("\n");
                }
            }
            Console.Write("\n");
        }
        unsafe static void byte_reverse(int n, byte* b)
        {
            //From Robin Becker's C code.
            /*  var bb = new byte[n];      //Testing with safe code
              for (int i = 0; i < n; ++i) bb[i] = b[i];
              byte_reverse(n, bb);
              for (int i = 0; i < n; ++i) b[i] = bb[i];
              return;*/
            var B = (byte*)b;
            var e = B + n;
            n >>= 1;
            while (n-- > 0)
            {
                swap(ref *--e, ref *B++);
            }
        }
        static void byte_reverseA<T>(int n, T[] b, int bstart = 0)
        {
            //Should be identical to byte_reverse<T>
            if (n > 1) Array.Reverse(b, bstart, n);
        }
        static void byte_reverse<T>(int n, T[] b, int bstart = 0)
        {
            //Safe version of Robin Becker's code, but can't use unless there's a safe way
            //to cast double[] to byte[]
            int N = n;
            n >>= 1;
            int ib = 0, ie = 0;
            while (n-- > 0)
            {
                swap(ref b[bstart + N + --ie], ref b[bstart + ib++]);
            }
        }

        public static void bound_reorganise(int f, int n, int nn, int m, double[] bb)
        {
            //	BEFORE bound_reorganise(1,n,temp_stocks,m,L);
            //	AFTER  bound_reorganise(0,n,temp_stocks,m,L);

            if (bb != null && n > nn)
            {
                int[] ns = new int[2];
                ns[1 - f] = n - nn + m;
                ns[f] = m;
                byte_reverse(ns[0], bb, nn);
                byte_reverse(ns[1], bb, nn);
                /*
                fixed (double* b = bb)
                    byte_reverse(ns[0] * sizeof(double), (byte*)(b + nn));
                fixed (double* b = bb)
                    byte_reverse(ns[1] * sizeof(double), (byte*)(b + nn));
                    */
            }
        }
    }
}
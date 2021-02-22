using System;
using System.Collections.Generic;
namespace Ordering
{
    class compare : IComparer<double>
    {
        public int Compare(double x, double y)
        {
            return (y - x) > 0 ? 0 : -1;//smallest have largest index
        }
    }
    class compareAbs : IComparer<double>
    {
        public int Compare(double x, double y)
        {
            return (Math.Abs(y) - Math.Abs(x)) > 0 ? 0 : -1;//smallest have largest index
        }
    }
    public class Order
    {
        public static void getorder(int n, double[] x, int[] order = null)
        {
            var cmp = new compare();
            if (order == null) Array.Sort(x, cmp);
            else
            {
                var xx = (double[])x.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xx, order, cmp);
            }
        }
        public static void getorderabs(int n, double[] x, int[] order = null)
        {
            var cmp = new compareAbs();
            if (order == null) Array.Sort(x, cmp);
            else
            {
                var xx = (double[])x.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xx, order, cmp);
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
                        var aa = array[k];
                        array[k] = array[j];
                        array[j] = aa;
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
                            var aa = array[j * (j + 3) / 2];
                            array[j * (j + 3) / 2] = array[k * (k + 3) / 2];
                            array[k * (k + 3) / 2] = aa;
                            for (int l = 0; l < n; ++l)
                            {
                                if (l < j && l < k)
                                {
                                    aa = array[j * (j + 1) / 2 + l];
                                    array[j * (j + 1) / 2 + l] = array[k * (k + 1) / 2 + l];
                                    array[k * (k + 1) / 2 + l] = aa;
                                }
                                else if (l > j && l > k)
                                {
                                    aa = array[l * (l + 1) / 2 + j];
                                    array[l * (l + 1) / 2 + j] = array[l * (l + 1) / 2 + k];
                                    array[l * (l + 1) / 2 + k] = aa;
                                }
                                else if (l > j && l < k)
                                {
                                    aa = array[l * (l + 1) / 2 + j];
                                    array[l * (l + 1) / 2 + j] = array[k * (k + 1) / 2 + l];
                                    array[k * (k + 1) / 2 + l] = aa;
                                }
                                else if (l < j && l > k)
                                {
                                    aa = array[j * (j + 1) / 2 + l];
                                    array[j * (j + 1) / 2 + l] = array[l * (l + 1) / 2 + k];
                                    array[l * (l + 1) / 2 + k] = aa;
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
                            var aa = array[j * n - j * (j - 1) / 2];
                            array[j * n - j * (j - 1) / 2] = array[k * n - k * (k - 1) / 2];
                            array[k * n - k * (k - 1) / 2] = aa;
                            marked[k] = true;
                            for (int l = 0; l < n; ++l)
                            {
                                if (l > j && l > k)
                                {
                                    aa = array[j * n - j * (j - 1) / 2 + l - j];
                                    array[j * n - j * (j - 1) / 2 + l - j] = array[k * n - k * (k - 1) / 2 + l - k];
                                    array[k * n - k * (k - 1) / 2 + l - k] = aa;
                                }
                                else if (l < j && l < k)
                                {
                                    aa = array[l * n - l * (l - 1) / 2 + j - l];
                                    array[l * n - l * (l - 1) / 2 + j - l] = array[l * n - l * (l - 1) / 2 + k - l];
                                    array[l * n - l * (l - 1) / 2 + k - l] = aa;
                                }
                                else if (l > j && l < k)
                                {
                                    aa = array[j * n - j * (j - 1) / 2 + l - j];
                                    array[j * n - j * (j - 1) / 2 + l - j] = array[l * n - l * (l - 1) / 2 + k - l];
                                    array[l * n - l * (l - 1) / 2 + k - l] = aa;
                                }
                                else if (l < j && l > k)
                                {
                                    aa = array[j * n - j * (j - 1) / 2 + l - j];
                                    array[j * n - j * (j - 1) / 2 + l - j] = array[k * n - k * (k - 1) / 2 + l - k];
                                    array[k * n - k * (k - 1) / 2 + l - k] = aa;
                                }
                            }
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
                                var aa = array[k * m + l];
                                array[k * m + l] = array[j * m + l];
                                array[j * m + l] = aa;
                            }
                        }
                        else
                        {
                            for (int l = 0, ll = 0; l < m; l++, ll += im)
                            {
                                var aa = array[k + ll];
                                array[k + ll] = array[j + ll];
                                array[j + ll] = aa;
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
        public unsafe static void byte_reverse(int n, byte* b)
        {
            byte* B = (byte*)b;
            byte* e = B + n;
            byte t;

            n >>= 1;
            while (n-- > 0)
            {
                t = *--e;
                *e = *B;
                *B++ = t;
            }
        }
        public unsafe static void bound_reorganise(int f, int n, int nn, int m, double[] bb)
        {
            //	BEFORE bound_reorganise(1,n,temp_stocks,m,L);
            //	AFTER  bound_reorganise(0,n,temp_stocks,m,L);

            if (bb != null && n > nn)
            {
                int[] ns = new int[2];
                ns[1 - f] = n - nn + m;
                ns[f] = m;
                fixed (double* b = bb)
                    byte_reverse(ns[0] * sizeof(double), (byte*)(b+nn));
                fixed (double* b = bb)
                    byte_reverse(ns[1] * sizeof(double), (byte*)(b+nn));
            }
        }
    }
}
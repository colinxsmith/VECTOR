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
    {
        public int Compare(compareitem x, compareitem y)
        {
            if ((y.x - x.x) > 0 && (y.bad == x.bad)) return 1;//smallest have largest index
            else if ((y.x - x.x) == 0 && (y.bad == x.bad)) return 0;//Should it be 1 or 0?
            else if ((y.bad == x.bad)) return -1;
            else if (y.bad > 0) return 1;
            else return -1;
        }
    }
    class compareAbs : IComparer<compareitem>
    {
        public int Compare(compareitem x, compareitem y)
        {
            if ((Math.Abs(y.x) - Math.Abs(x.x)) > 0 && (y.bad == x.bad)) return 1;//smallest have largest index
            else if ((Math.Abs(y.x) - Math.Abs(x.x)) == 0 && (y.bad == x.bad)) return 0;//Should it be 1 or 0?
            else if ((y.bad == x.bad)) return -1;
            else if (y.bad > 0) return 1;
            else return -1;
        }
    }
    public class Order
    {///<summary>Exchange the values of two references</summary>
        public static void swap<Q>(ref Q a, ref Q b)
        {
            var k2 = a;
            a = b;
            b = k2;
        }
        public static void getorder(int n, double[] x, int[] order = null, byte[] dropbad = null, double init = 0, short sign = 1, int xstart = 0)
        {
            var cmp = new compare();
            var xx = new compareitem[n];
            for (int i = 0; i < n; ++i)
            {
                xx[i] = new compareitem((x[i + xstart] - init) * sign, (dropbad == null) ? (byte)0 : dropbad[i]);
            }
            if (order == null) Array.Sort(xx, cmp);
            else
            {
                var xxx = (compareitem[])xx.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xxx, order, cmp);
            }
        }
        public static void getorderabs(int n, double[] x, int[] order = null, byte[] dropbad = null, int xstart = 0)
        {
            var cmp = new compareAbs();
            var xx = new compareitem[n];

            for (int i = 0; i < n; ++i)
            {
                xx[i] = new compareitem(x[i + xstart], dropbad == null ? (byte)0 : dropbad[i]);
            }
            if (order == null) Array.Sort(xx, cmp);
            else
            {
                var xxx = (compareitem[])xx.Clone();
                for (int i = 0; i < n; ++i) order[i] = i;
                Array.Sort(xxx, order, cmp);
            }
        }
        ///<summary>Re-order an array</summary>
        ///<param name="n"> length of array</param>
        ///<param name="order"> integer array defining the re-order</param>
        ///<param name="array"> array to re-order</param>
        ///<param name="astart"> start at index astart in array</param>
        public static void Reorder<T>(int n, int[] order, T[] array, int astart = 0)
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
                        swap(ref array[k + astart], ref array[j + astart]);
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }


        ///<summary>Re-order a symmetric matrix stored as an array</summary>
        ///<param name="n"> order of matrix, length is n*(n+1)/2</param>
        ///<param name="order"> integer array defining the re-order</param>
        ///<param name="array"> array to re-order</param>
        ///<param name="way"> either U or L for upper or lower triangle as defined by LAPACK</param>
        ///<param name="astart"> start at index astart in array</param>
        public static void ReorderSymm<T>(int n, int[] order, T[] array, char way = 'U', int astart = 0)
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
                                    swap(ref array[jj + j + astart], ref array[kk + k + astart]);
                                }
                                else if (l < j && l < k)
                                {
                                    swap(ref array[jj + l + astart], ref array[kk + l + astart]);
                                }
                                else if (l > j && l > k)
                                {
                                    swap(ref array[ll + j + astart], ref array[ll + k + astart]);
                                }
                                else if (l > j && l < k)
                                {
                                    swap(ref array[ll + j + astart], ref array[kk + l + astart]);
                                }
                                else if (l < j && l > k)
                                {
                                    swap(ref array[jj + l + astart], ref array[ll + k + astart]);
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
                                    swap(ref array[jj + astart], ref array[kk + astart]);
                                }
                                else if (l > j && l > k)
                                {
                                    swap(ref array[jj + l - j + astart], ref array[kk + l - k + astart]);
                                }
                                else if (l < j && l < k)
                                {
                                    swap(ref array[ll + j - l + astart], ref array[ll + k - l + astart]);
                                }
                                else if (l > j && l < k)
                                {
                                    swap(ref array[jj + l - j + astart], ref array[ll + k - l + astart]);
                                }
                                else if (l < j && l > k)
                                {
                                    swap(ref array[jj + l - j + astart], ref array[kk + l - k + astart]);
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
        ///<summary>Re-order an n by m matrix stored as a one variable array ai,j = a[i*m+j]</summary>
        ///<param name="n"> first dimension</param>
        ///<param name="m"> second dimension</param>
        ///<param name="order"> integer array defining the re-order</param>
        ///<param name="array"> 2d array as a single index array to re-order</param>
        ///<param name="columns"> transpose array if true (use this for constraint array A)</param>
        ///<param name="im"> jump if columns is false, otherwise 1</param>
        ///<param name="astart"> start at index astart in array</param>
        public static void Reorder_gen<T>(int n, int[] order, T[] array, int m = 1, int im = 1, bool columns = false, int astart = 0)
        {
            if (m == 1 && im == 1) { Reorder(n, order, array, astart); return; }
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
                                if (columns) swap(ref array[k * m + l + astart], ref array[j * m + l + astart]);
                                else swap(ref array[k + l + astart], ref array[j + l + astart]);
                            }
                        }
                        else
                        {
                            for (int l = 0; l < m; l++)
                            {
                                if (columns) swap(ref array[k * m + l + astart], ref array[j * m + l + astart]);
                                else swap(ref array[k + l * im + astart], ref array[j + l * im + astart]);
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
        ///<summary>
        ///Translation to c# from c of Robin Becker's subroutine for re-arranging the lower and
        ///upper bounds arrays. The scenario is that the optimiser expects the bounds to be
        ///[stock1,stock2,.......stockn,constraint1,constraint2,.....constraintm].
        ///When dropping for eg. basket constraint where basket size is nn we use this to re-order the
        ///upper and lower bounds so that the fixed stocks appear at the end of the array
        ///[stock1,stock2,.......stocknn,constraint1,constraint2,.....constraintm,stockn,stockn-1,....stocknn+1].
        ///IMPORTANT note how the order of the fixed stocks at the end is reversed, this doesn't affect the optimisation,
        ///but will give nightmares if the wrong order is used accidentally!
        ///</summary>
        ///<param name="f">organise forwards if 1, backwards if 0</param>
        ///<param name="n">number of true portfolio assets</param>
        ///<param name="nn">number of portfolio assets in optimisation</param>
        ///<param name="m">number of constraints</param>
        ///<param name="bb">array to re-order either upper or lower bounds</param>
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
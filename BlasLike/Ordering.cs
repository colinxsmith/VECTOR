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
            var marked = new bool[n];
            int j, k;
            for (i = 0; i < n; ++i) marked[i] = false;
            for (i = 0; i < n; i++)
            {
                if (!marked[i])
                {
                    for (j = i, k = order[j]; k != i; k = order[j = k])
                    {
                        T aa = array[k];
                        array[k] = array[j];
                        array[j] = aa;
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }

        public static void Reorder_gen<T>(int n, int[] order, T[] array, int m = 1, int im = 1)
        {
            if (m == 1 && im == 1) { Reorder(n, order, array); return; }
            if (array == null || order == null) return;
            var marked = new bool[n];
            for (var i = 0; i < n; i++)
            {
                if (!marked[i])
                {
                    for (int j = i, k = order[j]; k != i; k = order[j = k])
                    {
                        for (int l = 0, ll = 0; l < m; l++, ll += im)
                        {
                            var aa = array[k + ll];
                            array[k + ll] = array[j + ll];
                            array[j + ll] = aa;
                        }
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
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
    }
}
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
        public static void Reorder(int n, int[] order, double[] array)
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
                        double aa = array[k];
                        array[k] = array[j];
                        array[j] = aa;
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }

        public static void Reorder_gen(int n, int[] order, double[] array, int m = 1)
        {
            if (m == 1){ Reorder(n, order, array);return;}
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
                        for (int l = 0; l < m; ++l)
                        {
                            double aa = array[k + l];
                            array[k + l] = array[j + l];
                            array[j + l] = aa;
                        }
                        marked[k] = true;
                    }
                    marked[i] = true;
                }
            }
        }
        public static void Display(double[] arr1, string lab = "")
        {
            if (lab != "") Console.WriteLine(lab);
            for (int i = 0; i < arr1.Length; i++)
            {
                Console.Write(arr1[i] + " ");
            }
            Console.Write("\n");
        }
        public static void Display(int[] arr1, string lab = "")
        {
            if (lab != "") Console.WriteLine(lab);
            for (int i = 0; i < arr1.Length; i++)
            {
                Console.Write(arr1[i] + " ");
            }
            Console.Write("\n");
        }
    }
}
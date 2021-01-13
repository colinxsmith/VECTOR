using System;
using Blas;
namespace UseBlas
{
    class Program
    {
        static void Main(string[] args)
        {
            int n=3;
            double a=4;
            double[]x={1,2,3};
            double[]y={4,5,6};
            BlasLike.daxpyvec(n,a,x,y);
            for(int i=0;i<n;++i){
                Console.WriteLine($"y[{i}]={y[i]}");
            }
        }
    }
}

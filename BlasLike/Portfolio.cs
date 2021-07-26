using System;
using Blas;
using Solver;
using ActiveSet;
using System.Diagnostics;
namespace Portfolio
{
    public class Portfolio
    {
        public Portfolio(string file)
        {
            inFile = file;
            if (inFile != "")
                using (var OptData = new DataFile.InputSomeData())
                {
                    OptData.intFields = "n m";
                    OptData.doubleFields = "gamma alpha Q L U A buy sell initial bench";
                    OptData.stringFields = "names";
                    OptData.Read(inFile);
                    n = OptData.mapInt["n"][0];
                    m = OptData.mapInt["m"][0];
                    gamma = OptData.mapDouble["gamma"][0];
                    alpha = OptData.mapDouble["alpha"];
                    initial = OptData.mapDouble["initial"];
                    bench = OptData.mapDouble["bench"];
                    L = OptData.mapDouble["L"];
                    U = OptData.mapDouble["U"];
                    A = OptData.mapDouble["A"];
                    buy = OptData.mapDouble["buy"];
                    sell = OptData.mapDouble["sell"];
                    Q = OptData.mapDouble["Q"];
                    names = OptData.mapString["names"];
                    if (names != null) Array.Resize(ref names, n);
                    var back = makeQ();
                }
        }
        public string inFile = "";
        public int n;
        public int m;
        public double gamma;
        public double kappa;
        public double[] w = null;
        public double[] initial = null;
        public double[] bench = null;
        public double[] L = null;
        public double[] U = null;
        public double[] buy = null;
        public double[] sell = null;
        public double[] A = null;
        public double[] alpha = null;
        public double[] Q = null;
        public string[] names;
        public virtual int makeQ()
        {
            var nn = n * (n + 1) / 2;
            if (Q.Length == nn)
                return 0;
            else
                return -10;
        }
    }
    public class FPortfolio : Portfolio
    {
        public FPortfolio(string file) : base("")
        {
            inFile = file;
            using (var OptData = new DataFile.InputSomeData())
            {
                OptData.intFields = "n m nfac";
                OptData.doubleFields = "gamma alpha FC FL SV L U A buy sell initial bench";
                OptData.stringFields = "names";
                OptData.Read(inFile);
                n = OptData.mapInt["n"][0];
                m = OptData.mapInt["m"][0];
                nfac = OptData.mapInt["nfac"][0];
                gamma = OptData.mapDouble["gamma"][0];
                alpha = OptData.mapDouble["alpha"];
                initial = OptData.mapDouble["initial"];
                bench = OptData.mapDouble["bench"];
                L = OptData.mapDouble["L"];
                U = OptData.mapDouble["U"];
                A = OptData.mapDouble["A"];
                buy = OptData.mapDouble["buy"];
                sell = OptData.mapDouble["sell"];
                FC = OptData.mapDouble["FC"];
                FL = OptData.mapDouble["FL"];
                SV = OptData.mapDouble["SV"];
                names = OptData.mapString["names"];
                if (names != null) Array.Resize(ref names, n);
                var back = makeQ();
                var ok = ActiveOpt();
            }
        }
        void hessmull(int nn, int nrowh, int ncolh, int j, double[] QQ, double[] x, double[] hx)
        {
            Factorise.FacMul(n, nfac, Q, x, hx);
        }
        public int ActiveOpt()
        {
            var obj = 0.0;
            var iter = 10;
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            if (bench != null)
            {
                Factorise.FacMul(n, nfac, Q, bench, cextra);
                BlasLike.dnegvec(n, cextra);
            }
            BlasLike.daxpyvec(n, -gamma / (1 - gamma), c, cextra);
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var opt = new Optimise();
            opt.h = hessmull;
            return opt.QPopt(n, m, w, L, U, A, cextra, Q, ref obj, ref iter);
        }
        public int nfac;
        public double[] FL = null;
        public double[] SV = null;
        public double[] FC = null;
        public override int makeQ()
        {
            var nn = (nfac + 1) * n;
            Q = new double[nn];
            return Factorise.FMP(n, nfac, FC, SV, FL, Q);
        }
    }
}
using System;
using Blas;
using Solver;
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
                    OptData.intFields = "n m nfac";
                    OptData.doubleFields = "delta kappa gamma alpha Q L U A buy sell initial bench FC FL SV";
                    OptData.stringFields = "names";
                    OptData.Read(inFile);
                    n = OptData.mapInt["n"][0];
                    m = OptData.mapInt["m"][0];
                    gamma = OptData.mapDouble["gamma"][0];
                    kappa = OptData.mapDouble["kappa"][0];
                    delta = OptData.mapDouble["delta"][0];
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
                }
        }
        public virtual void Optimise()
        {
            var back = makeQ();
            var ip = InteriorOpt();
            Console.WriteLine($"Variance from IP:\t\t{Variance(w)}");
            var ok = ActiveOpt();
            ActiveSet.Optimise.printV("w from Active Set", w);
            Console.WriteLine($"Variance from Active Set:\t\t{Variance(w)}");
        }
        public string inFile = "";
        public int n;
        public int m;
        public double gamma;
        public double kappa;
        public double delta;
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
        public virtual void hessmull(int nn, int nrowh, int ncolh, int j, double[] QQ, double[] x, double[] hx)
        {
            Factorise.CovMul(n, Q, x, hx);
        }
        public virtual void hessmull(int nn, double[] QQ, double[] x, double[] hx)
        {
            Factorise.CovMul(n, Q, x, hx);
        }
        public double Variance(double[] w)
        {
            var Qx = new double[n];
            hessmull(n, Q, w, Qx);
            return BlasLike.ddotvec(n, w, Qx);
        }
        public int ActiveOpt()
        {
            var obj = 0.0;
            var iter = 10;
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            var opt = new ActiveSet.Optimise();
            opt.h = hessmull;
            if (bench != null)
            {
                hessmull(n, 0, 0, 0, Q, bench, cextra);
                BlasLike.dnegvec(n, cextra);
            }
            BlasLike.daxpyvec(n, -gamma / (1 - gamma), c, cextra);
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var back = opt.QPopt(n, m, w, L, U, A, cextra, Q, ref obj, ref iter);
            Console.WriteLine($"objective:\t\t{obj:F8}; {iter} iterations");
            return back;
        }

        public int InteriorOpt()
        {
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            var CTEST = new double[n];
            var b = new double[m];
            BlasLike.dcopyvec(m, U, b, n);
            var sign = new int[2 * n];
            var UL = new double[n];
            if (bench != null)
            {
                hessmull(n, Q, bench, cextra);
                BlasLike.dnegvec(n, cextra);
            }
            BlasLike.daxpyvec(n, -gamma / (1 - gamma), c, cextra);
            var zcount = 0;
            var signcount = 0;
            for (var i = 0; i < n; ++i)
            {
                //We do an LP to test for primal feasibility
                //but it's possible that c on its own gives an unbounded LP
                //i.e the dual is infeasible, so it's necessary to mess
                //about with c to get a dual feasible LP to test the primal
                //constraints.
                sign[i + n] = 1;
                if (L[i] == 0)
                {
                    sign[i] = 1;
                    UL[i] = L[i];
                    signcount++;
                }
                else if (U[i] == 0)
                {
                    sign[i] = -1;
                    UL[i] = U[i];
                    signcount++;
                }
                if (UL[i] == 0) zcount++;
                CTEST[i] = sign[i] * Math.Abs(cextra[i]);
            }
            Array.Resize(ref CTEST, 2 * n);
            Array.Resize(ref cextra, 2 * n);
            var AA = new double[2 * n * (n + m)];
            for (var con = 0; con < m; ++con)
            {
                BlasLike.dcopy(n, A, m, AA, n + m, con, con);
            }
            for (int i = 0, astart = m; i < n; ++i, astart++)
            {
                AA[astart + i * (n + m)] = 1;
                AA[astart + (i + n) * (n + m)] = 1;
            }
            if (zcount == n) UL = (double[])L.Clone();
            Array.Resize(ref UL, n);
            Array.Resize(ref UL, n*2);
            var bb = (double[])b.Clone();
            Array.Resize(ref bb, m + n);
            for (var i = 0; i < n; ++i)
            {
                bb[m + i] = U[i];
            }
            if (signcount != n) { CTEST = cextra; sign = null; }
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var HH = new double[n * (n + 1) / 2];
            Factorise.Fac2Cov(n, (int)(Q.Length / n) - 1, Q, HH);
            var ww = (double[])w.Clone();
            Array.Resize(ref ww, n * 2);
            // First do a homogenous LP do decide if the problem is feasible.
            // (homogenous QP only works if we're very lucky)
            var IOPT = new InteriorPoint.Optimise(n * 2, m + n, ww, AA, bb, CTEST);
            IOPT.alphamin = 1e-4;
            var back = IOPT.Opt("QP", null, null, true, UL, sign);
            if (back == 6) Console.WriteLine("INFEASIBLE");
            else
            {
                IOPT = new InteriorPoint.Optimise(n * 2, m + n, ww, AA, bb, cextra, n, HH);
                IOPT.h = hessmull;
                var testmul = new double[n];
                hessmull(n, Q, w, testmul);
                Console.WriteLine(BlasLike.ddotvec(n, ww, testmul));
                var kk = new Portfolio("");
                kk.hessmull(n, HH, ww, testmul);
                Console.WriteLine(BlasLike.ddotvec(n, ww, testmul));
                IOPT.alphamin = 1e-8;
                back = IOPT.Opt("QP", null, null, false, UL, sign);
                BlasLike.dcopyvec(n, ww, w);
            }
            return back;
        }

    }
    public class FPortfolio : Portfolio
    {
        public FPortfolio(string file) : base(file)
        {
            inFile = file;
            using (var OptData = new DataFile.InputSomeData())
            {
                OptData.intFields = "n m nfac";
                OptData.doubleFields = "delta kappa gamma alpha Q L U A buy sell initial bench FC FL SV";
                OptData.stringFields = "names";
                OptData.Read(inFile);
                nfac = OptData.mapInt["nfac"][0];
                FC = OptData.mapDouble["FC"];
                FL = OptData.mapDouble["FL"];
                SV = OptData.mapDouble["SV"];
                names = OptData.mapString["names"];
                if (names != null) Array.Resize(ref names, n);
            }
        }
        public override void hessmull(int nn, int nrowh, int ncolh, int j, double[] QQ, double[] x, double[] hx)
        {
            Factorise.FacMul(n, nfac, Q, x, hx);
        }
        public override void hessmull(int nn, double[] QQ, double[] x, double[] hx)
        {
            Factorise.FacMul(n, nfac, Q, x, hx);
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
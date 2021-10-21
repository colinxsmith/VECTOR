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
            if (true)
            {
                var AA = new double[n * (m + 1)];
                var LL = new double[n + m + 1];
                var UU = new double[n + m + 1];
                for (var i = 0; i < n + m; i++)
                {
                    if (i < n)
                    {
                        AA[i * (m + 1)] = 1;
                        AA[i * (m + 1) + 1] = alpha[i];
                    }
                    LL[i] = L[i];
                    UU[i] = U[i];
                }
                LL[n] = 1;
                UU[n] = 1;
                LL[n + 1] = -1e-2;//BlasLike.lm_max;
                UU[n + 1] = 0;
                m++;
                L = LL;
                U = UU;
                A = AA;
            }
            var Aw = new double[m];
            var ok = ActiveOpt();
            ActiveSet.Optimise.printV("w from Active Set", w);
            Console.WriteLine($"Variance from Active Set:\t\t{Variance(w)}");
            Factorise.dmxmulv(m, n, A, w, Aw);
            ActiveSet.Optimise.printV("Constraints", Aw);
            var ip = InteriorOpt();
            Console.WriteLine($"Variance from IP:\t\t{Variance(w)}");
            Factorise.dmxmulv(m, n, A, w, Aw);
            ActiveSet.Optimise.printV("Constraints", Aw);
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
            var dolarge = 0;
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            var CTEST = new double[n];
            var slackb = 0;
            var slackL = 0;
            var slackU = 0;
            for (var i = 0; i < m; ++i)
            {
                if (L[i + n] == -BlasLike.lm_max) slackU++;
                else if (U[i + n] == BlasLike.lm_max) slackL++;
                else if (U[i + n] != L[i + n]) slackb++;
            }
            var totalConstraintslack = slackU + slackL + 2 * slackb + slackL + slackU;
            var slackToConstraintL = slackL > 0 ? new int[slackL] : null;
            var slackToConstraintU = slackU > 0 ? new int[slackU] : null;
            var slackToConstraintBOTH = slackb > 0 ? new int[slackb] : null;
            var b = new double[m + slackb];
            BlasLike.dcopyvec(m, L, b, n);
            for (int i = 0, slack = 0, slackLL = 0, slackUU = 0; i < m; ++i)
            {
                if (U[i + n] != L[i + n])
                {
                    if (L[i + n] == -BlasLike.lm_max)
                    {
                        slackToConstraintU[slackUU++] = i;
                        b[i] = U[i + n];
                    }
                    else if (U[i + n] == BlasLike.lm_max)
                    {
                        slackToConstraintL[slackLL++] = i;
                    }
                    else
                    {
                        slackToConstraintBOTH[slack] = i;
                        b[slack++ + m] = U[i + n];
                    }
                }
            }
            var sign = new int[n * dolarge + n + totalConstraintslack];
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
                if (L[i] >= 0)
                {
                    sign[i] = 1;
                    if (dolarge == 1) sign[i + n] = 1;
                    UL[i] = L[i];
                    signcount++;
                }
                else if (U[i] <= 0)
                {
                    sign[i] = -1;
                    if (dolarge == 1) sign[i + n] = -1;
                    UL[i] = U[i];
                    signcount++;
                }
                else sign[i] = 1;
                if (UL[i] == 0) zcount++;
                CTEST[i] = sign[i] * Math.Abs(cextra[i]);
            }
            for (var i = 0; i < totalConstraintslack; ++i)
                sign[i + dolarge * n + n] = 1;
            Array.Resize(ref CTEST, dolarge * n + n + totalConstraintslack);
            Array.Resize(ref cextra, dolarge * n + n + totalConstraintslack);
            /*           var AA = new double[(dolarge * n + n + 2 * slackb) * (dolarge * n + m + slackb)];
                       for (int con = 0, slack = 0; con < m; ++con)
                       {
                           BlasLike.dcopy(n, A, m, AA, dolarge * n + m + slackb, con, con);
                           if (U[con + n] != L[con + n])
                           {
                               AA[con + (dolarge * n + n + slack + slackb) * (dolarge * n + m + slackb)] = -1;
                               AA[dolarge * n + m + slack + (dolarge * n + n + slack) * (dolarge * n + m + slackb)] = 1;
                               BlasLike.dcopy(n, A, m, AA, dolarge * n + m + slackb, con, dolarge * n + m + slack++);
                           }
                       }
                       for (int i = 0, astart = m; i < dolarge * n; ++i, astart++)
                       {
                           AA[astart + i * (dolarge * n + m + slackb)] = 1;
                           AA[astart + (i + n) * (dolarge * n + m + slackb)] = 1;
                       }*/
            if (zcount == n) UL = (double[])L.Clone();
            for (var i = 0; i < n; ++i) if (sign[i] == -1) UL[i] = U[i];
            Array.Resize(ref UL, n);
            Array.Resize(ref UL, dolarge * n + n + m + slackb);
            if (InteriorPoint.Optimise.lInfinity(UL) == 0) UL = null;
            var bb = (double[])b.Clone();
            Array.Resize(ref bb, m + dolarge * n + slackb);
            for (var i = 0; i < dolarge * n; ++i)
            {
                bb[m + i] = sign[i] == 1 ? U[i] : L[i];
            }
            for (var i = 0; i < slackb; ++i)
            {
                bb[m + i + dolarge * n] = b[m + i];
            }
            if (signcount != n) { CTEST = cextra; sign = null; }
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var HH = new double[n * (n + 1) / 2];
            Factorise.Fac2Cov(n, (int)(Q.Length / n) - 1, Q, HH);
            var ww = (double[])w.Clone();
            Array.Resize(ref ww, dolarge * n + n + totalConstraintslack);
            // First do a homogenous LP do decide if the problem is feasible.
            // (homogenous QP only works if we're very lucky)
            var IOPT = new InteriorPoint.Optimise(dolarge * n + n + totalConstraintslack, m + dolarge * n + slackb, ww, null, bb, cextra);
            IOPT.alphamin = 1e-10;
            IOPT.baseA = A;//We only need to pass the constraints without slack variables AA just use for testing
            IOPT.basen = n;
            IOPT.bases = dolarge * n;
            IOPT.basem = m;
            IOPT.slackToConstraintBOTH = slackToConstraintBOTH;
            IOPT.slackToConstraintL = slackToConstraintL;
            IOPT.slackToConstraintU = slackToConstraintU;
            var back =
            IOPT.Opt("QP", null, null, true, UL, sign);
            if (back < -10) Console.WriteLine($"Failed -- too many iterations");
            if (back < 0) Console.WriteLine($"Normal Matrix became ill-conditioned");
            if (back == 6) Console.WriteLine("INFEASIBLE");
            else
            {
                IOPT = new InteriorPoint.Optimise(dolarge * n + n + totalConstraintslack, m + dolarge * n + slackb, ww, null, bb, cextra, n, HH);
                IOPT.h = hessmull;
                IOPT.baseA = A;
                IOPT.basen = n;
                IOPT.bases = n * dolarge;
                IOPT.basem = m;
                IOPT.slackToConstraintBOTH = slackToConstraintBOTH;
                IOPT.slackToConstraintL = slackToConstraintL;
                IOPT.slackToConstraintU = slackToConstraintU;
                var testmul = new double[n];
                hessmull(n, Q, w, testmul);
                Console.WriteLine(BlasLike.ddotvec(n, ww, testmul));
                var kk = new Portfolio("");
                kk.hessmull(n, HH, ww, testmul);
                Console.WriteLine(BlasLike.ddotvec(n, ww, testmul));
                IOPT.alphamin = 1e-8;
                back = IOPT.Opt("QP", null, null, false, UL, sign);
                BlasLike.dcopyvec(n, ww, w);
                if (back < -10) Console.WriteLine($"Failed -- too many iterations");
                else if (back < 0) Console.WriteLine($"Normal Matrix became ill-conditioned");
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

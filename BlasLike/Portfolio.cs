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
        public void BuySellSetup(int n, int m, double[] A, double[] L, double[] U, double[] c, double[] Q, double[] initial, string[] names, bool useIP = true)
        {
            var delta = 0.98;
            var twoside = true; //twoside = false means treat sell side only
            var N = n + n;
            var M = m + n + (delta < 1.0 ? 1 : 0);
            if (twoside)
            {
                N += n;
                M += n;
            }
            var CC = new double[N];
            var AA = new double[N * M];
            var LL = new double[N + M];
            var UU = new double[N + M];
            var WW = new double[N];
            //Variables
            BlasLike.dcopyvec(n, L, LL);
            BlasLike.dcopyvec(n, U, UU);
            BlasLike.dsetvec(n, 0, LL, n);
            BlasLike.dsetvec(n, BlasLike.lm_max, UU, n);
            if (twoside)
            {
                BlasLike.dsetvec(n, 0, LL, n + n);
                BlasLike.dsetvec(n, BlasLike.lm_max, UU, n + n);
            }
            //Constraints
            BlasLike.dcopyvec(m, L, LL, n, N);
            BlasLike.dcopyvec(m, U, UU, n, N);
            for (var i = 0; i < m; i++)
            {
                BlasLike.dcopy(n, A, m, AA, M, i, i);
            }
            BlasLike.dsccopyvec(n, 1.0, initial, LL, 0, N + m);
            BlasLike.dsetvec(n, BlasLike.lm_max, UU, N + m);
            for (var i = m; i < m + n; ++i)
            {
                BlasLike.dset(1, 1.0, AA, M, i + M * (i - m));
                BlasLike.dset(1, 1.0, AA, M, i + M * (n + i - m));
            }
            if (twoside)
            {
                BlasLike.dsccopyvec(n, 1.0, initial, UU, 0, N + m + n);
                BlasLike.dsetvec(n, -BlasLike.lm_max, LL, N + m + n);
                for (var i = m + n; i < m + n + n; ++i)
                {
                    BlasLike.dset(1, 1.0, AA, M, i + M * (i - m - n));
                    BlasLike.dset(1, -1.0, AA, M, i + M * (n + n + i - m - n));
                }
            }
            BlasLike.dsccopyvec(n, 1, c, CC);
            if (delta < 1.0)
            {
                LL[N + M - 1] = 0;
                if (twoside)
                {
                    BlasLike.dset(n + n, 1.0, AA, M, n + n + m + M * n);
                    UU[N + M - 1] = delta * 2;
                }
                else
                {
                    BlasLike.dset(n, 1.0, AA, M, n + m + M * n);
                    UU[N + M - 1] = delta;
                }
            }
            this.L = LL;
            this.U = UU;
            this.A = AA;
            this.n = N;
            this.m = M;
            this.gamma = 0.5;
            BlasLike.dsetvec(n, 0, CC, n);
            if (twoside) BlasLike.dsetvec(n, 0, CC, n + n);
            this.alpha = CC;
            var back = InteriorOpt(1e-10, WW);
            var turnover = 0.0;
            for (var i = 0; i < n; ++i)
            {
                turnover += Math.Abs(WW[i] - initial[i]);
                var c1 = BlasLike.ddot(N, AA, M, WW, 1, i + m);
                if (twoside)
                {
                    var c2 = BlasLike.ddot(N, AA, M, WW, 1, n + i + m);
                    if (WW[i] <= initial[i]) Console.WriteLine($"{names[i]}\t{(WW[i] - initial[i]):F8}\t{WW[i + n]:F8} {c1:F8} {WW[i + n + n]:F8} {c2:F8} {initial[i]:F8}");
                    else Console.WriteLine($"{names[i]} {(WW[i] - initial[i]):F8}\t{WW[i + n]:F8} {c1:F8} {WW[i + n + n]:F8} {c2:F8} {initial[i]:F8}");
                }
                else
                {
                    if (WW[i] <= initial[i]) Console.WriteLine($"{names[i]}\t{(WW[i] - initial[i]):F8}\t{WW[i + n]:F8} {c1:F8}  {initial[i]:F8}");
                    else Console.WriteLine($"{names[i]} {(WW[i] - initial[i]):F8}\t{WW[i + n]:F8} {c1:F8}  {initial[i]:F8}");
                }
            }
            Console.WriteLine($"Turnover: {turnover * 0.5}");
        }
        public void GainLossSetUp(int n, int tlen, double[] DATA, string[] names, double R, double lambda, bool useIP = true)
        {
            var m = 1;
            var N = n + tlen;
            var M = m + tlen;
            var activeLp = 0;
            double[] ww = new double[n + tlen];
            double[] cc = new double[N];
            double[] LL = new double[N + M];
            double[] UU = new double[N + M];
            double[] AA = new double[N * M];
            double[] bb = new double[M];
            if ((activeLp == 0))
            {
                Q = new double[N * (N + 1) / 2];
                for (int i = 0, ij = 0; i < n; ++i, ij += i)
                {
                    Q[ij + i] = 0;
                }
            }


            double[] alpha = new double[n];
            for (var i = 0; i < n; i++)
            {
                alpha[i] = BlasLike.dsumvec(tlen, DATA, i * tlen) / tlen;
            }
            double maxret = BlasLike.lm_max;
            if (!useIP)
            {
                maxret = InteriorPoint.Optimise.lInfinity(DATA);
                maxret = Math.Max(maxret, R);
            }
            BlasLike.dsccopyvec(n, -1, alpha, cc);
            BlasLike.dsetvec(tlen, lambda / tlen, cc, n);

            BlasLike.dzerovec(n, LL);
            BlasLike.dzerovec(tlen, LL, n);
            BlasLike.dsetvec(m, 1, LL, N);
            BlasLike.dsetvec(tlen, R, LL, N + m);

            BlasLike.dsetvec(n, useIP ? BlasLike.lm_max : 1, UU);
            BlasLike.dsetvec(tlen, useIP ? BlasLike.lm_max : maxret, UU, n);
            BlasLike.dsetvec(m, 1, UU, N);
            BlasLike.dsetvec(tlen, useIP ? BlasLike.lm_max : maxret, UU, N + m);

            for (var i = 0; i < m; ++i)
            {
                BlasLike.dset(n, 1, AA, M, i);
            }
            for (var i = m; i < M; ++i)
            {
                BlasLike.dset(1, 1, AA, M, i + M * (n + i - m));
                BlasLike.dcopy(n, DATA, tlen, AA, M, i - m, i);
            }
            this.L = LL;
            this.U = UU;
            this.A = AA;
            this.n = N;
            this.m = M;
            this.gamma = 0.5;
            BlasLike.dnegvec(cc.Length, cc);
            this.alpha = cc;
            double alphamax = 0, alphamin = 0;
            BlasLike.dxminmax(this.alpha.Length, this.alpha, 1, ref alphamax, ref alphamin);
            Console.WriteLine($"c range ({alphamax},{alphamin}) ratio {(alphamin / alphamax):E8}");
            if (useIP)
            {
                var back = InteriorOpt(5e-11, ww);
                var loss = 0.0;
                var gain = 0.0;
                for (var i = 0; i < tlen; ++i)
                {
                    var GL = BlasLike.ddot(n, DATA, tlen, ww, 1, i) - R;
                    gain += Math.Max(0, GL);
                    loss += -Math.Min(0, GL);
                }
                var lossV = BlasLike.dsumvec(tlen, ww, n);
                Console.WriteLine($"Total GAIN = {gain:F8}");
                Console.WriteLine($"Total LOSS = \t\t\t\t{loss:F8}");
                Console.WriteLine($"Total LOSS (check from opt variables) = {lossV:F8}");
                for (var i = 0; i < n; ++i)
                {
                    Console.WriteLine($"{names[i]}\t{ww[i]:F8}");
                }
                for (var i = 0; i < M; ++i)
                {
                    var constraint = BlasLike.ddot(N, AA, M, ww, 1, i);
                    if (i < m) Console.WriteLine($"Constraint {i} = {constraint:F8}");
                    else Console.WriteLine($"Constraint {i} = {constraint:F8} G/L var {ww[n + i - 1]:F8}");
                }
            }
            else
            {
                BlasLike.dsetvec(n, 1.0 / n, ww);
                BlasLike.dsetvec(tlen, 1.0 / tlen, ww, n);
                var back = ActiveOpt(activeLp, ww);
                var gain = 0.0;
                var loss = 0.0;
                for (var i = 0; i < tlen; ++i)
                {
                    var GL = BlasLike.ddot(n, DATA, tlen, ww, 1, i) - R;
                    gain += Math.Max(0, GL);
                    loss += -Math.Min(0, GL);
                }
                var lossV = BlasLike.dsumvec(tlen, ww, n);
                Console.WriteLine($"Total GAIN = {gain:F8}");
                Console.WriteLine($"Total LOSS = \t\t\t\t{loss:F8}");
                Console.WriteLine($"Total LOSS (check from opt variables) = {lossV:F8}");
                for (var i = 0; i < n; ++i)
                {
                    Console.WriteLine($"{names[i]}\t{ww[i]:F8}");
                }
                for (var i = 0; i < M; ++i)
                {
                    var constraint = BlasLike.ddot(N, AA, M, ww, 1, i);
                    if (i < m) Console.WriteLine($"Constraint {i} = {constraint:F8}");
                    else Console.WriteLine($"Constraint {i} = {constraint:F8}  {UU[N + i]} G/L var {ww[n + i - 1]:F8}  {UU[n + i - 1]}");
                }
            }

        }
        public virtual void Optimise()
        {
            var back = makeQ();
            BlasLike.dscalvec(Q.Length, 1e5, Q);
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
                LL[n] = 0.9;
                UU[n] = 1;
                LL[n + 1] = -BlasLike.lm_max;
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
        public int ActiveOpt(int lp = 0, double[] www = null)
        {
            var obj = 0.0;
            var iter = 10;
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            var opt = new ActiveSet.Optimise();
            if (lp == 0) opt.h = hessmull;
            if (bench != null)
            {
                hessmull(n, 0, 0, 0, Q, bench, cextra);
                BlasLike.dnegvec(n, cextra);
            }
            BlasLike.daxpyvec(n, -gamma / (1 - gamma), c, cextra);
            if (www == null)
            {
                w = new double[n];
                BlasLike.dsetvec(n, 1.0 / n, w);
            }
            else
                w = www;
            var back = opt.QPopt(n, m, w, L, U, A, cextra, Q, ref obj, ref iter, lp);
            Console.WriteLine($"objective:\t\t{obj}; {iter} iterations");
            return back;
        }

        public int InteriorOpt(double conv = 1e-16, double[] wback = null)
        {
            for (var i = 0; i < n; ++i)
            {
                if (U[i] == 1 && L[i] == 0) U[i] = BlasLike.lm_max;
                if (L[i] == -1 && U[i] == 0) L[i] = -BlasLike.lm_max;
            }
            var slacklarge = 0;
            var c = (double[])alpha.Clone();
            var cextra = new double[n];
            var CTEST = new double[n];
            var slackb = 0;
            var slackL = 0;
            var slackU = 0;
            for (var i = 0; i < n; ++i)
            {
                if (U[i] != BlasLike.lm_max && U[i] != 0) slacklarge++;
                else if (L[i] != -BlasLike.lm_max && L[i] != 0) slacklarge++;
            }
            for (var i = 0; i < m; ++i)
            {
                if (L[i + n] == -BlasLike.lm_max) slackU++;
                else if (U[i + n] == BlasLike.lm_max) slackL++;
                else if (U[i + n] != L[i + n]) slackb++;
            }
            var totalConstraintslack = slackU + slackL + 2 * slackb; // + slackL + slackU;
            var slackToConstraintL = slackL > 0 ? new int[slackL] : null;
            var slackToConstraintU = slackU > 0 ? new int[slackU] : null;
            var slackToConstraintBOTH = slackb > 0 ? new int[slackb] : null;
            var slacklargeConstraint = slacklarge > 0 ? new int[slacklarge] : null;
            var b = new double[m + slackb];
            BlasLike.dcopyvec(m, L, b, n);
            for (int i = 0, slack = 0; i < n; i++)
            {
                if (U[i] != BlasLike.lm_max && U[i] != 0)
                {
                    slacklargeConstraint[slack++] = i;
                }
                else if (L[i] != -BlasLike.lm_max && L[i] != 0)
                {
                    slacklargeConstraint[slack++] = i;
                }
            }
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
            var sign = new int[slacklarge + n + totalConstraintslack];
            var UL = new double[n];
            if (bench != null)
            {
                hessmull(n, Q, bench, cextra);
                BlasLike.dnegvec(n, cextra);
            }
            BlasLike.daxpyvec(n, -gamma / (1 - gamma), c, cextra);
            var zcount = 0;
            var signfix = false;
            for (int i = 0, slack = 0; i < n; ++i)
            {
                //We do an LP to test for primal feasibility
                //but it's possible that c on its own gives an unbounded LP
                //i.e the dual is infeasible, so it's necessary to mess
                //about with c to get a dual feasible LP to test the primal
                //constraints.
                if (L[i] >= 0)
                {
                    sign[i] = 1;
                    if (U[i] != BlasLike.lm_max && slacklarge > 0 && slack < slacklarge && slacklargeConstraint[slack] == i) sign[slack++ + n] = 1;
                    UL[i] = L[i];
                }
                else if (U[i] <= 0)
                {
                    sign[i] = -1;
                    if (L[i] != -BlasLike.lm_max && slacklarge > 0 && slack < slacklarge && slacklargeConstraint[slack] == i) sign[slack++ + n] = -1;
                    UL[i] = U[i];
                    signfix = true;
                }
                else { sign[i] = 1; UL[i] = L[i]; sign[slack++ + n] = 1; }
                if (UL[i] == 0) zcount++;
                CTEST[i] = sign[i] * Math.Abs(cextra[i]);
            }
            for (var i = 0; i < totalConstraintslack; ++i)
                sign[i + slacklarge + n] = 1;
            Array.Resize(ref CTEST, slacklarge + n + totalConstraintslack);
            Array.Resize(ref cextra, slacklarge + n + totalConstraintslack);
            if (zcount == n) UL = (double[])L.Clone();
            for (var i = 0; i < n; ++i) if (sign[i] == -1) UL[i] = U[i];
            Array.Resize(ref UL, n);
            Array.Resize(ref UL, slacklarge + n + m + slackb);
            if (InteriorPoint.Optimise.lInfinity(UL) == 0) UL = null;
            var bb = (double[])b.Clone();
            Array.Resize(ref bb, m + slacklarge + slackb);
            for (int i = 0; i < slacklarge; ++i)
            {
                var ii = slacklargeConstraint[i];
                bb[m + i] = sign[ii] == 1 ? U[ii] : L[ii];
            }
            for (var i = 0; i < slackb; ++i)
            {
                bb[m + i + slacklarge] = b[m + i];
            }
            if (!signfix) { CTEST = cextra; sign = null; }
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var HH = new double[n * (n + 1) / 2];
            if (Q != null && Q.Length == (n * (n + 1) / 2)) HH = Q;
            else if (Q != null) Factorise.Fac2Cov(n, (int)(Q.Length / n) - 1, Q, HH);
            var ww = (double[])w.Clone();
            Array.Resize(ref ww, slacklarge + n + totalConstraintslack);
            // First do a homogenous LP do decide if the problem is feasible.
            // (homogenous QP only works if we're very lucky)
            var IOPT = new InteriorPoint.Optimise(slacklarge + n + totalConstraintslack, m + slacklarge + slackb, ww, null, bb, CTEST);
            IOPT.alphamin = 1e-10;
            IOPT.baseA = A;//We only need to pass the constraints without slack variables AA just use for testing
            IOPT.basen = n;
            IOPT.bases = slacklarge;
            IOPT.basem = m;
            IOPT.conv = conv;
            IOPT.compConv = Math.Max(conv, IOPT.compConv);
            IOPT.slacklargeConstraintToStock = slacklargeConstraint;
            IOPT.slackToConstraintBOTH = slackToConstraintBOTH;
            IOPT.slackToConstraintL = slackToConstraintL;
            IOPT.slackToConstraintU = slackToConstraintU;
            var back =
            IOPT.Opt("QP", null, null, true, UL, sign);
            if (back < -10) Console.WriteLine($"Failed -- too many iterations");
            if (back < 0) Console.WriteLine($"Normal Matrix became ill-conditioned");
            if (back == 6) Console.WriteLine("INFEASIBLE");
            else if (Q != null)
            {
                IOPT = new InteriorPoint.Optimise(slacklarge + n + totalConstraintslack, m + slacklarge + slackb, ww, null, bb, cextra, n, HH);
                IOPT.h = hessmull;
                IOPT.baseA = A;
                IOPT.basen = n;
                IOPT.bases = slacklarge;
                IOPT.basem = m;
                IOPT.slacklargeConstraintToStock = slacklargeConstraint;
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
            if (wback != null) BlasLike.dcopyvec(wback.Length, ww, wback);
            return back;
        }

    }
    public class FPortfolio : Portfolio
    {
        public FPortfolio(string file) : base(file)
        {
            inFile = file;
            if (inFile != "")
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

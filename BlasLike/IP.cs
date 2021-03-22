using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace InteriorPoint
{
    public class Optimise
    {
        bool homogenous = true;
        double mu;
        int maxiter = 100;
        int n;
        int m;
        double[] A = null;
        double[] db = null;
        double[] dc = null;
        double[] Ax = null;
        double[] Ay = null;
        double[] c = null;
        double[] cmod = null;
        double[] b = null;
        double[] x = null;
        double[] w1 = null;
        double[] dx = null;
        double[] y = null;
        double[] dy = null;
        double[] z = null;
        double[] dz = null;
        double[] H = null;
        char[] uplo = { 'U' };
        double[] HCOPY = null;
        double tau = 1;
        double dtau = 0;
        double dtau0 = 0;
        double kappa = 1;
        double dkappa = 0;
        double dkappa0 = 0;
        double[] rp;
        double[] rd;
        double[] rmu;
        double hrmu;
        double rkxy;
        int[] order;
        int[] horder;
        double[] M = null;
        double[] dx0 = null;
        double[] dz0 = null;
        Optimise(int n, int m, double[] x, double[] A, double[] b, double[] c, double[] H = null)
        {
            this.n = n;
            this.m = m;
            this.A = A;
            this.x = x;
            this.b = b;
            this.c = c;
            this.H = H;

            z = new double[n];
            dz = new double[n];
            dx = new double[n];
            dz0 = new double[n];
            dx0 = new double[n];
            y = new double[m];
            dy = new double[m];
            w1 = new double[n];
            Ax = new double[m];
            rp = new double[m];
            Ay = new double[n];
            rd = new double[n];
            rmu = new double[n];
            order = new int[m];
            horder = new int[n];
            M = new double[m * (m + 1) / 2];
            cmod = new double[n];
            db = new double[m];
            dc = new double[n];
        }
        void PrimalResidual()
        {
            Factorise.dmxmulv(m, n, A, x, Ax);
            BlasLike.dcopyvec(m, Ax, rp);
            BlasLike.dnegvec(m, rp);
            BlasLike.daxpyvec(m, tau, b, rp);
        }
        void SolvePrimary(double gamma, bool corrector = false)
        {
            var g1 = 1 - gamma;
            var fail = 0;
            for (var k = 0; k < order.Length; ++k) order.SetValue(0, k);
            if (!corrector)
            {
                BlasLike.dzerovec(M.Length, M);
                HCOPY = (double[])H.Clone();
                var lhs = new double[n];
                var res = new double[m];
                for (var i = 0; i < n; ++i)
                {
                    var ct = ActiveSet.Optimise.dprotdiv(ref z[i], ref x[i], ref fail);
                    if (fail != 0 && z[i] == 0) ct = BlasLike.lm_max;
                    HCOPY[i * (i + 3) / 2] += ct;
                }
                Factorise.Factor(uplo, n, HCOPY, horder);
                for (Int16 con = 0, ij = 0; con < m; ++con)
                {
                    BlasLike.dcopy(n, A, m, lhs, 1, con);
                    Factorise.Solve(uplo, n, 1, HCOPY, horder, lhs, n);
                    Factorise.dmxmulv(m, n, A, lhs, res);
                    for (var c1 = 0; c1 <= con; ++c1, ij++)
                    {
                        M[ij] = res[c1];
                    }
                }
                if (m != 1) Factorise.Factor(uplo, m, M, order);
                for (var i = 0; i < n; ++i)
                {
                    var top = (rmu[i] - g1 * mu);
                    var ct = ActiveSet.Optimise.dprotdiv(ref top, ref x[i], ref fail);
                    if (fail != 0 && top == 0) ct = BlasLike.lm_max;
                    w1[i] = rd[i] * g1 - ct;
                }
            }
            else
            {
                for (var i = 0; i < n; ++i)
                {
                    var top = (rmu[i] - g1 * mu - dx0[i] * dz0[i]);
                    var ct = ActiveSet.Optimise.dprotdiv(ref top, ref x[i], ref fail);
                    if (fail != 0 && top == 0) ct = BlasLike.lm_max;
                    w1[i] = rd[i] * g1 - ct;
                }
            }
            Factorise.Solve(uplo, n, 1, HCOPY, horder, w1, n);
            Factorise.dmxmulv(m, n, A, w1, dy);
            BlasLike.daxpyvec(m, g1, rp, dy);
            if (m == 1) dy[0] /= M[0];
            else Factorise.Solve(uplo, m, 1, M, order, dy, m);
            if (homogenous)
            {
                var cx = new double[n];
                for (var i = 0; i < n; ++i)
                {
                    var top = c[i] * x[i];
                    cx[i] = ActiveSet.Optimise.dprotdiv(ref top, ref z[i], ref fail);
                    if (fail != 0 && top == 0) cx[i] = BlasLike.lm_max;
                }
                Factorise.dmxmulv(m, n, A, cx, db);
                BlasLike.daddvec(m, db, b, db);
                if (m == 1) db[0] /= M[0];
                else
                {
                    Factorise.Solve(uplo, m, 1, M, order, db, m);
                }
                Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                for (var i = 0; i < n; ++i)
                {
                    var top = dx[i] * x[i];
                    var cc = ActiveSet.Optimise.dprotdiv(ref top, ref z[i], ref fail);
                    if (fail != 0 && top == 0) cc = BlasLike.lm_max;
                    dx[i] = cc - w1[i];
                }
                Factorise.dmxmulv(n, m, A, db, dc, 0, 0, 0, true);
                for (var i = 0; i < n; ++i)
                {
                    var top = dc[i] * x[i];
                    var cc = ActiveSet.Optimise.dprotdiv(ref top, ref z[i], ref fail);
                    if (fail != 0 && top == 0) cc = BlasLike.lm_max;
                    dc[i] = cc - cx[i];
                }
                var cdx = BlasLike.ddotvec(n, c, dx);
                var cdc = BlasLike.ddotvec(n, c, dc);
                var bdy = BlasLike.ddotvec(m, b, dy);
                var bdb = BlasLike.ddotvec(m, b, db);
                if (corrector) dtau = (cdx - bdy + rkxy * g1 + (hrmu - g1 * mu - dtau0 * dkappa0) / tau) / (bdb - cdc + kappa / tau);
                else dtau = (cdx - bdy + rkxy * g1 + (hrmu - g1 * mu) / tau) / (bdb - cdc + kappa / tau);
                BlasLike.daxpyvec(n, dtau, dc, dx);
                BlasLike.daxpyvec(m, dtau, db, dy);
                if (corrector) dkappa = (hrmu - g1 * mu - kappa * dtau - dtau0 * dkappa0) / tau;
                else dkappa = (hrmu - g1 * mu - kappa * dtau) / tau;
            }
            else
            {
                Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                Factorise.Solve(uplo, n, 1, HCOPY, horder, dx, n);
                BlasLike.dsubvec(n, dx, w1, dx);

                if (corrector)
                {
                    for (var i = 0; i < n; ++i)
                    {
                        var top = (rmu[i] - g1 * mu - dx[i] * z[i] - dx0[i] * dz0[i]);
                        dz[i] = ActiveSet.Optimise.dprotdiv(ref top, ref x[i], ref fail);
                        if (fail != 0 && top == 0) dz[i] = BlasLike.lm_max;
                    }
                }
                else
                {
                    for (var i = 0; i < n; ++i)
                    {
                        var top = (rmu[i] - g1 * mu - dx[i] * z[i]);
                        dz[i] = ActiveSet.Optimise.dprotdiv(ref top, ref x[i], ref fail);
                        if (fail != 0 && top == 0) dz[i] = BlasLike.lm_max;
                    }
                }
            }
        }

        void DualResudual()
        {
            Factorise.dmxmulv(n, m, A, y, Ay, 0, 0, 0, true);
            BlasLike.dcopyvec(n, Ay, rd);
            BlasLike.dnegvec(n, rd);
            BlasLike.daxpyvec(n, tau, c, rd);
            BlasLike.dsubvec(n, rd, z, rd);
        }
        void MuResidual()
        {
            for (int i = 0; i < n; ++i) rmu[i] = mu - x[i] * z[i];
            if (homogenous)
            {
                hrmu = mu - tau * kappa;
                rkxy = kappa + BlasLike.ddotvec(n, c, x) - BlasLike.ddotvec(m, b, y);
            }
        }
        double Complementarity()
        {
            return BlasLike.ddotvec(n, x, z);
        }
        void Mu()
        {
            if (homogenous) mu = (BlasLike.ddotvec(n, x, z) + tau * kappa) / (n + 1);
            else mu = BlasLike.ddotvec(n, x, z) / n;
        }
        public static double norm(double[] aa) => Math.Sqrt(BlasLike.ddotvec(aa.Length, aa, aa));
        public static int Opt(int n, int m, double[] w, double[] A, double[] b, double[] c, double[] H = null)
        {
            var opt = new Optimise(n, m, w, A, b, c, H);
            BlasLike.dcopyvec(c.Length, c, opt.cmod);//Can use cmod to include Hx in c for reporting
            BlasLike.dsetvec(n, 1, opt.x);
            BlasLike.dsetvec(n, 1, opt.z);
            opt.Mu();
            var mu0 = opt.mu;
            opt.PrimalResidual();
            opt.DualResudual();
            opt.MuResidual();
            var rp0 = norm(opt.rp);
            var rd0 = norm(opt.rd);
            var rp1 = rp0;
            var rd1 = rd0;
            var comp0 = opt.Complementarity();
            var compnow = comp0;
            var comp1 = comp0;
            var toosmall = 0.0;
            var i = 0;
            while (true)
            {
                rp1 = norm(opt.rp) / rp0;
                rd1 = norm(opt.rd) / rd0;
                comp1 = opt.Complementarity() / comp0;
                if (rp1 < BlasLike.lm_rooteps && rd1 < BlasLike.lm_rooteps && comp1 < BlasLike.lm_rooteps) break;
                if (i > opt.maxiter) break;
                opt.SolvePrimary(0);
            }
            return -1;
        }
    }
}
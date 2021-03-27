using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace InteriorPoint
{
    public class Optimise
    {
        bool usrH = false;
        bool homogenous = false;
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
        int nh;
        char[] uplo = { 'U' };
        double[] HCOPY = null;
        double tau = 1;
        double dtau = 1;
        double dtau0 = 1;
        double kappa = 1;
        double dkappa = 1;
        double dkappa0 = 1;
        double[] rp = null;
        double ddx;
        double ddz;
        double dd;
        double[] rd = null;
        double[] rmu = null;
        double[] cx = null;
        double hrmu;
        double rkxy;
        int[] order = null;
        int[] horder = null;
        double[] M = null;
        double[] dx0 = null;
        double[] dz0 = null;
        static double norm(double[] aa) => Math.Sqrt(BlasLike.ddotvec(aa.Length, aa, aa));
        double square(double a) => a * a;
        double gfunc(double a) => Math.Min(0.5, square(1 - a)) * (1 - a);
        double aob(double a, double b)
        {
            int fail = 21;
            var back = ActiveSet.Optimise.dprotdiv(ref a, ref b, ref fail);
            if (fail != 0 && b == 0) back = BlasLike.lm_max;
            return back;
        }
        Optimise(int n, int m, double[] x, double[] A, double[] b, double[] c, int nh = 0, double[] H = null)
        {
            this.n = n;
            this.m = m;
            this.A = A;
            this.x = x;
            this.b = b;
            this.c = c;
            this.H = H;
            this.nh = nh;

            z = new double[n];
            dz = new double[n];
            dx = new double[n];
            dz0 = new double[n];
            dx0 = new double[n];
            cx = new double[n];
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
        double Lowest()
        {
            var back = ddx;
            back = Math.Min(back, ddz);
            return Math.Min(back, dd);
        }
        void MaximumStep()
        {
            ddx = 1.0;
            ddz = 1.0;
            dd = 1.0;
            for (int i = 0; i < n; ++i)
            {
                if (dx[i] < 0) ddx = Math.Min(ddx, -aob(x[i], dx[i]));
                if (dz[i] < 0) ddz = Math.Min(ddz, -aob(z[i], dz[i]));
            }
            if (homogenous)
            {
                if (dtau < 0) dd = Math.Min(dd, -aob(tau, dtau));
                if (dkappa < 0) dd = Math.Min(dd, -aob(kappa, dkappa));
            }
        }
        void SolvePrimary(double gamma = 0.0, bool corrector = false)
        {
            var g1 = 1.0 - gamma;
            //            for (var k = 0; k < order.Length; ++k) order.SetValue(0, k);
            //            for (var k = 0; k < horder.Length; ++k) horder.SetValue(0, k);
            if (!corrector)
            {
                BlasLike.dzerovec(M.Length, M);
                if (usrH)
                {
                    var lhs = w1;//Highjack w1 and dy to save reallocation
                    var res = dy;
                    HCOPY = (double[])H.Clone();
                    for (int i = 0, ij = 0; i < nh; ++i, ij += i) HCOPY[i + ij] += aob(z[i], x[i]);
                    Factorise.Factor(uplo, nh, HCOPY, horder);
                    for (int con = 0, ij = 0; con < m; ++con, ij += con)
                    {
                        BlasLike.dcopy(n, A, m, lhs, 1, con);
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, lhs, nh);
                        for (int i = 0; i < (n - nh); ++i) lhs[i + nh] /= aob(z[i + nh], x[i + nh]);
                        Factorise.dmxmulv(m, n, A, lhs, res);
                        BlasLike.dcopyvec(con + 1, res, M, 0, ij);
                    }
                }
                else
                {
                    for (int i = 0, ij = 0; i < m; ++i, ij += i)
                    {
                        for (var k = 0; k < n; ++k)
                        {
                            if (A[k * m + i] != 0.0)
                            {
                                var xoz = aob(x[k], z[k]);
                                xoz *= A[k * m + i];
                                BlasLike.daxpyvec(i + 1, xoz, A, M, k * m, ij);
                            }
                        }
                    }
                }
                if (m != 1) Factorise.Factor(uplo, m, M, order);
                for (var i = 0; i < n; ++i) w1[i] = rd[i] * g1 - aob(rmu[i] - g1 * mu, x[i]);
            }
            else for (var i = 0; i < n; ++i) w1[i] = rd[i] * g1 - aob(rmu[i] - g1 * mu - dx0[i] * dz0[i], x[i]);
            if (usrH)
            {
                Factorise.Solve(uplo, nh, 1, HCOPY, horder, w1, nh);
                for (int i = 0; i < (n - nh); ++i) w1[i + nh] /= aob(z[i + nh], x[i + nh]);
            }
            else for (int i = 0; i < n; ++i) w1[i] *= aob(x[i], z[i]);
            Factorise.dmxmulv(m, n, A, w1, dy);
            BlasLike.daxpyvec(m, g1, rp, dy);
            if (m == 1) dy[0] /= M[0];
            else Factorise.Solve(uplo, m, 1, M, order, dy, m);
            if (homogenous)
            {
                BlasLike.dcopyvec(n, cmod, cx);
                if (usrH)
                {
                    Factorise.Solve(uplo, nh, 1, HCOPY, horder, cx, nh);
                    for (int i = 0; i < (n - nh); ++i) cx[i + nh] /= aob(z[i + nh], x[i + nh]);
                }
                else for (int i = 0; i < n; ++i) cx[i] *= aob(x[i], z[i]);
                Factorise.dmxmulv(m, n, A, cx, db);
                BlasLike.daddvec(m, db, b, db);
                if (m == 1) db[0] /= M[0];
                else Factorise.Solve(uplo, m, 1, M, order, db, m);
                Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                if (usrH)
                {
                    Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh);
                    for (int i = 0; i < (n - nh); ++i) dx[i + nh] /= aob(z[i + nh], x[i + nh]);
                }
                else for (int i = 0; i < n; ++i) dx[i] *= aob(x[i], z[i]);
                BlasLike.dsubvec(n, dx, w1, dx);
                Factorise.dmxmulv(n, m, A, db, dc, 0, 0, 0, true);
                if (usrH)
                {
                    Factorise.Solve(uplo, nh, 1, HCOPY, horder, dc, nh);
                    for (int i = 0; i < (n - nh); ++i) dc[i + nh] /= aob(z[i + nh], x[i + nh]);
                }
                else for (int i = 0; i < n; ++i) dc[i] *= aob(x[i], z[i]);
                BlasLike.dsubvec(n, dc, cx, dc);
                var cdx = BlasLike.ddotvec(n, cmod, dx);
                var cdc = BlasLike.ddotvec(n, cmod, dc);
                var bdy = BlasLike.ddotvec(m, b, dy);
                var bdb = BlasLike.ddotvec(m, b, db);
                //Check the equations in E. D. Andersen∗, C. Roos†, and T. Terlaky‡ they're wrong!!!!!!
                dtau = (cdx - bdy + rkxy * g1 + (hrmu - g1 * mu - (corrector ? dtau0 * dkappa0 : 0)) / tau) / (bdb - cdc + aob(kappa, tau));
                BlasLike.daxpyvec(n, dtau, dc, dx);
                BlasLike.daxpyvec(m, dtau, db, dy);
                for (int i = 0; i < n; ++i) dz[i] = aob((rmu[i] - g1 * mu - dx[i] * z[i] - (corrector ? dx0[i] * dz0[i] : 0)), x[i]);
                dkappa = (hrmu - g1 * mu - kappa * dtau - (corrector ? dtau0 * dkappa0 : 0)) / tau;
            }
            else
            {
                Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                if (usrH)
                {
                    Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh);
                    for (int i = 0; i < (n - nh); ++i) dx[i + nh] /= aob(z[i + nh], x[i + nh]);
                }
                else for (int i = 0; i < n; ++i) dx[i] *= aob(x[i], z[i]);
                BlasLike.dsubvec(n, dx, w1, dx);
                for (var i = 0; i < n; ++i) dz[i] = aob(rmu[i] - g1 * mu - dx[i] * z[i] - ((corrector) ? dx0[i] * dz0[i] : 0), x[i]);
            }
        }
        void PrimalResidual()
        {
            Factorise.dmxmulv(m, n, A, x, Ax);
            BlasLike.dcopyvec(m, Ax, rp);
            BlasLike.dnegvec(m, rp);
            BlasLike.daxpyvec(m, homogenous ? tau : 1, b, rp);
        }
        void DualResudual()
        {
            Factorise.dmxmulv(n, m, A, y, Ay, 0, 0, 0, true);
            BlasLike.dcopyvec(n, Ay, rd);
            BlasLike.dnegvec(n, rd);
            BlasLike.daxpyvec(n, homogenous ? tau : 1, cmod, rd);
            BlasLike.dsubvec(n, rd, z, rd);
        }
        void MuResidual()
        {
            for (int i = 0; i < x.Length; ++i) rmu[i] = mu - x[i] * z[i];
            if (homogenous)
            {
                hrmu = mu - tau * kappa;
                rkxy = kappa + BlasLike.ddotvec(x.Length, cmod, x) - BlasLike.ddotvec(y.Length, b, y);
            }
        }
        double Complementarity()
        {
            return BlasLike.ddotvec(n, x, z) + (homogenous ? tau * kappa : 0);
        }
        void Mu()
        {
            if (homogenous) mu = (BlasLike.ddotvec(x.Length, x, z) + tau * kappa) / (x.Length + 1);
            else mu = BlasLike.ddotvec(x.Length, x, z) / x.Length;
        }
        void update(double[] dx, double[] dy, double[] dz, double dtau, double dkappa, double step)
        {
            BlasLike.daxpyvec(n, step, dx, x);
            BlasLike.daxpyvec(n, step, dz, z);
            BlasLike.daxpyvec(m, step, dy, y);
            if (homogenous)
            {
                tau += dtau * step;
                kappa += dkappa * step;
            }
        }
        double Primal()
        {
            var lin = BlasLike.ddotvec(n, c, x);
            var linextra = BlasLike.ddotvec(n, cmod, x);
            return lin + 0.5 * (linextra - lin);
        }
        double Dual()
        {
            var lin = BlasLike.ddotvec(n, c, x);
            var linextra = BlasLike.ddotvec(n, cmod, x);
            var by = BlasLike.ddotvec(m, b, y);
            return by - 0.5 * (linextra - lin);
        }
        public static int Opt(int n, int m, double[] w, double[] A, double[] b, double[] c, int nh = 0, double[] H = null)
        {
            ActiveSet.Optimise.clocker(true);
            var opt = new Optimise(n, m, w, A, b, c, nh, H);
            opt.homogenous = true;
            opt.tau = 1;
            opt.kappa = 1;
            opt.usrH = nh > 0 && BlasLike.dsumvec(opt.H.Length, opt.H) != 0.0;
            BlasLike.dsetvec(n, 1, opt.x);
            BlasLike.dsetvec(n, 1, opt.z);
            var dxold = (double[])opt.dx.Clone();
            var dzold = (double[])opt.dz.Clone();
            var dyold = (double[])opt.dy.Clone();
            var dtauold = opt.dtau;
            var dkappaold = opt.kappa;
            opt.Mu();
            var mu0 = opt.mu;
            var i = 0;
            var extra = new double[nh];
            if (opt.usrH)
            {
                if (opt.homogenous)
                {
                    BlasLike.dcopyvec(opt.c.Length, opt.x, opt.cmod);
                    BlasLike.dscalvec(opt.c.Length, 1.0 / opt.tau, opt.cmod);
                    Factorise.dsmxmulv(nh, opt.H, opt.cmod, extra);
                }
                else Factorise.dsmxmulv(nh, opt.H, opt.x, extra);
                BlasLike.dzerovec(n - nh, opt.cmod, nh);
                BlasLike.daddvec(nh, opt.c, extra, opt.cmod);
            }
            else BlasLike.dcopyvec(opt.c.Length, opt.c, opt.cmod);
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
            var alpha1 = 0.0;
            var alpha2 = 0.0;
            var gamma = 0.0;
            while (true)
            {
                rp1 = norm(opt.rp) / rp0;
                rd1 = norm(opt.rd) / rd0;
                comp1 = opt.Complementarity() / comp0;
                if (rp1 < BlasLike.lm_eps && rd1 < BlasLike.lm_eps && comp1 < BlasLike.lm_eps) break;
                if (i > opt.maxiter) break;
                opt.SolvePrimary();
                BlasLike.dcopyvec(n, opt.dx, dxold);
                BlasLike.dcopyvec(n, opt.dz, dzold);
                BlasLike.dcopyvec(m, opt.dy, dyold);
                dtauold = opt.dtau;
                dkappaold = opt.dkappa;
                opt.MaximumStep();
                alpha1 = 0.99 * opt.Lowest();
                BlasLike.dsccopyvec(opt.n, alpha1, opt.dx, opt.dx0);
                BlasLike.dsccopyvec(opt.n, alpha1, opt.dz, opt.dz0);
                opt.dtau0 = alpha1 * opt.dtau;
                opt.dkappa0 = alpha1 * opt.dkappa;
                gamma = opt.gfunc(alpha1);
                opt.SolvePrimary(gamma, true);
                opt.MaximumStep();
                alpha2 = 0.99 * opt.Lowest();
                if (alpha1 > alpha2) opt.update(dxold, dyold, dzold, dtauold, dkappaold, alpha1);
                else opt.update(opt.dx, opt.dy, opt.dz, opt.dtau, opt.dkappa, alpha2);
                if (opt.usrH)
                {
                    if (opt.homogenous)
                    {
                        BlasLike.dcopyvec(opt.c.Length, opt.x, opt.cmod);
                        BlasLike.dscalvec(opt.c.Length, 1.0 / opt.tau, opt.cmod);
                        Factorise.dsmxmulv(nh, opt.H, opt.cmod, extra);
                    }
                    else Factorise.dsmxmulv(nh, opt.H, opt.x, extra);
                    BlasLike.dzerovec(n - nh, opt.cmod, nh);
                    BlasLike.daddvec(nh, opt.c, extra, opt.cmod);
                }
                opt.Mu();
                mu0 = opt.mu;
                opt.PrimalResidual();
                opt.DualResudual();
                opt.MuResidual();
                i++;
            }
            if (opt.homogenous)
            {
                Console.WriteLine($"tau = {opt.tau} kappa={opt.kappa}");
                if (opt.tau > opt.kappa)
                {
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                }
                else Console.WriteLine("INFEASIBLE");
            }
            Console.WriteLine($"{i} iterations out of {opt.maxiter}");
            Console.WriteLine($"Primal Utility:\t\t{opt.Primal()}");
            ActiveSet.Optimise.printV("x", opt.x);
            Console.WriteLine($"Dual Utility:\t\t{opt.Dual()}");
            ActiveSet.Optimise.printV("y", opt.y);
            ActiveSet.Optimise.printV("z", opt.z);
            Console.WriteLine($"Complementarity:\t{opt.Complementarity()}");
            Console.WriteLine($"Job took {ActiveSet.Optimise.clocker()} m secs");
            if (i >= opt.maxiter) return -1;
            else if (opt.homogenous && opt.tau < opt.kappa) return 6;
            else return 0;
        }
    }
}
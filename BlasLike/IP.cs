using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace InteriorPoint
{
    public enum conetype { QP, SOCP, SOCPR };
    public class Optimise
    {
        string optMode = "QP";
        int numberOfCones = 0;
        int[] cone = null;
        int[] typecone = null;
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
        double[] W = null;
        double[] W2 = null;
        double[] THETA = null;
        double[] dx = null;
        double[] y = null;
        double[] dy = null;
        double[] z = null;
        double[] dz = null;
        double[] H = null;
        double[] xbar = null;
        double[] zbar = null;
        double[] dxbar = null;
        double[] dzbar = null;
        int nh;
        char uplo = 'U';
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
        static double square(double a) => a * a;
        double gfunc(double a) => Math.Min(0.5, square(1 - a)) * (1 - a);
        double aob(double a, double b)
        {
            int fail = 21;
            var back = BlasLike.dprotdiv(ref a, ref b, ref fail);
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
            w1 = new double[2 * n];
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
        void MaximumStep(double gamma = 1.0)
        {
            double XX = BlasLike.ddotvec(n, x, x);
            double XdX = BlasLike.ddotvec(n, x, dx);
            double dXdX = BlasLike.ddotvec(n, dx, dx);
            double SS = BlasLike.ddotvec(n, z, z);
            double SdS = BlasLike.ddotvec(n, z, dz);
            double dSdS = BlasLike.ddotvec(n, dz, dz);
            double XS = BlasLike.ddotvec(n, x, z);
            double XdS = BlasLike.ddotvec(n, x, dz);
            double dXS = BlasLike.ddotvec(n, dx, z);
            double dXdS = BlasLike.ddotvec(n, dx, dz);
            double alpha = 1.0, desc;
            double lowest = 0.1, lowest1 = 1 - lowest;

            if (dXdX <= BlasLike.lm_eps)
            {
                if (XdX < -BlasLike.lm_eps)
                    alpha = Math.Min(alpha, lowest1 * (-XX / (XdX + XdX)));
            }
            else
            {
                if (Math.Abs(XdX) > BlasLike.lm_eps) desc = 1.0 - XX * dXdX / XdX / XdX;
                else desc = -XX * dXdX;
                if (desc > BlasLike.lm_eps * 8)
                {
                    if (Math.Abs(XdX) > BlasLike.lm_eps) desc = Math.Sqrt(desc) * Math.Abs(XdX);
                    else desc = Math.Sqrt(desc);
                    if (XdX + desc < 0)
                        alpha = Math.Min(alpha, lowest1 * ((-XdX - desc) / dXdX));
                }
                else if (desc >= -BlasLike.lm_eps && XdX < 0)
                    alpha = Math.Min(alpha, lowest1 * ((-XdX) / dXdX));
            }

            if (dSdS <= BlasLike.lm_eps)
            {
                if (SdS < -BlasLike.lm_eps)
                    alpha = Math.Min(alpha, lowest1 * (-SS / (SdS + SdS)));
            }
            else
            {
                if (Math.Abs(SdS) > BlasLike.lm_eps) desc = 1.0 - SS * dSdS / SdS / SdS;
                else desc = -SS * dSdS;
                if (desc > BlasLike.lm_eps * 8)
                {
                    if (Math.Abs(SdS) > BlasLike.lm_eps) desc = Math.Sqrt(desc) * Math.Abs(SdS);
                    else desc = Math.Sqrt(desc);
                    if (SdS + desc < 0)
                        alpha = Math.Min(alpha, lowest1 * ((-SdS - desc) / dSdS));
                }
                else if (desc >= -BlasLike.lm_eps * 8 && SdS < 0)
                    alpha = Math.Min(alpha, lowest1 * ((-SdS) / dSdS));
            }

            if (Math.Abs(dXdS) <= BlasLike.lm_eps)
            {
                if (dXS + XdS < -BlasLike.lm_eps)
                    alpha = Math.Min(alpha, lowest1 * (-XS / (XdS + dXS)));
            }
            else
            {
                if (Math.Abs(XdS + dXS) > BlasLike.lm_eps) desc = 1.0 - 4 * XS * dXdS / (XdS + dXS) / (XdS + dXS);
                else desc = -4 * XS * dXdS;
                if (desc > BlasLike.lm_eps * 8)
                {
                    if (Math.Abs(XdS + dXS) > BlasLike.lm_eps) desc = Math.Sqrt(desc) * Math.Abs(XdS + dXS);
                    else desc = Math.Sqrt(desc);
                    if ((XdS + dXS + desc) > 0 && dXdS < 0)
                        alpha = Math.Min(alpha, lowest1 * ((-XdS - dXS - desc) / 2.0 / dXdS));
                    else if ((XdS + dXS + desc) / dXdS < 0)
                        alpha = Math.Min(alpha, lowest1 * ((-XdS - dXS - desc) / 2.0 / dXdS));
                }
                else if (desc >= -BlasLike.lm_eps * 8 && (XdS + dXS) / dXdS < 0)
                    alpha = Math.Min(alpha, lowest1 * ((-XdS - dXS) / 2.0 / dXdS));
            }

            if (optMode == "SOCP")
            {
                ddx = alpha;
                ddz = alpha;
                dd = alpha;
                var vz1 = new double[cone.Length];
                var vz2 = new double[cone.Length];
                var vz3 = new double[cone.Length];
                var vx1 = new double[cone.Length];
                var vx2 = new double[cone.Length];
                var vx3 = new double[cone.Length];
                for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                {
                    var n = cone[icone];
                    if (typecone[icone] == (int)conetype.QP)
                    {
                        for (int i = cstart; i < n + cstart; ++i)
                        {
                            if (dx[i] < 0) ddx = Math.Min(ddx, -aob(x[i], dx[i]));
                            if (dz[i] < 0) ddz = Math.Min(ddz, -aob(z[i], dz[i]));
                        }
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        var Qz = new double[n];
                        var Qdz = new double[n];
                        var Qx = new double[n];
                        var Qdx = new double[n];
                        BlasLike.dcopyvec(n, z, Qz, cstart, 0);//Qz
                        Qmulvec(n, Qz);
                        BlasLike.dcopyvec(n, dz, Qdz, cstart, 0);//Qdz
                        Qmulvec(n, Qdz);
                        vz1[icone] = BlasLike.ddotvec(n, z, Qz, cstart, 0);//z.Qz
                        vz2[icone] = 2.0 * BlasLike.ddotvec(n, dz, Qz, cstart, 0);//dz.Qz
                        vz3[icone] = BlasLike.ddotvec(n, dz, Qdz, cstart, 0);//dz.Qdz
                        if (Qdz[n - 1] < 0) alpha = Math.Min(lowest1 * (-Qz[n - 1]) / Qdz[n - 1], alpha);
                        if (dz[n - 1 + cstart] < 0) alpha = Math.Min(lowest1 * (-z[n - 1 + cstart]) / dz[n - 1 + cstart], alpha);
                        BlasLike.dcopyvec(n, x, Qx, cstart, 0);//Qx
                        Qmulvec(n, Qx);
                        BlasLike.dcopyvec(n, dx, Qdx, cstart, 0);//Qdx
                        Qmulvec(n, Qdx);
                        vx1[icone] = BlasLike.ddotvec(n, x, Qx, cstart, 0);//x.Qx
                        vx2[icone] = 2.0 * BlasLike.ddotvec(n, dx, Qx, cstart, 0);//dx.Qx
                        vx3[icone] = BlasLike.ddotvec(n, dx, Qdx, cstart, 0);//dx.Qdx
                        if (Qdx[n - 1] < 0) alpha = Math.Min(lowest1 * (-Qx[n - 1]) / Qdx[n - 1], alpha);
                        if (dx[n - 1 + cstart] < 0) alpha = Math.Min(lowest1 * (-x[n - 1 + cstart]) / dx[n - 1 + cstart], alpha);
                        if (homogenous && dtau < 0) alpha = Math.Min(lowest1 * -tau / dtau, alpha);
                        if (homogenous && dkappa < 0) alpha = Math.Min(lowest1 * -kappa / dkappa, alpha);
                        double inner, r1, r2;
                        for (var i = 0; i < n; ++i)
                        {
                            if (vz1[i] + alpha * (vz2[i] + alpha * vz3[i]) < -BlasLike.lm_eps * 8)
                            {
                                if (Math.Abs(vz3[i]) <= BlasLike.lm_eps)
                                {
                                    if (vz2[i] < -BlasLike.lm_eps)
                                        alpha = Math.Min(lowest1 * -vz1[i] / vz2[i], alpha);
                                }
                                else if (Math.Abs(vz2[i]) > BlasLike.lm_eps && (inner = 1.0 - 4 * vz3[i] * vz1[i] / vz2[i] / vz2[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > BlasLike.lm_eps * 8 ? Math.Sqrt(inner) * Math.Abs(vz2[i]) : 0);
                                    r1 = (-vz2[i] - inner) / 2.0 / vz3[i]; r2 = (-vz2[i] + inner) / 2.0 / vz3[i];
                                    if (vz3[i] < -BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vz2[i] - inner) / 2.0 / vz3[i], alpha);
                                    }
                                    else if (vz3[i] > BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vz2[i] - inner) / 2.0 / vz3[i], alpha);
                                    }
                                }
                                else if (Math.Abs(vz2[i]) <= BlasLike.lm_eps && (inner = -4 * vz3[i] * vz1[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > BlasLike.lm_eps * 8 ? Math.Sqrt(inner) : 0);
                                    r1 = (-vz2[i] - inner) / 2.0 / vz3[i]; r2 = (-vz2[i] + inner) / 2.0 / vz3[i];
                                    if (vz3[i] < -BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vz2[i] - inner) / 2.0 / vz3[i], alpha);
                                    }
                                    else if (vz3[i] > BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vz2[i] - inner) / 2.0 / vz3[i], alpha);
                                    }
                                }
                                else
                                    Console.WriteLine("still negative");
                            }
                        }

                        for (var i = 0; i < n; ++i)
                        {
                            if (vx1[i] + alpha * (vx2[i] + alpha * vx3[i]) < -BlasLike.lm_eps * 8)
                            {
                                if (Math.Abs(vx3[i]) <= BlasLike.lm_eps)
                                {
                                    if (vx2[i] < -BlasLike.lm_eps)
                                        alpha = Math.Min(lowest1 * -vx1[i] / vx2[i], alpha);
                                }
                                else if (Math.Abs(vx2[i]) > BlasLike.lm_eps && (inner = 1.0 - 4 * vx3[i] * vx1[i] / vx2[i] / vx2[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > BlasLike.lm_eps * 8 ? Math.Sqrt(inner) * Math.Abs(vx2[i]) : 0);
                                    r1 = (-vx2[i] - inner) / 2.0 / vx3[i]; r2 = (-vx2[i] + inner) / 2.0 / vx3[i];
                                    if (vx3[i] < -BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vx2[i] - inner) / 2.0 / vx3[i], alpha);
                                    }
                                    else if (vx3[i] > BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vx2[i] - inner) / 2.0 / vx3[i], alpha);
                                    }
                                }
                                else if (Math.Abs(vx2[i]) <= BlasLike.lm_eps && (inner = -4 * vx3[i] * vx1[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > BlasLike.lm_eps * 8 ? Math.Sqrt(inner) : 0);
                                    r1 = (-vx2[i] - inner) / 2.0 / vx3[i]; r2 = (-vx2[i] + inner) / 2.0 / vx3[i];
                                    if (vx3[i] < -BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vx2[i] - inner) / 2.0 / vx3[i], alpha);
                                    }
                                    else if (vx3[i] > BlasLike.lm_eps)
                                    {
                                        alpha = Math.Min(lowest1 * (-vx2[i] - inner) / 2.0 / vx3[i], alpha);
                                    }
                                }
                                else
                                    Console.WriteLine("still negative");
                            }
                        }


                        double rhs, gamma1 = 1 - gamma, test1, test2 = 1, beta = 1e-6;
                        for (var l = 0; l < 1000; ++l)
                        {
                            rhs = beta * (1 - alpha * gamma1) * mu;
                            if (homogenous) test2 = (tau + alpha * dtau) * (kappa + alpha * dkappa);
                            test1 = 0;
                            for (var i = 0; i < n; ++i)
                            {
                                test1 += (vx1[i] + alpha * (vx2[i] + alpha * vx3[i])) * (vz1[i] + alpha * (vz2[i] + alpha * vz3[i]));
                            }

                            test1 = Math.Sqrt(test1);
                            if (test1 >= rhs && test2 >= rhs) break;
                            alpha *= lowest1;
                        }

                        ddx = ddz = dd = alpha;
                    }
                }
                if (homogenous)
                {
                    if (dtau < 0) dd = Math.Min(dd, -aob(tau, dtau));
                    if (dkappa < 0) dd = Math.Min(dd, -aob(kappa, dkappa));
                }
            }

            else if (optMode == "QP")
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
        }
        void SolvePrimary(double gamma = 0.0, bool corrector = false)
        {
            if (optMode == "QP")
            {
                var g1 = 1.0 - gamma;
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
            else if (optMode == "SOCP")
            {
                var g1 = 1.0 - gamma;
                if (!corrector)
                {
                    BlasLike.dzerovec(M.Length, M);
                    var lhs = w1;//Highjack w1 and dy to save reallocation
                    var res = dy;
                    HCOPY = (double[])H.Clone();
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            if (usrH)
                            {
                                for (int i = 0, ij = 0; i < nh; ++i, ij += i) HCOPY[i + ij] += 1.0 * W2[i + cstart];
                                Factorise.Factor(uplo, nh, HCOPY, horder);
                                for (int con = 0, ij = 0; con < m; ++con, ij += con)
                                {
                                    BlasLike.dcopy(n, A, m, lhs, 1, con + cstart * m, cstart);
                                    Factorise.Solve(uplo, nh, 1, HCOPY, horder, lhs, nh, 0, 0, cstart);
                                    for (int i = nh + cstart; i < n + cstart; ++i) lhs[i] /= W2[i];
                                    Factorise.dmxmulv(m, n, A, lhs, res, cstart * m, cstart, cstart);
                                    BlasLike.dcopyvec(con + 1, res, M, cstart, ij);
                                }
                            }
                            else
                            {
                                for (int i = 0, ij = 0; i < m; ++i, ij += i)
                                {
                                    for (var k = cstart; k < n + cstart; ++k)
                                    {
                                        if (A[k * m + i] != 0.0)
                                        {
                                            var xoz = 1.0 / W2[k];
                                            xoz *= A[k * m + i];
                                            BlasLike.daxpyvec(i + 1, xoz, A, M, k * m, ij);
                                        }
                                    }
                                }
                            }
                            for (var i = cstart; i < n + cstart; ++i)
                            {
                                w1[i] = rd[i] * g1 - aob(rmu[i] - g1 * mu, xbar[i]) * W[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            for (int i = 0, ij = 0; i < m; ++i, ij += i)
                            {
                                BlasLike.dcopy(n, A, m, lhs, 1, i + cstart * m, x.Length + cstart);
                                Qmulvec(n, lhs, x.Length + cstart);
                                W2trans(n, lhs, W, lhs, x.Length + cstart, cstart, cstart);
                                Qmulvec(n, lhs, cstart);
                                thetaScale(n, lhs, THETA[icone], true, true, cstart);
                                for (var k = cstart; k < n + cstart; ++k)
                                {
                                    if (lhs[k] != 0.0)
                                    {
                                        BlasLike.daxpyvec(i + 1, lhs[k], A, M, k * m, ij);
                                    }
                                }
                            }
                            rmu[n - 1 + cstart] -= g1 * mu;
                            applyXm1(n, xbar, rmu, lhs, cstart, cstart, x.Length + cstart);
                            rmu[n - 1 + cstart] += g1 * mu;
                            Wtrans(n, lhs, W, w1, x.Length + cstart, cstart, cstart);
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            for (var i = cstart; i < n + cstart; ++i)
                            {
                                w1[i] = rd[i] * g1 - w1[i];
                            }
                        }
                    }
                    if (m != 1) Factorise.Factor(uplo, m, M, order);
                }
                else
                {
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            for (var i = cstart; i < n + cstart; ++i)
                            {//-dx0[i]*W[i] *dz0[i]/W[i]
                                w1[i] = rd[i] * g1 - aob(rmu[i] - g1 * mu - dx0[i] * dz0[i], xbar[i]) * W[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            var lhs = w1;
                            rmu[n - 1 + cstart] -= g1 * mu;
                            applyXm1(n, xbar, rmu, lhs, cstart, cstart, x.Length + cstart);
                            rmu[n - 1 + cstart] += g1 * mu;
                            Wtrans(n, lhs, W, w1, x.Length + cstart, cstart, cstart);
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            for (var i = cstart; i < n + cstart; ++i)
                            {
                                w1[i] = rd[i] * g1 - w1[i];
                            }
                            BlasLike.dcopyvec(n, dz0, lhs, cstart, x.Length + cstart);
                            thetaScale(icone, lhs, THETA[icone], true, false, x.Length + cstart);
                            Qmulvec(icone, lhs, x.Length + cstart);
                            Wtrans(icone, lhs, W, dzbar, x.Length + cstart, cstart, cstart);
                            Qmulvec(icone, dzbar, cstart);
                            Wtrans(icone, dx0, W, dxbar, cstart, cstart, cstart);
                            thetaScale(icone, dxbar, THETA[icone], false, false, cstart);
                            applyX(icone, dxbar, dzbar, lhs, cstart, cstart, x.Length + cstart);
                            BlasLike.dsubvec(n, w1, lhs, w1, cstart, x.Length + cstart, cstart);
                        }
                    }
                }
                for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                {
                    var n = cone[icone];
                    if (typecone[icone] == (int)conetype.QP)
                    {
                        if (usrH)
                        {
                            Factorise.Solve(uplo, nh, 1, HCOPY, horder, w1, nh, 0, 0, cstart);
                            for (int i = nh + cstart; i < n + cstart; ++i) w1[i] /= W2[i];
                        }
                        else
                        {
                            for (int i = cstart; i < n + cstart; ++i) w1[i] /= W2[i];
                        }
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        Qmulvec(n, w1, cstart);
                        W2trans(n, w1, W, w1, cstart, cstart, x.Length + cstart);
                        BlasLike.dcopyvec(n, w1, w1, x.Length + cstart, cstart);
                        Qmulvec(n, w1, cstart);
                        thetaScale(n, w1, THETA[icone], true, true, cstart);
                    }
                }
                Factorise.dmxmulv(m, n, A, w1, dy);
                BlasLike.daxpyvec(m, g1, rp, dy);
                if (m == 1) dy[0] /= M[0];
                else Factorise.Solve(uplo, m, 1, M, order, dy, m);
                if (homogenous)
                {
                    BlasLike.dcopyvec(n, cmod, cx);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            if (usrH)
                            {
                                Factorise.Solve(uplo, nh, 1, HCOPY, horder, cx, nh, 0, 0, cstart);
                                for (int i = cstart + nh; i < n + cstart; ++i) cx[i] /= W2[i];
                            }
                            else
                            {
                                for (int i = cstart; i < n + cstart; ++i) cx[i] /= W2[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            Qmulvec(n, cx, cstart);
                            W2trans(n, cx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, cx, x.Length + cstart, cstart);
                            Qmulvec(n, cx, cstart);
                            thetaScale(n, cx, THETA[icone], true, true, cstart);
                        }
                    }
                    Factorise.dmxmulv(m, n, A, cx, db);
                    BlasLike.daddvec(m, db, b, db);
                    if (m == 1) db[0] /= M[0];
                    else Factorise.Solve(uplo, m, 1, M, order, db, m);
                    Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            if (usrH)
                            {
                                Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh, 0, 0, cstart);
                                for (int i = nh + cstart; i < n + cstart; ++i) dx[i] /= W2[i];
                            }
                            else
                            {
                                for (int i = cstart; i < n + cstart; ++i) dx[i] /= W2[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            Qmulvec(n, dx, cstart);
                            W2trans(n, dx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dx, x.Length + cstart, cstart);
                            Qmulvec(n, dx, cstart);
                            thetaScale(n, dx, THETA[icone], true, true, cstart);
                        }
                    }
                    BlasLike.dsubvec(n, dx, w1, dx);
                    Factorise.dmxmulv(n, m, A, db, dc, 0, 0, 0, true);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            if (usrH)
                            {
                                Factorise.Solve(uplo, nh, 1, HCOPY, horder, dc, nh, 0, 0, cstart);
                                for (int i = nh + cstart; i < n + cstart; ++i) dc[i] /= W2[i];
                            }
                            else
                            {
                                for (int i = cstart; i < n + cstart; ++i) dc[i] /= W2[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            Qmulvec(n, dc, cstart);
                            W2trans(n, dc, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dc, x.Length + cstart, cstart);
                            Qmulvec(n, dc, cstart);
                            thetaScale(n, dc, THETA[icone], true, true, cstart);
                        }
                    }
                    BlasLike.dsubvec(n, dc, cx, dc);
                    var cdx = BlasLike.ddotvec(n, cmod, dx);
                    var cdc = BlasLike.ddotvec(n, cmod, dc);
                    var bdy = BlasLike.ddotvec(m, b, dy);
                    var bdb = BlasLike.ddotvec(m, b, db);
                    //Check the equations in E. D. Andersen∗, C. Roos†, and T. Terlaky‡ they're wrong!!!!!!
                    dtau = (cdx - bdy + rkxy * g1 + (hrmu - g1 * mu - (corrector ? dtau0 * dkappa0 : 0)) / tau) / (bdb - cdc + aob(kappa, tau));
                    BlasLike.daxpyvec(n, dtau, dc, dx);
                    BlasLike.daxpyvec(m, dtau, db, dy);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            for (int i = cstart; i < n + cstart; ++i) dz[i] = aob((rmu[i] - g1 * mu - dx[i] * zbar[i] * W[i] - (corrector ? dx0[i] * dz0[i] : 0)), xbar[i]) * W[i];
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            double[] dxWzbarxbarm1 = new double[n];
                            Wtrans(n, dx, W, w1, cstart, cstart, x.Length + cstart);
                            thetaScale(n, w1, THETA[icone], false, false, x.Length + cstart);
                            applyX(n, zbar, w1, dxWzbarxbarm1, cstart, x.Length + cstart, 0);

                            rmu[n - 1 + cstart] -= g1 * mu;
                            BlasLike.dcopyvec(n, rmu, w1, cstart, cstart);
                            rmu[n - 1 + cstart] += g1 * mu;
                            if (corrector)
                            {
                                BlasLike.dcopyvec(n, dz0, w1, cstart, x.Length + cstart);//dz0
                                thetaScale(icone, w1, THETA[icone], true, false, x.Length + cstart);//theta-1dz0
                                Qmulvec(n, w1, x.Length + cstart);
                                Wtrans(n, w1, W, dzbar, x.Length + cstart, cstart, cstart);//Wtheta-1dz0
                                Qmulvec(n, dzbar, cstart);//dzbar
                                Wtrans(n, dx0, W, dxbar, cstart, cstart, cstart);
                                thetaScale(n, dxbar, THETA[icone], false, false, cstart);//dxbar
                                applyX(n, dxbar, dzbar, w1, cstart, cstart, x.Length + cstart);//dxbardzbar
                                BlasLike.dsubvec(n, w1, w1, w1, cstart, x.Length + cstart, cstart);
                            }
                            BlasLike.dsubvec(n, w1, dxWzbarxbarm1, w1, cstart, 0, cstart);
                            applyXm1(n, xbar, w1, w1, cstart, cstart, x.Length + cstart);//over xbar

                            Wtrans(n, w1, W, dz, x.Length + cstart, cstart, cstart);//times W
                            thetaScale(n, dz, THETA[icone], false, false, cstart);
                        }
                    }
                    dkappa = (hrmu - g1 * mu - kappa * dtau - (corrector ? dtau0 * dkappa0 : 0)) / tau;
                }
                else
                {
                    Factorise.dmxmulv(n, m, A, dy, dx, 0, 0, 0, true);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            if (usrH)
                            {
                                Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh, 0, 0, cstart);
                                for (int i = nh + cstart; i < n + cstart; ++i) dx[i] /= W2[i];
                            }
                            else
                            {
                                for (int i = cstart; i < n + cstart; ++i) dx[i] /= W2[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            Qmulvec(icone, dx, cstart);
                            W2trans(icone, dx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dx, x.Length + cstart, cstart);
                            Qmulvec(icone, dx, cstart);
                            thetaScale(icone, dx, THETA[icone], true, true, cstart);
                        }
                    }
                    BlasLike.dsubvec(n, dx, w1, dx);
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            for (var i = cstart; i < n + cstart; ++i) dz[i] = aob(rmu[i] - g1 * mu - dx[i] * zbar[i] * W[i] - ((corrector) ? dx0[i] * dz0[i] : 0), xbar[i]) * W[i];
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            double[] dxWzbarxbarm1 = new double[n];
                            Wtrans(n, dx, W, w1, cstart, cstart, x.Length + cstart);
                            thetaScale(n, w1, THETA[icone], false, false, x.Length + cstart);
                            applyX(n, zbar, w1, dxWzbarxbarm1, cstart, x.Length + cstart, 0);

                            rmu[n - 1 + cstart] -= g1 * mu;
                            BlasLike.dcopyvec(n, rmu, w1, cstart, cstart);
                            rmu[n - 1 + cstart] += g1 * mu;
                            if (corrector)
                            {
                                BlasLike.dcopyvec(n, dz0, w1, cstart, x.Length + cstart);//dz0
                                thetaScale(icone, w1, THETA[icone], true, false, x.Length + cstart);//theta-1dz0
                                Qmulvec(n, w1, x.Length + cstart);
                                Wtrans(n, w1, W, dzbar, x.Length + cstart, cstart, cstart);//Wtheta-1dz0
                                Qmulvec(n, dzbar, cstart);//dzbar
                                Wtrans(n, dx0, W, dxbar, cstart, cstart, cstart);
                                thetaScale(n, dxbar, THETA[icone], false, false, cstart);//dxbar
                                applyX(n, dxbar, dzbar, w1, cstart, cstart, x.Length + cstart);//dxbardzbar
                                BlasLike.dsubvec(n, w1, w1, w1, cstart, x.Length + cstart, cstart);
                            }
                            BlasLike.dsubvec(n, w1, dxWzbarxbarm1, w1, cstart, 0, cstart);
                            applyXm1(n, xbar, w1, w1, cstart, cstart, x.Length + cstart);//over xbar

                            Wtrans(n, w1, W, dz, x.Length + cstart, cstart, cstart);//times W
                            thetaScale(n, dz, THETA[icone], false, false, cstart);
                        }
                    }
                }
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
            if (optMode == "QP") for (int i = 0; i < x.Length; ++i) rmu[i] = mu - x[i] * z[i];
            else if (optMode == "SOCP")
            {
                for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                {
                    var n = cone[icone];
                    if (typecone[icone] == (int)conetype.QP)
                    {
                        for (int i = cstart; i < n + cstart; ++i)
                        {
                            W2[i] = aob(z[i], x[i]);
                            W[i] = Math.Sqrt(W2[i]);
                            THETA[icone] = -1e4;
                            xbar[i] = x[i] * W[i];
                            zbar[i] = z[i] / W[i];
                            rmu[i] = mu - xbar[i] * zbar[i];
                        }
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        var xQx = 0.0;
                        var zQz = 0.0;
                        if (n == 1) THETA[icone] = Math.Sqrt(aob(z[cstart], x[cstart]));
                        else
                        {
                            xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                            if (xQx <= 0)
                            {
                                if (x[n - 1 + cstart] > BlasLike.lm_rooteps)
                                {
                                    BlasLike.dscalvec(n - 1, .95 * x[n - 1 + cstart] / Math.Sqrt(x[n - 1 + cstart] * x[n - 1 + cstart] - xQx), x);
                                }
                                else
                                {
                                    BlasLike.dsetvec(n - 1, BlasLike.lm_rooteps, x, cstart);
                                    x[n - 1 + cstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)n);
                                    xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                }
                            }
                            zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                            if (zQz <= 0)
                            {
                                if (z[n - 1 + cstart] > BlasLike.lm_rooteps)
                                {
                                    BlasLike.dscalvec(n - 1, .95 * z[n - 1 + cstart] / Math.Sqrt(z[n - 1 + cstart] * z[n - 1 + cstart] - zQz), z);
                                }
                                else
                                {
                                    BlasLike.dsetvec(n - 1, BlasLike.lm_rooteps, z, cstart);
                                    z[n - 1 + cstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)n);
                                    zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                                }
                            }
                            THETA[icone] = Math.Sqrt(Math.Sqrt(zQz / xQx));
                            if (double.IsNaN(THETA[icone]))
                            {
                                if (double.IsNaN(zQz))
                                {
                                    BlasLike.dsetvec(n - 1, BlasLike.lm_rooteps, z, cstart);
                                    z[n - 1 + cstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)n);
                                    zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                                }
                                if (double.IsNaN(xQx))
                                {
                                    BlasLike.dsetvec(n - 1, BlasLike.lm_rooteps, x, cstart);
                                    x[n - 1 + cstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)n);
                                    xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                }
                                THETA[icone] = Math.Sqrt(Math.Sqrt(zQz / xQx));
                            }

                        }

                        if (!double.IsNaN(THETA[icone]) && THETA[icone] >= 0)
                        {
                            if (n == 1) W[cstart] = Math.Sqrt(aob(z[cstart], x[cstart]));
                            else
                            {
                                double zx = BlasLike.ddotvec(n, z, x, cstart, cstart);
                                double bot = Math.Sqrt((zx + Math.Sqrt(xQx * zQz)) * 2.0);
                                double z1 = THETA[icone] / bot;
                                double z2 = THETA[icone] * bot;
                                if (THETA[icone] == BlasLike.lm_eps) { z1 = 1.0; z2 = BlasLike.lm_eps; }
                                for (var i = cstart; i < n - 1 + cstart; ++i)
                                {
                                    W[i] = -z1 * x[i] + z[i] / z2;
                                }
                                W[n - 1 + cstart] = z1 * x[n - 1 + cstart] + z[n - 1 + cstart] / z2;
                            }
                        }
                        else
                        {
                            //DO something to stop main loop
                        }
                        double wcheck = BlasLike.ddotvec(n, W, W, cstart, cstart);
                        if (double.IsNaN(wcheck)) Console.WriteLine("BAD W");
                        Wtrans(n, x, W, xbar, cstart, cstart, cstart);
                        Wtrans(n, z, W, zbar, cstart, cstart, cstart);
                        thetaScale(n, xbar, THETA[icone], false, false, cstart);
                        thetaScale(n, zbar, THETA[icone], true, false, cstart);
                        Tmulvec(n, xbar, cstart);//Tmulvec does nothing for SOCP, needed for SOCPR
                        Tmulvec(n, zbar, cstart);
                        applyX(n, xbar, zbar, rmu, cstart, cstart, cstart);
                        BlasLike.dnegvec(n, rmu, cstart);
                        rmu[n - 1 + cstart] += mu;
                        Tmulvec(n, xbar, cstart);
                        Tmulvec(n, zbar, cstart);
                        Tmulvec(n, rmu, cstart);
                    }
                }
            }
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
            if (optMode == "QP")
            {
                if (homogenous) mu = (BlasLike.ddotvec(x.Length, x, z) + tau * kappa) / (x.Length + 1);
                else mu = BlasLike.ddotvec(x.Length, x, z) / x.Length;
            }
            else if (optMode == "SOCP")
            {
                int R1 = 0;
                for (int icone = 0; icone < numberOfCones; ++icone)
                {
                    R1 += Math.Min(1, cone[icone]);
                }
                if (homogenous) mu = (BlasLike.ddotvec(x.Length, x, z) + tau * kappa) / (R1 + 1);
                else mu = BlasLike.ddotvec(x.Length, x, z) / R1;
            }
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
        long timebase = -1233456;
        long timeaquired = -12233;
        long clocker(bool start = false)
        {
            long t = DateTimeOffset.UtcNow.ToUnixTimeMilliseconds();
            if (start)
                timeaquired = 0;
            else
                timeaquired += (t - timebase);
            timebase = t;
            return timeaquired;
        }

        public static int Opt(int n, int m, double[] w, double[] A, double[] b, double[] c, int nh = 0, double[] H = null, string mode = "QP", int[] cone = null, int[] typecone = null, bool homogenous = true)
        {
            var opt = new Optimise(n, m, w, A, b, c, nh, H);
            opt.optMode = mode;
            if (mode == "SOCP")
            {
                opt.cone = cone;
                opt.typecone = typecone;
                opt.numberOfCones = cone.Length;
                int ncheck = 0;
                foreach (int icc in opt.cone) ncheck += icc;
                Debug.Assert(ncheck == n);
                opt.xbar = new double[n];
                opt.zbar = new double[n];
                opt.dxbar = new double[n];
                opt.dzbar = new double[n];
                opt.W = new double[n];
                opt.W2 = new double[n];
                opt.THETA = new double[cone.Length];
            }
            opt.clocker(true);
            opt.homogenous = homogenous;
            opt.tau = 1;
            opt.kappa = 1;
            opt.usrH = nh > 0 && BlasLike.dsumvec(opt.H.Length, opt.H) != 0.0;
            if (mode == "QP")
            {
                BlasLike.dsetvec(n, 1, opt.x);
                BlasLike.dsetvec(n, 1, opt.z);
            }
            else if (mode == "SOCP")
            {
                for (int icone = 0, istart = 0; icone < opt.cone.Length; istart += opt.cone[icone], ++icone)
                {
                    if (opt.typecone[icone] == (int)conetype.QP)
                    {
                        var nn = opt.cone[icone];
                        BlasLike.dsetvec(nn, 1, opt.x, istart);
                        BlasLike.dsetvec(nn, 1, opt.z, istart);
                    }
                    else if (opt.typecone[icone] == (int)conetype.SOCP)
                    {
                        var nn = opt.cone[icone];
                        BlasLike.dzerovec(nn - 1, opt.x, istart);
                        BlasLike.dzerovec(nn - 1, opt.z, istart);
                        opt.x[nn - 1 + istart] = 1.0;
                        opt.z[nn - 1 + istart] = 1.0;
                    }
                }
            }
            var dxold = (double[])opt.dx.Clone();
            var dzold = (double[])opt.dz.Clone();
            var dyold = (double[])opt.dy.Clone();
            var dtauold = opt.dtau;
            var dkappaold = opt.kappa;
            opt.Mu();
            var mu0 = opt.mu;
            var i = 0;
            var extra = new double[nh];
            if (opt.optMode == "QP")
            {
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
            }
            else if (opt.optMode == "SOCP")
            {
                for (int icone = 0, istart = 0; icone < opt.cone.Length; istart += opt.cone[icone], ++icone)
                {
                    var nn = opt.cone[icone];
                    if (opt.typecone[icone] == (int)conetype.QP)
                    {
                        if (opt.usrH)
                        {
                            if (opt.homogenous)
                            {
                                BlasLike.dcopyvec(nn, opt.x, opt.cmod, istart, istart);
                                BlasLike.dscalvec(nn, 1.0 / opt.tau, opt.cmod, istart);
                                Factorise.dsmxmulv(nh, opt.H, opt.cmod, extra, istart);
                            }
                            else Factorise.dsmxmulv(nh, opt.H, opt.x, extra, istart);
                            BlasLike.dzerovec(nn - nh, opt.cmod, istart + nh);
                            BlasLike.daddvec(nh, opt.c, extra, opt.cmod, istart, 0, istart);
                        }
                        else BlasLike.dcopyvec(nn, opt.c, opt.cmod, istart, istart);
                    }
                    else if (opt.typecone[icone] == (int)conetype.SOCP)
                    {
                        BlasLike.dcopyvec(nn, opt.c, opt.cmod, istart, istart);
                    }
                }
            }
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
                opt.MaximumStep(gamma);
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
                    BlasLike.dzerovec(opt.c.Length - nh, opt.cmod, nh);
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
            Console.WriteLine($"Job took {opt.clocker()} m secs");
            if (i >= opt.maxiter) return -1;
            else if (opt.homogenous && opt.tau < opt.kappa) return 6;
            else return 0;
        }

        //Methods for Second Order Cone Programming operations
        static double root2m1 = Math.Sqrt(0.5);

        void Tmulvec(int ncone, double[] x, int xstart = 0, bool rotate = false)
        {
            int ncm1 = ncone - 1;
            if (rotate)
            {
                double a12 = (-x[ncone - 2 + xstart] + x[ncm1 + xstart]) * root2m1;
                double a11 = (x[ncone - 2 + xstart] + x[ncm1 + xstart]) * root2m1;
                x[ncm1 + xstart] = a11;
                x[ncone - 2 + xstart] = a12;
            }
        }
        void Qmulvec(int ncone, double[] x, int xstart = 0, bool rotate = false)
        {
            int ncm1 = ncone - 1;
            if (!rotate)
            {
                if (ncone > 1)
                {
                    BlasLike.dnegvec(ncm1, x, xstart);
                }
            }
            else
            {
                BlasLike.dnegvec(ncone - 2, x, xstart);
                Ordering.Order.swap(ref x[ncone - 2 + xstart], ref x[ncm1 + xstart]);
            }
        }
        void applyX(int ncone, double[] x, double[] A, double[] XA, int xstart = 0, int Astart = 0, int XAstart = 0)
        {
            int ncm1 = ncone - 1;
            double inner;
            if (ncone == 1)
            {
                XA[XAstart] = A[Astart] * x[xstart];
            }
            else
            {
                inner = BlasLike.ddotvec(ncm1, A, x, Astart, xstart);
                for (int j = 0; j < ncone; ++j)
                {
                    if (j != ncm1) { XA[j + XAstart] = x[xstart + ncm1] * A[j + Astart] + A[Astart + ncm1] * x[j + xstart]; }
                    else { XA[j + XAstart] = inner + x[xstart + ncm1] * A[Astart + ncm1]; }
                }
            }
        }
        void applyXm1(int ncone, double[] x, double[] A, double[] Xm1A, int xstart = 0, int Astart = 0, int Xm1Astart = 0)
        {
            int ncm1 = ncone - 1;
            double outer, inner;
            if (ncone == 1)
            {
                Xm1A[Xm1Astart] = A[Astart] / x[xstart];
            }
            else
            {
                outer = x[ncm1 + xstart] * x[ncm1 + xstart] - BlasLike.ddotvec(ncm1, x, x, xstart, xstart);
                if (outer <= 0)
                {
                    Console.WriteLine($"outer is not positive in applyXm1 {x[ncm1 + xstart]} {outer}");
                    if (x[ncm1] > BlasLike.lm_rooteps)
                    {
                        BlasLike.dscalvec(ncm1, .95 * x[ncm1 + xstart] / Math.Sqrt(x[ncm1 + xstart] * x[ncm1 + xstart] - outer), x, xstart);
                    }
                    else
                    {
                        BlasLike.dsetvec(ncm1, BlasLike.lm_rooteps, x, xstart);
                        x[ncm1 + xstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)ncone);
                    }
                    outer = x[ncm1 + xstart] * x[ncm1 + xstart] - BlasLike.ddotvec(ncm1, x, x, xstart, xstart);
                    Console.WriteLine($"fixed outer is now {outer}");
                }
                inner = BlasLike.ddotvec(ncm1, A, x, Astart, xstart);
                for (int j = 0; j < ncone; ++j)
                {
                    if (j != ncm1) { Xm1A[Xm1Astart + j] = x[j + xstart] * inner / x[ncm1 + xstart] / outer + A[j + Astart] / x[ncm1 + xstart] - A[ncm1 + Astart] * x[j + xstart] / outer; }
                    else { Xm1A[Xm1Astart + j] = (-inner + A[ncm1 + Astart] * x[ncm1 + xstart]) / outer; }
                }
            }
        }
        void thetaScale(int ncone, double[] A, double theta, bool recip = false, bool square = false, int Astart = 0)
        {
            if ((A != null || ncone > 1) && theta != 1)
            {
                if (square)
                    theta *= theta;
                if (recip)
                    theta = 1.0 / theta;
                BlasLike.dscalvec(ncone, theta, A, Astart);
            }
        }
        void Wtrans(int ncone, double[] A, double[] w, double[] WA, int Astart = 0, int wstart = 0, int WAStart = 0)
        {
            int ncm1 = ncone - 1;
            double wc, bot;
            if (ncone == 1)
            {
                WA[WAStart] = w[wstart] * A[Astart];
            }
            else
            {
                wc = BlasLike.ddotvec(ncm1, w, A, wstart, Astart);
                bot = 1.0 + w[ncm1 + wstart];
                if (bot < BlasLike.lm_eps)
                    Console.WriteLine($"bad cone in Wtrans {bot}");
                /*   for (int j = 0; j < ncm1; ++j)
                               {
                                   if (wc != 0) WA[j + WAStart] = w[j + wstart] * wc / bot;
                                   else WA[j + WAStart] = 0;
                                   WA[j + WAStart] += A[j + Astart] + w[j + wstart] * A[ncm1 + Astart];
                               }*/
                BlasLike.dsccopyvec(ncm1, wc / bot, w, WA, wstart, WAStart);
                BlasLike.daxpyvec(ncm1, A[ncm1 + Astart], w, WA, wstart, WAStart);
                BlasLike.daddvec(ncm1, WA, A, WA, WAStart, Astart, WAStart);

                WA[ncm1 + WAStart] = wc + w[ncm1 + wstart] * A[ncm1 + Astart];
            }
        }
        void WtransR(int ncone, double[] A, double[] w, double[] WA, int Astart = 0, int wstart = 0, int WAStart = 0)
        {
            int ncm1 = ncone - 1;
            double wc, bot;
            Tmulvec(ncone, A, Astart, true);
            wc = BlasLike.ddotvec(ncm1, w, A, wstart, Astart);
            bot = 1 + w[ncm1 + wstart];
            if (bot < BlasLike.lm_eps)
                Console.WriteLine($"bad cone in WtransR {bot}");
            /*  for (int j = 0; j < ncm1; ++j)
              {
                  if (wc != 0) WA[j + WAStart] = w[j + wstart] * wc / bot;
                  else WA[j + WAStart] = 0;
                  WA[j + WAStart] += A[j + Astart] + w[j + wstart] * A[ncm1 + Astart];
              }*/
            BlasLike.dsccopyvec(ncm1, wc / bot, w, WA, wstart, WAStart);
            BlasLike.daxpyvec(ncm1, A[ncm1 + Astart], w, WA, wstart, WAStart);
            BlasLike.daddvec(ncm1, WA, A, WA, WAStart, Astart, WAStart);
            WA[ncm1 + WAStart] = wc + w[ncm1 + wstart] * A[ncm1 + Astart];
            Tmulvec(ncone, A, Astart, true);
            Tmulvec(ncone, WA, WAStart, true);
        }
        void W2trans(int ncone, double[] A, double[] w, double[] W2A, int Astart = 0, int wstart = 0, int W2Astart = 0)
        {
            int ncm1 = ncone - 1;
            double wc;
            if (ncone == 1)
            {
                if (w[wstart] != 0 && A[Astart] != 0) W2A[W2Astart] = w[wstart] * w[wstart] * A[Astart];
                else W2A[W2Astart] = 0;
            }
            else
            {
                wc = BlasLike.ddotvec(ncone, w, A, wstart, Astart);
                if (wc != 0.0)
                {
                    /*   for (int j = 0; j < ncm1; ++j)
                       {
                           W2A[j + W2Astart] = 2 * w[j + wstart] * wc + A[j + Astart];
                       }*/
                    BlasLike.dsccopyvec(ncm1, 2 * wc, w, W2A, wstart, W2Astart);
                    BlasLike.daddvec(ncm1, W2A, A, W2A, W2Astart, Astart, W2Astart);
                    W2A[ncm1 + W2Astart] = 2 * w[ncm1 + wstart] * wc - A[ncm1 + Astart];
                }
                else
                {
                    BlasLike.dcopyvec(ncm1, A, W2A, Astart, W2Astart);
                    W2A[ncm1 + W2Astart] = -A[ncm1 + Astart];
                }
            }
        }
        void W2transR(int ncone, double[] A, double[] w, double[] W2A, int Astart = 0, int wstart = 0, int W2Astart = 0)
        {
            int ncm1 = ncone - 1;
            double wc;
            Tmulvec(ncone, A, Astart, true);
            wc = BlasLike.ddotvec(ncone, w, A, wstart, Astart);
            if (wc != 0.0)
            {
                /*   for (int j = 0; j < ncm1; ++j)
                           {
                               W2A[j + W2Astart] = 2 * w[j + wstart] * wc + A[j + Astart];
                           }*/
                BlasLike.dsccopyvec(ncm1, 2 * wc, w, W2A, wstart, W2Astart);
                BlasLike.daddvec(ncm1, W2A, A, W2A, W2Astart, Astart, W2Astart);
                W2A[ncm1 + W2Astart] = 2 * w[ncm1 + wstart] * wc - A[ncm1 + Astart];
            }
            else
            {
                BlasLike.dcopyvec(ncm1, A, W2A, Astart, W2Astart);
                W2A[ncm1 + W2Astart] = -A[ncm1 + Astart];
            }
            Tmulvec(ncone, A, Astart, true);
            Tmulvec(ncone, W2A, W2Astart, true);
        }
    }
}
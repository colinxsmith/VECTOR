using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace InteriorPoint
{
    public enum conetype { QP, SOCP, SOCPR };
    class BestResults
    {
        public double[] x;
        public double[] y;
        public double[] z;
        public double tau;
        public double kappa;
        public double norm;
        public BestResults(double[] x = null, double[] y = null, double[] z = null, double tau = -1, double kappa = -1, double rp = BlasLike.lm_max, double rd = BlasLike.lm_max, double comp = BlasLike.lm_max)
        {
            set(x, y, z, tau, kappa, rp, rd, comp);
        }
        public void set(double[] x, double[] y, double[] z, double tau, double kappa, double rp, double rd, double comp)
        {
            this.x = x != null ? (double[])x.Clone() : null;
            this.y = y == null ? null : (double[])y.Clone();
            this.z = z == null ? null : (double[])z.Clone();
            this.tau = tau;
            this.kappa = kappa;
            double[] pass = { rp, rd, 1e4 * comp };
            this.norm = Optimise.lInfinity(pass);
        }
        public void update(double[] x, double[] y, double[] z, double tau, double kappa, double rp, double rd, double comp)
        {
            double[] pass = { rp, rd, 1e4 * comp };
            double norm = Optimise.lInfinity(pass);
            if (norm < this.norm)
                set(x, y, z, tau, kappa, rp, rd, comp);
        }
    }
    public class Optimise
    {
        BestResults keep;
        bool copyKept = false;
        double alphamin = 1e-1;
        double conv = BlasLike.lm_eps;
        double compConv = BlasLike.lm_eps;
        int badindex = -1;
        string optMode = "QP";
        int numberOfCones = 0;
        int[] cone = null;
        int[] typecone = null;
        bool usrH = false;
        bool homogenous = false;
        double mu;
        int maxouter = 1000;
        int maxinner = 100;
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
        double laststep = 0;
        double[] lastdx = null;
        double[] lastdy = null;
        double[] lastdz = null;
        double lastdtau;
        double lastdkappa;
        double condition;
        double regularise;
        static double denomTest(double x) => x == 0 ? 1 : x;
        public static double lInfinity(double[] x)
        {
            var back = 0.0;
            foreach (var k in x)
            {
                back = Math.Max(back, Math.Abs(k));
            }
            return back;
        }
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
            keep = new BestResults();
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
        void MaximumStep(double gamma = 0)
        {
            if (optMode == "SOCP")
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
                double lowest = 1e-1, lowest1 = 1 - lowest;

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
                    else if ((desc >= -BlasLike.lm_eps) && (XdX < 0))
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
                    else if ((desc >= -BlasLike.lm_eps * 8) && (SdS < 0))
                        alpha = Math.Min(alpha, lowest1 * ((-SdS) / dSdS));
                }

                if (Math.Abs(dXdS) <= BlasLike.lm_eps)
                {
                    if (dXS + XdS < -BlasLike.lm_eps)
                        alpha = Math.Min(alpha, lowest1 * (-XS / (XdS + dXS)));
                }
                else
                {
                    if (Math.Abs(XdS + dXS) > BlasLike.lm_eps) desc = 1.0 - 4.0 * XS * dXdS / (XdS + dXS) / (XdS + dXS);
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
                    else if ((desc >= -BlasLike.lm_eps * 8) && ((XdS + dXS) / dXdS < 0))
                        alpha = Math.Min(alpha, lowest1 * ((-XdS - dXS) / 2.0 / dXdS));
                }

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
                        //      for (var i = icone; i < cone.Length; ++i)
                        {
                            var i = icone;
                            if (vz1[i] + alpha * (vz2[i] + alpha * vz3[i]) < -BlasLike.lm_eps * 8)
                            {
                                if (Math.Abs(vz3[i]) <= BlasLike.lm_eps)
                                {
                                    if (vz2[i] < -BlasLike.lm_eps)
                                        alpha = Math.Min(lowest1 * -vz1[i] / vz2[i], alpha);
                                }
                                else if (Math.Abs(vz2[i]) > BlasLike.lm_eps && (inner = 1.0 - 4.0 * vz3[i] * vz1[i] / vz2[i] / vz2[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > (BlasLike.lm_eps * 8) ? Math.Sqrt(inner) * Math.Abs(vz2[i]) : 0);
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
                                else if (Math.Abs(vz2[i]) <= BlasLike.lm_eps && (inner = -4.0 * vz3[i] * vz1[i]) > -BlasLike.lm_eps)
                                {
                                    inner = (inner > (BlasLike.lm_eps * 8) ? Math.Sqrt(inner) : 0);
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

                        //                       for (var i = 0; i < cone.Length; ++i)
                        {
                            var i = icone;
                            if (vx1[i] + alpha * (vx2[i] + alpha * vx3[i]) < -BlasLike.lm_eps * 8)
                            {
                                if (Math.Abs(vx3[i]) <= BlasLike.lm_eps * 8)
                                {
                                    if (vx2[i] < -BlasLike.lm_eps * 8)
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
                    }

                    double rhs, gamma1 = 1 - gamma, test1, test2 = 1, beta = 1e-8;
                    for (var l = 0; l < 1000; ++l)
                    {
                        rhs = beta * (1 - alpha * gamma1) * mu;
                        if (homogenous) test2 = (tau + alpha * dtau) * (kappa + alpha * dkappa);
                        test1 = 0;
                        for (var i = 0; i < cone.Length; ++i)
                        {
                            test1 += (vx1[i] + alpha * (vx2[i] + alpha * vx3[i])) * (vz1[i] + alpha * (vz2[i] + alpha * vz3[i]));
                        }

                        test1 = Math.Sqrt(test1);
                        if (test1 >= rhs && test2 >= rhs) break;
                        alpha *= lowest1;
                    }

                    ddx = ddz = dd = alpha;

                }
                if (homogenous)
                {
                    if (dtau < 0) alpha = Math.Min(alpha, -aob(tau, dtau) * lowest1);
                    if (dkappa < 0) alpha = Math.Min(alpha, -aob(kappa, dkappa) * lowest1);

                    ddx = ddz = dd = alpha;
                }
            }

            else if (optMode == "QP")
            {
                double lowest = 5e-2, lowest1 = 1 - lowest;
                ddx = 1.0;
                ddz = 1.0;
                dd = 1.0;
                for (int i = 0; i < n; ++i)
                {
                    if (dx[i] < 0) ddx = lowest1 * Math.Min(ddx, -aob(x[i], dx[i]));
                    if (dz[i] < 0) ddz = lowest1 * Math.Min(ddz, -aob(z[i], dz[i]));
                }
                if (homogenous)
                {
                    if (dtau < 0) dd = lowest1 * Math.Min(dd, -aob(tau, dtau));
                    if (dkappa < 0) dd = lowest1 * Math.Min(dd, -aob(kappa, dkappa));
                }
            }
        }
        void CreateNormalMatrix()
        {
            if (optMode == "QP")
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
            }
            else if (optMode == "SOCP")
            {
                BlasLike.dzerovec(M.Length, M);
                var lhs = w1;//Highjack w1 and dy to save reallocation
                var res = dy;
                HCOPY = nh > 0 ? (double[])H.Clone() : null;
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
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        for (int i = 0, ij = 0; i < m; ++i, ij += i)
                        {
                            BlasLike.dcopy(n, A, m, lhs, 1, i + cstart * m, x.Length + cstart);
                            W2m1trans(n, lhs, W, lhs, x.Length + cstart, cstart, cstart);
                            thetaScale(n, lhs, THETA[icone], true, true, cstart);
                            for (var k = cstart; k < n + cstart; ++k)
                            {
                                if (lhs[k] != 0.0)
                                {
                                    BlasLike.daxpyvec(i + 1, lhs[k], A, M, k * m, ij);
                                }
                            }
                        }
                    }
                }
            }
        }
        void SolvePrimaryDual(double gamma = 0.0, bool corrector = false)
        {
            if (optMode == "QP")
            {
                var g1 = 1.0 - gamma;
                if (!corrector)
                {
                    if (m != 1) badindex = Factorise.Factor(uplo, m, M, order);
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
                    var lhs = w1;//Highjack w1 and dy to save reallocation
                    var res = dy;
                    for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                    {
                        var n = cone[icone];
                        if (typecone[icone] == (int)conetype.QP)
                        {
                            for (var i = cstart; i < n + cstart; ++i)
                            {
                                w1[i] = rd[i] * g1 - aob(rmu[i] - g1 * mu, xbar[i]) * W[i];
                            }
                        }
                        else if (typecone[icone] == (int)conetype.SOCP)
                        {
                            BlasLike.dcopyvec(n, rmu, lhs, cstart, cstart);
                            lhs[n - 1 + cstart] -= g1 * mu;
                            applyXm1(n, xbar, lhs, lhs, cstart, cstart, x.Length + cstart);
                            Wtrans(n, lhs, W, w1, x.Length + cstart, cstart, cstart);
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            BlasLike.daxpyvec(n, -g1, rd, w1, cstart, cstart);
                            BlasLike.dnegvec(n, w1, cstart);
                        }
                    }
                    if (m != 1) badindex = Factorise.Factor(uplo, m, M, order);
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
                            BlasLike.dcopyvec(n, rmu, lhs, cstart, cstart);
                            lhs[n - 1 + cstart] -= g1 * mu;
                            applyX(n, dx0, dz0, lhs, cstart, cstart, x.Length + cstart);
                            BlasLike.dsubvec(n, lhs, lhs, lhs, cstart, x.Length + cstart, cstart);

                            applyXm1(n, xbar, lhs, lhs, cstart, cstart, x.Length + cstart);
                            Wtrans(n, lhs, W, w1, x.Length + cstart, cstart, cstart);
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            BlasLike.daxpyvec(n, -g1, rd, w1, cstart, cstart);
                            BlasLike.dnegvec(n, w1, cstart);
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
                        W2m1trans(n, w1, W, w1, cstart, cstart, x.Length + cstart);
                        BlasLike.dcopyvec(n, w1, w1, x.Length + cstart, cstart);
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
                            W2m1trans(n, cx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, cx, x.Length + cstart, cstart);
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
                            W2m1trans(n, dx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dx, x.Length + cstart, cstart);
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
                            W2m1trans(n, dc, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dc, x.Length + cstart, cstart);
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
                            double[] dxzbar = new double[n];
                            Wtrans(n, dx, W, w1, cstart, cstart, cstart);//times W
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            applyX(n, zbar, w1, dxzbar, cstart, cstart, 0);

                            BlasLike.dcopyvec(n, rmu, w1, cstart, cstart);
                            w1[n - 1 + cstart] -= g1 * mu;
                            if (corrector)
                            {
                                applyX(n, dx0, dz0, w1, cstart, cstart, x.Length + cstart);
                                BlasLike.dsubvec(n, w1, w1, w1, cstart, x.Length + cstart, cstart);
                            }
                            BlasLike.dsubvec(n, w1, dxzbar, w1, cstart, 0, cstart);
                            applyXm1(n, xbar, w1, w1, cstart, cstart, x.Length + cstart);//over xbar
                            Wtrans(n, w1, W, dz, x.Length + cstart, cstart, cstart);
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
                            W2m1trans(n, dx, W, w1, cstart, cstart, x.Length + cstart);
                            BlasLike.dcopyvec(n, w1, dx, x.Length + cstart, cstart);
                            thetaScale(n, dx, THETA[icone], true, true, cstart);
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
                            double[] dxzbar = new double[n];
                            Wtrans(n, dx, W, w1, cstart, cstart, cstart);
                            thetaScale(n, w1, THETA[icone], false, false, cstart);
                            applyX(n, w1, zbar, dxzbar, cstart, cstart, 0);

                            BlasLike.dcopyvec(n, rmu, w1, cstart, cstart);
                            w1[n - 1 + cstart] -= g1 * mu;
                            if (corrector)
                            {
                                applyX(n, dx0, dz0, w1, cstart, cstart, x.Length + cstart);
                                BlasLike.dsubvec(n, w1, w1, w1, cstart, x.Length + cstart, cstart);
                            }
                            BlasLike.dsubvec(n, w1, dxzbar, w1, cstart, 0, cstart);
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
                bool zcopy = false, xcopy = false;
                for (int icone = 0, cstart = 0; icone < cone.Length; cstart += cone[icone], icone++)
                {
                    var n = cone[icone];
                    if (typecone[icone] == (int)conetype.QP)
                    {
                        for (int i = cstart; i < n + cstart; ++i)
                        {
                            W2[i] = aob(z[i], x[i]);
                            W[i] = Math.Sqrt(W2[i]);
                            THETA[icone] = 1;
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
                            xcopy = false;
                            if (xQx <= BlasLike.lm_eps)
                            {
                                xcopy = true;
                                if (x[n - 1 + cstart] > BlasLike.lm_rooteps)
                                {
                                    BlasLike.dscalvec(n - 1, .99 * x[n - 1 + cstart] / Math.Sqrt(x[n - 1 + cstart] * x[n - 1 + cstart] - xQx), x);
                                    xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                }
                                else
                                {
                                    BlasLike.dsetvec(n - 1, BlasLike.lm_rooteps, x, cstart);
                                    x[n - 1 + cstart] = (1.0 + BlasLike.lm_rooteps) * BlasLike.lm_rooteps * Math.Sqrt((double)n);
                                    xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                }
                            }
                            zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                            zcopy = false;
                            if (zQz <= BlasLike.lm_eps)
                            {
                                zcopy = true;
                                if (z[n - 1 + cstart] > BlasLike.lm_eps * 0)
                                {
                                    BlasLike.dscalvec(n - 1, .99 * z[n - 1 + cstart] / Math.Sqrt(z[n - 1 + cstart] * z[n - 1 + cstart] - zQz), z);
                                    zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
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



                            if (!double.IsNaN(THETA[icone]) && THETA[icone] >= BlasLike.lm_eps)
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
                                Debug.Assert(false);
                                //DO something to stop main loop
                            }
                            double wcheck = BlasLike.ddotvec(n, W, W, cstart, cstart);
                            if (double.IsNaN(wcheck))
                            {
                                Console.WriteLine("BAD W"); Debug.Assert(false);
                            }
                            Wtrans(n, x, W, xbar, cstart, cstart, cstart); //xbar=thetaW.x
                            thetaScale(n, xbar, THETA[icone], false, false, cstart);
                            Wm1trans(n, z, W, zbar, cstart, cstart, cstart);//z=(thetaW)(thetaW)x
                            thetaScale(n, zbar, THETA[icone], true, false, cstart);//zbar=(Wtheta)m1.z=xbar
                            Tmulvec(n, xbar, cstart);//Tmulvec does nothing for SOCP, needed for SOCPR
                            Tmulvec(n, zbar, cstart);
                            /*       if (zcopy)
                                   {
                                       BlasLike.dcopyvec(n, xbar, zbar);
                                   }
                                   else if (xcopy)
                                   {
                                       BlasLike.dcopyvec(n, zbar, xbar);
                                   }*/
                            applyX(n, xbar, zbar, rmu, cstart, cstart, cstart);
                            Tmulvec(n, rmu, cstart);
                        }

                        BlasLike.dnegvec(n, rmu, cstart);
                        rmu[n - 1 + cstart] += mu;
                    }
                }
            }
            if (homogenous)
            {
                hrmu = mu - tau * kappa;
                rkxy = kappa + BlasLike.ddotvec(x.Length, cmod, x) - BlasLike.ddotvec(y.Length, b, y);
            }
            CreateNormalMatrix();
        }
        void ConditionEstimate()
        {
            if (m > 1)
            {
                var diags = new double[m];
                var order = new int[m];
                for (var i = 0; i < m; ++i)
                {
                    diags[i] = M[i * (i + 3) / 2];
                }
                Ordering.Order.getorder(m, diags, order, null, 0, 1);
                int o1 = Math.Max(order[0], order[1]), o2 = Math.Min(order[0], order[1]);
                double a1 = Math.Max(diags[o1], diags[o2]), a2 = Math.Min(diags[o2], diags[o1]), a12 = M[o1 * (o1 + 1) / 2 + o2];
                condition = a1 * (a1 - a12 * a12 / a2);//cond is a quick estimate of condition number using only 2 pivots.
                regularise = a1 * BlasLike.lm_rooteps;
            }
            else
            {
                condition = square(M[0]);
                regularise = M[0] * BlasLike.lm_rooteps;
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
        void update(double[] dx, double[] dy, double[] dz, double dtau, double dkappa, double step, double test = 0)
        {
            if (test == 0)
            {
                laststep = step;
                lastdx = (double[])dx.Clone();
                lastdy = (double[])dy.Clone();
                lastdz = (double[])dz.Clone();
                lastdtau = dtau;
                lastdkappa = dkappa;
            }
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
        double Gap()
        {
            return Math.Abs(Primal() - Dual() - (homogenous ? kappa : 0));
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
            var stepReduce = 1;
            opt.optMode = mode;
            if (mode == "SOCP")
            {
                opt.conv = (Math.Floor(1e-8 / BlasLike.lm_eps)) * BlasLike.lm_eps;
                opt.compConv = (Math.Floor(1e-11 / BlasLike.lm_eps)) * BlasLike.lm_eps;
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
            var ir = 1;
            var extra = new double[nh];
            if (opt.optMode == "QP")
            {
                opt.conv = (Math.Floor(1e-15 / BlasLike.lm_eps)) * BlasLike.lm_eps;
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
            var rp0 = lInfinity(opt.rp);
            var rd0 = lInfinity(opt.rd);
            var gap0 = denomTest(opt.Gap());
            rp0 = denomTest(rp0);
            rd0 = denomTest(rd0);
            var rp1 = rp0;
            var rd1 = rd0;
            var gap1 = gap0;
            var comp0 = opt.Complementarity();
            opt.keep.set(opt.x, opt.y, opt.z, opt.tau, opt.kappa, rp0, rd0, comp0);
            var compnow = comp0;
            var comp1 = comp0;
            var alpha1 = 0.0;
            var alpha2 = 0.0;
            var gamma = 0.0;
            var diff = new double[n];
            var gap = 1.0;
            while (true)
            {
                rp1 = lInfinity(opt.rp) / denomTest(rp0);
                rd1 = lInfinity(opt.rd) / denomTest(rd0);
                gap = opt.Gap();
                gap1 = gap / denomTest(gap0);
                comp1 = opt.Complementarity();
                opt.keep.update(opt.x, opt.y, opt.z, opt.tau, opt.kappa, lInfinity(opt.rp), lInfinity(opt.rd), comp1);
                if (/*Math.Abs(gap / opt.tau) < opt.conv && */ rp1 < opt.conv && rd1 < opt.conv && comp1 < opt.compConv /* * square(opt.tau)*/)
                    break;
                if (comp1 < opt.compConv && opt.tau < 1e-5 * opt.kappa) break;
                if (ir > opt.maxouter) break;
                if (i > opt.maxinner)
                {
                    ir++; i = 0;
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    opt.kappa /= opt.tau;
                    opt.tau = 1;
                    rp0 = denomTest(lInfinity(opt.rp));
                    rd0 = denomTest(lInfinity(opt.rd));
                    gap0 = denomTest(opt.Gap());
                };
                if (i > 2 && opt.optMode == "SOCP")
                {
                    BlasLike.dsubvec(n, opt.xbar, opt.zbar, diff);
                    double test = BlasLike.ddotvec(n, diff, diff);
                    if (test > BlasLike.lm_eps * 256)
                    {
                        Console.WriteLine($"xbar test = {test}");
                        Console.WriteLine($"rp1 = {rp1}");
                        Console.WriteLine($"rd1 = {rd1}");
                        Console.WriteLine($"comp1 = {comp1}");
                        //   break;
                    }
                }
                opt.SolvePrimaryDual();
                if (opt.badindex != 0) Console.WriteLine($"Normal matrix is unstable: badindex={opt.badindex}");
                BlasLike.dcopyvec(n, opt.dx, dxold);
                BlasLike.dcopyvec(n, opt.dz, dzold);
                BlasLike.dcopyvec(m, opt.dy, dyold);
                dtauold = opt.dtau;
                dkappaold = opt.dkappa;
                opt.MaximumStep();
                alpha1 = stepReduce * opt.Lowest();
                BlasLike.dsccopyvec(opt.n, 1.0, opt.dx, opt.dx0);//was alpha1
                BlasLike.dsccopyvec(opt.n, 1.0, opt.dz, opt.dz0);//was alpha1
                opt.dtau0 = 1.0 * opt.dtau;//was alpha1
                opt.dkappa0 = 1.0 * opt.dkappa;//was alpha1
                gamma = opt.gfunc(alpha1);
                opt.SolvePrimaryDual(gamma, true);
                opt.MaximumStep(gamma);
                alpha2 = stepReduce * opt.Lowest();
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
                var t1 = 0.0;
                if (((t1 = Math.Max(alpha1, alpha2)) < opt.alphamin))
                {
                    var scl = 1.0;
                    opt.update(opt.lastdx, opt.lastdy, opt.lastdz, opt.lastdtau, opt.lastdkappa, -opt.laststep, 1);
                    opt.update(opt.lastdx, opt.lastdy, opt.lastdz, opt.lastdtau, opt.lastdkappa, 0.9 * opt.laststep, 1);
                    BlasLike.dscalvec(opt.y.Length, scl / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, scl / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, scl / opt.tau, opt.z);
                    gap = opt.Primal() - opt.Dual();
                    if (gap < 0)
                    {
                        var dgap = BlasLike.ddotvec(opt.y.Length, opt.y, opt.b);
                        BlasLike.dsetvec(opt.y.Length, 1.0, opt.y);
                        dgap -= BlasLike.ddotvec(opt.y.Length, opt.y, opt.b);
                        var step = -gap / dgap / opt.tau;
                        if (dgap != 0) BlasLike.dsetvec(opt.y.Length, (1.0 + BlasLike.lm_rooteps) * step, opt.y);
                    }

                    opt.Mu();
                    mu0 = opt.mu;
                    opt.PrimalResidual();
                    opt.DualResudual();
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    opt.kappa = opt.mu; //*=  scl / opt.tau;
                    opt.tau = scl;
                    i = 0; ir++;
                    //    opt.conv *= 1.01;
                    rp0 = denomTest(lInfinity(opt.rp));
                    rd0 = denomTest(lInfinity(opt.rd));
                    gap0 = denomTest(opt.Gap());
                }
                opt.Mu();
                mu0 = opt.mu;
                opt.PrimalResidual();
                opt.DualResudual();
                opt.MuResidual();
                opt.ConditionEstimate();
                if (false && opt.condition > BlasLike.lm_reps)
                {
                    opt.update(opt.lastdx, opt.lastdy, opt.lastdz, opt.lastdtau, opt.lastdkappa, -opt.laststep, 1);
                    opt.update(opt.lastdx, opt.lastdy, opt.lastdz, opt.lastdtau, opt.lastdkappa, 0.5 * opt.laststep, 1);
                    /*    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                                        BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                                        BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                                        opt.kappa /= opt.tau;
                                        opt.tau = 1;*/
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    i = 0;
                    rp0 = denomTest(lInfinity(opt.rp));
                    rd0 = denomTest(lInfinity(opt.rd));
                    gap0 = denomTest(opt.Gap());
                }
                if (false && opt.condition > BlasLike.lm_reps)
                {
                    for (int ii = 0, id = 0; ii < m; ++ii, id += ii)
                    {
                        if (opt.M[id + ii] < opt.regularise)
                            opt.M[id + ii] += opt.regularise;
                    }
                    i = 0;
                    rp0 = denomTest(lInfinity(opt.rp));
                    rd0 = denomTest(lInfinity(opt.rd));
                    gap0 = denomTest(opt.Gap());
                }
                gap = opt.Primal() - opt.Dual();
                if (opt.tau < 1e-5 && opt.kappa < 1e-5)
                {
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    opt.kappa /= opt.tau;
                    opt.tau = 1;
                    ir++; i = 0;
                    opt.Mu();
                    mu0 = opt.mu;
                    opt.PrimalResidual();
                    opt.DualResudual();
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    rp0 = denomTest(lInfinity(opt.rp));
                    rd0 = denomTest(lInfinity(opt.rd));
                    gap0 = denomTest(opt.Gap());
                }
                i++;
            }
            if (opt.copyKept)
            {
                BlasLike.dcopyvec(opt.x.Length, opt.keep.x, opt.x);
                BlasLike.dcopyvec(opt.y.Length, opt.keep.y, opt.y);
                BlasLike.dcopyvec(opt.z.Length, opt.keep.z, opt.z);
                opt.tau = opt.keep.tau;
                opt.kappa = opt.keep.kappa;
            }
            var infease = !(opt.homogenous && (opt.tau > 1e2 * opt.kappa));
            if (opt.homogenous)
            {
                Console.WriteLine($"tau = {opt.tau} kappa={opt.kappa}");
                if (opt.tau != 0)
                {
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                }
                if (infease) Console.WriteLine("INFEASIBLE");
            }
            Console.WriteLine($"{ir} outer iterations out of {opt.maxouter}");
            Console.WriteLine($"{i} iterations out of {opt.maxinner}");
            Console.WriteLine($"Relative Primal Residual\t\t {rp1}");
            Console.WriteLine($"Relative Dual Residual\t\t\t {rd1}");
            Console.WriteLine($"Relative Complementarity\t\t {comp1}");
            Console.WriteLine($"Primal Utility:\t\t{opt.Primal()}");
            ActiveSet.Optimise.printV("x", opt.x);
            Console.WriteLine($"Dual Utility:\t\t{opt.Dual()}");
            ActiveSet.Optimise.printV("y", opt.y);
            ActiveSet.Optimise.printV("z", opt.z);
            Console.WriteLine($"Complementarity:\t{BlasLike.ddotvec(opt.n, opt.x, opt.z)}");
            Console.WriteLine($"Gap:\t\t\t{opt.Primal() - opt.Dual()}");
            Console.WriteLine($"Job took {opt.clocker()} m secs");
            Console.WriteLine($"Last conv {opt.conv}");

            if (i >= opt.maxinner || ir >= opt.maxouter) return -1;
            else if (opt.homogenous && infease) return 6;
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
        void Wm1trans(int ncone, double[] A, double[] w, double[] WA, int Astart = 0, int wstart = 0, int WAStart = 0)
        {
            Qmulvec(ncone, A, Astart);//A is changed to QA
            Wtrans(ncone, A, w, WA, Astart, wstart, WAStart); //WQA
            Qmulvec(ncone, WA, WAStart);// QWQA = W-1A since QWQ = inverse(W)
            Qmulvec(ncone, A, Astart);//Change QA back to A, since QQ=I (avoid having to copy A at start)
        }
        void W2m1trans(int ncone, double[] A, double[] w, double[] WA, int Astart = 0, int wstart = 0, int WAStart = 0)
        {
            Qmulvec(ncone, A, Astart);//A is changed to QA
            W2trans(ncone, A, w, WA, Astart, wstart, WAStart); //WWQA
            Qmulvec(ncone, WA, WAStart);// QWWQA = WW-1A since QWWQ = inverse(WW)
            Qmulvec(ncone, A, Astart);//Change QA back to A, since QQ=I (avoid having to copy A at start)
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
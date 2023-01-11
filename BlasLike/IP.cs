using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace InteriorPoint
{
    public delegate void hessmull(int n, double[] hess, double[] wrk, double[] hx,int hstart=0);

    public enum conetype { QP, SOCP, SOCPR };
    class BestResults
    {
        public double[] x;
        public double[] y;
        public double[] z;
        public double tau;
        public double kappa;
        public double norm;
        public BestResults(double[] x = null, double[] y = null, double[] z = null, double tau = -1, double kappa = -1, double rp = BlasLike.lm_max, double rd = BlasLike.lm_max, double comp = BlasLike.lm_max, double norm = BlasLike.lm_max)
        {
            set(x, y, z, tau, kappa, rp, rd, comp, norm);
        }
        public void set(double[] x, double[] y, double[] z, double tau, double kappa, double rp, double rd, double comp, double norm)
        {
            this.x = x != null ? (double[])x.Clone() : null;
            this.y = y == null ? null : (double[])y.Clone();
            this.z = z == null ? null : (double[])z.Clone();
            this.tau = tau;
            this.kappa = kappa;
            this.norm = norm;
        }
        public void update(double[] x, double[] y, double[] z, double tau, double kappa, double rp, double rd, double comp)
        {
            double[] pass = { rp, rd, comp };
            double norm = Optimise.lInfinity(pass);
            if (norm < this.norm)
                set(x, y, z, tau, kappa, rp, rd, comp, norm);
        }
    }
    public class Optimise
    {
        public double[] baseA = null;
        double alpha1;
        double alpha2;
        int innerIteration;
        public int basen = 0;
        public int basem = 0;
        public int bases = 0;
        public hessmull h = null;
        BestResults keep;
        public bool copyKept = true;
        public double alphamin = 1e-1;
        public double alphamin2 = 1e-5;
        public bool copyBest = true;
        public double conv = BlasLike.lm_eps * 16;
        public double compConv = BlasLike.lm_eps * 16;
        int badindex = -1;
        string optMode = "QP";
        int numberOfCones = 0;
        int[] cone = null;
        int[] typecone = null;
        bool usrH = false;
        bool specialDebug = false;
        bool homogenous = false;
        double mu;
        public int maxouter = 5;
        public int maxinner = 50;
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
        public double[] y = null;
        double[] dy = null;
        public double[] z = null;
        double[] dz = null;
        ///<summary>H is the symmetric second order matrix. 
        ///For an optimisation that uses a factor model, this is 
        ///supplied separately and not what's used in hessmull.</summary>
        public double[] H = null;
        double[] xbar = null;
        double[] zbar = null;
        int[] sign = null;
        double[] dxbar = null;
        double[] dzbar = null;
        public int nh;
        char uplo = 'U';
        double[] HCOPY = null;
        //    double[] HHCopy = null;
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
        public int[] slackToConstraintBOTH = null;
        public int[] slackToConstraintBOTH_inverse = null;
        int slackmboth = 0;
        int slackmL = 0;
        int slackmU = 0;
        int slacklarge;
        public int[] slacklargeConstraintToStock = null;
        public int[] slacklargeConstraintToStock_inverse = null;
        public int[] slackToConstraintL = null;
        public int[] slackToConstraintL_inverse = null;
        public int[] slackToConstraintU = null;
        public int[] slackToConstraintU_inverse = null;
        static double denomTest(double x) => x * x <= 1 ? 1 : x;
        public double[] conetop(double[] resid)
        {
            if (true || resid.Length != cone.Length) return resid;
            else
            {
                var back = new double[cone.Length];
                int cstart, i;
                for (i = 0, cstart = 0; i < cone.Length; cstart += cone[i], i++)
                {
                    if (cone[i] == 1)
                        back[i] = resid[cone[i] + cstart - 1];
                    else
                    {
                        back[i] = resid[cone[i] + cstart - 1];
                    }
                }
                return back;
            }
        }
        public static double LNorm(double[] x)
        {
            return Math.Sqrt(BlasLike.ddotvec(x.Length, x, x));
        }
        public static double lInfinity(double[] x)
        {
            var back = 0.0;
            foreach (var k in x)
            {
                back = Math.Max(back, Math.Abs(k));
            }
            return back;
        }
        public void qphess1(int n, double[] hess, double[] wrk, double[] hx,int hstart=0)
        {
            Solver.Factorise.dsmxmulv(n, hess, wrk, hx);
        }
        static double norm(double[] aa) => Math.Sqrt(BlasLike.ddotvec(aa.Length, aa, aa));
        static public double square(double a) => a * a;
        double corrector_gamma = 0.7;
        double gfunc(double a) => Math.Min(corrector_gamma, square(1 - a)) * (1 - a);//Try 0.7 instead of 0.5
        double aob(double a, double b)
        {
            int fail = 21;
            var back = BlasLike.dprotdiv(ref a, ref b, ref fail);
            if (fail != 0 && b == 0) back = BlasLike.lm_max;
            return back;
        }
        public Optimise(int n, int m, double[] x, double[] A, double[] b, double[] c, int nh = 0, double[] H = null)
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
            basen = n;
            bases = 0;
            basem = m;

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
                double alpha = 1.0;
                double lowest = 1e-2, lowest1 = 1 - lowest;
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
                        alpha = Math.Min(ddx, alpha);
                        alpha = Math.Min(ddz, alpha);
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
                        BlasLike.dcopyvec(n, x, Qx, cstart, 0);//Qx
                        Qmulvec(n, Qx);
                        BlasLike.dcopyvec(n, dx, Qdx, cstart, 0);//Qdx
                        Qmulvec(n, Qdx);
                        vx1[icone] = BlasLike.ddotvec(n, x, Qx, cstart, 0);//x.Qx
                        vx2[icone] = 2.0 * BlasLike.ddotvec(n, dx, Qx, cstart, 0);//dx.Qx
                        vx3[icone] = BlasLike.ddotvec(n, dx, Qdx, cstart, 0);//dx.Qdx
                        if (Qdz[n - 1] < 0) alpha = Math.Min(-aob(Qz[n - 1], Qdz[n - 1]), alpha);
                        if (dz[n - 1 + cstart] < 0) alpha = Math.Min(-aob(z[n - 1 + cstart], dz[n - 1 + cstart]), alpha);
                        if (Qdx[n - 1] < 0) alpha = Math.Min(-aob(Qx[n - 1], Qdx[n - 1]), alpha);
                        if (dx[n - 1 + cstart] < 0) alpha = Math.Min(-aob((x[n - 1 + cstart]), dx[n - 1 + cstart]), alpha);
                        double inner;
                        var i = icone;
                        if ((vz1[i] + alpha * (vz2[i] + alpha * vz3[i])) < 0)
                        {
                            if (Math.Abs(vz3[i]) <= 0)
                            {
                                ColourConsole.WriteEmbeddedColourLine($"[yellow]No[/yellow] [cyan]vz3[/cyan] [magenta]{vz3[i]}[/magenta]");
                                if (vz2[i] < 0)
                                    alpha = Math.Min(-aob(vz1[i], vz2[i]), alpha);
                            }
                            else if (Math.Abs(vz2[i]) > 0 && (inner = 1.0 - 4.0 * aob(vz3[i], vz2[i]) * aob(vz1[i], vz2[i])) > 0)
                            {
                                double r1, r2;
                                inner = (inner > 0 ? Math.Sqrt(inner) * Math.Abs(vz2[i]) : 0);
                                r1 = aob((-vz2[i] - inner), vz3[i]) / 2.0; r2 = aob((-vz2[i] + inner), vz3[i]) / 2.0;
                                if (vz3[i] > 0)
                                {
                                    if (r1 < 0) r1 = r2;
                                    alpha = Math.Min(Math.Max(r1, 0), alpha);
                                }
                                else
                                {
                                    if (r2 < 0) r2 = r1;
                                    alpha = Math.Min(Math.Max(r2, 0), alpha);
                                }
                            }
                            else if (Math.Abs(vz2[i]) <= 0 && (inner = -4.0 * vz3[i] * vz1[i]) >= 0)
                            {
                                inner = (inner > 0 ? Math.Sqrt(inner) : 0);
                                double r1, r2;
                                r1 = aob((-vz2[i] - inner), vz3[i]) / 2.0; r2 = aob((-vz2[i] + inner), vz3[i]) / 2.0;
                                if (vz3[i] > 0)
                                {
                                    if (r1 < 0) r1 = r2;
                                    alpha = Math.Min(Math.Max(r1, 0), alpha);
                                }
                                else
                                {
                                    if (r2 < 0) r2 = r1;
                                    alpha = Math.Min(Math.Max(r2, 0), alpha);
                                }
                            }
                            else if ((inner = Math.Abs(1.0 - 4.0 * aob(vz3[i], vz2[i]) * aob(vz1[i], vz2[i]))) < BlasLike.lm_eps)
                            {
                                var r2 = aob((-vz2[i]), vz3[i]) / 2.0;
                                alpha = Math.Min(Math.Max(r2, 0), alpha);
                            }
                            else
                                ColourConsole.WriteEmbeddedColourLine($"[red]ROUNDING cone z{i}[/red][magenta] step {alpha}[/magenta] [green]{vz1[i]}[/green] [yellow]{vz2[i]}[/yellow] [cyan]{vz3[i]}[/cyan]");
                        }
                        if ((vx1[i] + alpha * (vx2[i] + alpha * vx3[i])) < 0)
                        {
                            if (Math.Abs(vx3[i]) <= 0)
                            {
                                ColourConsole.WriteEmbeddedColourLine($"[yellow]No[/yellow] [cyan]vx3[/cyan] [magenta]{vx3[i]}[/magenta]");
                                if (Math.Abs(vx2[i]) > 0)
                                    alpha = Math.Min(-aob(vx1[i], vx2[i]), alpha);
                            }
                            else if (Math.Abs(vx2[i]) > 0 && (inner = 1.0 - 4 * aob(vx3[i], vx2[i]) * aob(vx1[i], vx2[i])) >= 0)
                            {
                                double r1, r2;
                                inner = (inner > 0 ? Math.Sqrt(inner) * Math.Abs(vx2[i]) : 0);
                                r1 = aob((-vx2[i] - inner), vx3[i]) / 2.0; r2 = aob((-vx2[i] + inner), vx3[i]) / 2.0;
                                if (vx3[i] > 0)
                                {
                                    if (r1 < 0) r1 = r2;
                                    alpha = Math.Min(Math.Max(r1, 0), alpha);
                                }
                                else
                                {
                                    if (r2 < 0) r2 = r1;
                                    alpha = Math.Min(Math.Max(r2, 0), alpha);
                                }
                            }
                            else if (Math.Abs(vx2[i]) <= 0 && (inner = -4 * vx3[i] * vx1[i]) >= 0)
                            {
                                double r1, r2;
                                inner = (inner > 0 ? Math.Sqrt(inner) : 0);
                                r1 = aob((-vx2[i] - inner), vx3[i]) / 2.0; r2 = aob((-vx2[i] + inner), vx3[i]) / 2.0;
                                if (vx3[i] > 0)
                                {
                                    if (r1 < 0) r1 = r2;
                                    alpha = Math.Min(Math.Max(r1, 0), alpha);
                                }
                                else
                                {
                                    if (r2 < 0) r2 = r1;
                                    alpha = Math.Min(Math.Max(r2, 0), alpha);
                                }
                            }
                            else if ((inner = Math.Abs(1.0 - 4.0 * aob(vz3[i], vz2[i]) * aob(vz1[i], vz2[i]))) < BlasLike.lm_eps)
                            {
                                var r2 = aob(-vx2[i], vx3[i]) / 2.0;
                                alpha = Math.Min(Math.Max(r2, 0), alpha);
                            }
                            else
                                ColourConsole.WriteEmbeddedColourLine($"[red]ROUNDING cone x{i}[/red][magenta] step {alpha}[/magenta] [green]{vx1[i]}[/green] [yellow]{vx2[i]}[/yellow] [cyan]{vx3[i]}[/cyan]");
                        }
                        if (alpha == 0) break;
                    }
                }
                if (homogenous)
                {
                    if (dtau < 0) alpha = Math.Min(alpha, -aob(tau, dtau));
                    if (dkappa < 0) alpha = Math.Min(alpha, -aob(kappa, dkappa));
                }
                double gamma1 = 1 - gamma, test1, test2 = 1, beta = 1e-7;
                bool bad = true;
                var rhs = beta * (1.0 - alpha * gamma1) * mu;
                double ratio, ratiolimL = BlasLike.lm_eps16384, ratiolimU = 1.0 / ratiolimL;

                for (int i = 0, cstart = 0; i < cone.Length; cstart += cone[i], i++)
                {
                    var n = cone[i];
                    test1 = 0; bad = true;
                    while (bad == true && alpha > BlasLike.lm_eps)
                    {
                        rhs = beta * (1.0 - alpha * gamma1) * mu;
                        if (homogenous && kappa > BlasLike.lm_eps && dkappa != 0) test2 = (tau + alpha * dtau) * (kappa + alpha * dkappa);
                        if (typecone[i] == (int)conetype.SOCP)
                        {
                            var roundlim = BlasLike.lm_eps;
                            if (n > 1)
                            {
                                ratio = ((vz1[i] + alpha * (vz2[i] + alpha * vz3[i])) / (vx1[i] + alpha * (vx2[i] + alpha * vx3[i])));
                                if ((test1 = (vx1[i] + alpha * (vx2[i] + alpha * vx3[i])) * (vz1[i] + alpha * (vz2[i] + alpha * vz3[i]))) > square(roundlim)
                                && (ratio > square(ratiolimL) && ratio < square(ratiolimU)))
                                {
                                    if (test1 > square(rhs) && test2 > rhs) { bad = false; }
                                }
                            }
                            else
                            {
                                ratio = ((z[cstart] + alpha * dz[cstart]) / (x[cstart] + alpha * dx[cstart]));
                                if ((test1 = (x[cstart] + alpha * dx[cstart]) * (z[cstart] + alpha * dz[cstart])) > roundlim
                                && (ratio > ratiolimL && ratio < ratiolimU))
                                {
                                    if (test1 > rhs && test2 > rhs) { bad = false; }
                                }
                            }
                        }
                        else if (typecone[i] == (int)conetype.QP)
                        {
                            for (var ii = cstart; ii < n + cstart; ++ii)
                            {
                                if ((test1 = (x[ii] + alpha * dx[ii]) * (z[ii] + alpha * dz[ii])) > 0)
                                {
                                    if (test1 > rhs && test2 > rhs) { bad = false; continue; }
                                }
                                bad = true;
                                break;
                            }
                        }
                        if (!bad) break;
                        alpha *= lowest1;
                    }
                }
                ddx = ddz = dd = alpha;
            }

            else if (optMode == "QP")
            {
                double lowest = 1e-3, lowest1 = 1 - lowest;
                ddx = 1.0;
                ddz = 1.0;
                dd = 1.0;
                for (int i = 0; i < n; ++i)
                {
                    if (sign == null)
                    {
                        if (dx[i] < 0) ddx = Math.Min(ddx, -aob(x[i], dx[i]));
                        if (dz[i] < 0) ddz = Math.Min(ddz, -aob(z[i], dz[i]));
                    }
                    else
                    {
                        if (dx[i] * sign[i] < 0)
                            ddx = Math.Min(ddx, -aob(x[i], dx[i]));
                        if (dz[i] * sign[i] < 0)
                            ddz = Math.Min(ddz, -aob(z[i], dz[i]));
                    }
                }
                if (homogenous)
                {
                    if (dtau < 0) dd = Math.Min(dd, -aob(tau, dtau));
                    if (dkappa < 0) dd = Math.Min(dd, -aob(kappa, dkappa));
                }
                ddx *= lowest1;
                ddz *= lowest1;
                dd *= lowest1;
            }
        }
        void CreateNormalMatrix()
        {
            var clock1 = clocker();
            if (optMode == "QP")
            {
                BlasLike.dzerovec(M.Length, M);
                /*    for (var k = 0; k < Math.Min(20, x.Length); ++k)
                    {
                        Console.Write(x[k]);
                        Console.Write(",");
                    }
                    Console.Write("\n");
                    for (var k = 0; k < Math.Min(20, z.Length); ++k)
                    {
                        Console.Write(z[k]);
                        Console.Write(",");
                    }
                    Console.Write("\n");*/
                if (usrH)
                {
                    var lhs = w1;//Highjack w1 to save reallocation
                    HCOPY = (double[])H.Clone();
                    for (int i = 0, ij = 0; i < nh; ++i, ij += i) HCOPY[i + ij] += aob(z[i], x[i]);
                    //HHCopy = (double[])HCOPY.Clone();
                    Factorise.Factor(uplo, nh, HCOPY, horder);
                    for (int con = 0, ij = 0; con < m; ++con, ij += con)
                    {
                        if (A != null)
                            BlasLike.dcopy(n, A, m, lhs, 1, con);
                        else
                        {
                            if (con < basem)
                            {
                                BlasLike.dcopy(basen, baseA, basem, lhs, 1, con);
                            }
                            else if (con < basem + bases)
                            {
                                BlasLike.dzerovec(basen, lhs);
                                var qq = slacklargeConstraintToStock[con - basem];
                                lhs[qq] = 1;
                            }
                            else if (con < basem + bases + slackmboth)
                            {
                                var conn = slackToConstraintBOTH[con - basem - bases];
                                BlasLike.dcopy(basen, baseA, basem, lhs, 1, conn);
                            }
                        }
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, lhs, nh);
                        //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, lhs);
                        if (basen != nh)
                        {
                            for (int i = nh; i < basen; ++i)
                            {
                                if (lhs[i] != 0)
                                    lhs[i] /= aob(z[i], x[i]);
                            }
                        }
                        if (A != null)
                        {
                            for (var ii = 0; ii < con + 1; ++ii)
                            {
                                M[ij + ii] = BlasLike.ddot(n, A, m, lhs, 1, ii);
                            }
                        }
                        else
                        {
                            for (var ii = 0; ii <= con; ++ii)
                            {
                                if (ii < basem)
                                {
                                    M[ij + ii] = BlasLike.ddot(basen, baseA, basem, lhs, 1, ii);
                                    if (ii == con && slackmboth > 0)
                                    {
                                        var iii = slackToConstraintBOTH_inverse[ii];
                                        if (iii != -1)
                                        {
                                            var qq = basen + bases + slackmboth + iii;
                                            M[ij + ii] += aob(x[qq], z[qq]);
                                        }
                                    }
                                    if (ii == con && slackmL > 0)
                                    {
                                        var iii = slackToConstraintL_inverse[ii];
                                        if (iii != -1)
                                        {
                                            var qq = basen + bases + slackmboth * 2 + iii;
                                            M[ij + ii] += aob(x[qq], z[qq]);
                                        }
                                    }
                                    if (ii == con && slackmU > 0)
                                    {
                                        var iii = slackToConstraintU_inverse[ii];
                                        if (iii != -1)
                                        {
                                            var qq = basen + bases + slackmboth * 2 + slackmL + iii;
                                            M[ij + ii] += aob(x[qq], z[qq]);
                                        }
                                    }
                                }
                                else if (ii < basem + bases)
                                {
                                    var qq = slacklargeConstraintToStock[ii - basem];
                                    M[ij + ii] = lhs[qq];
                                    if (ii == con)
                                    {
                                        qq = basen + ii - basem;
                                        M[ij + ii] += aob(x[qq], z[qq]);
                                    }
                                }
                                else if (ii < basem + bases + slackmboth)
                                {
                                    var iii = slackToConstraintBOTH[ii - basem - bases];
                                    M[ij + ii] = BlasLike.ddot(basen, baseA, basem, lhs, 1, iii);
                                    if (ii == con)
                                    {
                                        var qq = (ii - basem - bases) + basen + bases;
                                        M[ij + ii] += aob(x[qq], z[qq]);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (A != null)
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
                    else
                    {
                        for (int i = 0, ij = 0; i < m; ++i, ij += i)
                        {
                            if (i < basem)
                            {
                                for (var k = 0; k < basen; ++k)
                                {
                                    if (baseA[k * basem + i] != 0.0)
                                    {
                                        var xoz = aob(x[k], z[k]);
                                        xoz *= baseA[k * basem + i];
                                        BlasLike.daxpyvec(i + 1, xoz, baseA, M, k * basem, ij);
                                    }
                                }
                                if (slackmboth > 0)
                                {
                                    var ii = slackToConstraintBOTH_inverse[i];
                                    if (ii != -1)
                                    {
                                        var kk = ii + basen + bases + slackmboth;
                                        var xozz = -aob(x[kk], z[kk]);
                                        M[ij + i] += -xozz;
                                    }
                                }
                                if (slackmL > 0)
                                {
                                    var ii = slackToConstraintL_inverse[i];
                                    if (ii != -1)
                                    {
                                        var kk = ii + basen + bases + slackmboth * 2;
                                        var xozz = -aob(x[kk], z[kk]);
                                        M[ij + i] += -xozz;
                                    }
                                }
                                if (slackmU > 0)
                                {
                                    var ii = slackToConstraintU_inverse[i];
                                    if (ii != -1)
                                    {
                                        var kk = ii + basen + bases + slackmboth * 2 + slackmL;
                                        var xozz = aob(x[kk], z[kk]);
                                        M[ij + i] += xozz;
                                    }
                                }
                            }
                            else if (i < basem + bases)
                            {
                                var k = slacklargeConstraintToStock[i - basem];
                                var xoz = aob(x[k], z[k]);
                                for (var j = 0; j < basem; ++j)
                                {
                                    M[ij + j] += baseA[k * basem + j] * xoz;
                                }
                                M[ij + i] += xoz;
                                k = i - basem + basen;
                                xoz = aob(x[k], z[k]);
                                M[ij + i] += xoz;
                            }
                            else if (i < basem + bases + slackmboth)
                            {
                                var ii = slackToConstraintBOTH[i - basem - bases];
                                for (var k = 0; k < basen; ++k)
                                {
                                    if (baseA[k * basem + ii] != 0.0)
                                    {
                                        var xoz = aob(x[k], z[k]);
                                        var cons = 0;
                                        xoz *= baseA[k * basem + ii];
                                        BlasLike.daxpyvec(basem, xoz, baseA, M, k * basem, ij);
                                        if (bases > 0 && (cons = slacklargeConstraintToStock_inverse[k]) != -1)
                                            M[ij + cons + basem] += xoz;
                                        for (var j = 0; j <= i - basem - bases; j++)
                                        {
                                            var jj = slackToConstraintBOTH[j];
                                            M[ij + basem + bases + j] += xoz * baseA[k * basem + jj];
                                        }
                                    }
                                }
                                var kk = i - basem + basen;
                                M[ij + i] += aob(x[kk], z[kk]);
                            }
                        }
                    }
                }
                if (specialDebug)
                {
                    ColourConsole.WriteInfo($"usrH {usrH} NORMAL MATRIX");
                    for (var k = 0; k < M.Length; ++k)
                    {
                        Console.Write($"{M[k]:F10}");
                        Console.Write(",");
                        if (k % 20 == 19)
                            Console.Write("\n");
                    }
                    Console.Write("\n");
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
                                AmultSparse(lhs, res, cstart * m, cstart, cstart);
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
                            BlasLike.dcopy(n, A, m, lhs, 1, i + cstart * m, n);
                            if (n == 1 && lhs[n] == 0) continue;
                            if (n == 2 && lhs[n] == 0 && lhs[n + 1] == 0) continue;
                            W2m1trans(n, lhs, W, lhs, n, cstart);
                            thetaScale(n, lhs, THETA[icone], true, true);
                            for (var k = 0; k < n; ++k)
                            {
                                if (lhs[k] != 0.0)
                                {
                                    BlasLike.daxpyvec(i + 1, lhs[k], A, M, (k + cstart) * m, ij);
                                }
                            }
                        }
                    }
                }
            }
            ColourConsole.WriteInfo($"Time in CreateNormalMatrix {clocker() - clock1} m secs");
        }
        void adaptedResiduals(double[] rp, double[] rd, double[] rm, double delta, double[] rpnew, double[] rdnew, double[] rmnew)
        {
            //Try to use approximate residuals when step length is small
            var rpr = lInfinity(rp);
            var rdr = lInfinity(conetop(rd));
            var rmr = lInfinity(rm);
            if (rpr == 0) rpr = 1;
            if (rdr == 0) rdr = 1;
            if (rmr == 0) rmr = 1;
            BlasLike.dscalvec(rp.Length, -delta, rpnew);
            BlasLike.dscalvec(rd.Length, -delta, rdnew);
            BlasLike.dscalvec(rm.Length, -delta, rmnew);
            BlasLike.daddvec(rp.Length, rp, rpnew, rpnew);
            BlasLike.daddvec(rd.Length, rd, rdnew, rdnew);
            BlasLike.daddvec(rm.Length, rm, rmnew, rmnew);
        }
        void SolvePrimaryDual(double gamma = 0.0, bool corrector = false)
        {//BlasLike.dcopyvec(xbar.Length,zbar,xbar);
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
                    //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, w1);
                    for (int i = nh; i < n; ++i)
                        if (w1[i] != 0)
                            w1[i] /= aob(z[i], x[i]);
                }
                else for (int i = 0; i < n; ++i) w1[i] *= aob(x[i], z[i]);
                if (specialDebug)
                {
                    ColourConsole.WriteInfo($"{usrH} w1");
                    for (var k = 0; k < n; ++k)
                    {
                        Console.Write($"{w1[k]:F10}");
                        Console.Write(",");
                        if (k % 20 == 19)
                            Console.Write("\n");
                    }
                    Console.Write("\n");
                    ColourConsole.WriteInfo($"{usrH} x");
                    for (var k = 0; k < n; ++k)
                    {
                        Console.Write($"{x[k]:F10}");
                        Console.Write(",");
                        if (k % 20 == 19)
                            Console.Write("\n");
                    }
                    Console.Write("\n");
                    ColourConsole.WriteInfo($"{usrH} z");
                    for (var k = 0; k < n; ++k)
                    {
                        Console.Write($"{z[k]:F10}");
                        Console.Write(",");
                        if (k % 20 == 19)
                            Console.Write("\n");
                    }
                    Console.Write("\n");
                    ColourConsole.WriteInfo($"{usrH} dy");
                    for (var k = 0; k < Math.Min(20, dy.Length); ++k)
                    {
                        Console.Write($"{dy[k]:F10}");
                        Console.Write(",");
                    }
                    Console.Write("\n");
                }
                AmultSparse(w1, dy);
                if (specialDebug)
                {
                    ColourConsole.WriteInfo($"{usrH} dy after Amult");
                    for (var k = 0; k < Math.Min(20, dy.Length); ++k)
                    {
                        Console.Write($"{dy[k]:F10}");
                        Console.Write(",");
                    }
                    Console.Write("\n");
                }
                BlasLike.daxpyvec(m, g1, rp, dy);
                if (specialDebug)
                {
                    ColourConsole.WriteInfo($"{usrH} dy After axpy");
                    for (var k = 0; k < Math.Min(20, dy.Length); ++k)
                    {
                        Console.Write($"{dy[k]:F10}");
                        Console.Write(",");
                    }
                    Console.Write("\n");
                }
                if (m == 1) dy[0] /= M[0];
                else Factorise.Solve(uplo, m, 1, M, order, dy, m);
                if (homogenous)
                {
                    BlasLike.dcopyvec(n, cmod, cx);
                    if (usrH)
                    {
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, cx, nh);
                        //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, cx);
                        for (int i = 0; i < (n - nh); ++i) cx[i + nh] /= aob(z[i + nh], x[i + nh]);
                    }
                    else for (int i = 0; i < n; ++i) cx[i] *= aob(x[i], z[i]);
                    AmultSparse(cx, db);
                    BlasLike.daddvec(m, db, b, db);
                    if (m == 1) db[0] /= M[0];
                    else Factorise.Solve(uplo, m, 1, M, order, db, m);
                    AmultSparseT(dy, dx);
                    if (usrH)
                    {
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh);
                        //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, dx);
                        for (int i = 0; i < (n - nh); ++i) dx[i + nh] /= aob(z[i + nh], x[i + nh]);
                    }
                    else for (int i = 0; i < n; ++i) dx[i] *= aob(x[i], z[i]);
                    BlasLike.dsubvec(n, dx, w1, dx);
                    AmultSparseT(db, dc);
                    if (usrH)
                    {
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, dc, nh);
                        //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, dc);
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
                    if (specialDebug)
                    {
                        ColourConsole.WriteInfo($"{usrH} dy after Solve");
                        for (var k = 0; k < Math.Min(20, dy.Length); ++k)
                        {
                            Console.Write($"{dy[k]:F10}");
                            Console.Write(",");
                        }
                        Console.Write("\n");
                    }
                    AmultSparseT(dy, dx);
                    if (usrH)
                    {
                        Factorise.Solve(uplo, nh, 1, HCOPY, horder, dx, nh);
                        //Factorise.SolveRefine(nh, HHCopy, HCOPY, horder, dx);
                        for (int i = 0; i < (n - nh); ++i) dx[i + nh] /= aob(z[i + nh], x[i + nh]);
                    }
                    else for (int i = 0; i < n; ++i) dx[i] *= aob(x[i], z[i]);
                    /* ColourConsole.WriteInfo($"{usrH} dx");
                     for (var k = 0; k < Math.Min(20, dx.Length); ++k)
                     {
                         Console.Write(dx[k]);
                         Console.Write(",");
                     }
                     Console.Write("\n");*/
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
                        /*   ColourConsole.WriteInfo($"{usrH} w1");
                           for (var k = 0; k < Math.Min(20, w1.Length); ++k)
                           {
                               Console.Write(w1[k]);
                               Console.Write(",");
                           }
                           Console.Write("\n");*/
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        W2m1trans(n, w1, W, w1, cstart, cstart, x.Length + cstart);
                        BlasLike.dcopyvec(n, w1, w1, x.Length + cstart, cstart);
                        thetaScale(n, w1, THETA[icone], true, true, cstart);
                    }
                }
                AmultSparse(w1, dy);
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
                    AmultSparse(cx, db);
                    BlasLike.daddvec(m, db, b, db);
                    if (m == 1) db[0] /= M[0];
                    else Factorise.Solve(uplo, m, 1, M, order, db, m);
                    AmultSparseT(dy, dx);
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
                    AmultSparseT(db, dc);
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
                    AmultSparseT(dy, dx);
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
                            /*   ColourConsole.WriteInfo($"{usrH} dx");
                               for (var k = 0; k < Math.Min(20, dx.Length); ++k)
                               {
                                   Console.Write(dx[k]);
                                   Console.Write(",");
                               }
                               Console.Write("\n");*/
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
        void PrimalResidual(double usetau = 0.0)
        {
            var tau = usetau == 0.0 ? this.tau : usetau;
            AmultSparse(x, Ax);
            BlasLike.dcopyvec(m, Ax, rp);
            BlasLike.dnegvec(m, rp);
            BlasLike.daxpyvec(m, homogenous ? tau : 1, b, rp);
        }
        void DualResudual(double usetau = 0.0)
        {
            var tau = usetau == 0.0 ? this.tau : usetau;
            AmultSparseT(y, Ay);
            BlasLike.dcopyvec(n, Ay, rd);
            BlasLike.dnegvec(n, rd);
            BlasLike.daxpyvec(n, (homogenous ? tau : 1), cmod, rd);
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
                            THETA[icone] = 1;
                            xbar[i] = x[i] * W[i];
                            zbar[i] = z[i] / W[i];
                            rmu[i] = mu - xbar[i] * zbar[i];
                        }
                    }
                    else if (typecone[icone] == (int)conetype.SOCP)
                    {
                        var xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                        var zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);

                        if (n > 1) THETA[icone] = Math.Sqrt(Math.Sqrt(zQz / xQx));
                        else THETA[icone] = Math.Sqrt(z[n - 1 + cstart] / x[n - 1 + cstart]);
                        if (false && (double.IsNaN(THETA[icone]) || THETA[icone] == 0 || xQx == 0))
                        {
                            if (double.IsNaN(zQz) || zQz <= 0)
                            {
                                double z1 = z[n - 1 + cstart] * z[n - 1 + cstart];
                                double inner = z1 - zQz;
                                double fac = zQz <= 0 ? 1 - BlasLike.lm_eps16 : Math.Sqrt((z1 - BlasLike.lm_eps) / inner);
                                BlasLike.dscalvec(n - 1, fac, x, cstart);
                                BlasLike.dscalvec(n - 1, fac, z, cstart);
                                xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                            }
                            if (double.IsNaN(xQx) || xQx <= 0)
                            {
                                double x1 = x[n - 1 + cstart] * x[n - 1 + cstart];
                                double inner = x1 - xQx;
                                double fac = xQx <= 0 ? 1 - BlasLike.lm_eps16 : Math.Sqrt((x1 - BlasLike.lm_eps) / inner);
                                BlasLike.dscalvec(n - 1, fac, x, cstart);
                                BlasLike.dscalvec(n - 1, fac, z, cstart);
                                xQx = x[n - 1 + cstart] * x[n - 1 + cstart] - BlasLike.ddotvec(n - 1, x, x, cstart, cstart);
                                zQz = z[n - 1 + cstart] * z[n - 1 + cstart] - BlasLike.ddotvec(n - 1, z, z, cstart, cstart);
                            }
                            THETA[icone] = Math.Sqrt(Math.Sqrt(zQz / xQx));
                        }



                        if (!double.IsNaN(THETA[icone]) && THETA[icone] != 0)
                        {
                            double zx = BlasLike.ddotvec(n, z, x, cstart, cstart);
                            double bot = Math.Sqrt((zx + Math.Sqrt(xQx * zQz)) * 2.0);
                            double z1 = THETA[icone] / bot;
                            double z2 = THETA[icone] * bot;
                            //  if (THETA[icone] == BlasLike.lm_eps) { z1 = 1.0; z2 = BlasLike.lm_eps; }
                            for (var i = cstart; i < n - 1 + cstart; ++i)
                            {
                                W[i] = -z1 * x[i] + z[i] / z2;
                            }
                            W[n - 1 + cstart] = z1 * x[n - 1 + cstart] + z[n - 1 + cstart] / z2;
                        }
                        else
                        {
                            ColourConsole.WriteEmbeddedColourLine($"[yellow]THETA for cone {icone} is [/yellow][red]{THETA[icone]}[/red] [magenta]({zQz} / {xQx})[/magenta]");
                            Debug.Assert(false);
                            //DO something to stop main loop

                        }
                        double wcheck = BlasLike.ddotvec(n, W, W, cstart, cstart);
                        if (double.IsNaN(wcheck))
                        {
                            ColourConsole.WriteInfo("BAD W"); Debug.Assert(false);
                        }
                        Wtrans(n, x, W, xbar, cstart, cstart, cstart); //xbar=thetaW.x
                        thetaScale(n, xbar, THETA[icone], false, false, cstart);
                        Wm1trans(n, z, W, zbar, cstart, cstart, cstart);//z=(thetaW)(thetaW)x
                        thetaScale(n, zbar, THETA[icone], true, false, cstart);//zbar=(Wtheta)m1.z=xbar
                        Tmulvec(n, xbar, cstart);//Tmulvec does nothing for SOCP, needed for SOCPR
                        Tmulvec(n, zbar, cstart);
                        //            BlasLike.dcopyvec(n, xbar, zbar, cstart, cstart);
                        applyX(n, xbar, zbar, rmu, cstart, cstart, cstart);
                        Tmulvec(n, rmu, cstart);
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
        void AmultSparseT(double[] y, double[] x, int astart = 0, int ystart = 0, int xstart = 0)
        {
            if (A != null)
                Factorise.dmxmulv(n, m, A, y, x, astart, ystart, xstart, true);
            else
            {
                for (var k = 0; k < slackmboth; ++k)
                {
                    var km = slackToConstraintBOTH[k];
                    y[ystart + km] += y[ystart + basem + bases + k];
                }
                Factorise.dmxmulv(basen, basem, baseA, y, x, astart, ystart, xstart, true);
                for (var k = 0; k < slackmboth; ++k)
                {
                    var km = slackToConstraintBOTH[k];
                    y[ystart + km] -= y[ystart + basem + bases + k];
                }
                for (var ik = 0; ik < bases; ++ik)
                {
                    var k = slacklargeConstraintToStock[ik];
                    x[xstart + basen + ik] = y[ystart + basem + ik];
                    x[xstart + k] += x[xstart + basen + ik];
                }
                for (var k = 0; k < slackmboth; ++k)
                {
                    var km = slackToConstraintBOTH[k];
                    x[xstart + basen + bases + k] = y[ystart + basem + bases + k];
                    x[xstart + basen + bases + slackmboth + k] = -y[ystart + km];
                }
                for (var k = 0; k < slackmL; ++k)
                {
                    var km = slackToConstraintL[k];
                    x[xstart + basen + bases + slackmboth * 2 + k] = -y[ystart + km];
                }
                for (var k = 0; k < slackmU; ++k)
                {
                    var km = slackToConstraintU[k];
                    x[xstart + basen + bases + slackmboth * 2 + slackmL + k] = y[ystart + km];
                }
            }
        }
        void AmultSparse(double[] x, double[] y, int astart = 0, int xstart = 0, int ystart = 0)
        {
            if (A != null)
                Factorise.dmxmulv(m, n, A, x, y, astart, xstart, ystart);
            else
            {
                Factorise.dmxmulv(basem, basen, baseA, x, y, astart, xstart, ystart);
                for (int km, k = 0; k < slackmboth; ++k)
                {
                    km = slackToConstraintBOTH[k];
                    y[k + ystart + basem + bases] = y[km + ystart] + x[xstart + k + basen + bases];
                    y[ystart + km] -= x[xstart + k + basen + bases + slackmboth];
                }
                for (int km, k = 0; k < slackmL; ++k)
                {
                    km = slackToConstraintL[k];
                    y[ystart + km] -= x[xstart + k + basen + bases + slackmboth * 2];
                }
                for (int km, k = 0; k < slackmU; ++k)
                {
                    km = slackToConstraintU[k];
                    y[ystart + km] += x[xstart + k + basen + bases + slackmboth * 2 + slackmL];
                }
                for (var ik = 0; ik < bases; ++ik)
                {
                    var k = slacklargeConstraintToStock[ik];
                    y[ystart + ik + basem] = x[xstart + k] + x[xstart + ik + basen];
                }
            }
        }
        void ConditionEstimate()
        {
            if (m > 1)
            {
                var diags = new double[m];
                var order = new int[m];
                for (int i = 0, ij = 0; i < m; ++i, ij += i + 1)
                {
                    //    Debug.Assert(i * (i + 3) / 2 == ij);
                    diags[i] = M[ij];
                }
                try
                {
                    Ordering.Order.getorder(m, diags, order, null, 0.0, 1);
                    int o1 = Math.Max(order[0], order[m - 1]), o2 = Math.Min(order[0], order[m - 1]);
                    double a1 = Math.Max(diags[o1], diags[o2]), a2 = Math.Min(diags[o2], diags[o1]), a12 = M[o1 * (o1 + 1) / 2 + o2];
                    condition = a1 / (a2 - a12 * a12 / a1);//cond is a quick estimate of condition number using only 2 pivots.
                    regularise = a1 * BlasLike.lm_eps;
                }
                catch
                {
                    ColourConsole.WriteInfo("######################################### Fix-up for getorder exception in condition estimate ###################################");
                    double a1 = 0, a2 = 0;
                    BlasLike.dxminmax(m, diags, 1, ref a1, ref a2);
                    condition = a1 / a2;
                    regularise = a1 * BlasLike.lm_eps;
                }
            }
            else
            {
                condition = square(M[0]);
                regularise = M[0] * BlasLike.lm_eps;
            }
        }
        double Complementarity(bool scale = false)
        {
            double tau = this.tau;
            double kappa = this.kappa;
            if (scale) { tau = 1; kappa = this.kappa / this.tau; }
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
        public void ConeReset(double correction = 0)
        {
            if (correction == 0) correction = BlasLike.lm_eps;
            for (int icc = 0, conestart = 0; icc < numberOfCones; ++icc)
            {
                var reduce = false;
                if (typecone[icc] == (int)InteriorPoint.conetype.SOCP && cone[icc] > 1)
                {
                    var din = BlasLike.ddotvec(cone[icc] - 1, x, x, conestart, conestart);
                    var dtop = square(x[conestart + cone[icc] - 1]);
                    var scale = dtop / din;

                    din = BlasLike.ddotvec(cone[icc] - 1, z, z, conestart, conestart);
                    dtop = square(z[conestart + cone[icc] - 1]);
                    scale = Math.Max(scale, dtop / din);

                    if (scale < (1 + correction))
                        reduce = true;

                    if (reduce)
                    {
                        var reduction = 0.95;//Mess with the cone drastically
                        ColourConsole.WriteEmbeddedColourLine($"\t\t\t[red]REDUCTION[/red]\t\t\t\t[green]{reduction}[/green]");
                        //    BlasLike.dscalvec(cone[icc] - 1, reduction, z, conestart);
                        //    BlasLike.dscalvec(cone[icc] - 1, reduction, x, conestart);
                        x[conestart + cone[icc] - 1] /= reduction;
                        z[conestart + cone[icc] - 1] /= reduction;
                    }
                }
                conestart += cone[icc];
            }
        }
        public bool coneQcheck(double step = 0, double[] dx = null, double[] dz = null)
        {
            double testxQx, testzQz;
            var back = false;
            var eps = BlasLike.lm_eps;
            int nn, cstart, ic;
            double limL = BlasLike.lm_eps, limU = 1.0 / limL, tt;
            limL = square(limL);
            limU = square(limU);
            for (ic = 0, cstart = 0; ic < numberOfCones; cstart += nn, ++ic)
            {
                nn = cone[ic];
                if (nn > 1)
                {
                    testxQx = x[nn - 1 + cstart] * x[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, x, x, cstart, cstart);
                    testzQz = z[nn - 1 + cstart] * z[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, z, z, cstart, cstart);
                    if (step > 0)
                    {
                        var dxQx = 2.0 * (x[nn - 1 + cstart] * dx[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, dx, x, cstart, cstart));
                        var dzQz = 2.0 * (z[nn - 1 + cstart] * dz[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, dz, z, cstart, cstart));
                        var dxQdx = (dx[nn - 1 + cstart] * dx[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, dx, dx, cstart, cstart));
                        var dzQdz = (dz[nn - 1 + cstart] * dz[nn - 1 + cstart] - BlasLike.ddotvec(nn - 1, dz, dz, cstart, cstart));
                        testxQx += step * (dxQx + step * dxQdx);
                        testzQz += step * (dzQz + step * dzQdz);
                    }
                    tt = testzQz / testxQx;
                    if (tt < limL || tt > limU || testxQx <= eps || testzQz <= eps)
                    {
                        back = true;
                    }
                    if (back)
                    {
                        ColourConsole.WriteEmbeddedColourLine($"[cyan]xQx[/cyan][magenta]{ic} [/magenta] [green]{testxQx}[/green]");
                        ColourConsole.WriteEmbeddedColourLine($"[cyan]zQz[/cyan][magenta]{ic} [/magenta] [green]{testzQz}[/green]");
                    }
                }
            }
            return back;
        }
        public int Opt(string mode = "QP", int[] cone = null, int[] typecone = null, bool homogenous = true, double[] L = null, int[] sign = null)
        {
            var smallestStep = 1e-3;
            var rd1ONE = false;
            var opt = this;
            this.sign = sign;
            double[] bl = null;
            double[] QL = null;
            double zL = 0;
            ///<summary>stepReduce is the factor by which the step length to the boundary is reduced</summary>
            var stepReduce = 1 - 5e-2;
            opt.optMode = mode;
            if (mode == "SOCP")
            {
                opt.conv = (Math.Floor(1e-8 / BlasLike.lm_eps)) * BlasLike.lm_eps;
                opt.compConv = (Math.Floor(1e-8 / BlasLike.lm_eps)) * BlasLike.lm_eps;
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
            opt.usrH = (h == null && nh > 0 && (opt.H != null && BlasLike.dsumvec(opt.H.Length, opt.H) != 0.0)) || opt.h != null;
            if (mode == "QP")
            {
                if (slacklargeConstraintToStock != null)
                {
                    slacklarge = slacklargeConstraintToStock.Length;
                    Debug.Assert(slacklarge == bases);
                    slacklargeConstraintToStock_inverse = new int[basen];
                    for (var ib = 0; ib < basen; ++ib) slacklargeConstraintToStock_inverse[ib] = -1;
                    for (var ib = 0; ib < slacklarge; ++ib)
                    {
                        slacklargeConstraintToStock_inverse[slacklargeConstraintToStock[ib]] = ib;
                    }
                }
                if (slackToConstraintL != null)
                {
                    slackmL = slackToConstraintL.Length;
                    slackToConstraintL_inverse = new int[basem];
                    for (var ib = 0; ib < basem; ++ib) slackToConstraintL_inverse[ib] = -1;
                    for (var ib = 0; ib < slackToConstraintL.Length; ++ib)
                    {
                        slackToConstraintL_inverse[slackToConstraintL[ib]] = ib;
                    }
                }
                if (slackToConstraintU != null)
                {
                    slackmU = slackToConstraintU.Length;
                    slackToConstraintU_inverse = new int[basem];
                    for (var ib = 0; ib < basem; ++ib) slackToConstraintU_inverse[ib] = -1;
                    for (var ib = 0; ib < slackToConstraintU.Length; ++ib)
                    {
                        slackToConstraintU_inverse[slackToConstraintU[ib]] = ib;
                    }
                }
                if (slackToConstraintBOTH != null)
                {
                    slackmboth = slackToConstraintBOTH.Length;
                    slackToConstraintBOTH_inverse = new int[basem];
                    for (var ib = 0; ib < basem; ++ib) slackToConstraintBOTH_inverse[ib] = -1;
                    for (var ib = 0; ib < slackToConstraintBOTH.Length; ++ib)
                    {
                        slackToConstraintBOTH_inverse[slackToConstraintBOTH[ib]] = ib;
                    }
                }
                if (h == null) h = qphess1;
                if (sign == null)
                {
                    BlasLike.dsetvec(n, 1.0, opt.x);
                    BlasLike.dsetvec(n, 1.0, opt.z);
                }
                else
                {
                    for (var ii = 0; ii < n; ++ii)
                    {
                        opt.z[ii] = opt.x[ii] = (double)sign[ii];
                    }
                }
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
            innerIteration = 0;
            var ir = 1;
            var extra = new double[nh];
            if (opt.optMode == "QP")
            {
                if (L != null)
                {
                    bl = new double[m];
                    opt.AmultSparse(L, bl);
                    if (opt.H != null)
                    {
                        QL = new double[n];
                        h(n, opt.H, L, QL);
                        BlasLike.daddvec(n, opt.c, QL, opt.c);
                    }
                    BlasLike.dsubvec(m, opt.b, bl, opt.b);
                }
                if (opt.usrH)
                {
                    if (opt.homogenous)
                    {
                        BlasLike.dcopyvec(opt.c.Length, opt.x, opt.cmod);
                        BlasLike.dscalvec(opt.c.Length, 1.0 / opt.tau, opt.cmod);
                        h(nh, opt.H, opt.cmod, extra);
                    }
                    else h(nh, opt.H, opt.x, extra);
                    BlasLike.dzerovec(n - nh, opt.cmod, nh);
                    BlasLike.daddvec(nh, opt.c, extra, opt.cmod);
                    BlasLike.dcopyvec(n - nh, opt.c, opt.cmod, nh, nh);
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
                            BlasLike.dcopyvec(nn - nh, opt.c, opt.cmod, istart + nh, istart + nh);
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
            var rd0 = lInfinity(conetop(rd));
            var gap0 = denomTest(opt.Gap());
            rp0 = denomTest(rp0);
            rd0 = denomTest(rd0);
            var rp1 = BlasLike.lm_max;
            var rd1 = BlasLike.lm_max;
            var rp1prev = BlasLike.lm_max;
            var rd1prev = BlasLike.lm_max;
            var r_increase = 0;
            var gap1 = gap0;
            var comp0 = opt.Complementarity();
            opt.keep.set(opt.x, opt.y, opt.z, opt.tau, opt.kappa, rp0, rd0, comp0, Math.Max(Math.Max(rp0, rd0), comp0));
            var compnow = comp0;
            var comp1 = comp0;
            alpha1 = 1.0;
            alpha2 = 1.0;
            var gamma = 0.0;
            var diff = new double[n];
            var gap = 1.0;
            var iup = 0;
            while (true)
            {
                rd1ONE = (rd1 >= 1 || rp1 >= 1 || Double.IsNaN(rd1) || Double.IsNaN(rp1));
                rp1prev = rp1;
                rd1prev = rd1;
                rp1 = lInfinity(opt.rp) / denomTest(rp0);
                rd1 = lInfinity(conetop(rd)) / denomTest(rd0);
                if (rd1 > rd1prev || rp1 > rp1prev) r_increase++;
                else r_increase = 0;
                //      if (rd1ONE && (rd1 >= 1 || rp1 >= 1 || Double.IsNaN(rd1) || Double.IsNaN(rp1)))
                //          break;
                ColourConsole.WriteEmbeddedColourLine($"[darkgreen]{innerIteration,4}[/darkgreen] [magenta]rp1 {rp1:E10}[/magenta]\t[cyan]rd1 {rd1:E10}[/cyan]");
                gap = opt.Gap();
                gap1 = gap / denomTest(gap0);
                comp1 = opt.Complementarity();
                opt.keep.update(opt.x, opt.y, opt.z, opt.tau, opt.kappa, rp1, rd1, comp1);
                if (rp1 <= opt.conv * tau && rd1 <= opt.conv * tau && comp1 <= opt.compConv * tau)
                { copyBest = false; break; }
                if (condition > BlasLike.lm_reps)
                {
                    double[] pass = { rp1, rd1, comp1 };
                    if (lInfinity(pass) > opt.keep.norm)
                    {
                        iup++;
                    }
                    //                if (iup > 100)
                    //                    break;
                }
                //   if (comp1 < opt.compConv && opt.tau < 1e-5 * opt.kappa) break;
                if (ir > opt.maxouter)
                    break;
                if (false && (rp1 > rp1prev || rd1 > rd1prev) && innerIteration > 1 && r_increase > 10)
                {
                    ir++; innerIteration = 0; r_increase = 0;
                    /*       BlasLike.dcopyvec(opt.x.Length, opt.keep.x, opt.x);
                           BlasLike.dcopyvec(opt.y.Length, opt.keep.y, opt.y);
                           BlasLike.dcopyvec(opt.z.Length, opt.keep.z, opt.z);
                           opt.tau = opt.keep.tau;
                           opt.kappa = opt.keep.kappa;*/

                    //   if (alpha1 > alpha2) opt.update(dxold, dyold, dzold, dtauold, dkappaold, -alpha1, 1);
                    //   else opt.update(opt.dx, opt.dy, opt.dz, opt.dtau, opt.dkappa, -alpha2, 1);

                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    opt.kappa /= opt.tau * 0.5;
                    opt.tau = 1;
                    opt.Mu();
                    mu0 = opt.mu;
                    opt.PrimalResidual();
                    opt.DualResudual();
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    rp0 = lInfinity(opt.rp);
                    rd0 = lInfinity(conetop(rd));
                    gap0 = denomTest(opt.Gap());
                }

                else if (innerIteration > opt.maxinner)
                {
                    ir++; innerIteration = 0;
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    opt.kappa /= opt.tau;
                    opt.tau = 1;
                    opt.Mu();
                    mu0 = opt.mu;
                    opt.PrimalResidual();
                    opt.DualResudual();
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    rp0 = lInfinity(opt.rp);
                    rd0 = lInfinity(conetop(rd));
                    gap0 = denomTest(opt.Gap());
                }
                if (innerIteration > 2 && opt.optMode == "SOCP")
                {
                    BlasLike.dsubvec(n, opt.xbar, opt.zbar, diff);
                    double test = Math.Sqrt(BlasLike.ddotvec(n, diff, diff));
                    if (test > BlasLike.lm_eps)
                    {
                        ColourConsole.WriteInfo($"xbar test = {test}");
                        ColourConsole.WriteInfo($"rp1 = {rp1}");
                        ColourConsole.WriteInfo($"rd1 = {rd1}");
                        ColourConsole.WriteInfo($"comp1 = {comp1}");
                        //                 BlasLike.dcopyvec(xbar.Length, xbar, zbar);
                        //       break;
                    }
                }
                opt.SolvePrimaryDual();
                if (opt.badindex != 0 && m > 1) ColourConsole.WriteInfo($"Normal matrix is unstable: badindex={opt.badindex}");
                BlasLike.dcopyvec(n, opt.dx, dxold);
                BlasLike.dcopyvec(n, opt.dz, dzold);
                BlasLike.dcopyvec(m, opt.dy, dyold);
                dtauold = opt.dtau;
                dkappaold = opt.dkappa;
                opt.MaximumStep(0);
                alpha1 = stepReduce * opt.Lowest();
                BlasLike.dsccopyvec(opt.n, 1.0, opt.dx, opt.dx0);//was alpha1
                BlasLike.dsccopyvec(opt.n, 1.0, opt.dz, opt.dz0);//was alpha1
                opt.dtau0 = opt.dtau;//was *alpha1
                opt.dkappa0 = opt.dkappa;//was *alpha1
                gamma = opt.gfunc(alpha1);
                opt.SolvePrimaryDual(gamma, true);
                opt.MaximumStep(gamma);
                alpha2 = stepReduce * opt.Lowest();
                if (false && alpha2 < alphamin2 && alpha1 < alphamin2)
                {
                    double gammafirst = alpha2 > alpha1 ? gamma : 0;
                    bool clarify = alpha2 > alpha1;
                    if (!clarify) alphamin2 = Math.Max(alphamin2 * 0.1, BlasLike.lm_rooteps);
                    double alphak = Math.Max(alpha1, alpha2), alphanew, gammaold, alphaold = alphak, delt = 1e-12;
                    double gammabest = gammafirst;
                    if (clarify)
                    {
                        BlasLike.dcopyvec(n, opt.dx, dxold);
                        BlasLike.dcopyvec(n, opt.dz, dzold);
                        BlasLike.dcopyvec(m, opt.dy, dyold);
                        dtauold = opt.dtau;
                        dkappaold = opt.dkappa;
                        alpha1 = alpha2;
                    }
                    gamma = 0;
                    opt.SolvePrimaryDual(gamma, true);
                    opt.MaximumStep(gamma);
                    alphanew = stepReduce * opt.Lowest();
                    if (alphanew > alphak && alphanew >= alphaold)
                    {
                        gammabest = gamma;
                    }

                    gammaold = gamma;
                    alphaold = alphanew;
                    gamma = 0.999999999;
                    opt.SolvePrimaryDual(gamma, true);
                    opt.MaximumStep(gamma);
                    alphanew = stepReduce * opt.Lowest();
                    if (alphanew > alphak && alphanew >= alphaold)
                    {
                        gammabest = gamma;
                    }
                    delt = 1e-1;
                    while (gamma > delt && alphanew < 0.9)
                    {
                        alphaold = Math.Max(alphanew, alphaold);
                        gamma -= delt;
                        opt.SolvePrimaryDual(gamma, true);
                        opt.MaximumStep(gamma);
                        alphanew = stepReduce * opt.Lowest();
                        if (alphanew > alphak && alphanew >= alphaold)
                        {
                            gammabest = gamma;
                        }
                        if (alphanew >= alphaold)
                        {
                            gammaold = gamma;
                            if (alphanew == alphaold) break;
                        }
                    }
                    if (gamma != gammaold)
                    {
                        gamma = gammaold;
                        opt.SolvePrimaryDual(gamma, true);
                        opt.MaximumStep(gamma);
                        alphanew = stepReduce * opt.Lowest();
                        if (alphanew > alphak && alphanew >= alphaold)
                        {
                            gammabest = gamma;
                        }
                    }
                    delt = gamma * 1e-1;
                    while (gamma > delt && alphanew < 0.9)
                    {
                        alphaold = Math.Max(alphanew, alphaold);
                        gamma -= delt;
                        opt.SolvePrimaryDual(gamma, true);
                        opt.MaximumStep(gamma);
                        alphanew = stepReduce * opt.Lowest();
                        if (alphanew > alphak && alphanew >= alphaold)
                        {
                            gammabest = gamma;
                        }
                        if (alphanew >= alphaold)
                        {
                            gammaold = gamma;
                        }
                    }
                    if (gamma != gammaold)
                    {
                        gamma = gammaold;
                        opt.SolvePrimaryDual(gamma, true);
                        opt.MaximumStep(gamma);
                        alphanew = stepReduce * opt.Lowest();
                        if (alphanew > alphak && alphanew >= alphaold)
                        {
                            gammabest = gamma;
                        }
                    }
                    if ((alphanew = Math.Max(alphaold, alphanew)) >= alphak)
                    {
                        if (gamma != gammabest)
                        {
                            opt.SolvePrimaryDual(gammabest, true);
                            opt.MaximumStep(gammabest);
                            alphanew = stepReduce * opt.Lowest();
                        }
                        alpha2 = alphanew;
                        if (alphanew > alphak) ColourConsole.WriteEmbeddedColourLine($"[green]Found better step size[/green] [cyan]{alphanew}[/cyan] [red]{alphak} ({gamma}[/red],[magenta]{gammafirst}[/magenta][red])[/red]");
                        corrector_gamma = gammabest;
                    }

                    if (clarify && alphanew < alphak)
                    {
                        ColourConsole.WriteEmbeddedColourLine($"[red]CLARIFY[/red] [yellow]alpha2 not improved {alpha2} {alphanew}[/yellow]");
                    }

                }
                if (alpha1 >= alpha2) coneQcheck(alpha1, dxold, dzold);
                else coneQcheck(alpha2, dx, dz);  //Only need this for testing
                if (alpha1 > alpha2)
                {
                    opt.update(dxold, dyold, dzold, dtauold, dkappaold, alpha1);

                    while (coneQcheck() && alpha1 > alpha2)
                    {
                        opt.update(dxold, dyold, dzold, dtauold, dkappaold, -alpha1);
                        alpha1 *= 0.5;
                        opt.update(dxold, dyold, dzold, dtauold, dkappaold, alpha1);
                    }
                    if (alpha1 <= alpha2)
                        opt.update(dxold, dyold, dzold, dtauold, dkappaold, -alpha1);
                }
                if (alpha1 < alpha2)
                {
                    opt.update(opt.dx, opt.dy, opt.dz, opt.dtau, opt.dkappa, alpha2);

                    while (coneQcheck() && alpha2 > alpha1)
                    {
                        opt.update(opt.dx, opt.dy, opt.dz, opt.dtau, opt.dkappa, -alpha2);
                        alpha2 *= 0.5;
                        opt.update(opt.dx, opt.dy, opt.dz, opt.dtau, opt.dkappa, alpha2);
                    }
                }
                if (opt.usrH)
                {
                    if (opt.homogenous)
                    {
                        BlasLike.dcopyvec(opt.c.Length, opt.x, opt.cmod);
                        BlasLike.dscalvec(opt.c.Length, 1.0 / opt.tau, opt.cmod);
                        h(nh, opt.H, opt.cmod, extra);
                    }
                    else h(nh, opt.H, opt.x, extra);
                    BlasLike.dzerovec(opt.c.Length - nh, opt.cmod, nh);
                    BlasLike.daddvec(nh, opt.c, extra, opt.cmod);
                    BlasLike.dcopyvec(n - nh, opt.c, opt.cmod, nh, nh);
                }
                var t1 = 0.0;
                if ((homogenous && (t1 = Math.Max(alpha1, alpha2)) < opt.alphamin))
                {
                    ColourConsole.WriteEmbeddedColourLine($"[red]Small step[/red] [green]{t1}[/green] [magenta] Condition {condition}[/magenta]");
                    if (alpha1 < smallestStep && alpha2 < smallestStep)
                    {
                        ColourConsole.WriteEmbeddedColourLine($"\t\t\t[red]BREAK[/red] [cyan]due to zero step length[/cyan]");
                        break;
                    }
                    //       var scl = 1.0;
                    /*       if(Math.Abs(tau-1)>BlasLike.lm_rooteps){
                           BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                           BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                           BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);

                               opt.kappa /= opt.tau;
                               gap = opt.Primal() - opt.Dual();
                               opt.tau = scl;}
                           */
                    if ((condition <= BlasLike.lm_reps / 2))
                    {
                        ColourConsole.WriteEmbeddedColourLine($"\t\t\t[red]Modify cone[/red]");
                        opt.ConeReset(5e-1);
                    }
                    if (gap < 0)
                    {
                        var dgap = BlasLike.ddotvec(opt.y.Length, opt.y, opt.b);
                        BlasLike.dsetvec(opt.y.Length, 1.0, opt.y);
                        dgap -= BlasLike.ddotvec(opt.y.Length, opt.y, opt.b);
                        var step = -gap / dgap / opt.tau;
                        if (dgap != 0) BlasLike.dsetvec(opt.y.Length, (1.0 + BlasLike.lm_rooteps) * step, opt.y);
                    }
                }
                opt.Mu();
                mu0 = opt.mu;
                opt.PrimalResidual();
                opt.DualResudual();
                opt.MuResidual();
                opt.ConditionEstimate();
                if (opt.condition > BlasLike.lm_reps)
                {
                    var mult = rp1 / rd1;
                    if (mult < 1) mult = 1.0 / mult;
                    if (mult > 2)
                        for (int ii = 0, id = 0; ii < m; ++ii, id += ii)
                        {
                            opt.M[id + ii] *= 1e2;
                        }
                }
                //         gap = opt.Primal() - opt.Dual();
                if (homogenous && opt.tau < 1e-8 && opt.kappa < 1e-8)
                {
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    opt.kappa /= opt.tau * 0.5;
                    opt.tau = 1;
                    ir++; innerIteration = 0;
                    opt.Mu();
                    mu0 = opt.mu;
                    opt.PrimalResidual();
                    opt.DualResudual();
                    opt.MuResidual();
                    opt.ConditionEstimate();
                    rp0 = lInfinity(opt.rp);
                    rd0 = lInfinity(conetop(rd));
                    gap0 = denomTest(opt.Gap());
                }
                innerIteration++;
            }
            if (/*ir > 1 &&*/ opt.copyKept && copyBest)
            {
                BlasLike.dcopyvec(opt.x.Length, opt.keep.x, opt.x);
                BlasLike.dcopyvec(opt.y.Length, opt.keep.y, opt.y);
                BlasLike.dcopyvec(opt.z.Length, opt.keep.z, opt.z);
                opt.tau = opt.keep.tau;
                opt.kappa = opt.keep.kappa;
            }
            var infease = !(opt.homogenous && (opt.tau > 1e1 * opt.kappa));
            if (opt.homogenous)
            {
                ColourConsole.WriteEmbeddedColourLine($"[green]tau = {opt.tau}[/green] [cyan]kappa={opt.kappa}[/cyan]");
                if (opt.tau != 0 && opt.tau != 1)
                {
                    BlasLike.dscalvec(opt.x.Length, 1.0 / opt.tau, opt.x);
                    BlasLike.dscalvec(opt.z.Length, 1.0 / opt.tau, opt.z);
                    BlasLike.dscalvec(opt.y.Length, 1.0 / opt.tau, opt.y);
                }
                if (infease) ColourConsole.WriteError("INFEASIBLE");
            }
            if (optMode == "QP")
            {
                if (L != null)
                {
                    BlasLike.daddvec(n, opt.x, L, opt.x);
                    if (opt.H != null) BlasLike.dsubvec(n, opt.c, QL, opt.c);
                    BlasLike.daddvec(m, opt.b, bl, opt.b);
                    zL = BlasLike.ddotvec(n, L, opt.z);
                }
            }

            /*   if (innerIteration >= opt.maxinner || ir >= opt.maxouter) return -100;
               else*/
            if (opt.homogenous && infease) return 6;
            else
            {
                Mu();
                mu0 = opt.mu;
                PrimalResidual(1.0);
                DualResudual(1.0);
                comp1 = Complementarity(true);
                rp1 = lInfinity(rp);
                rd1 = lInfinity(conetop(rd));
                ColourConsole.WriteInfo($"{ir} outer iterations out of {opt.maxouter}");
                ColourConsole.WriteInfo($"{innerIteration} iterations out of {opt.maxinner}");
                ColourConsole.WriteInfo($"Relative Primal\tResidual\t\t {rp1}");
                ColourConsole.WriteInfo($"Relative Dual\tResidual\t\t {rd1}");

                ColourConsole.WriteInfo($"Relative Complementarity\t\t {comp1}");
                ColourConsole.WriteEmbeddedColourLine($"[green]Primal Utility:[/green]\t\t\t\t [cyan]{opt.Primal()}[/cyan]");
                ColourConsole.WriteEmbeddedColourLine($"[green]Dual Utility:[/green]\t\t\t\t [cyan]{opt.Dual() + zL}[/cyan]");
                ActiveSet.Optimise.printV("x", opt.x);
                ActiveSet.Optimise.printV("y", opt.y);
                ActiveSet.Optimise.printV("z", opt.z);
                ColourConsole.WriteInfo($"Complementarity:\t{BlasLike.ddotvec(opt.n, opt.x, opt.z) - zL}");
                ColourConsole.WriteInfo($"Gap:\t\t\t{opt.Primal() - opt.Dual() - zL}");
                ColourConsole.WriteInfo($"Job took {opt.clocker()} m secs");
                ColourConsole.WriteInfo($"Last conv {opt.conv}");
                return -iup;
            }
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
                    ColourConsole.WriteInfo($"outer is not positive in applyXm1 {x[ncm1 + xstart]} {outer}");
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
                    ColourConsole.WriteInfo($"fixed outer is now {outer}");
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
                    ColourConsole.WriteInfo($"bad cone in Wtrans {bot}");
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
                ColourConsole.WriteInfo($"bad cone in WtransR {bot}");
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

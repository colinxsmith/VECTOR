using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace ActiveSet
{
    public delegate void hessmull(int n, int nrowh, int ncolh, int j, double[] hess, double[] wrk, double[] hx,int hstart=0,int xstart=0);
    public class Optimise
    {
        static void shifter<S>(ref S a, ref S b, ref S c, ref S d)
        {
            a = b; b = c; c = d;
        }///<summary>Prototype for one dimensional function. Used as an argument for PathMin and Solve1D</summary>
         ///<param name="x">One dimensional variable</param>
        public delegate double OneD(double x, object info);
        ///<summary>Find x such that f(x)=0 such that x is between x=gammabot and x=gammatop and return the value of x which gives f(x)=0</summary>
        ///<param name="OneDimensionalFunction"> Function to be minimised</param>
        ///<param name="gammabot"> Lower value of optimisation domain</param>
        ///<param name="gammatop"> Upper value of optimisation domain</param>
        ///<param name="tol"> tolerance</param>
        public static double Solve1D(OneD OneDimensionalFunction, double gammabot = 0,
                                       double gammatop = 1.0, double tol = 0, object info = null)
        {
            if (tol == 0) { tol = BlasLike.lm_rooteps;  }
            if (gammatop == 1.0) { gammatop = 1-BlasLike.lm_eps;  }
            int iter, itmax = 200;
            short signk = 1;
            double c = 0, d = 0, e = 0, min1, min2, fc, p, q, r, s, tol1, xm;
            double gamma_opt, a = gammatop, b = gammabot;
            double fa = OneDimensionalFunction(a, info);
            if (Math.Abs(fa) < BlasLike.lm_rooteps)
                return a;
            double fb = OneDimensionalFunction(b, info);
            if (Math.Abs(fb) < BlasLike.lm_rooteps)
                return b;
            double scalelimit = Math.Max(((Math.Abs(fa) + Math.Abs(fb)) * .5), 1.0);
            if (tol == 0)
                tol = BlasLike.lm_eps;
            if (fa * fb > 0)
            {
                return 1e+10;
            }
            if (fa < fb)
            {
                double kk = fa;
                signk = -1;
                fa = fb; fb = kk;
            }
            fc = fb;
            gamma_opt = signk == 1 ? gammabot : gammatop;
            for (iter = 0; iter < itmax; iter++)
            {
                if (fb * fc > 0)
                {
                    c = a;
                    fc = fa;
                    e = d = b - a;
                }
                if (Math.Abs(fc) < Math.Abs(fb))
                {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                tol1 = 2 * BlasLike.lm_eps * Math.Abs(b) + 0.5 * tol;
                xm = 0.5 * (c - b);
                if (Math.Abs(xm) <= tol1 || Math.Abs(fb) < 1e-12)
                {
                    if (Math.Abs(fb) > 1e-02 * scalelimit)
                    {
                        ColourConsole.WriteEmbeddedColourLine($"[red]Problem with the 1d function shape (many valued or not continuous) error[/red] [yellow]g1={b} f({b})={fb} g2={c} f({c})={fc}[/yellow]");
                    }
                    //It used to be: return gamma_opt
			if (signk == 1)return b;
			else return (gammatop - b) + gammabot;
                }
                if (Math.Abs(e) >= tol1 && Math.Abs(fa) > Math.Abs(fb))
                {
                    s = fb / fa;
                    if (a == c)
                    {
                        p = 2 * xm * s;
                        q = 1 - s;
                    }
                    else
                    {
                        q = fa / fc;
                        r = fb / fc;
                        p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
                        q = (q - 1) * (r - 1) * (s - 1);
                    }
                    if (p > 0) q = -q;
                    p = Math.Abs(p);
                    min1 = 3 * xm * q - Math.Abs(tol1 * q);
                    min2 = Math.Abs(e * q);
                    if (2 * p < (min1 < min2 ? min1 : min2))
                    {
                        e = d;
                        d = p / q;
                    }
                    else
                    {
                        d = xm;
                        e = d;
                    }
                }
                else
                {
                    d = xm;
                    e = d;
                }
                a = b;
                fa = fb;
                b += (Math.Abs(d) > tol1) ? d : (xm > 0 ? Math.Abs(tol1) : -Math.Abs(tol1));
                if (signk == 1)
                {
                    fb = OneDimensionalFunction(b,info);
                    gamma_opt = b;
                }
                else
                {
                    double new_b = (gammatop - b) + gammabot;
                    fb = OneDimensionalFunction(new_b,info);
                    gamma_opt = new_b;
                }
            }
            return gamma_opt;
        }

        ///<summary>Minimise a function f(x) with 1 variable between x=gammabot and x=gammatop and return the value of x which gives this minimum</summary>
        ///<param name="OneDimensionalFunction"> Function to be minimised</param>
        ///<param name="gammabot"> Lower value of optimisation domain</param>
        ///<param name="gammatop"> Upper value of optimisation domain</param>
        ///<param name="tol"> tolerance</param>
        ///<param name="stopifpos"> stop if the function becomes positive if 1</param>
        public static double PathMin(OneD OneDimensionalFunction, double gammabot,
                                       double gammatop, double tol, int stopifpos, object info = null)
        {
            //	Based on routine in "Numerical Recipes in C" page 301.
            double ax = gammabot, bx = 0.5, cx = gammatop, ret_val;
            double a, b, d = 1e34, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, e = 0.0, basek;
            int error = 0;
            int count = 0, maxcount = 100;

            bx = (ax + cx) * .5;
            a = ((ax < cx) ? ax : cx);
            b = ((ax > cx) ? ax : cx);

            x = w = v = bx;

            fw = OneDimensionalFunction(x, info);
            basek = Math.Abs(fw); if (fw < 0) { fw = fv = fx = -1; } else { fw = fv = fx = 1; }
            if (basek == 0)
            {
                fw = fv = fx = basek;
                basek = 1;
            }
            ret_val = x;
            while (error == 0 && count < maxcount)
            {
                xm = 0.5 * (a + b);
                tol2 = 2.0 * (tol1 = tol * Math.Abs(x) + BlasLike.lm_eps);
                if (Math.Abs(x - xm) <= (tol2 - 0.5 * (b - a)))
                {
                    if (x > (gammatop + gammabot) * 0.5)
                        return Math.Min(Math.Max(gammatop, gammabot), ret_val + tol2);//If x is close to gammatop then x = gammatop is probably correct (same pieces throughout)
                    else
                        return Math.Max(Math.Min(gammatop, gammabot), ret_val - tol2);//If x is close to gammabot then x = gammabot
                }
                if (Math.Abs(e) > tol1)
                {
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0)
                        p = -p;
                    q = Math.Abs(q);
                    etemp = e;
                    e = d;
                    if (Math.Abs(p) >= Math.Abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                    {
                        d = (e = (x >= xm ? a - x : b - x)) * BlasLike.lm_golden_ratio;
                    }
                    else
                    {
                        d = p / q;
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                            d = BlasLike.dsign(tol1, xm - x);
                    }
                }
                else
                {
                    d = (e = (x >= xm ? a - x : b - x)) * BlasLike.lm_golden_ratio;
                }
                u = (Math.Abs(d) >= tol1 ? x + d : x + BlasLike.dsign(tol1, d));
                fu = OneDimensionalFunction(u, info) / basek;
                ret_val = u;
                if ((stopifpos * fu) > 0) break;
                if (fu <= fx)
                {
                    if (u >= x)
                        a = x;
                    else
                        b = x;
                    shifter(ref v, ref w, ref x, ref u);
                    shifter(ref fv, ref fw, ref fx, ref fu);
                }
                else
                {
                    if (u < x)
                        a = u;
                    else
                        b = u;
                    if (fu <= fw || w == x)
                    {
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u;
                        fv = fu;
                    }
                }
                count++;
            }
            if (count >= maxcount)
                ColourConsole.WriteError("PathMin did not converge");
            return ret_val;
        }
        ///<summary>
        ///Return the next lowest multiple of double base
        ///</summary>
        ///<param name="eps">Typically a convergence number</param>
        public static double small_round(double eps)
        {
            double reps = BlasLike.lm_eps * Math.Floor(eps / BlasLike.lm_eps);
            return (Math.Abs(eps - reps) > BlasLike.lm_eps * 8 ? eps : reps);//In case eps is too big
        }

        public hessmull h = null;
        double[] xnorm;
        int NROWA;
        bool UNITQ;
        double[] Q;
        double[] W;
        double[] A;
        double[] L;
        double[] U;
        double[] c;
        double AccuracyModify = -1.0;
        int msg;
        bool scldqp;
        int[] ISTATE;
        double[] LAMBDA;
        double[] FEATOL;
        int[] KACTV;
        int[] KFREE;
        int istart;
        double[] parm = new double[4];
        int nq;
        double asize;
        int NROWRT;
        int NCOLRT;
        double dtmax;
        double dtmin;
        double[] ANORM;
        double[] WRK;
        double[] LWRK;
        double[] RT;
        double[] ZY;
        double[] PX;
        double[] QTG;
        double[] RLAM;
        double[] AP;
        long timebase = -1233456;
        long timeaquired = -12233;
        void wrexit(string name, int inform, int iter)
        {
            ColourConsole.WriteInfo($"{name} Inform={inform} Iter={iter}");
        }
        void wdinky(string name, double ztgnrm, double dinky)
        {
            ColourConsole.WriteInfo($"\n//{name}//         ZTGNRM         DINKY\n//{name}//{ztgnrm}{dinky}");
            //lm_wmsg("\n//%s//         ZTGNRM         DINKY\n//%s//%14.5lg%14.5lg",
            //name,name,ztgnrm,dinky);
        }
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
        void dlpcore(bool lp, int minsum, bool orthog, int vertex, ref int inform, ref int iter,
               int itmax, byte lcrash, int n, int nclin, ref int nctotl, ref int nactiv,
               ref int nfree, ref int numinf, ref double obj, double[] xnorm)
        {
            /*
            lpcore finds a feasible point for the general linear constraints
            and bounds. the sum of the infeasibilities is minimized using
            a linear programming algorithm which may perform non-simplex
            steps. at each iteration the direction of search is defined as
            the projection of the steepest-descent direction. this
            projection is computed using an orthogonal factorization of the
            matrix of constraints in the working set
            if  lp = 1,  lpcore will solve the linear programming problem
            defined by the objective cvec, the constraint matrix  a	 and the
            bounds	bl, bu

            values of istate(j)...
            - 2	    - 1		0	    1	       2	 3
            a*x < bl   a*x > bu	a*x free   a*x = bl   a*x = bu	 bl = bu

            if  vertex = 1,	 the initial point  x  will be made into a
            vertex by regarding some of the free variables	x(j)  as being
            on an temporary bound.	some of these variables may remain on
            their temporary bounds.	 if so, their state will be
            istate(j) = 4
        */
            string lprob = "LP";

            int a_dim1, a_offset;

            int iadd, jadd;
            double alfa;
            int jdel, kdel, ndel;
            int ifix = 123;
            bool prnt;
            bool added;
            double palfa = 44;
            int ifail = 19;
            double bigdx;
            int isdel = 0;
            double objlp, condt;
            double anorm, dinky;
            int ncnln;
            bool stall;
            double wgfix = 1e9;
            int ncolz = 89;
            double pnorm = 1e23;
            bool nullr;
            //    int nrowj;
            int nclin0, kb;
            double bigalf, epspt9;
            double feamax = 1, feamin = 2;
            bool delete_;
            int nfixed;
            int jbigst = 8, kbigst = 45;
            double tolact;
            bool modfyg;
            double condmx, atphit = 45, cslast = 46, rdlast = 1e12, objsiz, snlast = 76, suminf = 12,
                 trulam = 99;
            int idummy, msglvl;
            int jsmlst = 1, ksmlst = 7;
            double smllst = 1e34;
            int mstall;
            double ztgnrm;
            int nstall;
            bool firstv;
            int hitlow = 6;
            bool unitpg;
            double bnd;
            double gtp = 1e45;

            //  --w;
            // --iw;
            // --x;
            //   --featol;
            // --cvec;
            // --clamda;
            //          --bu;
            // --bl;
            //         --ax;
            a_dim1 = NROWA;
            a_offset = a_dim1 + 1;
            //       a -= a_offset;
            //  --kfree;
            // --kactiv;
            //  --istate;
            /*     SPECIFY MACHINE-DEPENDENT PARAMETERS. */
            /*INITIALIZE */
            ncnln = 0;
            nclin0 = Math.Max(nclin, 1);
            //    nrowj = 1;
            inform = 0;
            iter = 0;
            jadd = 0;
            jdel = 0;
            ndel = 0;
            nstall = 0;
            numinf = 1;
            msglvl = msg;
            msg = 0;
            if (iter >= istart)
            {
                msg = msglvl;
            }
            var bigbnd = parm[0];
            bigdx = parm[1];
            tolact = parm[2];
            epspt9 = parm[3];
            alfa = 0;
            condmx = BlasLike.lm_max;
            objlp = 0;
            added = true;
            firstv = false;
            modfyg = true;
            nullr = true;
            unitpg = false;
            BlasLike.dxminmax(nctotl, FEATOL, 1, ref feamax, ref feamin);
            /*
            ---------------------------------------------------------------------
            given an initial point	x, compute the following....
            (1)	the initial working set
            (2)	the tq factorization of the matrix of constraints in the
                working set
            (3)	the value and gradient of the sum of infeasibilities at the
                point x. if x is feasible and the solution of an lp is
                required, the linear objective function and gradient is
                computed
                the array rlamda is used as temporary work space
            ---------------------------------------------------------------------
            */
            dlpcrsh(orthog, vertex, lcrash, n, nclin, nctotl,
     NROWRT, NCOLRT, ref nactiv, ref ncolz, ref nfree);
            dlpgrad(lp, n, nctotl, feamin, ref numinf, ref suminf);
            dzyprod(6, n, nactiv, ncolz, nfree, nq, KACTV, KFREE
                , QTG, WRK);
            obj = suminf;
            if (lp) objlp = BlasLike.ddotvec(n, c, W);
            if (lp && numinf == 0) obj = objlp;
            if (numinf == 0 && !(lp)) goto L320;
            /* .......................START OF THE MAIN LOOP........................ 
            */
            /*     DEFINE SMALL QUANTITIES THAT REFLECT THE MAGNITUDE OF  C,  X, */
            /*     AND THE NORM OF THE CONSTRAINTS IN THE WORKING SET. */
            L20:
            objsiz = (1 + Math.Abs(obj)) / (1 + xnorm[0]);
            if (numinf == 0)
            {
                objsiz = (BlasLike.lm_eps + Math.Abs(obj)) / (BlasLike.lm_eps + xnorm[0]);
            }
            anorm = 0;
            if (nactiv > 0)
            {
                anorm = Math.Abs(dtmax);
            }
            dinky = epspt9 * Math.Max(anorm, objsiz);
            /*     COMPUTE THE NORMS OF THE PROJECTED GRADIENT AND THE GRADIENT WITH 
            */
            /*     RESPECT TO THE FREE VARIABLES. */
            ztgnrm = 0;
            if (ncolz > 0) ztgnrm = dnrm2vec(ncolz, QTG);

            if (msg >= 80) wdinky("LPCORE", ztgnrm, dinky);
            delete_ = ztgnrm <= dinky;
            /*     PRINT THE DETAILS OF THIS ITERATION. */
            prnt = added || ndel > 1;
            if (!prnt) goto L40;

            condt = BlasLike.dprotdiv(ref dtmax, ref dtmin, ref ifail);
            if (ifail != 0 && dtmax == 0) condt = BlasLike.lm_max;
            dlpprt(lp, NROWRT, n, nclin, nfree, isdel, nactiv,
                ncolz, iter, jadd, jdel, alfa, condt, numinf,
                suminf, objlp);
            added = false;
            jadd = 0;
            jdel = 0;
        L40:
            if (numinf == 0 && !lp)
            {
                goto L320;
            }
            if (!delete_)
            {
                goto L100;
            }
            /* --------------------------------------------------------------------- 
            */
            /*     THE PROJECTED GRADIENT IS NEGLIGIBLE. */
            /*     WE HAVE TO DELETE A CONSTRAINT BEFORE A MOVE CAN BE MADE. */
            /* --------------------------------------------------------------------- 
            */
            dgetlamd(lprob, n, nactiv, ncolz, nfree,
                ref jsmlst, ref ksmlst, ref smllst);
            /* --------------------------------------------------------------------- 
            */
            /*     TEST FOR CONVERGENCE.  IF THE LEAST (ADJUSTED) MULTIPLIER IS */
            /*     GREATER THAN A SMALL NEGATIVE QUANTITY, AN ADEQUATE  LP  SOLUTION 
            */
            /*     HAS BEEN FOUND. */
            /* --------------------------------------------------------------------- 
            */
            if (smllst >= -dinky)
            {
                jsmlst = 0;
            }
            if (jsmlst == 0)
            {
                goto L60;
            }
            if (vertex != 0 && ncolz >= 1)
            {
                goto L60;
            }
            /*     PREPARE TO DELETE THE CONSTRAINT WITH INDEX  JSMLST. */
            jdel = jsmlst;
            kdel = ksmlst;
            isdel = ISTATE[jdel - 1];
            ISTATE[jdel - 1] = 0;
            goto L80;
        /* --------------------------------------------------------------------- 
        */
        /*     IF STILL INFEASIBLE, WE CAN REDUCE THE SUM OF INFEASIBILITIES */
        /*     IF THERE IS A MULTIPLIER GREATER THAN ONE. */
        /* --------------------------------------------------------------------- 
        */
        /*     INSTEAD OF LOOKING FOR THE LAST VIOLATED CONSTRAINT IN BNDALF, */
        /*     WE MUST NOW LOOK FOR THE FIRST VIOLATED CONSTRAINT ALONG  P. */
        /*     THIS WILL ENSURE THAT THE WEIGHTED SUM OF INFEASIBILITIES */
        /*     DECREASES. */
        L60:
            if (numinf == 0 || minsum != 0)
            {
                goto L280;
            }
            /*     FIND THE BIGGEST MULTIPLIER LARGER THAN UNITY. */
            /*     FOR THE PURPOSES OF THE TEST,  THE  J-TH  MULTIPLIER IS SCALED */
            /*     BY  FEATOL(J)/FEAMIN.  THIS FORCES CONSTRAINTS WITH LARGER  FEATOL 
            */
            /*     VALUES TO BE DELETED FIRST. */
            dlpbgst(n, nactiv, nfree, ref jbigst, ref kbigst, dinky, feamin,
                ref trulam);
            if (jbigst == 0)
            {
                goto L280;
            }
            jdel = jbigst;
            kdel = kbigst;
            isdel = ISTATE[jbigst - 1];

            ISTATE[jbigst - 1] = (trulam > 0) ? -2 : -1;
            firstv = true;
        /* --------------------------------------------------------------------- 
        */
        /*     UPDATE THE  TQ  FACTORIZATION OF THE MATRIX OF CONSTRAINTS IN THE 
        */
        /*     WORKING SET. */
        /* --------------------------------------------------------------------- 
        */
        L80:
            ++ndel;
            ddelcon(modfyg, orthog, jdel, kdel, nactiv, ncolz, nfree, n,
                 NROWRT);
            ++ncolz;
            if (jdel <= (int)n)
            {
                ++(nfree);
            }
            if (jdel > (int)n)
            {
                --nactiv;
            }
            goto L20;
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE SEARCH DIRECTION,  P = - Z*(PROJECTED GRADIENT). */
        /* --------------------------------------------------------------------- 
        */
        L100:
            if (iter >= itmax)
            {
                goto L400;
            }
            ++(iter);
            if (iter >= istart)
            {
                msg = msglvl;
            }
            dfindp(nullr, unitpg, n, nclin,
                NROWRT, ncolz, ncolz, ref nfree,
                delete_, ref gtp, ref pnorm, ref rdlast);
            /* --------------------------------------------------------------------- 
            */
            /*     FIND THE CONSTRAINT WE BUMP INTO ALONG  P. */
            /*     UPDATE  X  AND  AX  IF THE STEP  ALFA  IS NONZERO. */
            /* --------------------------------------------------------------------- 
            */
            /*     ALFA  IS INITIALIZED TO  BIGALF.  IF IT REMAINS THAT WAY AFTER */
            /*     THE CALL TO BNDALF, IT WILL BE REGARDED AS INFINITE. */

            bigalf = BlasLike.dprotdiv(ref bigdx, ref pnorm, ref ifail);
            if (ifail != 0 && bigdx == 0) bigalf = BlasLike.lm_max;
            inform = dbndalf(firstv, ref hitlow, ref jadd, n,
    nctotl, numinf, ref alfa, ref palfa, ref atphit, ref bigalf,
    pnorm);
            if (inform != 0 || jadd == 0)
            {
                goto L300;
            }
            /*     TEST IF  ALFA*PNORM  IS NEGLIGIBLE. */
            stall = (Math.Abs(alfa * pnorm)) <= epspt9 * xnorm[0];
            if (!stall)
            {
                goto L120;
            }
            /*     TAKE A ZERO STEP. */
            /*     IF A NON-ORTHOGONAL  TQ  FACTORIZATION IS BEING RECURRED AND  X */
            /*     IS NOT YET FEASIBLE,  THE GRADIENT OF THE SUM OF INFEASIBILITIES */

            /*     MUST BE RECOMPUTED. */
            alfa = 0;
            ++nstall;
            mstall = 2000;
            if (nstall <= mstall && orthog)
            {
                goto L160;
            }
            if (nstall <= mstall && !orthog)
            {
                goto L140;
            }
            goto L380;
        /*     CHANGE  X  TO  X + ALFA*P.  UPDATE  AX  ALSO. */
        L120:
            nstall = 0;
            BlasLike.daxpyvec(n, alfa, PX, W);

            if (nclin > 0)
                BlasLike.daxpyvec(nclin, alfa, AP, LWRK);
            xnorm[0] = dnrm2vec(n, W);
            if (lp) objlp = BlasLike.ddotvec(n, c, W);
            /*     IF  X  IS NOT YET FEASIBLE,  COMPUTE  OBJ  AND  GRAD  AS THE VALUE 
            */
            /*     AND GRADIENT OF THE SUM OF INFEASIBILITIES (IF  X  IS FEASIBLE, */
            /*     THE VECTOR  QTG  IS UPDATED AND  GRAD  NEED NOT BE COMPUTED). */
            L140:
            if (numinf == 0) goto L160;
            dlpgrad(lp, n, nctotl, feamin, ref numinf, ref suminf);
            if (!orthog && jadd <= (int)n)
            {
                wgfix = QTG[jadd - 1];
            }
            dzyprod(6, n, nactiv, ncolz, nfree, nq, KACTV, KFREE
                , QTG, WRK);
            obj = suminf;
        /* --------------------------------------------------------------------- 
        */
        /*     ADD A CONSTRAINT TO THE WORKING SET. */
        /* --------------------------------------------------------------------- 
        */
        /*     UPDATE  ISTATE. */
        L160:
            if (lp && numinf == 0)
            {
                obj = objlp;
            }
            if (hitlow != 0)
            {
                ISTATE[jadd - 1] = 1;
            }
            if (hitlow == 0)
            {
                ISTATE[jadd - 1] = 2;
            }
            if (L[jadd - 1] == U[jadd - 1])
            {
                ISTATE[jadd - 1] = 3;
            }
            /*     IF A BOUND IS TO BE ADDED, MOVE  X  EXACTLY ONTO IT, EXCEPT WHEN */

            /*     A NEGATIVE STEP WAS TAKEN.  (BNDALF  MAY HAVE HAD TO MOVE TO SOME 
            */
            /*     OTHER CLOSER CONSTRAINT.) */
            iadd = jadd - n;
            if (jadd > (int)n) goto L200;
            bnd = (hitlow != 0 ? L : U)[jadd - 1];
            if (alfa >= 0) W[jadd - 1] = bnd;
            for (ifix = 1; ifix <= nfree; ++ifix)
            {
                if (KFREE[ifix - 1] == jadd)
                {
                    goto L200;
                }
                /* L180: */
            }
        /*     UPDATE THE  TQ  FACTORS OF THE MATRIX OF CONSTRAINTS IN THE */
        /*     WORKING SET.  USE THE ARRAY  P  AS WORK SPACE. */
        L200:
            added = true;
            ndel = 0;
            inform = daddcon(modfyg, false, orthog, ifix, iadd, jadd,
                nactiv, ncolz, ncolz, nfree, n, NROWA, KFREE, condmx, cslast, snlast,
                 WRK, PX);
            --ncolz;
            nfixed = n - nfree;
            if (nfixed == 0)
            {
                goto L240;
            }
            kb = nactiv + nfixed;
            for (idummy = 1; idummy <= nfixed; ++idummy)
            {
                KACTV[kb] = KACTV[kb - 1];
                --kb;
                /* L220: */
            }
        L240:
            if (jadd <= (int)n)
            {
                /*
                ADD A BOUND.  IF STABILIZED ELIMINATIONS ARE BEING USED TO UPDATE 
                THE  TQ  FACTORIZATION,  RECOMPUTE THE COMPONENT OF THE GRADIENT
                CORRESPONDING TO THE NEWLY FIXED VARIABLE
                */
                --(nfree);
                KACTV[nactiv] = jadd;
                if (!orthog)
                {
                    if (numinf > 0) QTG[nfree - 1] = wgfix;
                    else if (lp) QTG[nfree - 1] = c[jadd - 1];
                }
            }
            else
            {
                /*ADD A GENERAL LINEAR CONSTRAINT*/
                ++nactiv;
                KACTV[nactiv - 1] = iadd;
            }
            goto L20;
        /* .........................END OF MAIN LOOP............................ 
        */
        /*     NO CONSTRAINTS TO DROP. */
        L280:
            if (numinf > 0)
            {
                goto L340;
            }
            goto L320;
        /*     ERROR IN  BNDALF  --  PROBABLY UNBOUNDED LP. */
        L300:
            if (numinf == 0)
            {
                goto L360;
            }
            goto L340;
        /*     FEASIBLE SOLUTION FOUND, OR OPTIMAL LP SOLUTION. */
        L320:
            inform = 0;
            goto L420;
        /*     THE LINEAR CONSTRAINTS AND BOUNDS APPEAR TO BE INFEASIBLE. */
        L340:
            inform = 1;
            goto L420;
        /*     UNBOUNDED LP. */
        L360:
            inform = 2;
            goto L420;
        /*     TOO MANY ITERATIONS WITHOUT CHANGING  X. */
        L380:
            inform = 3;
            goto L420;
        /*     TOO MANY ITERATIONS. */
        L400:
            inform = 4;
        /* --------------------------------------------------------------------- 
        */
        /*     PRINT FULL SOLUTION.  IF NECESSARY, RECOMPUTE THE MULTIPLIERS. */
        /* --------------------------------------------------------------------- 
        */
        L420:
            msg = msglvl;
            if (msg >= 1)
            {
                wrexit("LP", inform, iter);
            }
            if (inform > 0) dgetlamd(lprob, n, nactiv, ncolz, nfree,
                ref jsmlst, ref ksmlst, ref smllst);
            if (!lp && inform == 0) BlasLike.dzerovec(n, RLAM);
            dprtsol(nfree, n, nclin, ncnln, nctotl,
                nactiv);
        }
        double dnrm2vec(int n, double[] x, int xstart = 0)
        {
            if (n == 1) return (x[xstart] < 0.0 ? -x[xstart] : x[xstart]);
            else
            {
                double scale = 0.0, ssq = 1.0;
                BlasLike.dsssqvec(n, x, ref scale, ref ssq, xstart);
                return sc_norm(scale, ssq);
            }
        }
        double sc_norm(double scale, double ssq)
        {
            ssq = Math.Sqrt(ssq);
            return (scale < BlasLike.lm_max / ssq ? scale * ssq : BlasLike.lm_max);
        }
        void dlpprt(bool lp, int Nrowrt, int n, int nclin, int nfree, int isdel, int nactiv,
               int ncolz, int iter, int jadd, int jdel, double alfa, double condt,
               int numinf, double suminf, double objlp)
        {

            char[] lstate = { ' ', 'L', 'U', 'E', 'T' };

            char ladd, ldel;
            int inct, k, laddi, ldeli;

            //    --wrk2;
            //    --wrk1;
            //    --x;
            //rt -= Nrowrt + 1;
            var rt_offset = Nrowrt + 1;
            //  a -= NROWA + 1;
            var a_offset = NROWA + 1;
            //       --kfree;
            // --istate;

            if (msg < 5)
            {
                return;
            }
            ldeli = 0;
            laddi = 0;
            if (jdel > 0)
            {
                ldeli = isdel;
            }
            if (jadd > 0)
            {
                laddi = ISTATE[jadd];
            }
            ldel = lstate[ldeli];
            ladd = lstate[laddi];
            if (msg < 15)
            {
                /*
                ------------------------------------------------------------------
                print heading (possibly) and terse line
                ------------------------------------------------------------------
                */
                if (iter == 0)
                {//Haven't bother to program lpprnt
                 //    lpprnt(lp, iter, jdel, ldel, jadd, ladd, alfa, condt, numinf, suminf, objlp);
                    return;
                }
            }
            /*
            ---------------------------------------------------------------------
            print terse line,  x,  istate  and  kfree
            ---------------------------------------------------------------------
            */
            itrwri("LP", iter);
            //   lpprnt(lp, iter, jdel, ldel, jadd, ladd, alfa, condt, numinf, suminf, objlp);
            lm_mdvwri("\nLP VARIABLES", n, W);
            lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", n, ISTATE);
            if (nclin > 0)
                lm_mdvwri("\nSTATUS OF THE LP GENERAL CONSTRAINTS", nclin,
                    ISTATE, n);
            if (nfree > 0) lm_mdvwri("\nLIST OF FREE LP VARIABLES", nfree, KFREE);
            /*
            ---------------------------------------------------------------------
            compute and print  ax.	use  work = ap	to avoid side effects
            ---------------------------------------------------------------------
            */
            if (msg < 20)
            {
                return;
            }
            if (nclin > 0)
            {
                for (k = 1; k <= nclin; ++k)
                {
                    AP[k - 1] = BlasLike.ddot(n, A, NROWA, W, 1, k + NROWA - a_offset);
                }
                lm_mdvwri("\nVALUES OF LP GENERAL LINEAR CONSTRAINTS", nclin, AP);
            }
            /*
            ---------------------------------------------------------------------
            print the diagonals of	t
            ---------------------------------------------------------------------
            */
            if (msg >= 30)
            {
                inct = Nrowrt - 1;
                if (nactiv != 0)
                {
                    BlasLike.dcopy(nactiv, RT, inct, WRK, 1, nactiv + (ncolz + 1) * Nrowrt - rt_offset);
                    lm_mdvwri("\nDIAGONALS OF LP WORKING SET FACTOR  T", nactiv, WRK);
                }
            }
        }
        void lm_mdvwri<T>(string nn, int na, T[] wrk, int wstart = 0)
        {
            if (wrk == null) return;
            ColourConsole.WriteInfo(nn);
            for (int i = 0; i < na; ++i)
            {
                Console.Write($"{wrk[i + wstart]} ");
                if (i % 6 == 5) Console.Write("\n");
            }
            Console.Write("\n");
        }
        void itrwri(string name, int n)
        {
            ColourConsole.WriteInfo($"{name} iteration {n}");
        }
        void dlpcrsh(bool orthog, int vertex, byte lcrash, int n, int nclin, int nctotl,
         int Nrowrt, int Ncolrt, ref int nactiv, ref int ncolz,
       ref int nfree)
        {

            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1,
                i__2;
            int c__1 = 1;
            int n__ = n;
            var bigbnd = parm[0];


            int iadd = 34, jadd;
            double amin;
            int imin = 1000000000, jmin = 1000000000, ifix = -111111111;
            double resl, resu;
            int nact1;
            int i, j, k;
            int idiag;
            double b1, b2;
            double rnorm;
            bool nolow, noupp;
            int ncolz1, kb;
            int was_is, nfixed;
            double colmin, toobig;
            int nartif;
            double condmx, cslast = 78;
            int inform = 10;
            double resmin, colsiz, snlast = 0;
            int idummy;
            double rowmax;
            double bnd, res;
            double tolact = parm[2];


            //      --wrk2;
            //         --wrk1;
            //        --p;
            zy_dim1 = nq;
            zy_offset = zy_dim1 + 1;
            //        zy -= zy_offset;
            i__1 = Nrowrt * Ncolrt;
            BlasLike.dzerovec(i__1, RT);
            rt_dim1 = Nrowrt;
            rt_offset = rt_dim1 + 1;
            //         rt -= rt_offset;
            // --qtg;
            // --x;
            // --bu;
            //--bl;
            // --ax;
            //   --anorm;
            a_dim1 = NROWA;
            a_offset = a_dim1 + 1;
            //        a -= a_offset;
            //  --kfree;
            //    --kactiv;
            //     --istate;

            /*
            set the maximum allowable condition estimator of the constraints
            in the working set.  note that the conservative value used in
            lpcrsh is smaller than that used when a constraint is added to
            the working set during a typical iteration
            */
            condmx = 1 / BlasLike.lm_rooteps;
            if (msg >= 80)
            {
                //lm_wmsg("\n//LPCRSH//  LCRASH NCLIN NCTOTL\n//LPCRSH//%7d%7ld%7ld",
                //    lcrash, CL(*nclin), CL(*nctotl));
                ColourConsole.WriteInfo($"{lcrash},{nclin}, {nctotl}");
                lm_mdvwri("\nLP VARIABLES BEFORE CRASH...", n, W);
                lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", nctotl, ISTATE);
            }
            nfixed = 0;
            nactiv = 0;
            nartif = 0;
            /*     IF A COLD START IS BEING MADE, INITIALIZE  ISTATE. */
            /*     IF  BL(J) = BU(J),  SET  ISTATE(J)=3  FOR ALL VARIABLES AND LINEAR 
            */
            /*     CONSTRAINTS. */
            if (lcrash > 0) goto L60;
            i__1 = nctotl;
            for (j = 1; j <= i__1; ++j)
                ISTATE[j - 1] = L[j - 1] == U[j - 1] ? 3 : 0;
            /* L40: */
            /*     INITIALIZE  NFIXED,  NACTIV  AND  KACTIV. */
            /*     ENSURE THAT THE NUMBER OF BOUNDS AND GENERAL CONSTRAINTS IN THE */
            /*     WORKING SET DOES NOT EXCEED  N. */
            L60:
            i__1 = nctotl;
            for (j = 1; j <= i__1; ++j)
            {
                if (nfixed + nactiv == (int)n || ISTATE[j - 1] == 4) ISTATE[j - 1] = 0;
                if (ISTATE[j - 1] > 0)
                {
                    if (j <= (int)n)
                    {
                        ++nfixed;
                        if (ISTATE[j - 1] == 1) W[j - 1] = L[j - 1];
                        if (ISTATE[j - 1] >= 2) W[j - 1] = U[j - 1];
                    }
                    else
                    {
                        ++nactiv;
                        if (lcrash < 2) KACTV[nactiv - 1] = j - n;
                    }
                }
            }
            nfree = n - nfixed;
            ncolz = nfree - nactiv;
            if (msg >= 80)
            {
                lm_mdvwri("\nLP VARIABLES AFTER CRASH INIT...", n, W);
                lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", nctotl, ISTATE);
            }
            /*if a hot start is required, the tq factorization is already known*/
            if (lcrash > 1)
            {
                goto L600;
            }
            dtmax = 1;
            dtmin = 1;
            //    *unitq = 1;
            UNITQ = true;
            /*     COMPUTE THE 2-NORMS OF THE CONSTRAINT ROWS. */
            asize = 1;
            if (nclin == 0)
            {
                goto L140;
            }
            i__1 = nclin;
            for (j = 1; j <= i__1; ++j)
            {
                //        anorm[j] = dnrm2(n, &a[j + a_dim1], *nrowa);
                ANORM[j - 1] = dnrm2(n, A, NROWA, j + a_dim1 - a_offset);
                /* L120: */
            }
            amin = 1e10;
            BlasLike.dxminmax(nclin, ANORM, 1, ref asize, ref amin);
        L140:
            if (lcrash > 0) goto L320;
            /*
            ---------------------------------------------------------------------
            if a cold start is required, an attempt is made to add as many
            constraints as possible to the working set
            ---------------------------------------------------------------------
            */
            if (nfixed + nactiv == (int)n)
            {
                goto L460;
            }
            /* see if any variables are outside their bounds */
            for (j = 1; j <= (int)n; ++j)
            {
                if (ISTATE[j - 1] != 0) continue;
                b1 = L[j - 1];
                b2 = U[j - 1];
                nolow = b1 <= -bigbnd;
                noupp = b2 >= bigbnd;

                was_is = 0;
                if (nolow)
                {
                    goto L160;
                }
                if (W[j - 1] - b1 <= (1 + Math.Abs(b1)) * tolact)
                {
                    was_is = 1;
                }
            L160:
                if (noupp)
                {
                    goto L180;
                }
                if (b2 - W[j - 1] <= (1 + Math.Abs(b2)) * tolact)
                {

                    was_is = 2;
                }
            L180:
                if (was_is == 0)
                {
                    continue;
                }
                /*        SET VARIABLE EQUAL TO ITS BOUND. */
                ISTATE[j - 1] = was_is;
                if (was_is == 1) W[j - 1] = b1;
                if (was_is == 2) W[j - 1] = b2;
                ++nfixed;
                if (nfixed + nactiv == (int)n) goto L460;
            }
            /* --------------------------------------------------------------------- 
            */
            /*     THE FOLLOWING LOOP FINDS THE LINEAR CONSTRAINT THAT IS CLOSEST */
            /*     TO BEING SATISFIED.  IF IT IS SUFFICIENTLY CLOSE, IT WILL BE ADDED 
            */
            /*     TO THE WORKING SET AND THE PROCESS WILL BE REPEATED. */
            /* --------------------------------------------------------------------- 
            */
            /*     FIRST, COMPUTE  AX  FOR INEQUALITY LINEAR CONSTRAINTS. */
            if (nclin == 0) goto L320;
            i__1 = nclin;
            for (i = 1; i <= i__1; ++i)
            {
                j = n + i;
                if (ISTATE[j - 1] > 0)
                {
                    goto L220;
                }
                LWRK[i - 1] = BlasLike.ddot(n, A, NROWA, W, 1, i + a_dim1 - a_offset);
            L220:
                ;
            }
            toobig = tolact + tolact;
            i__1 = n;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                resmin = toobig;
                was_is = 0;
                i__2 = nclin;
                for (i = 1; i <= i__2; ++i)
                {
                    j = n + i;
                    if (ISTATE[j - 1] > 0)
                    {
                        goto L280;
                    }
                    b1 = L[j - 1];
                    b2 = U[j - 1];
                    nolow = b1 <= -bigbnd;
                    noupp = b2 >= bigbnd;
                    resl = toobig;
                    resu = toobig;
                    if (nolow)
                    {
                        goto L240;
                    }
                    resl = Math.Abs(LWRK[i - 1] - b1) / (1 + Math.Abs(b1));
                L240:
                    if (noupp)
                    {
                        goto L260;
                    }
                    resu = Math.Abs(LWRK[i - 1] - b2) / (1 + Math.Abs(b2));
                L260:
                    res = Math.Min(resl, resu);
                    if (res >= tolact)
                    {
                        goto L280;
                    }
                    if (res >= resmin)
                    {
                        goto L280;
                    }
                    resmin = res;
                    imin = i;
                    was_is = 1;
                    if (resl > resu)
                    {
                        was_is = 2;
                    }
                L280:
                    ;
                }
                if (was_is == 0)
                {
                    goto L320;
                }
                ++nactiv;
                KACTV[nactiv - 1] = imin;
                j = n + imin;
                ISTATE[j - 1] = was_is;
                if (nfixed + nactiv == (int)n)
                {
                    goto L460;
                }
                /* L300: */
            }
        /* --------------------------------------------------------------------- 
        */
        /*     IF NECESSARY, ADD TEMPORARY BOUNDS TO MAKE A VERTEX. */
        /* --------------------------------------------------------------------- 
        */
        L320:
            ncolz = n - nfixed - nactiv;
            if (vertex == 0 || ncolz == 0)
            {
                goto L460;
            }
            /*     COMPUTE LENGTHS OF COLUMNS OF SELECTED LINEAR CONSTRAINTS */
            /*     (JUST THE ONES CORRESPONDING TO FREE VARIABLES). */
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (ISTATE[j - 1] != 0)
                {
                    goto L380;
                }
                colsiz = 0;
                if (nclin == 0)
                {
                    goto L360;
                }
                i__2 = nclin;
                for (k = 1; k <= i__2; ++k)
                {
                    i = n + k;
                    if (ISTATE[i - 1] > 0)
                    {
                        colsiz += Math.Abs(A[k + j * a_dim1 - a_offset]);
                    }
                    /* L340: */
                }
            L360:
                WRK[j - 1] = colsiz;
            L380:
                ;
            }
            /*     FIND THE  NARTIF  SMALLEST SUCH COLUMNS. */
            /*     THIS IS AN EXPENSIVE LOOP.  LATER WE CAN REPLACE IT */
            /*     BY A 4-PASS PROCESS (SAY), ACCEPTING THE FIRST COL THAT */
            /*     IS WITHIN  T  OF  COLMIN, WHERE  T = 0.0, 0.001, 0.01, 0.1 (SAY). 
            */
            i__1 = ncolz;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                colmin = BlasLike.lm_max;
                i__2 = n;
                for (j = 1; j <= i__2; ++j)
                {
                    if (ISTATE[j - 1] != 0)
                    {
                        goto L400;
                    }
                    if (nclin == 0)
                    {
                        goto L420;
                    }
                    colsiz = WRK[j - 1];
                    if (colmin <= colsiz)
                    {
                        goto L400;
                    }
                    colmin = colsiz;
                    jmin = j;
                L400:
                    ;
                }
                j = jmin;
            L420:
                ISTATE[j - 1] = 4;
                ++nartif;
                /* L440: */
            }
        /* --------------------------------------------------------------------- 
        */
        /*     A TRIAL WORKING SET HAS NOW BEEN SELECTED. */
        /* --------------------------------------------------------------------- 
        */
        /*     SET  KFREE  TO POINT TO THE FREE VARIABLES. */
        L460:
            nfree = 0;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (ISTATE[j - 1] != 0)
                {
                    goto L480;
                }
                ++(nfree);
                KFREE[nfree - 1] = j;
            L480:
                ;
            }
            /*     COMPUTE THE TQ FACTORIZATION FOR THE SELECTED LINEAR CONSTRAINTS. 
            */
            /*     FIRST, THE COLUMNS CORRESPONDING TO SIMPLE BOUNDS IN THE WORKING */

            /*     SET ARE REMOVED. THE RESULTING  NACTIV BY NFREE  MATRIX (NACTIV */
            /*     LE NFREE) IS FACTORIZED BY ADDING THE CONSTRAINTS ONE AT A TIME */
            /*     AND UPDATING USING PLANE ROTATIONS OR STABILIZED ELIMINATIONS. */
            /*     THE  NACTIV BY NACTIV  TRIANGULAR MATRIX  T  AND THE NFREE BY */
            /*     NFREE MATRIX  Q  ARE STORED IN THE ARRAYS  RT  AND  ZY. */
            ncolz = nfree;
            if (nactiv == 0)
            {
                goto L540;
            }
            nact1 = nactiv;
            nactiv = 0;
            dtqadd(orthog, ref inform, c__1, nact1, ref nactiv, ref ncolz, ref nfree, n__,
                      condmx);
            /*     IF A VERTEX IS REQUIRED BUT  TQADD  WAS UNABLE TO ADD ALL OF THE */

            /*     SELECTED GENERAL CONSTRAINTS, ADD MORE TEMPORARY BOUNDS. */
            if (vertex == 0 || ncolz == 0)
            {
                goto L540;
            }
            ncolz1 = ncolz;
            i__1 = ncolz1;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                rowmax = 0;
                i__2 = nfree;
                for (i = 1; i <= i__2; ++i)
                {
                    //       rnorm = dnrm2(*ncolz, &zy[i + zy_dim1], *Nq);
                    rnorm = dnrm2(ncolz, ZY, nq, i + zy_dim1 - zy_offset);
                    if (rowmax >= rnorm)
                    {
                        goto L500;
                    }
                    rowmax = rnorm;
                    ifix = i;
                L500:
                    ;
                }
                jadd = KFREE[ifix - 1];
                inform = daddcon(false, false, orthog, ifix, iadd, jadd,
                    nactiv, ncolz, ncolz, nfree, n, NROWA,
                    KFREE, condmx, cslast, snlast,
                     WRK, RLAM);
                --(nfree);
                --(ncolz);
                ++nartif;
                ISTATE[jadd - 1] = 4;
                /* L520: */
            }
        /*     SET ELEMENTS  NACTIV + 1, ......, NACTIV + NFIXED  OF  KACTIV TO */

        /*     POINT TO THE FIXED VARIABLES. */
        L540:
            kb = nactiv;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (ISTATE[j - 1] == 0) continue;
                ++kb;
                KACTV[kb - 1] = j;
            }
            /* --------------------------------------------------------------------- 
            */
            /*     THE TQ FACTORIZATION HAS BEEN COMPUTED.  FIND THE POINT CLOSEST TO 
            */
            /*     THE USER-SUPPLIED  X  THAT LIES ON THE INITIAL WORKING SET. */
            /* --------------------------------------------------------------------- 
            */
            /*     SET WRK1 = RESIDUALS FOR CONSTRAINTS IN THE WORKING SET. */
            if (nactiv == 0)
            {
                goto L600;
            }
            i__1 = nactiv;
            for (i = 1; i <= i__1; ++i)
            {
                k = KACTV[i - 1];
                j = n + k;
                bnd = L[j - 1];
                if (ISTATE[j - 1] > 1)
                {
                    bnd = U[j - 1];
                }
                WRK[i - 1] = bnd - BlasLike.ddot(n, A, NROWA, W, 1, k + a_dim1 - a_offset);
                /* L580: */
            }
            /*     SOLVE FOR P, THE SMALLEST CORRECTION TO X THAT GIVES A POINT */
            /*     ON THE CONSTRAINTS IN THE WORKING SET. */
            /*     FIRST SOLVE  T*WRK1 = RESIDUALS, THEN GET  P = Y*WRK1. */
            idiag = 1;
            {
                drtmxsolve(2, nactiv, RT, Nrowrt, WRK,
                   ref idiag, (ncolz + 1) * rt_dim1 + 1 - rt_offset);
                BlasLike.dzerovec(n, PX);
                BlasLike.dcopyvec(nactiv, WRK, PX, 0, ncolz);
                dzyprod(2, n, nactiv, ncolz, nfree, nq, KACTV, KFREE,
                    PX, WRK);
                BlasLike.daxpyvec(n, 1, PX, W);
            }
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE 2-NORM OF  X. */
        /*     INITIALIZE  AX  FOR ALL GENERAL CONSTRAINTS. */
        /* --------------------------------------------------------------------- 
        */
        L600:
            xnorm[0] = dnrm2vec(n, W);
            if (nclin == 0)
            {
                goto L640;
            }
            BlasLike.dzerovec(nclin, LWRK);
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (W[j - 1] != 0)
                {
                    //   BlasLike.daxpyvec(*nclin, W[j-1], &a[j * a_dim1 + 1],  &ax[1]);
                    BlasLike.daxpyvec(nclin, W[j - 1], A, LWRK, j * a_dim1 + 1 - a_offset);
                }
                /* L620: */
            }
        /*     A POINT THAT SATISFIES THE INITIAL WORKING SET HAS BEEN FOUND. */
        L640:
            ncolz = nfree - nactiv;
            nfixed = n - nfree;
            if (msg >= 80)
            {
                //  lm_wmsg(
                //"\nLPCRSH. WORKING SET SELECTED ...\nBOUNDS = %ld TEMPORARY BOUNDS = %ld GENERAL LINEAR = %ld",
                //  CL(nfixed), CL(nartif), CL(*nactiv));
                ColourConsole.WriteInfo($"{nfixed},{nartif}, {nactiv}");
                lm_mdvwri("\nLP VARIABLES AFTER  CRASH...", n, W);
            }
        }
        void dlpgrad(bool lp, int n, int nctotl, double feamin,
       ref int numinf, ref double suminf)
        {

            /*
                if numinf > 0,	lpgrad	finds the number and weighted sum of
                infeasibilities for the bounds and linear constraints. an
                appropriate gradient vector is returned in  grad
                if numinf = 0,	and if an lp problem is being solved,  grad  will
                be loaded with the true linear objective
                positive values of  istate(j)  will not be altered.  these mean
                the following..

                    1	   2	     3
                    a*x = bl   a*x = bu   bl = bu

                other values of	 istate(j)  will be reset as follows..

                    a*x < bl   a*x > bu	a*x free
                    - 2	    - 1		0
            */
            int j, k;
            double s, feasj;
            bool nolow, noupp;
            double weight, atx = 1e9;
            var bigbnd = parm[0];

            // --x;
            // --grad;
            // --featol;
            //      --cvec;
            // --bu;
            //       --bl;
            //        a -= NROWA + 1;
            var a_offset = NROWA + 1;
            //      --istate;

            if (numinf != 0)
            {
                numinf = 0;
                suminf = 0;
                BlasLike.dzerovec(n, QTG);
                for (j = 1; j <= (int)nctotl; ++j)
                {
                    /* do nothing if the variable or constraint is at a bound */
                    if (ISTATE[j - 1] > 0) continue;
                    feasj = FEATOL[j - 1];
                    nolow = L[j - 1] <= -bigbnd;
                    noupp = U[j - 1] >= bigbnd;
                    k = j - n;
                    if (j <= (int)n)
                    {
                        atx = W[j - 1];
                    }
                    if (j > (int)n)
                    {
                        atx = BlasLike.ddot(n, A, NROWA, W, 1, k + NROWA - a_offset);
                    }
                    ISTATE[j - 1] = 0;
                    /* see if the lower bound is violated */
                    if (nolow)
                    {
                        goto L20;
                    }
                    s = L[j - 1] - atx;
                    if (s <= feasj)
                    {
                        goto L20;
                    }
                    ISTATE[j - 1] = -2;
                    weight = -feamin / feasj;
                    goto L40;
                /*see if the upper bound is violated*/
                L20:
                    if (noupp)
                    {
                        continue;
                    }
                    s = atx - U[j - 1];
                    if (s <= feasj)
                    {
                        continue;
                    }
                    ISTATE[j - 1] = -1;
                    weight = feamin / feasj;
                /* add the infeasibility */
                L40:
                    ++(numinf);
                    suminf += Math.Abs(weight) * s;
                    if (j <= (int)n)
                    {
                        QTG[j - 1] = weight;
                    }
                    if (j > (int)n) BlasLike.daxpy(n, weight, A, NROWA, QTG, 1, k + NROWA - a_offset);
                }
            }
            /* if feasible, install true objective */
            if (lp && numinf == 0) BlasLike.dcopyvec(n, c, QTG);
        }
        void dzyprod(short mode, int n, int nactiv, int ncolz, int nfree, int nq, int[] kactiv, int[] kfree, double[] v, double[] wrk)
        {
            /*
                dzyprod transforms the vector  v  in various ways using the 
                matrix  Q = ( Z  Y )  defined by the input parameters
                MODE	RESULT 
                ----	------ 
                1	V = Z*V 
                2	V = Y*V 
                3	V = Q*V       (NOT YET USED) 
                            ON INPUT,  V  IS ASSUMED TO BE ORDERED AS  ( V(FREE)  V(FIXED) ). 
                            ON OUTPUT, V  IS A FULL N-VECTOR. 

                4	V = Z(T)*V 
                5	V = Y(T)*V 
                6	V = Q(T)*V 
                            ON INPUT,  V  IS A FULL N-VECTOR. 
                            ON OUTPUT, V  IS ORDERED AS  ( V(FREE)  V(FIXED) ). 

                7	V = Y(T)*V 
                8	V = Q(T)*V 
                            ON INPUT,  V  IS A FULL N-VECTOR. 
                            ON OUTPUT, V  IS AS IN MODES 5 AND 6 EXCEPT THAT V(FIXED) IS NOT SET. 

                BEWARE THAT  NCOLZ  WILL SOMETIMES BE  NCOLR. 
                ALSO, MODES  1, 4, 7 AND 8  DO NOT INVOLVE  V(FIXED).
                NACTIV  AND  THE ARRAY  KACTIV  ARE NOT USED FOR THOSE CASES. 
            */
            int lenv;
            int j, k, l;
            int j1, j2, ka, kw, nfixed;
            // --wrk;
            //   zy -= nq + 1;
            int zy_offset = nq + 1;

            //        --v;
            //   --kfree;
            //        --kactiv;

            nfixed = n - nfree;
            j1 = 1;
            j2 = nfree;
            if (mode == 1 || mode == 4) j2 = ncolz;
            else if (mode == 2 || mode == 5 || mode == 7) j1 = ncolz + 1;
            lenv = j2 - j1 + 1;
            if (mode < 4)
            {
                /*MODE = 1, 2  OR  3*/
                if (nfree > 0) BlasLike.dzerovec(nfree, wrk);
                /*COPY  V(FIXED)  INTO THE END OF  WRK. */
                if (mode != 1 && nfixed != 0)
                    BlasLike.dcopyvec(nfixed, v, wrk, nfree, nfree);
                /*SET  WRK  =  RELEVANT PART OF  ZY * V. */
                if (lenv > 0)
                {
                    if (UNITQ) BlasLike.dcopyvec(lenv, v, wrk, j1 - 1, j1 - 1);
                    else for (j = j1; j <= j2; ++j)
                            if (v[j - 1] != 0)
                                BlasLike.daxpy(nfree, v[j - 1], ZY, 1, wrk, 1, j * nq + 1 - zy_offset);
                }
                /*EXPAND  WRK  INTO  V  AS A FULL N-VECTOR. */
                BlasLike.dzerovec(n, v);
                if (nfree > 0)
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k - 1];
                        v[j - 1] = wrk[k - 1];
                    }
                /*COPY  WRK(FIXED)  INTO THE APPROPRIATE PARTS OF  V*/
                if (mode == 1 || nfixed == 0) return;
                for (l = 1; l <= nfixed; ++l)
                {
                    kw = nfree + l;
                    ka = nactiv + l;
                    j = kactiv[ka - 1];
                    v[j - 1] = wrk[kw - 1];
                }
            }
            else
            {
                /*
                MODE = 4, 5, 6, 7  OR  8
                PUT THE FIXED COMPONENTS OF  V  INTO THE END OF  WRK
                */
                if (mode != 4 && mode <= 6 && nfixed != 0)
                    for (l = 1; l <= nfixed; ++l)
                    {
                        kw = nfree + l;
                        ka = nactiv + l;
                        j = kactiv[ka - 1];
                        wrk[kw - 1] = v[j - 1];
                    }
                /*PUT THE FREE  COMPONENTS OF  V  INTO THE BEGINNING OF  WRK. */
                if (nfree != 0)
                {
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k - 1];
                        wrk[k - 1] = v[j - 1];
                    }
                    /*SET  V  =  RELEVANT PART OF  ZY(T) * WRK*/
                    if (lenv > 0)
                    {
                        if (UNITQ) BlasLike.dcopyvec(lenv, wrk, v, j1 - 1, j1 - 1);
                        else for (j = j1; j <= j2; ++j)
                                v[j - 1] = BlasLike.ddotvec(nfree, ZY, wrk, j * nq + 1 - zy_offset);
                    }
                }
                /*COPY THE FIXED COMPONENTS OF WRK INTO THE END OF  V*/
                if (mode != 4 && mode <= 6 && nfixed != 0)
                    BlasLike.dcopyvec(nfixed, wrk, v, nfree, nfree);
            }
        }

        void dzyprod(short mode, int n, int nactiv, int ncolz, int nfree, int nq, int[] kactiv, int[] kfree, double[] v, double[] wrk, int kc = 0, int kfr = 0, int vst = 0, int wst = 0)
        {
            /*
                dzyprod transforms the vector  v  in various ways using the 
                matrix  Q = ( Z  Y )  defined by the input parameters
                MODE	RESULT 
                ----	------ 
                1	V = Z*V 
                2	V = Y*V 
                3	V = Q*V       (NOT YET USED) 
                            ON INPUT,  V  IS ASSUMED TO BE ORDERED AS  ( V(FREE)  V(FIXED) ). 
                            ON OUTPUT, V  IS A FULL N-VECTOR. 

                4	V = Z(T)*V 
                5	V = Y(T)*V 
                6	V = Q(T)*V 
                            ON INPUT,  V  IS A FULL N-VECTOR. 
                            ON OUTPUT, V  IS ORDERED AS  ( V(FREE)  V(FIXED) ). 

                7	V = Y(T)*V 
                8	V = Q(T)*V 
                            ON INPUT,  V  IS A FULL N-VECTOR. 
                            ON OUTPUT, V  IS AS IN MODES 5 AND 6 EXCEPT THAT V(FIXED) IS NOT SET. 

                BEWARE THAT  NCOLZ  WILL SOMETIMES BE  NCOLR. 
                ALSO, MODES  1, 4, 7 AND 8  DO NOT INVOLVE  V(FIXED).
                NACTIV  AND  THE ARRAY  KACTIV  ARE NOT USED FOR THOSE CASES. 
            */
            int lenv;
            int j, k, l;
            int j1, j2, ka, kw, nfixed;
            //           --wrk;
            wst--;
            //   zy -= nq + 1;
            int zy_offset = nq + 1;

            //          --v;
            vst--;
            //       --kfree;
            kfr--;
            //     --kactiv;
            kc--;

            nfixed = n - nfree;
            j1 = 1;
            j2 = nfree;
            if (mode == 1 || mode == 4) j2 = ncolz;
            else if (mode == 2 || mode == 5 || mode == 7) j1 = ncolz + 1;
            lenv = j2 - j1 + 1;
            if (mode < 4)
            {
                /*MODE = 1, 2  OR  3*/
                if (nfree > 0) BlasLike.dzerovec(nfree, wrk, 1 + wst);
                /*COPY  V(FIXED)  INTO THE END OF  WRK. */
                if (mode != 1 && nfixed != 0)
                    //   BlasLike.dcopyvec(nfixed, &v[nfree + 1], &wrk[nfree + 1]);
                    BlasLike.dcopyvec(nfixed, v, wrk, nfree + 1 + vst, nfree + 1 + wst);
                /*SET  WRK  =  RELEVANT PART OF  ZY * V. */
                if (lenv > 0)
                {
                    //  if (UNITQ) BlasLike.dcopyvec(lenv, &v[j1], &wrk[j1]);
                    if (UNITQ) BlasLike.dcopyvec(lenv, v, wrk, j1 + vst, j1 + wst);
                    else for (j = j1; j <= j2; ++j)
                            if (v[j + vst] != 0)
                                // BlasLike.daxpy(nfree, v[j], &pZY[j * nq + 1 - zy_offset], 1, &wrk[1], 1);
                                BlasLike.daxpyvec(nfree, v[j + vst], ZY, wrk, j * nq + 1 - zy_offset, 1 + wst);
                }
                /*EXPAND  WRK  INTO  V  AS A FULL N-VECTOR. */
                BlasLike.dzerovec(n, v, 1 + vst);
                if (nfree > 0)
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k + kfr];
                        v[j + vst] = wrk[k + wst];
                    }
                /*COPY  WRK(FIXED)  INTO THE APPROPRIATE PARTS OF  V*/
                if (mode == 1 || nfixed == 0) return;
                for (l = 1; l <= nfixed; ++l)
                {
                    kw = nfree + l;
                    ka = nactiv + l;
                    j = kactiv[ka + kc];
                    v[j + vst] = wrk[kw + wst];
                }
            }
            else
            {
                /*
                MODE = 4, 5, 6, 7  OR  8
                PUT THE FIXED COMPONENTS OF  V  INTO THE END OF  WRK
                */
                if (mode != 4 && mode <= 6 && nfixed != 0)
                    for (l = 1; l <= nfixed; ++l)
                    {
                        kw = nfree + l;
                        ka = nactiv + l;
                        j = kactiv[ka + kc];
                        wrk[kw + wst] = v[j + vst];
                    }
                /*PUT THE FREE  COMPONENTS OF  V  INTO THE BEGINNING OF  WRK. */
                if (nfree != 0)
                {
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k + kfr];
                        wrk[k + wst] = v[j + vst];
                    }
                    /*SET  V  =  RELEVANT PART OF  ZY(T) * WRK*/
                    if (lenv > 0)
                    {
                        if (UNITQ) BlasLike.dcopyvec(lenv, wrk, v, j1 + wst, j1 + vst);
                        else for (j = j1; j <= j2; ++j)
                                v[j + vst] = BlasLike.ddotvec(nfree, ZY, wrk, j * nq + 1 - zy_offset, 1 + vst);
                    }
                }
                /*COPY THE FIXED COMPONENTS OF WRK INTO THE END OF  V*/
                if (mode != 4 && mode <= 6 && nfixed != 0)
                    BlasLike.dcopyvec(nfixed, wrk, v, nfree + 1 + wst, nfree + 1 + vst);
            }
        }
        void dgetlamd(string lprob, int n, int nactiv, int ncolz, int nfree, ref int jsmlst, ref int ksmlst, ref double smllst)
        {/*
	dgetlamd first computes the lagrange multiplier estimates for the
	given working set.  it then determines the values and indices of
	certain significant multipliers.  in this process, the multipliers
	for inequalities at their upper bounds are adjusted so that a
	negative multiplier for an inequality constraint indicates
	non-optimality.  in the following, the term minimum refers to the
	ordering of numbers on the real line, and not to their magnitude.
	smllst  is the minimum among the inequality constraints of the
	(adjusted) multipliers scaled by the 2-norm of the
	associated constraint row
	jsmlst  is the index of the constraint corresponding to  smllst
	ksmlst  marks its position in  kactiv
	on exit,  elements  1  thru  nactiv  of   rlamda  contain the
	(unadjusted) multipliers for the general constraints.  elements
	nactiv  onwards of  rlamda  contain the (unadjusted) multipliers
	for the simple bounds
*/
            double blam, rlam, anormj;
            int nlam;
            int i, j, k, l, jgfxd, ka, kb, was_is = -111111111, nfixed;
            int idiag;

            //       rt -= Nrowrt + 1;
            var rt_offset = NROWRT + 1;
            // --rlamda;
            // --qtg;
            // --anorm;
            //   a -= NROWA + 1;
            var a_offset = NROWA + 1;
            // --kactiv;
            //         --istate;

            /*
            first, compute the lagrange multipliers for the general
            constraints in the working set, by solving
            t(transpose)*rlamda = y(t)*grad
            */
            nfixed = n - nfree;
            nlam = nfixed + nactiv;
            if (nactiv > 0)
            {
                //   BlasLike.dcopyvec(nactiv, &qtg[ncolz + 1], pRLAM);
                BlasLike.dcopyvec(nactiv, QTG, RLAM, ncolz);
                idiag = 1;
                drtmxsolve(-2, nactiv, RT, NROWRT,
                    RLAM, ref idiag, (ncolz + 1) * NROWRT + 1 - rt_offset);
            }
            /*
            now set elements nactiv, nactiv+1,... of rlamda equal to the
            multipliers for the bound constraints in the working set
            */
            for (l = 1; l <= nfixed; ++l)
            {
                kb = nactiv + l;
                j = KACTV[kb - 1];
                jgfxd = nfree + l;
                blam = QTG[jgfxd - 1];
                for (ka = 1; ka <= nactiv; ++ka)
                {
                    i = KACTV[ka - 1];
                    blam -= A[i + j * NROWA - a_offset] * RLAM[ka - 1];
                }
                RLAM[kb - 1] = blam;
            }

            /*find  allmax and smllst*/
            smllst = BlasLike.lm_max;
            jsmlst = 0;
            ksmlst = 0;
            for (k = 1; k <= nlam; ++k)
            {
                j = KACTV[k - 1];
                if (k > nactiv) anormj = 1;
                else
                {
                    anormj = ANORM[j - 1];
                    j += n;
                }
                was_is = ISTATE[j - 1];
                /*
                change the sign of the estimate if the constraint is in the
                working set (or violated) at its upper bound
                if (was_is == 3) rlam = Math.Abs(rlam);
                */
                if (was_is != 3)
                {
                    rlam = RLAM[k - 1] * anormj;
                    /*not a fixed variable or an equality constraint*/
                    if (was_is == 2) rlam = -rlam;
                    else if (was_is == 4) rlam = -Math.Abs(rlam);
                    /*find the smallest multiplier for the inequalities*/
                    if (rlam < smllst)
                    {
                        smllst = rlam;
                        jsmlst = j;
                        ksmlst = k;
                    }
                }
            }

            /*if required, print the multipliers*/
            if (msg < 20) return;
            if (nactiv > 0) w_lam(lprob, "CONSTRAINTS...", nactiv, KACTV, RLAM);
            l = nactiv + 1;
            if (l <= nlam) w_lam(lprob, "BOUND CONSTRAINTS...", nlam - nactiv, KACTV, RLAM);
            if (msg >= 80)
                lm_wmsg("\n//dgetlam//  JSMLST     SMLLST     KSMLST\n//dgetlam//%8ld%11.2lg%11ld",
                    jsmlst, smllst, ksmlst);
        }
        void w_lam(string msg1, string msg2, int n, int[] a, double[] r)
        {
            int i, j;

            lm_wmsg("\nMULTIPLIERS FOR THE ", msg1, msg2);
            for (i = j = 0; i < n; i++)
            {
                if (j == 4)
                {
                    Console.Write("\n");
                    j = 1;
                }
                else j++;
                //lm_printf((char*)"%5ld%11.2le",CL(a[i]),r[i]);
                Console.Write($"{a[i]} {r[i]}");
            }
            Console.Write("\n");
        }
        void lm_wmsg<T>(string mess, T n1)
        {
            ColourConsole.WriteInfo($"{mess} {n1}");
        }
        void lm_wmsg(string mess, int n1, double n2, int n3)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3}");
        }
        void lm_wmsg(string mess, double n1, double n2, double n3)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3}");
        }
        void lm_wmsg<T>(string mess, T n1, T n2, double n3, double n4)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4}");
        }
        void lm_wmsg<T>(string mess, int n1, int n2, int n3, int n4, int n5, T n6)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6}");
        }
        void lm_wmsg(string mess, string n1, string n2, string n3, int n4, double n5, double n6)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6}");
        }
        void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, bool n6, int n7, int n8, int n9, int n10)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9} {n10}");
        }
        void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, bool n6, double n7, int n8, double n9)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9}");
        }
        void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, int n6, double n7, int n8, double n9)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9}");
        }
        void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, double n6, double n7, double n8, double n9, double n10, double n11, double n12, double n13, double n14)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9} {n10} {n11} {n12} {n13} {n14}");
        }
        void lm_wmsg(string mess, string n1, string n2, int n3)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2} {n3}");
        }
        void lm_wmsg<T>(string mess, string n1, T n2)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2}");
        }
        void lm_wmsg(string mess, double n1, double n2)
        {
            ColourConsole.WriteInfo($"{mess} {n1} {n2}");
        }
        void dlpbgst(int n, int nactiv, int nfree, ref int jbigst, ref int kbigst,
         double dinky, double feamin, ref double trulam)
        {
            double rlam, biggst;
            int nlam, j, k, was_is, nfixed;

            // --rlamda;
            // --featol;
            //      --kactiv;
            //   --istate;

            jbigst = 0;
            nfixed = n - nfree;
            nlam = nfixed + nactiv;
            if (nlam == 0) return;
            biggst = 1 + dinky;
            for (k = 1; k <= nlam; ++k)
            {
                j = KACTV[k - 1];
                if (k <= nactiv) j += n;

                was_is = ISTATE[j - 1];
                if (was_is >= 1)
                {
                    rlam = RLAM[k - 1];
                    if (was_is == 2) rlam = -rlam;
                    if (was_is == 3) rlam = Math.Abs(rlam);
                    rlam *= FEATOL[j - 1] / feamin;
                    if (biggst < rlam)
                    {
                        biggst = rlam;
                        trulam = RLAM[k - 1];
                        jbigst = j;
                        kbigst = k;
                    }
                }
            }
            if (msg >= 80) lm_wmsg("\n//LPBGST// JBIGST         BIGGST\n//LPBGST//%7ld%15.4lg", jbigst, biggst);
        }
        void ddelcon(bool modfyg, bool orthog, int jdel, int kdel, int nactiv, int ncolz, int nfree, int n, int Nrowrt)
        {/*
	ddelcon updates the factorization of the matrix of
	constraints in the working set,  A(free)*(Z Y) = (0 T)
	if there are no general constraints in the working set and the
	matrix  Q = (Z Y)  is the identity, Q will not be touched
*/
            int i, j, k, l, ldiag;
            int nfree1, nactp1, nactv1, ka;
            double store;
            double cs = 1e3, sn = 1e2;
            int ibegin, ifreed;
            int nfreei, nactpi, istore;


            //         zy -= Nq + 1;
            int zy_offset = nq + 1;
            //  rt -= Nrowrt + 1;
            var rt_offset = Nrowrt + 1;
            //    --qtg;
            //   a -= NROWA + 1;
            var a_offset = NROWA + 1;
            //  --kfree;
            //  --kactiv;

            if (jdel <= n)
            {
                /*A SIMPLE BOUND IS BEING DELETED FROM THE WORKING SET. */
                ifreed = kdel - nactiv;
                if (msg >= 80)
                    lm_wmsg("BOUND DELETED", nactiv, ncolz, nfree, ifreed,
                        jdel, UNITQ);
                nactv1 = nactiv;
                nfree1 = nfree + 1;
                ibegin = 1;
                KFREE[nfree1 - 1] = jdel;

                /*
                ADD THE GRADIENT CORRESPONDING TO THE NEWLY-FREED VARIABLE TO THE 
                END OF  Q(FREE)(T)G(FREE).  THIS IS DONE BY INTERCHANGING THE
                APPROPRIATE ELEMENTS OF  QTG  AND  KACTIV
                */
                if (modfyg && ifreed != 1)
                {
                    nfreei = nfree + ifreed;
                    nactp1 = nactiv + 1;
                    nactpi = nactiv + ifreed;
                    store = QTG[nfree1 - 1];
                    QTG[nfree1 - 1] = QTG[nfreei - 1];
                    QTG[nfreei - 1] = store;
                    istore = KACTV[nactp1 - 1];
                    KACTV[nactp1 - 1] = KACTV[nactpi - 1];
                    KACTV[nactpi - 1] = istore;
                }

                /*COPY THE INCOMING COLUMN OF A INTO THE END OF  T. */
                if (!UNITQ)
                {
                    for (ka = 1; ka <= nactiv; ++ka)
                    {
                        i = KACTV[ka - 1];
                        RT[ka + nfree1 * Nrowrt - rt_offset] = A[i + jdel * NROWA - a_offset];
                    }
                    /*EXPAND  Q  BY ADDING A UNIT ROW AND COLUMN. */

                    {
                        BlasLike.dzero(nfree, ZY, nq, nfree1 + nq - zy_offset);
                        BlasLike.dzerovec(nfree, ZY, nfree1 * nq + 1 - zy_offset);
                    }
                    ZY[nfree1 + nfree1 * nq - zy_offset] = 1;
                }
            }
            else
            {
                /*A GENERAL CONSTRAINT IS BEING DELETED FROM THE WORKING SET*/
                if (msg >= 80)
                    lm_wmsg("CONSTRAINT DELETED", nactiv, ncolz, nfree, kdel,
                        jdel, UNITQ);
                nactv1 = nactiv - 1;
                nfree1 = nfree;
                ibegin = kdel;
                /*DELETE A ROW OF  T  AND MOVE THE ONES BELOW IT UP. */
                for (i = kdel; i <= nactv1; ++i)
                {
                    j = i + 1;
                    KACTV[i - 1] = KACTV[j - 1];
                    ldiag = nfree - i;
                    BlasLike.dcopy(j, RT, Nrowrt, RT, Nrowrt, j + ldiag * Nrowrt - rt_offset, i + ldiag * Nrowrt - rt_offset);
                }
            }

            /*
            ELIMINATE THE SUPER-DIAGONAL ELEMENTS OF  T
            USING A BACKWARD SWEEP OF 2*2 TRANFORMATIONS
            */
            if (ibegin <= (int)nactv1)
            {
                k = nfree1 - ibegin;
                l = nactv1 - ibegin;
                for (i = ibegin; i <= nactv1; ++i)
                {
                    delmgen(orthog, ref RT[i + (k + 1) * Nrowrt - rt_offset], ref RT[+i + k * Nrowrt - rt_offset], ref cs, ref sn);
                    if (l > 0) delm(orthog, l, RT, 1, RT, 1, cs, sn, i + 1 + (k + 1) * Nrowrt - rt_offset, i + 1 + k * Nrowrt - rt_offset);
                    if (nactv1 > 0) delm(orthog, nfree1, ZY, 1, ZY, 1, cs, sn, (k + 1) * nq + 1 - zy_offset, k * nq + 1 - zy_offset);
                    if (modfyg) delm(orthog, 1, QTG, 1, QTG, 1, cs, sn, k, k - 1);
                    --k;
                    --l;
                }
            }

            /*COMPRESS THE ELEMENTS OF  KACTIV  CORRESPONDING TO FIXED VARIABLES*/
            i = n - nfree1;
            if (i != 0)
            {
                j = nactv1 + 1;
                for (k = 1; k <= i; ++k)
                {
                    KACTV[j - 1] = KACTV[j];
                    ++j;
                }
            }
            /*ESTIMATE THE CONDITION NUMBER OF  T. */
            if (nactv1 > 0)
            {
                BlasLike.dxminmax(nactv1, RT, Nrowrt - 1, ref dtmax, ref dtmin, nactv1 + (ncolz + 2) * Nrowrt - rt_offset);
            }
        }
        void dfindp(bool nullr, bool unitpg, int n, int nclin, int Nrowrt, int ncolr, int ncolz, ref int nfree, bool negligible, ref double gtp, ref double pnorm, ref double rdlast)
        {
            /*
                findp computes the following quantities for  lpcore,  qpcore  and
                    lccore ..
                (1)	the search direction  p  (and its 2-norm)
                (2)	the vector	v  such that  r(t)v = - z(t)g(free).  this vector
                    is required by	lccore	only
                (3)	the vector	ap,  where  a  is the matrix of linear
                    constraints. and, if  nullr  is false,
                (4)	the	 (ncolr)-th diagonal element of the cholesky factor of the
                    projected hessian
            */
            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset;
            int j;
            short idiag;

            // --work;
            zy_dim1 = nq;
            zy_offset = zy_dim1 + 1;
            //    zy -= zy_offset;
            //  --v;
            rt_dim1 = Nrowrt;
            rt_offset = rt_dim1 + 1;
            //    rt -= rt_offset;
            // --qtg;
            //  --p;
            //  --ap;
            a_dim1 = NROWA;
            a_offset = a_dim1 + 1;
            //   a -= a_offset;
            // --kfree;
            // --istate;
            BlasLike.dcopyvec(ncolr, QTG, PX);
            BlasLike.dscalvec(ncolr, -1, PX);
            if (nullr) goto L60;
            rdlast = RT[ncolr + ncolr * rt_dim1 - rt_offset];
            /*     *** */
            /*     CORRECTION INSERTED BY MHW, 22 OCT 1985. */
            /*     THIS ENSURES A NON-ZERO SEARCH DIRECTION. */
            /*     *** */
            if (ncolr < ncolz && negligible) PX[ncolr - 1] = rdlast;
            /* --------------------------------------------------------------------- 
            */
            /*     SOLVE THE SYSTEM   R(T)R (PZ) = - Z(T)G(FREE). */
            /* --------------------------------------------------------------------- 
            */
            if (unitpg)
            {
                goto L20;
            }
            /*     PERFORM THE FORWARD SUBSTITUTION  R(T)V = - Z(T)G(FREE). */
            idiag = dtmxsolve((short)-1, ncolr, RT, Nrowrt, PX, (short)1);
            goto L40;
        /*     THE PROJECTED GRADIENT IS A MULTIPLE OF THE UNIT VECTOR, THE */
        /*     FORWARD SUBSTITUTION MAY BE AVOIDED. */
        L20:
            if (negligible) PX[ncolr - 1] = -1;
            else PX[ncolr - 1] /= rdlast;
            /*     PERFORM THE BACKWARD SUBSTITUTION   R(PZ) = P. */
            L40:
            BlasLike.dcopyvec(ncolr, PX, WRK);
            idiag = dtmxsolve((short)1, ncolr, RT, Nrowrt, PX, (short)1);
        /* --------------------------------------------------------------------- 
        */
        /*     THE VECTOR  (PZ)  HAS BEEN COMPUTED. */
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE DIRECTIONAL DERIVATIVE  G(T)P = (GZ)(T)(PZ). */
        L60:
            gtp = BlasLike.ddotvec(ncolr, QTG, PX);
            /* --------------------------------------------------------------------- 
            */
            /*     COMPUTE  P = Z * PZ. */
            /* --------------------------------------------------------------------- 
            */
            /*     NACTIV  AND  KACTIV  ARE NOT USED IN  ZYPROD.  N  AND  KFREE */
            /*     SERVE AS ARGUMENTS FOR  NACTIV  AND  KACTIV. */
            dzyprod(1, n, n, ncolr, nfree, nq, KFREE, KFREE, PX,
                  WRK);
            pnorm = dnrm2vec(nfree, WRK);
            if (msg >= 80) lm_mdvwri("\n//FINDP//   P ... ", n, PX);
            /* --------------------------------------------------------------------- 
            */
            /*     COMPUTE  AP. */
            /* --------------------------------------------------------------------- 
            */
            if (nclin > 0)
            {
                BlasLike.dzerovec(nclin, AP);
                for (j = 1; j <= (int)n; ++j)
                {
                    if (ISTATE[j - 1] <= 0)
                        //     BlasLike.daxpyvec(*nclin, PX[j - 1], pA + j * a_dim1 + 1 - a_offset, pAP);
                        BlasLike.daxpyvec(nclin, PX[j - 1], A, AP, j * a_dim1 + 1 - a_offset);
                    /* L100: */
                }
                if (msg >= 80) lm_mdvwri("\n//FINDP//  AP ... ", nclin, AP);
            }
            return;
        }
        short dbndalf(bool firstv, ref int hitlow, ref int jadd, int n, int nctotl, int numinf, ref double alfa, ref double palfa, ref double atphit, ref double bigalf, double pnorm)
        {
            /*
                dbndalf finds a step  alfa  such that the point  x + alfa*p
                reaches one of the linear constraints (including bounds).
                Two possible steps are defined as follows...
                alfa1	is the maximum step that can be taken without violating
                    one of the linear constraints that is currently satisfied.
                alfa2	reaches a linear constraint that is currently violated.
                    usually this will be the furthest such constraint along  p,
                    but if  firstv = 1 it will be the first one along  p.
                    this is used only by dlpcore when the problem has been
                    determined to be infeasible, and we are now minimizing the
                    sum of infeasibilities. (alfa2 is not defined if numinf==0)

                alfa will usually be the minimum of  alfa1  and  alfa2.

                alfa could be negative (since we allow inactive constraints
                to be violated by as much as  featol).  in such cases, a
                third possible step is computed, to find the nearest satisfied
                constraint (perturbed by  featol) along the direction  - p.
                alfa  will be reset to this step if it is shorter.  this is the
                only case for which the final step  alfa  does not move  x
                exactly onto a constraint (the one denoted by  jadd).
                constraints in the working set are ignored  (istate(j) ge 1).

                jadd	denotes which linear constraint is reached.
                hitlow	indicates whether it is the lower or upper bound that
                    has restricted  alfa.

                values of istate(j)....
                    -2         -1         0           1           2          3
                    a*x < bl   a*x > bu   a*x free   a*x == bl   a*x == bu   bl == bu

                the values -2 and -1 do not occur once lpcore finds a feasible point.
            */
            //#define FMT2 \
            //            "\n//BNDALF//  NEGATIVE STEP.\n//BNDALF//           ALFA          PALFA\n//BNDALF//%15.5lg%15.5lg"

            short inform;
            double d__1;
            int jadd1 = 9;
            double alfa1, alfa2;
            int jadd2 = 8;
            bool hlow1, hlow2, lastv;
            int i, j;
            double palfa1 = 1e34, palfa2 = 1e5, apmax1, apmax2;
            int jsave1, jsave2;
            double epspt9;
            int js;
            double absatp;
            double rownrm, atp, res, atx, atp1, atp2;
            // --x;
            //          --p;
            // --featol;
            //          --bu;
            //          --bl;
            // --ax;
            // --ap;
            // --anorm;
            //      --istate;

            epspt9 = parm[3];
            inform = 0;

            /*
            FIRST PASS -- FIND STEPS TO PERTURBED CONSTRAINTS, SO THAT
            PALFA1  WILL BE SLIGHTLY LARGER THAN THE TRUE STEP, AND
            PALFA2  WILL BE SLIGHTLY SMALLER THAN IT SHOULD BE.  IN DEGENERATE 
            CASES, THIS STRATEGY GIVES US SOME FREEDOM IN THE SECOND PASS.
            THE GENERAL IDEA FOLLOWS THAT DESCRIBED BY P.M.J. HARRIS, P.21 OF 
            MATHEMATICAL PROGRAMMING 5, 1 (1973), 1--28.
            */
            dbdpert(firstv, false, bigalf, pnorm, ref jadd1, ref jadd2, ref palfa1, ref
                palfa2, n, nctotl);
            jsave1 = jadd1;
            jsave2 = jadd2;
            var bigbnd = parm[0];
            /*
            SECOND PASS -- RECOMPUTE STEP-LENGTHS WITHOUT PERTURBATION
            AMONGST CONSTRAINTS THAT ARE CLOSE TO THE PERTURBED STEPS
            CHOOSE THE ONE (OF EACH TYPE) THAT MAKES THE LARGEST ANGLE
            WITH THE SEARCH DIRECTION
            */
            if (msg == 99) ColourConsole.WriteInfo(
        "BNDALF ENTERED\n    J  JS         FEATOL         AX             AP     JADD1        ALFA1     JADD2        ALFA2");
            alfa1 = bigalf;
            alfa2 = 0.0;
            if (firstv) alfa2 = bigalf;
            apmax1 = 0.0;
            apmax2 = 0.0;
            atp1 = 0.0;
            atp2 = 0.0;
            hlow1 = false;
            hlow2 = false;
            lastv = !firstv;
            for (j = 1; j <= nctotl; ++j)
            {
                js = ISTATE[j - 1];
                if (js > 0) continue;
                if (j > n)
                {
                    /*GENERAL LINEAR CONSTRAINT. */
                    i = j - n;
                    atx = LWRK[i - 1];
                    atp = AP[i - 1];
                    /*			lm_wmsg((char*)"atx %e  atp %e  istate[%d]=%d",atx,atp,j,js); */
                    rownrm = ANORM[i - 1] + 1.0;
                }
                else
                {
                    /*BOUND CONSTRAINT. */
                    atx = W[j - 1];
                    atp = PX[j - 1];
                    rownrm = 1.0;
                }
                if (Math.Abs(atp) <= epspt9 * rownrm * pnorm) res = -1.0;
                else if (atp > BlasLike.lm_eps)
                {
                    /*ATX IS INCREASING. */
                    /*TEST FOR SMALLER ALFA1 IF UPPER BOUND IS SATISFIED. */
                    if (js != -1)
                    {
                        if (U[j - 1] < bigbnd)
                        {
                            res = U[j - 1] - atx;
                            if (((int)j == jsave1 || palfa1 * atp >= res)
                                && apmax1 * rownrm * pnorm < atp)
                            {
                                apmax1 = atp / (rownrm * pnorm);
                                alfa1 = res / atp;
                                jadd1 = j;
                                atp1 = atp;
                                hlow1 = false;
                            }
                        }
                        /*TEST FOR BIGGER ALFA2 IF LOWER BOUND IS VIOLATED. */
                        if (js == -2)
                        {
                            res = L[j - 1] - atx;
                            if ((firstv || (int)j == jsave2 || palfa2 * atp <= res)
                                && (lastv || (int)j == jsave2 || palfa2 * atp >= res)
                                && (apmax2 * rownrm * pnorm < atp))
                            {
                                apmax2 = atp / (rownrm * pnorm);
                                if (atp >= 1.0) alfa2 = res / atp;
                                else if (res < bigalf * atp) alfa2 = res / atp;
                                else alfa2 = bigalf;
                                jadd2 = j;
                                atp2 = atp;
                                hlow2 = true;
                            }
                        }
                    }
                }
                else if (js != -2)
                {
                    /*ATX IS DECREASING. */
                    /*TEST FOR SMALLER ALFA1 IF LOWER BOUND IS SATISFIED*/
                    absatp = -atp;
                    if (L[j - 1] > -bigbnd)
                    {
                        res = atx - L[j - 1];
                        if (((int)j == jsave1 || palfa1 * absatp >= res)
                            && (apmax1 * rownrm * pnorm < absatp))
                        {
                            apmax1 = absatp / (rownrm * pnorm);
                            alfa1 = res / absatp;
                            jadd1 = j;
                            atp1 = atp;
                            hlow1 = true;
                        }
                    }
                    /*TEST FOR BIGGER ALFA2 IF UPPER BOUND IS VIOLATED*/
                    if (js == -1)
                    {
                        res = atx - U[j - 1];
                        if ((firstv || (int)j == jsave2 || palfa2 * absatp <= res)
                            && (lastv || (int)j == jsave2 || palfa2 * absatp >= res)
                            && (apmax2 * rownrm * pnorm < absatp))
                        {
                            apmax2 = absatp / (rownrm * pnorm);
                            if (absatp >= 1.0) alfa2 = res / absatp;
                            else if (res < bigalf * absatp) alfa2 = res / absatp;
                            else alfa2 = bigalf;
                            jadd2 = j;
                            atp2 = atp;
                            hlow2 = false;
                        }
                    }
                }
                if (msg == 99)
                    lm_wmsg("%5ld%4ld%15.5lg%15.5lg%15.5lg%6ld%17.7lg%6ld%17.7lg",
                        j, js, FEATOL[j - 1], atx, atp, jadd1, alfa1,
                        jadd2, alfa2);
            }

            /*IF FEASIBLE, ONLY ALFA1 WILL HAVE BEEN SET. */
            alfa = alfa1;
            palfa = palfa1;
            jadd = jadd1;
            atphit = atp1;
            hitlow = hlow1 ? 1 : 0;
            if (numinf != 0 && jadd2 != 0 && (alfa2 < alfa1 || (alfa2 <= palfa1 && apmax2 >= apmax1)))
            {
                /*
                INFEASIBLE -- SEE IF WE STEP TO THE FURTHEST VIOLATED CONSTRAINT. 
                BE PREPARED TO STEP IN THE RANGE  (ALFA1, PALFA1)  IF THE VIOLATED 
                CONSTRAINT HAS A LARGER VALUE OF  AP
                */
                alfa = alfa2;
                jadd = jadd2;
                atphit = atp2;
                hitlow = hlow2 ? 1 : 0;
            }
            else if (alfa < -BlasLike.lm_eps)
            {
                /*
                NEGATIVE STEP
                JADD  WILL RETAIN ITS CURRENT VALUE, BUT WE MAY SHORTEN  ALFA
                TO BE  - PALFA1,  THE STEP TO THE NEAREST PERTURBED SATISFIED
                CONSTRAINT ALONG THE DIRECTION -P
                */
                dbdpert(firstv, true, bigalf, pnorm, ref jadd1, ref jadd2, ref palfa1, ref
                    palfa2, n, nctotl);
                if (msg >= 80) lm_wmsg("NEGATIVE STEP", alfa, palfa1);
                d__1 = Math.Abs(alfa);
                alfa = -Math.Min(d__1, palfa1);
            }

            /*
            TEST FOR UNDEFINED OR INFINITE STEP.  THIS SHOULD MEAN THAT THE
            SOLUTION IS UNBOUNDED
            */
            if (jadd == 0)
            {
                alfa = bigalf;
                palfa = bigalf;
                inform = 2;
            }
            if (alfa >= bigalf) inform = 3;
            if (msg >= 80 && inform > 0) lm_wmsg(
        "\n//BNDALF//  UNBOUNDED STEP.\n//BNDALF//  JADD          ALFA\n//BNDALF//  %4ld%15.5lg",
                jadd, alfa);
            return inform;
        }
        short daddcon(bool modfyg, bool modfyr, bool orthog, int ifix, int iadd, int jadd, int nactiv, int ncolr, int ncolz, int nfree, int n, int nrowart, int[] kfree, double condmx, double cslast, double snlast, double[] wrrk1, double[] wrk2, int kfr = 0, int wk1 = 0, int wk2 = 0)
        {

            /*
                daddcon updates the factorization of the matrix of
                constraints in the working set,  a(free) * (z y) = (0 t)
                if the int argument modfyr  is true, the cholesky
                factorization of the projected hessian, r(t)*r, is updated also
                there are three separate cases to consider (although each case
                shares code with another)..

                (1)	a free variable becomes fixed on one of its bounds when there
                    are already some general constraints in the working set
                (2)	a free variable becomes fixed on one of its bounds when there
                    are only bound constraints in the working set
                (3)	a general constraint (corresponding to row  iadd  of  a) is
                    added to the working set
                in cases (1) and (2), we assume that  kfree(ifix) = jadd
                in all cases,  jadd  is the index of the constraint being added
                if there are no general constraints in the working set,  the
                matrix  q = (z y)  is the identity and will not be touched
                if modfyr  is true and  ncolz is greater than one on entry,
                cslast and snlast contain the last of the sequence of givens
                rotations used to reduce the intermediate upper-hessenberg matrix
                to upper-triangular form.  these elements are needed by qpcore
                if  modfyg  is true on entry, the column transformations are
                applied to the vector  q(t)grad,  stored in  qtg
            */
            const double zero = 0.0;

            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1;

            double beta = 41, cond;
            int nelm, inct, nact1;
            double d;
            int i, j, k;
            int ldiag;
            int ifail = -27;
            double delta;
            double dtnew = -9999999999;
            int iswap = 10;
            int lrowr, nfree1 = 0, ncolz1;
            double cs = 5, condbd, sn = 6;
            double tdtmin, tdtmax;
            int itrans = 12, kp1;

            //       --wrk2;
            wk2--;
            //        --wrk1;
            wk1--;
            zy_dim1 = nq;
            zy_offset = zy_dim1 + 1;
            //    zy -= zy_offset;
            rt_dim1 = NROWRT;
            rt_offset = rt_dim1 + 1;
            //    rt -= rt_offset;
            // --qtg;
            a_dim1 = nrowart; //Either called as NROWA or NROWRT
            a_offset = a_dim1 + 1;
            // a -= a_offset;
            //      --kfree;
            kfr--;

            /*
            if the condition estimator of the updated factors is greater than
            condbd,  a warning message is printed
            */
            condbd = Math.Pow(BlasLike.lm_eps, -0.9);
            ncolz1 = ncolz - 1;
            if (jadd > (int)n)
            {
                goto L60;
            }
            /*
                 a simple bound has entered the working set.  iadd  is not used
            */
            if (msg >= 80)
                lm_wmsg(
            "\n//ADDCON//  SIMPLE BOUND ADDED.\n//ADDCON// NACTIV NCOLZ NFREE  IFIX  JADD UNITQ\n//ADDCON//%7ld%6ld%6ld%6ld%6ld%6ld",
                    nactiv, ncolz, nfree, ifix,
                    jadd, UNITQ);
            /*
                 set  wrk1 = appropriate row of  q
                 reorder the elements of  kfree (this requires reordering the
                 corresponding rows of  q)
            */
            nfree1 = nfree - 1;
            nact1 = nactiv;
            if (UNITQ)
            {
                goto L20;
            }
            /*
                 q  is stored explicitly.  interchange components  ifix  and  nfree
                 of  kfree  and swap the corresponding rows of  q
            */
            // BlasLike.dcopy(*nfree, &pZY[*ifix + zy_dim1 - zy_offset], nq, &wrk1[1], 1);
            BlasLike.dcopy(nfree, ZY, nq, WRK, 1, ifix + zy_dim1 - zy_offset);
            if (ifix == nfree) goto L180;
            kfree[ifix + kfr] = kfree[nfree + kfr];
            // BlasLike.dcopy(*nfree, &pZY[*nfree + zy_dim1 - zy_offset], nq, &pZY[*ifix + zy_dim1 - zy_offset], nq);
            BlasLike.dcopy(nfree, ZY, nq, ZY, nq, nfree + zy_dim1 - zy_offset, ifix + zy_dim1 - zy_offset);
            goto L180;
        /*
             q  is not stored, but  kfree  defines an ordering of the columns
             of the identity matrix that implicitly define  z
             reorder  kfree  so that variables  ifix+1,...,nfree  are moved one
             position to the left
        */
        L20:
            BlasLike.dzerovec(nfree, WRK);
            WRK[ifix - 1] = 1.0;
            if (ifix == nfree)
            {
                goto L180;
            }
            i__1 = nfree1;
            for (i = ifix; i <= i__1; ++i)
            {
                kfree[i + kfr] = kfree[i + 1 + kfr];
                /* L40: */
            }
            goto L180;
        /*
             a general constraint has entered the working set
             ifix is not used
        */
        L60:
            if (msg >= 80)
                lm_wmsg("\n//ADDCON//  GENERAL CONSTRAINT ADDED.\n//ADDCON// NACTIV NCOLZ NFREE  IADD  JADD UNITQ\n//ADDCON//%7ld%6ld%6ld%6ld%6ld%6ld",
                    nactiv, ncolz, nfree, iadd, jadd, UNITQ);
            nact1 = nactiv + 1;
            /*
                 transform the incoming row of  a  by  q(t)
            */
            //     BlasLike.dcopy(n, &a[*iadd + a_dim1], NROWA, &wrk1[1], 1);
            BlasLike.dcopy(n, A, NROWA, WRK, 1, iadd + a_dim1 - a_offset);
            dzyprod(8, n, nactiv, ncolz, nfree, nq, kfree, kfree,
                WRK, wrk2, 1 + kfr, 1 + kfr, 0, 1 + wk2);
            if (!UNITQ)
            {
                goto L100;
            }
            /*
                 this is the first general constraint to be added  --  set  q = i.
            */
            i__1 = nfree;
            for (j = 1; j <= i__1; ++j)
            {
                // BlasLike.dzerovec(*nfree, &pZY[j * zy_dim1 + 1 - zy_offset]);
                BlasLike.dzerovec(nfree, ZY, j * zy_dim1 + 1 - zy_offset);
                ZY[j + j * zy_dim1 - zy_offset] = 1.0;
            }
            //      *unitq = 0;
            UNITQ = false;
        /*
             check that the incoming row is not dependent upon those
             already in the working set
        */
        L100:
            dtnew = dnrm2vec(ncolz, WRK);
            if (nact1 > 1)
            {
                goto L140;
            }
            /*
                 this is the only general constraint in the working set
            */
            cond = BlasLike.dprotdiv(ref asize, ref dtnew, ref ifail);
            if (ifail != 0 && asize == 0)
            {
                cond = BlasLike.lm_max;
            }
            if (cond >= condmx)
            {
                goto L480;
            }
            if (cond >= condbd && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n", jadd);
            dtmax = dtnew;
            dtmin = dtnew;
            goto L180;
        /*
             there are already some general constraints in the working set
             update the estimate of the condition number
        */
        L140:
            tdtmax = Math.Max(dtnew, dtmax);
            tdtmin = Math.Min(dtnew, dtmin);
            cond = BlasLike.dprotdiv(ref tdtmax, ref tdtmin, ref ifail);
            if (ifail != 0 && tdtmax == 0)
            {
                cond = BlasLike.lm_max;
            }
            if (cond >= condmx)
            {
                goto L480;
            }
            if (cond >= condbd && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n", jadd);
            dtmax = tdtmax;
            dtmin = tdtmin;
        /*
             use one or more column transformations to reduce the first  ncolz1
             elements of  wrk1  to zero.  this affects  zy,  except if (unitq).
             the transformations may also be applied to  qtg  and  r
        */
        L180:
            if (ncolz1 == 0)
            {
                goto L360;
            }
            if (modfyr || UNITQ)
            {
                goto L320;
            }
            /*
            there is no  r.  use a single elimination or householder matrix
            */
            if (orthog)
            {
                goto L240;
            }
            /*
                 elimination we use  elm( ..., zero, zero )   to perform an interchange
            */
            BlasLike.detagen(ncolz1, ref WRK[ncolz - 1], WRK, 1, ref iswap, ref itrans);

            if (iswap > 0)
                delm(orthog, nfree, ZY, 1, ZY, 1, zero, zero, ncolz * zy_dim1 + 1 - zy_offset, iswap * zy_dim1 + 1 - zy_offset);
            if (itrans == 0)
            {
                goto L220;
            }
            i__1 = ncolz1;
            for (j = 1; j <= i__1; ++j)
            {
                d = WRK[j - 1];
                if (d == 0)
                {
                    goto L200;
                }
                //   BlasLike.daxpy(*nfree), d, &zy[*ncolz * zy_dim1 + 1], 1, &zy[j * zy_dim1 + 1], 1);
                BlasLike.daxpy(nfree, d, ZY, 1, ZY, 1, ncolz * zy_dim1 + 1 - zy_offset, j * zy_dim1 + 1 - zy_offset);
            L200:
                ;
            }
        L220:
            if (!modfyg)
            {
                goto L360;
            }
            if (iswap > 0)
                // delm(orthog, 1, &qtg[*ncolz], 1, &qtg[iswap], 1, zero, zero);
                delm(orthog, 1, QTG, 1, QTG, 1, zero, zero, ncolz - 1, iswap - 1);
            if (itrans > 0)
            {
                // BlasLike.daxpy(ncolz1, qtg[*ncolz], &wrk1[1], 1, &qtg[1], 1);
                BlasLike.daxpy(ncolz1, QTG[ncolz - 1], WRK, 1, QTG, 1);
            }
            goto L360;

        /*
        orthogonal transformation
        we use a householder reflection,   i  -  1/beta  v v(t)
        there are two ways of applying the reflection.  the update to  z
        is done via   w  =  z * v,   z  =  z  -  1/beta  w v(t),
        where  v = wrk1 (from householder), and  w = wrk2 (workspace)
        the update to  qtg  is the more usual  d =  - qtg(t)*v / beta,
        qtg  =  qtg  +  d * v
        note that  delta  has to be stored after the reflection is used
        */

        L240:
            delta = WRK[ncolz - 1];
            var c__1 = 1;
            dhhrflctgen(ncolz1, ref delta, WRK, c__1, BlasLike.lm_eps, ref beta);
            if (beta != 0)
            {
                WRK[ncolz - 1] = beta;
            }
            if (beta <= 0)
            {
                goto L360;
            }
            BlasLike.dzerovec(nfree, wrk2, 1 + wk2);
            i__1 = ncolz;
            for (j = 1; j <= i__1; ++j)
            {
                d = WRK[j - 1];
                if (d == 0)
                {
                    goto L260;
                }
                //     BlasLike.daxpy(*nfree, d, &zy[j * zy_dim1 + 1], 1, &wrk2[1], 1);
                //    BlasLike.daxpy(nfree, d, pZY + j * zy_dim1 + 1 - zy_offset, 1, &wrk2[1], 1);
                BlasLike.daxpy(nfree, d, ZY, 1, wrk2, 1, j * zy_dim1 + 1 - zy_offset, 1 + wk2);
            L260:
                ;
            }
            i__1 = ncolz;
            for (j = 1; j <= i__1; ++j)
            {
                d = WRK[j - 1];
                if (d == 0)
                {
                    goto L280;
                }
                d = -d / beta;
                //  BlasLike.daxpy(nfree, d, wrk2[1], 1, &pZY[j * zy_dim1 + 1 - zy_offset], 1);
                BlasLike.daxpy(nfree, d, wrk2, 1, ZY, 1, 1 + wk2, j * zy_dim1 + 1 - zy_offset);
            L280:
                ;
            }
            if (!modfyg)
            {
                goto L300;
            }
            d = BlasLike.ddotvec(ncolz, WRK, QTG);
            d = -d / beta;
            BlasLike.daxpyvec(ncolz, d, WRK, QTG);
        L300:
            WRK[ncolz - 1] = delta;
            goto L360;

        /*r  has to be modified.  use a sequence of 2*2 transformations*/
        L320:
            lrowr = ncolr;
            i__1 = ncolz1;
            for (k = 1; k <= i__1; ++k)
            {
                /*
                compute the transformation that reduces wrk1(k) to zero,
                then apply it to the relevant columns of  z  and  grad(t)q.
                */
                kp1 = k + 1;
                //     delmgen(orthog, WRK[kp1-1], WRK[k-1], cs, sn);
                delmgen(orthog, ref WRK[kp1 - 1], ref WRK[k - 1], ref cs, ref sn);
                if (!UNITQ)
                {
                    delm(orthog, nfree, ZY, 1, ZY, 1, cs, sn, kp1 * zy_dim1 + 1 - zy_offset, k * zy_dim1 + 1 - zy_offset);
                }
                if (modfyg)
                    //  delm(orthog, 1, &qtg[kp1], 1, &qtg[k], 1, cs, sn);
                    delm(orthog, 1, QTG, 1, QTG, 1, cs, sn, kp1 - 1, k - 1);
                /*
                apply the same transformation to the cols of  r  if relevant
                this generates a subdiagonal element in  r  which must be
                eliminated by a row rotation.  the last such row rotation
                is needed by  qpcore
                */
                if (!(modfyr && k < ncolr))
                {
                    goto L340;
                }
                RT[kp1 + k * rt_dim1 - rt_offset] = 0;
                {
                    delm(orthog, kp1, RT, 1, RT, 1, cs, sn, kp1 * rt_dim1 + 1 - rt_offset, k * rt_dim1 + 1 - rt_offset);
                    BlasLike.drotg(ref RT[k + k * rt_dim1 - rt_offset], ref RT[kp1 + k * rt_dim1 - rt_offset], ref cslast, ref snlast);
                }
                RT[kp1 + k * rt_dim1 - rt_offset] = 0;
                --lrowr;
                dsymplanerotate(lrowr, RT,
                        NROWRT, RT,
                        NROWRT, cslast, snlast, k + kp1 * rt_dim1 - rt_offset, kp1 + kp1 * rt_dim1 - rt_offset);
            L340:
                ;
            }
        /* if adding a general constraint, insert the new row of  t  and exit */
        L360:
            if (jadd > (int)n)
            {
                //     BlasLike.dcopy(nact1, &wrk1[ncolz], 1, &pRT[nact1 + ncolz * rt_dim1 - rt_offset], NROWRT);
                BlasLike.dcopy(nact1, WRK, 1, RT, NROWRT, ncolz + wk1, nact1 + ncolz * rt_dim1 - rt_offset);
                return 0;
            }
            /*
                 we are adding a bound.  continue reducing the elements of  wrk1
                 to zero.  this affects  y,  t  and  qtg
                 first, set the super-diagonal elements of t to zero
            */
            if (nactiv == 0)
            {
                goto L440;
            }
            //BlasLike.dzero(*nactiv, &rt[*nactiv + *ncolz * rt_dim1], nrowrt - 1);
            BlasLike.dzero(nactiv, RT, NROWRT - 1, nactiv + ncolz * rt_dim1 - rt_offset);
            nelm = 1;
            ldiag = nactiv;
            i__1 = nfree1;
            for (k = ncolz; k <= i__1; ++k)
            {
                // delmgen(orthog, &wrk1[k + 1], &wrk1[k], &cs, &sn);
                delmgen(orthog, ref WRK[k + 1 - 1], ref WRK[k - 1], ref cs, ref sn);
                delm(orthog, nfree, ZY, 1, ZY, 1, cs, sn, (k + 1) * zy_dim1 + 1 - zy_offset, k * zy_dim1 + 1 - zy_offset);
                delm(orthog, nelm, RT, 1,
                    RT, 1, cs, sn, ldiag + (k + 1) * rt_dim1 - rt_offset, ldiag + k * rt_dim1 - rt_offset);
                if (modfyg)
                    //   delm(orthog, 1, &qtg[k + 1], 1, &qtg[k], 1, cs, sn);
                    delm(orthog, 1, QTG, 1, QTG, 1, cs, sn, k, k - 1);
                ++nelm;
                --ldiag;
            }
            /*
            the diagonals of  t  have been altered.  recompute the largest and
            smallest values
            */
            inct = NROWRT - 1;
            BlasLike.dxminmax(nactiv, RT, inct, ref dtmax, ref dtmin, nactiv + (ncolz1 + 1) * rt_dim1 - rt_offset);
            if (dtmin / dtmax * condmx < 1.0) goto L480;
            if (dtmin / dtmax * condbd < 1.0 && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n",
                jadd);
            /*
            the last row of  zy  has been transformed to a multiple of the
            unit vector  e(nfree).  if orthogonal transformations have been
            used throughout, the last column of  zy  is the same.   we can
            therefore resurrect the gradient element of the newly fixed
            variable
            */
            L440:
            if (orthog && modfyg)
            {
                QTG[nfree - 1] /= WRK[nfree - 1];
            }
            /*the factorization has been successfully updated*/
            return 0;

        /*THE PROPOSED WORKING SET APPEARS TO BE LINEARLY DEPENDENT*/
        L480:
            if (msg >= 80)
            {
                ColourConsole.WriteInfo("\n//ADDCON//  DEPENDENT CONSTRAINT REJECTED");
                if (jadd <= (int)n)
                    lm_wmsg(
                    "\n//ADDCON//     ASIZE     DTMAX     DTMIN\n//ADDCON//%10.2le%10.2le%10.2le",
                        asize, dtmax, dtmin);
                else lm_wmsg(
             "\n//ADDCON//     ASIZE     DTMAX     DTMIN     DTNEW\n//ADDCON//%10.2le%10.2le%10.2le",
                 asize, dtmax, dtmin, dtnew);
            }
            return 1;
        }
        void dprtsol(int nfree, int n, int nclin, int ncnln, int nctotl, int nactiv)
        {

            /*
                dprtsol	expands the lagrange multipliers into  clamda
                if  msg >= 10  or  msg == 1,  prtsol  then prints  x, a*x,
                c(x), their bounds,  the multipliers, and the residuals
                (distance to the nearest bound)

                dprtsol	is called by  lpcore, qpcore, lccore and npcore	 just
                before they exit
            */
            var bigbnd = parm[0];
            string id = "VLN";
            string lstate = "--++FRLLULEQTB";


            int nlam, nplin, nfixed, j, k, ip, was_is;
            double v, b1, b2, res, res2, wlam;
            string ls;
            char[] id3 = new char[1];

            // --x;
            //       --rlamda;
            // --clamda;
            // --c;
            //        --bu;
            //        --bl;
            //       a -= NROWA + 1;
            var a_offset = NROWA + 1;
            // --kactiv;
            //--istate;

            nplin = n + nclin;
            /*EXPAND BOUND, LINEAR AND NONLINEAR MULTIPLIERS INTO CLAMDA*/
            BlasLike.dzerovec(nctotl, LAMBDA);
            nfixed = n - nfree;
            nlam = nactiv + nfixed;
            for (k = 1; k <= nlam; ++k)
            {
                j = KACTV[k - 1];
                if (k <= nactiv) j += n;
                LAMBDA[j - 1] = RLAM[k - 1];
            }
            if (msg < 10 && msg != 1) return;
            ColourConsole.WriteInfo("\n\nVARBL STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
            id3[0] = id[0];
            for (j = 1; j <= nctotl; ++j)
            {
                b1 = L[j - 1];
                b2 = U[j - 1];
                wlam = LAMBDA[j - 1];
                was_is = ISTATE[j - 1];
                ls = lstate.Substring(((was_is + 2) << 1)); //IS THIS RIGHT
                                                            //		ls = lstate + ((was_is + 2) << 1);
                if (j <= n)
                {
                    /* SECTION 1 -- THE VARIABLES  X. */
                    /* ------------------------------ */
                    k = j;
                    v = W[j - 1];
                }
                else if (j <= nplin)
                {
                    /* SECTION 2 -- THE LINEAR CONSTRAINTS  A*X. */
                    /* ----------------------------------------- */
                    if (j == n + 1)
                    {
                        ColourConsole.WriteInfo("\n\nLNCON STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
                        id3[0] = id[1];
                    }
                    k = j - n;
                    v = BlasLike.ddot(n, A, NROWA, W, 1, k + NROWA - a_offset);
                }
                else
                {
                    /*        SECTION 3 -- THE NONLINEAR CONSTRAINTS  C(X). */
                    /*        --------------------------------------------- */
                    if (ncnln <= 0) continue;
                    if (j == nplin + 1)
                    {
                        ColourConsole.WriteInfo("\n\nNLCON STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
                        id3[0] = id[2];
                    }
                    k = j - nplin;
                    v = c[k];
                }
                /* PRINT A LINE FOR THE J-TH VARIABLE OR CONSTRAINT. */
                /* ------------------------------------------------- */
                res = v - b1;
                res2 = b2 - v;
                if (Math.Abs(res) > Math.Abs(res2)) res = res2;
                ip = 1;
                if (b1 <= -bigbnd) ip = 2;
                if (b2 >= bigbnd) ip += 2;
                ColourConsole.WriteInfo($"{id3[0]},{k},{ls[0]},{ls[1]},{v}");
                if (ip != 0) ColourConsole.WriteInfo($"{b1}");
                else ColourConsole.WriteInfo("     NONE      ");
                if (ip < 3) ColourConsole.WriteInfo($"{b2}");
                else ColourConsole.WriteInfo("     NONE      ");
                ColourConsole.WriteInfo($"{wlam},{res}");
            }
        }
        double dnrm2(int n, double[] x, int incx, int xstart = 0)
        {

            if (n == 1) return (x[xstart] < 0.0 ? -x[xstart] : x[xstart]);
            else
            {
                double scale = 0.0;
                double ssq = 1.0;
                if (incx == 1) BlasLike.dsssqvec(n, x, ref scale, ref ssq, xstart);
                else BlasLike.dsssq(n, x, incx, ref scale, ref ssq, xstart);
                return sc_norm(scale, ssq);
            }
        }
        void dtqadd(bool orthog, ref int inform, int k1, int k2, ref int nactiv, ref int ncolz, ref int nfree, int n, double condmx)
        {
            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1;

            int iadd, jadd, ifix = 27, i, k, l, iswap;
            double cslast = 12, snlast = 8;

            /*     TQADD INCLUDES GENERAL LINEAR CONSTRAINTS  K1  THRU  K2  AS NEW */
            /*     COLUMNS OF THE TQ FACTORIZATION STORED IN  RT, ZY. */
            // --wrk2;
            //  --wrk1;
            zy_dim1 = nq;
            zy_offset = zy_dim1 + 1;
            //   zy -= zy_offset;
            rt_dim1 = NROWRT;
            rt_offset = rt_dim1 + 1;
            //   rt -= rt_offset;
            // --qtg;
            a_dim1 = NROWA;
            a_offset = a_dim1 + 1;
            //  a -= a_offset;
            //  --kfree;
            //  --kactiv;
            // --istate;

            i__1 = k2;
            for (k = k1; k <= i__1; ++k)
            {
                iadd = KACTV[k - 1];
                jadd = n + iadd;
                if (nactiv == nfree)
                {
                    goto L20;
                }
                inform = daddcon(false, false, orthog, ifix, iadd, jadd,
                    nactiv, ncolz, ncolz, nfree, n, NROWA, KFREE, condmx, cslast, snlast,
                     WRK, RLAM);
                if (inform > 0)
                {
                    goto L20;
                }
                ++nactiv;
                --(ncolz);
                goto L40;
            L20:
                ISTATE[jadd - 1] = 0;
                KACTV[k - 1] = -KACTV[k - 1];
            L40:
                ;
            }
            if (nactiv == k2)
            {
                return;
            }
            /*     SOME OF THE CONSTRAINTS WERE CLASSED AS DEPENDENT AND NOT INCLUDED 
            */
            /*     IN THE FACTORIZATION.  MOVE ACCEPTED INDICES TO THE FRONT OF */
            /*     KACTIV  AND SHIFT REJECTED INDICES (WITH NEGATIVE VALUES) TO */
            /*     THE END. */
            l = k1 - 1;
            i__1 = k2;
            for (k = k1; k <= i__1; ++k)
            {
                i = KACTV[k - 1];
                if (i < 0)
                {
                    goto L60;
                }
                ++l;
                if (l == k)
                {
                    goto L60;
                }
                iswap = KACTV[l - 1];
                KACTV[l - 1] = i;
                KACTV[k - 1] = iswap;
            L60:
                ;
            }
        }
        void drtmxsolve(int job, int n, double[] t, int nrt, double[] b, ref int idiag, int tstart = 0, int bstart = 0)
        {

            /*
                purpose
                =======
                drtmxsolve solves systems of reverse triangular equations

                description
                ===========
                drtmxsolve solves the equations

                    T*x = b ,   or   ( T' )*x = b ,

                where T is an n by n upper or lower reverse triangular matrix
                an upper reverse triangular matrix has the form illustrated by

                T = ( x  x  x  x  x )
                    ( x  x  x  x  0 )
                    ( x  x  x  0  0 )
                    ( x  x  0  0  0 )
                    ( x  0  0  0  0 )

                 and a lower reverse triangular matrix has the form illustrated by

                 T = ( 0  0  0  0  x )
                     ( 0  0  0  x  x )
                     ( 0  0  x  x  x )
                     ( 0  x  x  x  x )
                     ( x  x  x  x  x )

                 the type of triangular system solved is controlled by the
                 parameter job as described in the parameter section below

                 parameters
                 ==========
                 job   - integer
                    on entry, job must contain one of the values -2, -1, 0, 1, 2
                    to specify the type of reverse triangular system to be solved
                    as follows
                    job =  0 	T is assumed to be reverse diagonal and the
                            equations T*x = b are solved
                    job =  1	T is assumed to be upper reverse triangular and the
                            equations T*x = b are solved
                    job = -1	T is assumed to be upper reverse triangular and the
                            equations ( T' )*x = b are solved
                    job =  2	T is assumed to be lower reverse triangular and the
                            equations T*x = b are solved
                    job = -2	T is assumed to be lower reverse triangular and the
                            equations ( T' )*x = b are solved
                 n     - integer
                    on entry, n specifies the order of the matrix T. n must be at
                    least unity
                 T     - real array of dimension ( nrt, nct ). nct must be at least n
                    before entry T must contain the triangular elements. only
                    those elements contained in the reverse triangular part of T
                    are referenced by this routine
                    unchanged on exit

                 nrt   - integer
                    on entry, nrt specifies the first dimension of T as declared
                    in the calling program. nrt must be at least n
                 b     - real array of dimension ( n )
                    before entry, b must contain the right hand side of the
                    equations to be solved
                    on successful exit, b contains the solution vector x
                 idiag - integer
                    before entry, idiag must be assigned a value. for users
                    unfamiliar with this parameter the recommended value is zero
                    on successful exit idiag will be zero. a positive value
                    of idiag denotes an error as follows
                    idiag = 1     one of the input parameters n, or nra, or job
                            has been incorrectly specified
                    idiag > 1	the element T( j, n - j + 1 ), where
                            j = idiag - 1, is either zero or is too small
                            to avoid overflow in computing an element of x.
                            note that this element is a reverse diagonal element of T

                 further comments
                 ================
                 if T is part of a matrix a partitioned as

                    A =	( A1  A2 )
                            ( A3  T  )

                 where A1 is an m by k matrix ( m>=0, k>=0), then this routine
                 may be called with the parameter T as a( m + 1, k + 1 ) and nrt as
                 the first dimension of A as declared in the calling (sub) program.
            */
            int t_dim1, t_offset, i__1;

            int fail = -900;
            double temp;
            int j, k;
            int jlast;

            //        --b;
            bstart--;
            t_dim1 = nrt;
            t_offset = t_dim1 + 1;
            //        t -= t_offset;
            tstart -= t_offset;

            if (n >= 1 && nrt >= n && Math.Abs(job) <= 2)
            {
                goto L20;
            }
            idiag = (int)lm_check_fail((short)(idiag), (short)1, "drtmxsolve");
            return;
        L20:
            if (job != 0)
            {
                goto L60;
            }
            k = n;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                b[j + bstart] = BlasLike.dprotdiv(ref b[j + bstart], ref t[j + k * t_dim1 + tstart], ref fail);
                if (fail != 0)
                {
                    goto L280;
                }
                --k;
                /* L40: */
            }
            goto L140;
        L60:
            if (job != 1)
            {
                goto L100;
            }
            j = n;
            i__1 = n;
            for (k = 1; k <= i__1; ++k)
            {
                b[j + bstart] = BlasLike.dprotdiv(ref b[j + bstart], ref t[j + k * t_dim1 + tstart], ref fail);
                if (fail != 0)
                {
                    goto L280;
                }
                if (j > 1)
                {
                    // BlasLike.daxpy(j - 1, -b[j+bstart], &t[k * t_dim1 + 1+tstart], 1, &b[1+bstart], 1);
                    BlasLike.daxpyvec(j - 1, -b[j + bstart], t, b, k * t_dim1 + 1 + tstart, 1 + bstart);
                }
                --j;
                /* L80: */
            }
            goto L140;
        L100:
            if (job != 2)
            {
                goto L140;
            }
            k = n;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                b[j + bstart] = BlasLike.dprotdiv(ref b[j + bstart], ref t[j + k * t_dim1 + tstart], ref fail);
                if (fail != 0)
                {
                    goto L280;
                }
                if (j < (int)n)
                {
                    //   BlasLike.daxpy(n - j, -b[j+bstart], &t[j + 1 + k * t_dim1+tstart], 1, &b[j + 1+bstart], 1);
                    BlasLike.daxpyvec(n - j, -b[j + bstart], t, b, j + 1 + k * t_dim1 + tstart, j + 1 + bstart);
                }
                --k;
                /* L120: */
            }
        L140:
            if (n == 1)
            {
                goto L180;
            }
            jlast = n / 2;
            k = n;
            i__1 = jlast;
            for (j = 1; j <= i__1; ++j)
            {
                temp = b[j + bstart];
                b[j + bstart] = b[k + bstart];
                b[k + bstart] = temp;
                --k;
                /* L160: */
            }
        L180:
            if (job != -1)
            {
                goto L220;
            }
            k = n;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (j > 1) b[j + bstart] -= BlasLike.ddotvec((int)(j - 1), t, b, k * t_dim1 + 1 + tstart, 1 + bstart);
                b[j + bstart] = BlasLike.dprotdiv(ref b[j + bstart], ref t[j + k * t_dim1 + bstart], ref fail);
                if (fail != 0) goto L280;
                --k;
                /* L200: */
            }
            goto L260;
        L220:
            if (job != -2) goto L260;
            j = n;
            i__1 = n;
            for (k = 1; k <= i__1; ++k)
            {
                if (j < (int)n) b[j + bstart] -= BlasLike.ddotvec((int)(n - j), t, b, j + 1 + k * t_dim1 + tstart, j + 1 + bstart);
                b[j + bstart] = BlasLike.dprotdiv(ref b[j + bstart], ref t[j + k * t_dim1 + tstart], ref fail);
                if (fail != 0)
                {
                    goto L280;
                }
                --j;
                /* L240: */
            }
        L260:
            idiag = 0;
            return;
        L280:
            idiag = (int)lm_check_fail((short)(idiag), (short)(j + 1), "drtmxsolve");
        }
        short lm_check_fail(short ifail, short ierror, string srname)
        {
            if (ierror != 0)
            {
                /*Abnormal exit from calling routine */
                if (ifail <= 0)
                {
                    /*Noisy exit */
                    Console.WriteLine($"Mathematics routine {srname}: exited with ifail={ierror}",
                    srname, ierror);
                    if (ifail != 0) { ColourConsole.WriteInfo("Hard failure"); return -50; }
                }
                //#if 0
                else
                {   /*Soft failure*/
                    Console.WriteLine($"Mathematics routine {srname}: exited with ifail={ierror}",
                    srname, ierror);
                    ColourConsole.WriteInfo(" ** Soft failure - control returned");
                }
                //#endif
            }
            return ierror;
        }
        void delmgen(bool orthog, ref double x, ref double y, ref double cs, ref double sn)
        {/*
	If orthog delmgen generates a plane rotation.  Otherwise,
	delmgen  generates an elimination transformation  E  such that
	(X Y)*E  =  (X  0)   OR   (Y  0),  depending on the relative
	sizes of  X  &  Y.
*/
            if (orthog) BlasLike.drotg(ref x, ref y, ref cs, ref sn);
            else
            {
                cs = 1;
                sn = 0;
                if (y != 0)
                {
                    if (Math.Abs(x) < Math.Abs(y))
                    {
                        cs = 0;
                        sn = -x / y;
                        x = y;
                    }
                    else sn = -y / x;
                }
            }
            y = 0;
        }
        void delm(bool orthog, int n, double[] x, int incx, double[] y, int incy, double cs, double sn, int xstart = 0, int ystart = 0)
        {/*
	If  orthog  is true, delm  applies a plane rotation.  otherwise,
	elm computes the transformation (x y)*e  and returns the result
	in  (x y),  where the 2 by 2 matrix  e  is defined by  cs  and  sn 

	as follows...
	e  = 	( 1  sn )	if  cs>0 else	e  =	(     1 )
		(     1 )				( 1  sn )
*/

            if (!orthog)
            {
                if (cs <= 0) BlasLike.dswap(n, x, incx, y, incy, xstart, ystart);
                if (sn != 0) BlasLike.daxpy(n, sn, x, incx, y, incy, xstart, ystart);
            }
            else dsymplanerotate(n, x, incx, y, incy, cs, sn, xstart, ystart);
        }
        void dsymplanerotate(int n, double[] x, int incx, double[] y, int incy, double c, double s, int xstart = 0, int ystart = 0)
        {/*
	dsymplanerotate performs the symmetric plane rotation
	( x  y ) = ( x  y )*( c   s )   s != 0
                    	    ( s  -c )

	If s is supplied as zero then x and y are unaltered
*/
            int i__1, i__2;

            double temp1;
            int i, ix, iy;
            //        --y;
            //        --x;
            ystart -= 1;
            xstart -= 1;

            if (n > 0 && s != 0.0)
            {
                if (c == 0 && s == 1)
                {
                    if (incx == incy && incx > 0)
                    {
                        i__1 = (n - 1) * incx + 1;
                        i__2 = incx;
                        for (ix = 1; (i__2 < 0 ? ix >= i__1 : ix <= i__1); ix += i__2)
                        {
                            temp1 = x[ix + xstart];
                            x[ix + xstart] = y[ix + ystart];
                            y[ix + ystart] = temp1;
                        }
                    }
                    else
                    {
                        if (incy >= 0)
                        {
                            iy = 1;
                        }
                        else
                        {
                            iy = 1 - (n - 1) * incy;
                        }
                        if (incx > 0)
                        {
                            i__2 = (n - 1) * incx + 1;
                            i__1 = incx;
                            for (ix = 1; (i__1 < 0 ? ix >= i__2 : ix <= i__2); ix += i__1)
                            {
                                temp1 = x[ix + xstart];
                                x[ix + xstart] = y[iy + ystart];
                                y[iy + ystart] = temp1;
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[ix + xstart];
                                x[ix + xstart] = y[iy + ystart];
                                y[iy + ystart] = temp1;
                                ix += incx;
                                iy += incy;
                            }
                        }
                    }
                }
                else if (c == 0 && s == -1)
                {
                    if (incx == incy && incx > 0)
                    {
                        i__1 = (n - 1) * incx + 1;
                        i__2 = incx;
                        for (ix = 1; (i__2 < 0 ? ix >= i__1 : ix <= i__1); ix += i__2)
                        {
                            temp1 = -x[ix + xstart];
                            x[ix + xstart] = -y[ix + ystart];
                            y[ix + ystart] = temp1;
                        }
                    }
                    else
                    {
                        if (incy >= 0)
                        {
                            iy = 1;
                        }
                        else
                        {
                            iy = 1 - (n - 1) * incy;
                        }
                        if (incx > 0)
                        {
                            i__2 = (n - 1) * incx + 1;
                            i__1 = incx;
                            for (ix = 1; (i__1 < 0 ? ix >= i__2 : ix <= i__2); ix += i__1)
                            {
                                temp1 = -x[ix + xstart];
                                x[ix + xstart] = -y[iy + ystart];
                                y[iy + ystart] = temp1;
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = -x[ix + xstart];
                                x[ix + xstart] = -y[iy + ystart];
                                y[iy + ystart] = temp1;
                                ix += incx;
                                iy += incy;
                            }
                        }
                    }
                }
                else
                {
                    if (incx == incy && incx > 0)
                    {
                        i__1 = (n - 1) * incx + 1;
                        i__2 = incx;
                        for (ix = 1; (i__2 < 0 ? ix >= i__1 : ix <= i__1); ix += i__2)
                        {
                            temp1 = x[ix + xstart];
                            x[ix + xstart] = c * temp1 + s * y[ix + ystart];
                            y[ix + ystart] = s * temp1 - c * y[ix + ystart];
                        }
                    }
                    else
                    {
                        if (incy >= 0)
                        {
                            iy = 1;
                        }
                        else
                        {
                            iy = 1 - (n - 1) * incy;
                        }
                        if (incx > 0)
                        {
                            i__2 = (n - 1) * incx + 1;
                            i__1 = incx;
                            for (ix = 1; (i__1 < 0 ? ix >= i__2 : ix <= i__2); ix += i__1)
                            {
                                temp1 = x[ix + xstart];
                                x[ix + xstart] = c * temp1 + s * y[iy + ystart];
                                y[iy + ystart] = s * temp1 - c * y[iy + ystart];
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[ix + xstart];
                                x[ix + xstart] = c * temp1 + s * y[iy + ystart];
                                y[iy + ystart] = s * temp1 - c * y[iy + ystart];
                                ix += incx;
                                iy += incy;
                            }
                        }
                    }
                }
            }
        }
        short dtmxsolve(short job, int n, double[] t, int nrt, double[] b, short idiag, int tstart = 0, int bstart = 0)
        {
            /*
                purpose
                =======
                    dtmxsolve solves systems of triangular equations
                    important.	in double precision implementations the real
                    ---------	declarations should be interpreted to mean
                            double precision
                description
                ===========
                    dtmxsolve solves the equations
                    t*x = b ,   or	 ( t**t )*x = b ,
                    where t is an n by n upper or lower triangular matrix
                    the type of triangular system solved is controlled by the
                    parameter job as described in the parameter section below
                parameters
                ==========
                    job   - short
                        on entry, job must contain one of the values -2, -1, 0, 1, 2
                        to specify the type of triangular system to be solved as
                        follows
                        job =  0	t is assumed to be diagonal and the
                                equations t*x = b are solved
                        job =  1	t is assumed to be upper triangular and the
                                equations t*x = b are solved
                        job = -1	t is assumed to be upper triangular and the
                                equations ( t**t )*x = b are solved
                        job =  2	t is assumed to be lower triangular and the
                                equations t*x = b are solved
                        job = -2	t is assumed to be lower triangular and the
                                equations ( t**t )*x = b are solved
                                unchanged on exit
                    n     - dimen
                        on entry, n specifies the order of the matrix t. n must be at
                        least unity
                        unchanged on exit
                    t     - real matrix of dimension ( nrt, nct ). nct must be at least n
                        before entry t must contain the triangular elements. only
                        those elements contained in the triangular part of t are
                        referenced by this routine
                        unchanged on exit
                    nrt   - dimen
                        on entry, nrt specifies the first dimension of t as declared
                        in the calling (sub) program. nrt must be at least n
                    b     - real vector of dimension ( n )
                        before entry, b must contain the right hand side of the
                        equations to be solved
                        on successful exit, b contains the solution vector x
                    idiag - short
                        before entry, idiag must be assigned a value. for users
                        unfamiliar with this parameter (described in chapter p01)
                        the recommended value is zero

                    return value
                        0	success
                        1	invalid arguments (n, nrt, job)
                        > 1	(return - 1)th diagonal element of t is
                            either zero or is too small to avoid overflow
                            in computing an element of x

                    further comments
                    ================
                        if t is part of a matrix a partitioned as
                        a =	( a1  a2 ) ,
                            ( a3  t	 )
                        where a1 is an m by k matrix ( m>=0, k>=0), then this routine
                        may be called with the parameter t as a( m + 1, k + 1 ) and nrt as
                        the first dimension of a as declared in the calling (sub) program.
            */

            int fail = -99;
            int k;

            //       --b;
            //       t -= nrt + 1;
            bstart--;
            tstart -= nrt + 1;
            if (n < 1 || nrt < n || Math.Abs(job) > 2) return lm_check_fail(idiag, (short)1, "dtmxsolve");

            if (job == 1 || job == -2)
            {
                for (k = n; k >= 1; --k) if ((fail = (b[k + bstart] != 0.0 ? 1 : 0)) != 0) break;
            }
            else
            {
                for (k = 1; k <= n; ++k) if ((fail = (b[k + bstart] != 0.0 ? 1 : 0)) != 0) break;
            }

            if (fail != 0)
                switch (job)
                {
                    case 0:
                        for (; k <= n; ++k)
                        {
                            b[k + bstart] = BlasLike.dprotdiv(ref b[k + bstart], ref t[k + k * nrt + tstart], ref fail);
                            if (fail != 0) break;
                        }
                        break;
                    case 1:
                        for (; k >= 1; --k)
                        {
                            b[k + bstart] = BlasLike.dprotdiv(ref b[k + bstart], ref t[k + k * nrt + tstart], ref fail);
                            if (fail != 0) break;
                            if (k > 1) BlasLike.daxpyvec(k - 1, -b[k + bstart], t, b, k * nrt + 1 + tstart, 1 + bstart);
                        }
                        break;
                    case 2:
                        for (; k <= n; ++k)
                        {
                            b[k + bstart] = BlasLike.dprotdiv(ref b[k + bstart], ref t[k + k * nrt + tstart], ref fail);
                            if (fail != 0) break;
                            if (k < n) BlasLike.daxpyvec(n - k, -b[k + bstart], t, b, k + 1 + k * nrt + tstart, k + 1 + bstart);
                        }
                        break;
                    case -1:
                        for (; k <= n; ++k)
                        {
                            if (k > 1) b[k + bstart] -= BlasLike.ddotvec((k - 1), t, b, k * nrt + 1 + tstart, 1 + bstart);
                            b[k + bstart] = BlasLike.dprotdiv(ref b[k + bstart], ref t[k + k * nrt + tstart], ref fail);
                            if (fail != 0) break;
                        }
                        break;
                    case -2:
                        for (; k >= 1; --k)
                        {
                            if (k < n) b[k + bstart] -= BlasLike.ddotvec((n - k), t, b, k + 1 + k * nrt + tstart, k + 1 + bstart);
                            b[k + bstart] = BlasLike.dprotdiv(ref b[k + bstart], ref t[k + k * nrt + tstart], ref fail);
                            if (fail != 0) break;
                        }
                        break;
                }

            return (short)(fail != 0 ? lm_check_fail(idiag, (short)(k + 1), "dtmxsolve") : 0);
        }
        void dbdpert(bool firstv, bool negstp, double bigalf, double pnorm, ref int jadd1, ref int jadd2, ref double palfa1, ref double palfa2, int n, int nctotl)
        {

            /*
                dbdpert  finds steps  palfa1, palfa2  such that
                the point  X + palfa1*P  reaches a linear constraint that is
                currently not in the working set but is satisfied,
                the point  X + palfa2*P  reaches a linear constraint that is
                currently not in the working set but is violated.
                The constraints are perturbed by an amount  featol, so that
                palfa1  is slightly larger than it should be, and
                palfa2  is slightly smaller than it should be.  This gives
                some leeway later when the exact steps are computed by dbndalf.
                Constraints in the working set are ignored  (istate(J) >= 1).
                If  negstp  is true, the search direction will be taken to be -P.

                VALUES OF ISTATE(J)....
                - 2         - 1         0           1          2         3
                A*X < BL   A*X > BU   A*X free   A*X = BL   A*X = BU   BL = BU

                The values -2 and -1 do not occur once dlpcore finds a feasible point.
            */
            var bigbnd = parm[0];
            int i, j;
            bool lastv;
            double epspt9;
            int js;
            double absatp, rownrm, atp, res, atx;
            //  --x;
            // --p;
            // --featol;
            // --bu;
            //            --bl;
            //  --ax;
            // --ap;
            //        --anorm;
            // --istate;

            epspt9 = parm[3];
            if (msg == 99) ColourConsole.WriteInfo(
        "\n   J  JS         FEATOL         AX             AP     JADD1       PALFA1     JADD2       PALFA2\n");
            lastv = !firstv;
            jadd1 = 0;
            jadd2 = 0;
            palfa1 = bigalf;
            palfa2 = 0.0;
            if (firstv) palfa2 = bigalf;
            for (j = 1; j <= nctotl; ++j)
            {
                js = ISTATE[j - 1];
                if (js > 0) continue;
                if (j > n)
                {
                    /*GENERAL LINEAR CONSTRAINT. */
                    i = j - n;
                    atx = LWRK[i - 1];
                    atp = AP[i - 1];
                    rownrm = 1.0 + ANORM[i - 1];
                }
                else
                {
                    /*BOUND CONSTRAINT. */
                    atx = W[j - 1];
                    atp = PX[j - 1];
                    rownrm = 1.0;
                }
                if (negstp) atp = -atp;
                if (Math.Abs(atp) <= epspt9 * rownrm * pnorm) res = -1.0;
                else if (atp > BlasLike.lm_eps)
                {
                    /*AX IS INCREASING*/
                    /*TEST FOR SMALLER PALFA1 IF UPPER BOUND IS SATISFIED. */
                    if (js != -1)
                    {
                        if (U[j - 1] < bigbnd)
                        {
                            res = U[j - 1] - atx + FEATOL[j - 1];
                            if (bigalf * atp > Math.Abs(res) && palfa1 * atp > res)
                            {
                                palfa1 = res / atp;
                                jadd1 = j;
                            }
                        }
                        /*TEST FOR DIFFERENT PALFA2 IF LOWER BOUND IS VIOLATED. */
                        if (js == -2)
                        {
                            res = L[j - 1] - atx - FEATOL[j - 1];
                            if (bigalf * atp > Math.Abs(res)
                                && (firstv || palfa2 * atp < res)
                                && (lastv || palfa2 * atp > res))
                            {
                                palfa2 = res / atp;
                                jadd2 = j;
                            }
                        }
                    }
                }
                else if (js != -2)
                {
                    /*AX IS DECREASING TEST FOR SMALLER PALFA1 IF LOWER BOUND IS SATISFIED*/
                    absatp = -atp;
                    if (L[j - 1] > -bigbnd)
                    {
                        res = atx - L[j - 1] + FEATOL[j - 1];
                        if (bigalf * absatp > Math.Abs(res) && palfa1 * absatp > res)
                        {
                            palfa1 = res / absatp;
                            jadd1 = j;
                        }
                    }
                    /*TEST FOR DIFFERENT PALFA2 IF UPPER BOUND IS VIOLATED*/
                    if (js == -1)
                    {
                        res = atx - U[j - 1] - FEATOL[j - 1];
                        if (bigalf * absatp > Math.Abs(res)
                            && (firstv || palfa2 * absatp < res)
                            && (lastv || palfa2 * absatp > res))
                        {
                            palfa2 = res / absatp;
                            jadd2 = j;
                        }
                    }
                }
                if (msg == 99)
                    lm_wmsg("%5ld%4ld%15.5lg%15.5lg%15.5lg%6ld%17.7lg%6ld%17.7lg",
                    j, js, FEATOL[j - 1], atx, atp,
                    jadd1, palfa1, jadd2, palfa2);
            }
        }
        void dhhrflctgen(int n, ref double alpha, double[] x, int incx, double tol, ref double z1, int xstart = 0)
        {

            /*
                dhhrflctgen generates details of a householder reflection, p, such that
                    p*( alpha ) = ( beta ),   p'*p = i
                      (   x   )   (   0  )

                p is given in the form
                    p = i - ( 1/z( 1 ) )*z*z',

                where z is an ( n + 1 ) element vector
                z( 1 ) is returned in z1. if the elements of x are all zero, or if
                the elements of x are all less than tol*abs( alpha ) in absolute
                value, then z1 is returned as zero and p can be taken to be the
                unit matrix. otherwise z1 always lies in the range ( 1.0, 2.0 )
                if tol is not in the range ( 0.0, 1.0 ) then the value 0.0 is used in
                place of tol
                the remaining elements of z are overwritten on x and beta is
                overwritten on alpha
            */
            double beta;
            var work = new double[1];
            double scale, tl, ssq;

            //      --x;
            xstart--;

            if (n < 1)
            {
                z1 = 0;
            }
            else
            {
                if (tol <= 0 || tol > 1)
                {
                    tl = 0;
                }
                else
                {
                    tl = Math.Abs(alpha) * tol;
                }
                ssq = 1;
                scale = 0;
                BlasLike.dsssq(n, x, incx, ref scale, ref ssq, 1 + xstart);
                if (scale == 0 || scale < tl)
                {
                    z1 = 0;
                }
                else
                {
                    if (alpha != 0)
                    {
                        work[0] = alpha;
                        BlasLike.dsssqvec(1, work, ref scale, ref ssq);
                        beta = -BlasLike.dsign(sc_norm(scale, ssq), alpha);
                        z1 = (beta - alpha) / beta;
                    }
                    else
                    {
                        beta = -sc_norm(scale, ssq);
                        z1 = 1;
                    }
                    BlasLike.dscal(n, -1.0 / beta, x, incx, 1 + xstart);
                    alpha = beta;
                }
            }
        }
        short dqpsol(int itmax, short msglvl, int n, int nclin, int nctotl, int nrowa, int nrowh, int ncolh, int cold, int lp, int orthog, ref int iter, ref double obj, int leniw, int lenw, short ifail)
        {
            parm[0] = 1e10;
            parm[1] = 1e20;
            parm[2] = 1e-6;
            parm[3] = Math.Pow(BlasLike.lm_eps, 0.9);
            /*
                dqpsol solves quadratic programming (QP) problems of the form 
                minimize     c'*x  +  1/2 x'*H*X 
                subject to		(  x  ) 
                        BL  <=  (     ) <=  BU 
                            ( A*x ) 

                where ' denotes the transpose of a column vector. 
                The symmetric matrix  H  may be positive-definite, positive 
                semi-definite, or indefinite. 
                n  is the number of variables (dimension of  x). 
                nclin  is the number of general linear constraints (rows of  a). 
                (nclin may be zero.) 

                The matrix   H  is defined by the subroutine  qphess, which 
                must compute the matrix-vector product  H*x  for any vector  x. 
                the vector  c  is entered in the one-dimensional array  cvec. 
                the first  n  components of  bl  and   bu  are lower and upper 
                bounds on the variables.  the next  nclin  components are 
                lower and upper bounds on the general linear constraints. 
                the matrix  a  of coefficients in the general linear constraints 
                is entered as the two-dimensional array a (of dimension 
                nrowa  by  n). if nclin = 0,  a  is not accessed. 
                the vector  x  must contain an initial estimate of the solution, 
                and will contain the computed solution on output. 
            */
            int itmx;
            string l = (lp != 0 ? ((lp & 2) != 0 ? "FP" : "LP") : "QP");

            //double bigdx;
            int nfree = -2, ncnln;
            xnorm = new double[1];// epspt9;
            int maxact, minact;
            byte lcrash;
            //double tolact;
            int minfxd, inform = -2, mxfree, nactiv = -3, numinf = 23;
            var litotl = 0;
            int minsum;
            int mxcolz;
            int vertex;
            int lax;

            //#define NCLIN &nclin_
            int nclin_ = nclin;
            //#define NCTOTL &nctotl_
            int nctotl_ = nctotl;
            //#define NROWA &nrowa_
            NROWA = nrowa;
            int nrowa_ = nrowa;
            int iter_ = 0;
            //#define NROWH &nrowh_
            int nrowh_ = nrowh;
            //#define NCOLH &ncolh_
            int ncolh_ = ncolh;
            // --w;
            // --iw;
            //--clamda;
            //   --istate;
            // --x;
            //       --featol;
            //  --cvec;
            //         --bu;
            //        --bl;

            /*IF ITMAX IS NOT POSITIVE ON ENTRY SET IT TO 50*/
            itmx = itmax;
            if (itmx <= 0) itmx = 50;

            /*
            IF THERE IS NO FEASIBLE POINT FOR THE LINEAR CONSTRAINTS AND
            BOUNDS, COMPUTE THE MINIMUM SUM OF INFEASIBILITIES
            IT IS NOT NECESSARY TO START THE QP PHASE AT A VERTEX
            */
            minsum = (lp & 2) != 0 ? 0 : 1;
            vertex = 0;
            /*
            ANY CHANGE IN X THAT IS GREATER THAN  BIGDX  WILL BE REGARDED
            AS AN INFINITE STEP
            */
            //bigdx = 1e20;

            /*
            DURING SELECTION OF THE INITIAL WORKING SET (BY CRASH),
            CONSTRAINTS WITH RESIDUALS LESS THAN  TOLACT  WILL BE MADE ACTIVE. 
            */
            //tolact = .01;
            //epspt9 = pow(lm_eps, 0.9);
            /*	parm[0] = *bigbnd;
                parm[1] = bigdx;
                parm[2] = tolact;
                parm[3] = epspt9;*/

            /*
            assign the dimensions of arrays in the parameter list of qpcore
            economies of storage are possible if the minimum number of active
            constraints and the minimum number of fixed variables are known in
            advance.  the expert user should alter  minact  and  minfxd
            accordingly
            if a linear program is being solved and the matrix of general
            constraints is fat,  i.e.,  nclin < n,  a non-zero value is
            known for  minfxd.  note that in this case,  vertex  must be 1
            */
            minact = 0;
            minfxd = 0;
            /*	Vadim Moroz's data showed that this can lead to errors!  Colin 19-7-2000
                if (lp && nclin < n)
                {
                    minfxd = n - nclin - 1;
                    vertex = 1;
                }*/
            mxfree = n - minfxd;
            maxact = Math.Min(n, nclin);
            maxact = Math.Max(1, maxact);
            mxcolz = n - (minfxd + minact);
            nq = Math.Max(1, mxfree);
            NROWRT = Math.Max(mxcolz, maxact);
            NCOLRT = Math.Max(1, mxfree);
            ncnln = 0;

            /*allocate certain arrays that are not done in  alloc*/
            litotl = 0;
            lax = 1;
            var lwtotl = lax + nrowa - 1;
            /*allocate remaining work arrays*/
            dalloc(2, n, nclin, ncnln, nctotl, ref litotl, ref lwtotl);
            /*set the message level for  lpdump, qpdump, chkdat  and  lpcore*/
            msg = 0;
            if (msglvl >= 5) msg = 5;
            if (lp != 0 || msglvl >= 15) msg = msglvl;
            /*
            *** the following statement must be executed if  istart   ***
            *** is not set in the calling routine.                    ***
            */
            istart = 0;
            lcrash = 1;
            if (cold != 0) lcrash = 0;
            /*check input parameters and storage limits*/
            if (msglvl == 99)
                dlpdump(n, nclin, nctotl, nrowa, lcrash, lp, minsum,
                vertex);
            if (msglvl == 99)
                dqpdump(n, nrowh, ncolh);
            if (badboundcount() > 0) return 10;
            iter = 0;


            /*
            no scaling is provided by this version of  dqpsol
            give a fake value for the start of the scale array
            */
            scldqp = false;


            /*
            ---------------------------------------------------------------------
            call  lpcore  to obtain a feasible point, or solve a linear
            problem
            ---------------------------------------------------------------------
            */
            dlpcore((lp & 1) != 0, minsum, orthog != 0, vertex, ref inform, ref iter_,
            itmx, lcrash, n, nclin, ref nctotl, ref nactiv, ref nfree, ref numinf,
                 ref obj, xnorm);
            iter = (short)iter_;
            if (lp != 0)
            {
                /*THE PROBLEM WAS AN LP, NOT A QP*/
                if (inform > 2) inform += 4;
                if (inform == 1) inform = 6;
            }
            else if (inform == 0)
            {
                /*
                ---------------------------------------------------------------------
                call  qpcore  to solve a quadratic problem
                ---------------------------------------------------------------------
                */
                msg = msglvl;
                /*
                *** the following statement must be executed if  istart   ***
                *** is not set in the calling routine.                    ***
                */
                istart = 0;

                dqpcore(orthog, ref inform, ref iter_, itmx, n, nclin, nctotl,
                         nrowh, ncolh, nactiv, nfree,
                         ref obj, xnorm);
                iter = (short)iter_;
            }
            else
            {
                /*
                trouble in  lpcore
                inform cannot be given the value  2  when finding a feasible
                point, so it is necessary to decrement all the values of  inform
                that are greater than  2
                */
                if (inform > 2) --inform;
                inform += 5;
            }
            /*print messages if required*/
            if (msglvl > 0)
            {
                switch (inform)
                {
                    case 0:
                        lm_wmsg("\nEXIT QPSOL- OPTIMAL %.2s SOLUTION.", l);
                        break;
                    case 1:
                        ColourConsole.WriteInfo("WEAK LOCAL MINIMUM.");
                        break;
                    case 2:
                        lm_wmsg("\nEXIT dqpsol- %.2s SOLUTION IS UNBOUNDED.", l);
                        break;
                    case 3:
                        ColourConsole.WriteInfo("ZERO MULTIPLIERS.");
                        break;
                    case 4:
                        ColourConsole.WriteInfo("TOO MANY ITERATIONS WITHOUT CHANGING X.");
                        break;
                    case 5:
                        ColourConsole.WriteInfo("TOO MANY ITERATIONS.");
                        break;
                    case 6:
                        ColourConsole.WriteInfo("CANNOT SATISFY THE LINEAR CONSTRAINTS.");
                        break;
                    case 7:
                        ColourConsole.WriteInfo("TOO MANY ITERATIONS WITHOUT CHANGING X IN THE LP PHASE.");
                        break;
                    case 8:
                        ColourConsole.WriteInfo("TOO MANY ITERATIONS DURING THE LP PHASE.");
                        break;
                }
                if (numinf == 0) lm_wmsg("\n FINAL %.2s OBJECTIVE VALUE =%20.9lg", l, obj);
                else lm_wmsg("\n FINAL SUM OF INFEASIBILITIES =%20.9lg", obj);
            }

            return (short)(inform == 0 ? 0 : lm_check_fail((short)ifail, (short)inform, "QPSOL"));
        }
        void dlpdump(int n, int nclin, int nctotl, int nrowa, int lcrash, int lp, int minsum, int vertex)
        {
            int j, k;
            double atx;

            ColourConsole.WriteInfo("\n\n\n\n\n\nOUTPUT FROM LPDUMP\n ******************");
            lm_wmsg("\nLCRASH =%d LP=%d MINSUM=%d VERTEX=%d", lcrash, lp, minsum, vertex);

            /*PRINT  A  BY ROWS AND COMPUTE  AX = A*X. */
            for (k = 0; k < nclin; ++k)
            {
                lm_wmsg("\nROW%6ld OF A ...", k + 1);
                lm_gdvwri(n, A, nrowa, k);
                LWRK[k] = BlasLike.ddot(n, A, nrowa, W, 1, k);
            }

            /*PRINT  BL, BU  AND  X OR AX. */
            ColourConsole.WriteInfo("\n              J      BL(J)          BU(J)           X(J)");
            for (j = 0; j < nctotl; ++j)
            {
                if (j < n)
                {
                    k = j;
                    atx = W[j];
                }
                else
                {
                    k = j - n;
                    atx = LWRK[k];
                    if (k != 0)
                        ColourConsole.WriteInfo("\n              I    BL(N+I)        BU(N+I)         A(I)*X");
                }
                lm_wmsg("", k + 1, L[j], U[j], atx);
            }
            if ((lp & 1) != 0) lm_mdvwri("\nCVEC ...", n, c);
            if (lcrash != 0) lm_mdvwri("\nISTATE ...", nctotl, ISTATE);
        }
        void lm_gdvwri(int n, double[] x, int inc, int xstart = 0)
        {
            int i, ii;
            for (i = 0, ii = 0; i < n; i++)
            {
                Console.Write($"{x[xstart + ii]} ");
                ii += inc;
                if (i % 6 == 5) Console.Write("\n");
            }
            Console.Write("\n");
        }
        void dalloc(byte nalg, int n, int nclin, int ncnln, int nctotl, ref int litotl, ref int lwtotl)
        {
            int lqtg, lwrk, lrlam, lkfree, lkactv, lanorm, lap, lrt, lpx, lzy;
            lkactv = litotl + 1;
            KACTV = new int[n];
            lkfree = lkactv + n;
            KFREE = new int[n];
            litotl = lkfree + n - 1;
            lanorm = lwtotl + 1;
            ANORM = new double[nclin];
            lap = lanorm + nclin;
            AP = new double[nclin];
            lpx = lap + nclin;
            PX = new double[n];
            lqtg = lpx + n;
            QTG = new double[n];
            lrlam = lqtg + n;
            RLAM = new double[n];
            lrt = lrlam + n;
            RT = new double[NROWRT * NCOLRT];
            lzy = lrt + NROWRT * NCOLRT;
            ZY = new double[nq * nq];
            lwrk = lzy + nq * nq;
            WRK = new double[n];
            lwtotl = lwrk + n - 1;
        }
        void dqpdump(int n, int nrowh, int ncolh)
        {
            int j, i;

            ColourConsole.WriteInfo("\n\n\n\n\n\nOUTPUT FROM QPDUMP\n******************");
            lm_mdvwri("\nCVEC ...", n, c);

            /*PRINT  HESS  UNLESS IT APPEARS TO BE IMPLICIT. */
            lm_wmsg("\nNROWH =%6ld NCOLH =%6ld", nrowh, ncolh);
            if (Q == null) return;
            if (nrowh > 1 || ncolh > 1)
            {
                if (ncolh == 1) lm_mdvwri("\nHESS ...", nrowh, Q);
                else
                {
                    i = Math.Min(ncolh, n);
                    for (j = 1; j <= i; ++j)
                    {
                        lm_wmsg("\nCOLUMN%6ld OF  HESS ...", j);
                        lm_gdvwri(i, Q, nrowh, j - 1);
                    }
                }
            }
            /*CALL  QPHESS  TO COMPUTE EACH COLUMN OF THE HESSIAN. */
            ColourConsole.WriteInfo("\n\n THE FOLLOWING IS RETURNED BY  QPHESS.");
            BlasLike.dzerovec(n, WRK);
            for (j = 1; j <= n; ++j)
            {
                WRK[i = j - 1] = 1.0;
                qphess(n, nrowh, ncolh, j, Q, WRK, PX);
                lm_wmsg("\nCOLUMN%6ld FROM  QPHESS ...", j);
                lm_gdvwri(n, PX, 1);
                WRK[i] = 0.0;
            }
        }
        void qphess(int n, int nrowh, int ncolh, int j, double[] hess, double[] wrk, double[] hx,int hstart=0,int xstart=0)
        {
            h(n, nrowh, ncolh, j, hess, wrk, hx,hstart:hstart,xstart:xstart);
        }
        public void qphess1(int n, int nrowh, int ncolh, int j, double[] hess, double[] wrk, double[] hx,int hstart=0,int xstart=0)
        {
            Solver.Factorise.dsmxmulv(n, hess, wrk, hx,ystart:hstart,xstart:xstart);
        }
        void dqpcore(int orthog, ref int inform, ref int iter, int itmax, int n, int nclin, int nctotl, int nrowh, int ncolh, int nactiv, int nfree, ref double objqp, double[] xnorm)
        {

            /*
                dqpcore, a subroutine for indefinite quadratic programming.
                it is assumed that a previous call to either  dlpcore  or  dqpcore

                has defined an initial working set of linear constraints and
                bounds. istate, kactiv	and  kfree  will have been set
                accordingly, and the arrays  rt	 and  zy  will contain the TQ
                factorization of the matrix whose rows are the gradients of the
                active linear constraints with the columns corresponding to the
                active bounds removed.	the TQ factorization of the resulting
                (nactiv by nfree) matrix is  a(free)*q = (0 t),	 where	q  is
                (nfree by nfree) and  t	 is reverse-triangular.
                values of istate(j) for the linear constraints.......

                istate(j)
                ---------
                0		constraint	j	is not in the working set
                1		constraint	j	is in the working set at its lower bound
                2		constraint	j	is in the working set at its upper bound
                3		constraint	j	is in the working set as an equality
                4		variable	j	is temporarily fixed at the value x(j)
                                    the corresponding artificial bound is included in the
                                    working set (the  TQ  factorization is adjusted accordingly)

                constraint  j  may be violated by as much as  featol(j)
            */
            string lprob = "QP";
            int mstall = 2000;

            int i;

            int iadd, jadd;
            double alfa;
            int jdel, kdel;
            double emax;
            int ifix = 9;
            double palfa = -12;
            int ifail = 12;
            double condh, bigdx;
            int isdel = 0;
            double condt, drmin, anorm, dinky, drmax, hsize;
            int ncnln, ncolr, ncolz;
            double pnorm = 1e24;
            //	int nrowj;
            int nhess;
            bool nullr, stall, uncon;
            int nclin0, kb;
            double bigalf, epspt9;
            double gfixed = 89, alfhit = 4;
            bool refine;
            int jdsave, nfixed;
            bool posdef;
            int modfyg, negligible;
            double condmx, atphit = 234, cslast = 5, gfnorm, rdlast = 12, objsiz, snlast = 6;
            int idummy, issave = 0, msglvl;
            int jsmlst = 9, ksmlst = 40;
            double smllst = 6e23;
            int nstall, numinf;
            double ztgnrm;
            int firstv, hitlow = 99;
            bool nocurv, renewr = false, unitpg, zerolm;
            int modfyr;
            double bnd;
            double gtp = 0;
            //         --w;
            //  --iw;
            //  --x;
            // --scale;
            //  --featol;
            // --cvec;
            //  --clamda;
            //--bu;
            //        --bl;
            //        --ax;
            // --kfree;
            //  --kactiv;
            //    --istate;

            /*INITIALIZE */
            iter = 0;
            jadd = 0;
            jdel = 0;
            jdsave = 0;
            nclin0 = Math.Max(nclin, 1);
            ncnln = 0;
            ncolz = nfree - nactiv;
            //	nrowj = 1;
            nstall = 0;
            nhess = 0;
            numinf = 0;
            msglvl = msg;
            msg = 0;
            if (istart == 0) msg = msglvl;
            var bigbnd = parm[0];
            bigdx = parm[1];
            epspt9 = parm[3];
            alfa = 0;
            condmx = BlasLike.lm_max;
            drmax = 1;
            drmin = 1;
            emax = 0;
            hsize = 1;
            firstv = 0;
            modfyr = 1;
            modfyg = 1;
            nocurv = false;
            nullr = false;
            posdef = true;
            refine = false;
            stall = true;
            uncon = false;
            unitpg = false;
            zerolm = false;

            /*
            given the  TQ  factorization of the matrix of constraints in the
            working set, compute the following quantities....
            (1)	the cholesky factor  r,	 of  z(t)hz  (if  z(t)hz  is not
                positive definite, find a positive-definite  (ncolr)-th	 order
                principal submatrix of	z(t)h z,
            (2)	the  qp	 objective function,
            (3)	the vector  q(free)(t)g(free),
            (4)	the vector  g(fixed).
                use the array  rlam  as temporary work space.
            */
            int nhess_conv = nhess;
            ncolr = dqpcrsh(n, ncolz, nfree, ref nhess_conv,
                nq, nrowh, ncolh, NROWRT, ref hsize,
                LWRK);
            dqpgrad(1, n, nactiv, nfree, ref nhess_conv, nq,
                nrowh, ncolh, jadd, alfa, ref objqp, ref gfixed,
                gtp, RLAM);
            nhess = nhess_conv;

            /*
            during the main loop, one of three things will happen
                (i)	the convergence criterion will be satisfied and the
                    algorithm will terminate.
                (ii)	a linear constraint will be deleted.
                (iii)	a direction of search will be computed and a constraint may
                    be added to the working set (note that a zero step may be taken
                    along the search direction).

            these computations occur in sections i, ii, and iii of the main loop
            */
            while (true)
            {
                /*
                ******* section i.  test for convergence *******
                compute the norms of the projected gradient and the gradient with
                respect to the free variables.
                */
                ztgnrm = 0;

                if (ncolr > 0) ztgnrm = dnrm2vec(ncolr, QTG);
                gfnorm = ztgnrm;
                if (nfree > 0 && nactiv > 0) gfnorm = dnrm2vec(nfree, QTG);

                /*
                define small quantities that reflect the magnitude of  c,  x,  h
                and the matrix of constraints in the working set
                */
                objsiz = (BlasLike.lm_eps + Math.Abs(objqp)) / (BlasLike.lm_eps + xnorm[0]);
                anorm = 0;
                if (nactiv > 0) anorm = Math.Abs(dtmax);
                /*Computing MAX */
                dinky = Math.Max(anorm, objsiz);
                dinky = Math.Max(epspt9, AccuracyModify) * Math.Max(dinky, gfnorm);

                if (msg >= 80) wdinky("QPCORE", ztgnrm, dinky);

                /*
                print the details of this iteration
                use the largest and smallest diagonals of r to estimate the
                condition number of the projected hessian matrix
                */
                condt = BlasLike.dprotdiv(ref dtmax, ref dtmin, ref ifail);
                if (ifail != 0 && dtmax == 0) condt = BlasLike.lm_max;
                if (ncolr > 0) BlasLike.dxminmax(ncolr, RT, NROWRT + 1, ref drmax, ref drmin);
                condh = BlasLike.dprotdiv(ref drmax, ref drmin, ref ifail);
                if (ifail != 0 && drmax == 0) condh = BlasLike.lm_max;
                if (condh >= BlasLike.lm_rootmax) condh = BlasLike.lm_max;
                if (condh < BlasLike.lm_rootmax) condh *= condh;
                dqpprt(orthog, isdel, iter, jadd, jdel, nactiv, ncolz, nfree,
                    n, nclin, nhess,
                    alfa, condh, condt, objqp, gfnorm,
                    ztgnrm, emax, AP);
                jadd = 0;
                jdel = 0;
                negligible = (ztgnrm <= dinky) ? 1 : 0;
                if (posdef)
                {
                    if (negligible != 0 || (uncon && (ztgnrm <= Math.Sqrt(dinky) || refine)))
                    {
                        /*
                        the projected gradient is negligible and the projected hessian
                        is positive definite.  if  r  is not complete it must be
                        expanded.  otherwise, if the current point is not optimal,
                        a constraint must be deleted from the working set
                        */
                        unitpg = negligible != 0;
                        alfa = 0;
                        uncon = false;
                        refine = false;
                        jdel = -(ncolr + 1);
                        if (ncolr >= ncolz)
                        {
                            dgetlamd(lprob, n, nactiv, ncolz, nfree,
                                ref jsmlst, ref ksmlst, ref smllst);

                            /*
                            test for convergence.  if the least (adjusted) multiplier is
                            greater than the small positive quantity  dinky,  an adequate
                            solution has been found
                            */
                            if (smllst > dinky)
                            {
                                //AddLog((char*)"smllst %20.8e dinky %20.8e\n",smllst,dinky);
                                /*OPTIMAL QP SOLUTION FOUND*/
                                i = 0;
                                break;
                            }

                            /*
                            ***** section ii.  delete a constraint from the working set *****
                            delete the constraint with the least (adjusted) multiplier
                            first check if there are any tiny multipliers
                            */
                            if (smllst > -dinky) zerolm = true;
                            jdel = jsmlst;
                            jdsave = jsmlst;
                            kdel = ksmlst;
                            isdel = ISTATE[jdel - 1];
                            issave = isdel;
                            ISTATE[jdel - 1] = 0;

                            /*update the TQ factorization of the matrix of constraints in the working set*/


                            ddelcon(modfyg != 0, orthog != 0, jdel, kdel,
                                nactiv, ncolz, nfree, n,
                                 NROWRT);
                            ++ncolz;
                            if (jdel <= (int)n) ++(nfree);
                            else --nactiv;
                        }

                        /*
                        the projected hessian is expanded by a row and column. compute
                        the elements of the new column of the cholesky factor r
                        use the vector p as temporary work space
                        */
                        renewr = true;
                        ++ncolr;
                        dqpcolr(ref nocurv, ref posdef, ref renewr, n, ref ncolr,
                            ref nfree, ref nrowh, ref ncolh, ref nhess, KFREE, ref
                            cslast, ref snlast, ref drmax, ref emax, ref hsize, ref rdlast, PX);
                        /*REPEAT THE MAIN LOOP*/
                        continue;
                    }
                    else if ((refine = uncon)) unitpg = false;
                }
                /*
                ******* section iii.  compute the search direction *******
                first, check for a weak local minimum. exit if the norm of the
                projected gradient is small and the curvature along p is not
                significant.  also, check for too many iterations and update the
                iteration count.  the iteration counter is only updated when a
                search direction is computed
                */
                if ((negligible != 0 && ncolr == ncolz && nocurv) || (zerolm && nocurv))
                {
                    //	AddLog("%d\tobjsiz ztgnrm gfnorm dinky%-.10e %-.10e %-.10e %-.10e\n",
                    //		*iter,objsiz,ztgnrm,gfnorm,dinky);
                    //	AddLog("\tnocurv %ld\n",nocurv);
                    /*WEAK LOCAL MINIMUM*/
                    i = 1;
                    break;
                }
                if (iter >= itmax)
                {
                    /*too many iterations*/
                    i = 5;
                    break;
                }

                ++iter;
                if (iter >= istart) msg = msglvl;
                //ColourConsole.WriteInfo($"iter {iter} objective {objqp:E16}");
                dfindp(nullr, unitpg, n, nclin,
                    NROWRT, ncolr, ncolz, ref nfree,
                    negligible != 0, ref gtp, ref pnorm, ref rdlast);

                /*
                if a constraint has just been deleted and the projected gradient
                is small (this can only occur here when the projected hessian is
                indefinite), the sign of  p  may be incorrect because of rounding
                errors in the computation of  ztg.  fix the sign of  p	by
                forcing it to satisfy the constraint that was just deleted
                */
                if ((jdsave > 0 && negligible != 0) || zerolm) dqpchkp(n, nclin, issave, jdsave);

                /*
                find the constraint we bump into along p
                update x and a*x if the step alfa is nonzero
                alfhit	is initialized to  bigalf.  if it remains that way after
                the call to bndalf, it will be regarded as infinite
                */
                bigalf = BlasLike.dprotdiv(ref bigdx, ref pnorm, ref ifail);
                if (ifail != 0 && bigdx == 0) bigalf = BlasLike.lm_max;
                dbndalf(firstv != 0, ref hitlow, ref jadd, n,
    nctotl, numinf, ref alfhit, ref palfa, ref atphit, ref bigalf,
     pnorm);


                /*
                if the projected hessian is positive definite, the step
                alfa = 1.0 will be the step to the minimum of the quadratic
                function on the current subspace
                */
                alfa = 1;

                /*
                if the step to the minimum on the subspace is less than the
                distance to the nearest constraint,  the constraint is not added
                to the working set
                */
                uncon = palfa > 1 && posdef;
                if (uncon) jadd = 0;
                else alfa = alfhit;
                /*check for an unbounded solution*/
                if (alfa >= bigalf)
                {
                    /*unbounded QP*/
                    i = 2;
                    break;
                }
                /*test if the change in	 x  is negligible*/
                stall = Math.Abs(alfa * pnorm) <= epspt9 * xnorm[0];
                if (stall)
                {
                    /*
                    take a zero step
                    exit if more than  50  iterations occur without changing  x
                    if such an exit is made when there are some near-zero
                    multipliers, the user should call a separate routine that
                    checks the solution
                    */
                    alfa = 0;
                    ++nstall;
                    if (nstall > mstall)
                    {
                        /*too many iterations without changing x*/
                        i = zerolm ? 3 : 4;
                        break;
                    }
                }
                else
                {
                    /*
                    compute the new value of the qp objective function.  if its
                    value has not increased,  update  objqp,  q(free)(t)g(free)  and
                    g(fixed). an increase in the objective can occur only after a
                    move along a direction of negative curvature from a point with
                    tiny multipliers. use the array	 rlam  as temporary storage
                    */

                    if (dqpgrad(2, n, nactiv, nfree, ref nhess_conv, nq,
                        nrowh, ncolh, jadd, alfa, ref objqp, ref
                        gfixed, gtp,
                           RLAM) < -BlasLike.lm_eps)
                    {
                        nhess = nhess_conv;
                        /*weak local minimum*/
                        //		AddLog("PLACE 2\n");
                        i = 1;
                        break;
                    }

                    nhess = nhess_conv;
                    /*
                    change	x  to  x + alfa*p.  update  ax	also
                    we no longer need to remember jdsave, the last constraint deleted
                    */
                    nstall = 0;
                    jdsave = 0;
                    zerolm = false;
                    BlasLike.daxpy(n, alfa, PX, 1, W, 1);
                    if (nclin > 0) BlasLike.daxpy(nclin, alfa, AP, 1, LWRK, 1);
                    xnorm[0] = dnrm2vec(n, W);
                }

                /*if an unconstrained step was taken, repeat the main loop*/
                if (uncon) continue;

                /*add a constraint to the working set. update  istate*/
                if (hitlow != 0) ISTATE[jadd - 1] = 1;
                else ISTATE[jadd - 1] = 2;
                if (L[jadd - 1] == U[jadd - 1]) ISTATE[jadd - 1] = 3;

                /*
                if a bound is to be added, move x exactly onto it, except when
                a negative step was taken.  (bndalf  may have had to move to some
                other closer constraint.)
                */
                iadd = jadd - n;
                if (jadd <= n)
                {
                    if (hitlow != 0) bnd = L[jadd - 1];
                    else bnd = U[jadd - 1];
                    if (alfa >= 0) W[jadd - 1] = bnd;
                    i = nfree;
                    for (ifix = 1; ifix <= (int)i; ++ifix) if (KFREE[ifix - 1] == jadd) break;
                }
                /*
                update the TQ factors of the matrix of constraints in the
                working set.  use the array  p	as temporary work space
                */

                daddcon(modfyg != 0, modfyr != 0, orthog != 0, ifix, iadd, jadd,
                    nactiv, ncolr, ncolz, nfree, n, NROWRT,
                    KFREE, condmx, cslast, snlast,
                     WRK, PX);
                --ncolr;
                --ncolz;
                nfixed = n - nfree;
                if (nfixed != 0)
                {
                    kb = nactiv + nfixed;
                    i = nfixed;
                    for (idummy = 1; idummy <= nfixed; ++idummy)
                    {
                        KACTV[kb] = KACTV[kb - 1];
                        --kb;
                    }
                }
                if (jadd <= n)
                {
                    /*
                    add a bound.  if stabilized eliminations are being used to update
                    the  TQ	 factorization,	 recompute the component of the gradient
                    corresponding to the newly fixed variable
                    use the array  p  as temporary work space
                    */
                    --nfree;
                    KACTV[nactiv] = jadd;
                    if (orthog == 0)
                    {

                        dqpgrad(3, n, nactiv, nfree, ref nhess_conv, nq,
                            nrowh, ncolh, jadd, alfa, ref objqp,
                            ref QTG[nfree], gtp, PX);
                        nhess = nhess_conv;
                    }
                }
                else
                {
                    /*add a general linear constraint*/
                    ++nactiv;
                    KACTV[nactiv - 1] = iadd;
                }

                /*
                repeat the main loop if the projected hessian that was used to
                compute this search direction was positive definite
                */
                if (ncolr == 0) posdef = true;
                if (ncolr == 0) emax = 0;

                if (!posdef)
                {
                    /*
                    the projected hessian was not sufficiently positive definite
                    before the constraint was added.  either compute the true value
                    of the last diagonal of	 r  or	recompute the whole of its last
                    column. use the array  rlam  as temporary work space
                    */
                    dqpcolr(ref nocurv, ref posdef, ref renewr, n, ref ncolr,
                        ref nfree, ref nrowh, ref ncolh, ref nhess, KFREE, ref cslast, ref snlast, ref drmax, ref emax, ref hsize, ref rdlast,
                        RLAM);
                }
                /*.........................END OF MAIN LOOP............................*/
            }


            inform = i;
            /*PRINT FULL SOLUTION*/
            msg = msglvl;
            if (msg >= 1) wrexit(lprob, i, iter);
            if (i > 0) dgetlamd(lprob, n, nactiv, ncolz, nfree,
                      ref jsmlst, ref ksmlst, ref smllst);
            dprtsol(nfree, n, nclin, ncnln, nctotl,
                nactiv);
        }
        int dqpcrsh(int n, int ncolz, int nfree, ref int nhess, int Nq, int nrowh, int ncolh, int Nrowrt, ref double hsize, double[] scale)
        {

            /*
                 dqpcrsh  computes the cholesky factor  r  of the projected hessian
                 z(t) h z,  given  z  and its dimensions  nfree by ncolz
                 if the projected hessian is indefinite, a smaller cholesky
                 factorization  r1(t) r1 = z1(t) h z1  is returned, where  z1  is
                 composed of  ncolr  columns of  z.  column interchanges are
                 used to maximize  ncolr.  these are applied to  z
            */
            int zy_offset;
            double dmin_, dmax_, d, t;
            int i, j, k, kmax, jthcol, ncolr;


            if (ncolz == 0) return 0;
            //    --wrk;
            zy_offset = Nq + 1;
            //   zy -= zy_offset;
            //   --scale;
            // rt -= Nrowrt + 1;
            var rt_offset = Nrowrt + 1;
            //  --kfree;
            ncolr = 0;
            /*
            compute  z(t) h z  and store the upper-triangular symmetric part
            in the first  ncolz  columns of  rt
            */
            for (k = 1; k <= ncolz; ++k)
            {

                BlasLike.dzerovec(n, WRK);
                if (!UNITQ)
                {
                    /* expand the column of  z  into an  n-vector */
                    for (i = 1; i <= nfree; ++i)
                    {
                        j = KFREE[i - 1];
                        WRK[j - 1] = ZY[i + k * nq - zy_offset];
                    }
                    if (scldqp) Factorise.ddmxmulv(n, scale, 1, WRK, 1);
                    jthcol = 0;
                }
                else
                {
                    /*
                        only bounds are in the working set.  the  k-th column of  z is
                        just a column of the identity matrix
                    */
                    jthcol = KFREE[k - 1];
                    WRK[jthcol - 1] = 1;
                }
                /*set  rt(*,k)  =  top of   h * (column of  z)*/
                qphess(n, nrowh, ncolh, jthcol, Q, WRK, RLAM);
                ++nhess;
                if (UNITQ && scldqp) BlasLike.dscalvec(n, scale[jthcol - 1], RLAM);
                if (scldqp) Factorise.ddmxmulv(n, scale, 1, RLAM, 1);
                dzyprod(4, n, nfree, ncolz, nfree, Nq, KFREE, KFREE,
                    RLAM, WRK);
                //         BlasLike.dcopyvec(ncolz, hz1, &rt[k * Nrowrt + 1]);
                //              BlasLike.dcopy(ncolz, RLAM, 1, RT, 1, 0, (k-1) * Nrowrt);
                BlasLike.dcopyvec(ncolz, RLAM, RT, 0, (k - 1) * Nrowrt);
                /*update an estimate of the size of the projected hessian*/
                //t = Math.Abs(rt[k + k * Nrowrt]);
                t = Math.Abs(RT[k + (k - 1) * Nrowrt - 1]);
                if (t > hsize) hsize = t;
            }
            /*
                 form the cholesky factorization  r(t) r  =  z(t) h z  as far as
                 possible, using symmetric row and column interchanges
            */
            dmin_ = BlasLike.lm_eps * hsize;
            for (j = 1; j <= ncolz; ++j)
            {
                /*FIND THE MAXIMUM REMAINING DIAGONAL*/
                kmax = j;
                //   dmax_ = rt[j + j * Nrowrt];
                dmax_ = RT[j + (j - 1) * Nrowrt - 1];
                for (k = j; k <= ncolz; ++k)
                {
                    //d = rt[k + k * Nrowrt];
                    d = RT[k + (k - 1) * Nrowrt - 1];
                    if (dmax_ < d)
                    {
                        dmax_ = d;
                        kmax = k;
                    }
                }
                /*see if the diagonal is big enough*/
                if (dmax_ <= dmin_) break;
                ncolr = j;
                /*permute the columns of z*/
                if (kmax != j)
                {
                    if (!UNITQ)
                    {
                        BlasLike.dcopyvec(nfree, ZY, WRK, kmax * Nq + 1 - zy_offset);
                        BlasLike.dcopyvec(nfree, ZY, ZY, j * Nq + 1 - zy_offset, kmax * Nq + 1 - zy_offset);
                        BlasLike.dcopyvec(nfree, WRK, ZY, 0, j * Nq + 1 - zy_offset);
                    }
                    else
                    {
                        Ordering.Order.swap(ref KFREE[kmax - 1], ref KFREE[j - 1]);
                    }
                    /*interchange rows and columns of the projected hessian*/

                    BlasLike.dswapvec(j, RT, RT, 1 + kmax * Nrowrt - rt_offset, 1 + j * Nrowrt - rt_offset);
                    BlasLike.dswap(kmax - j + 1, RT, 1, RT, Nrowrt, j + kmax * Nrowrt - rt_offset, j + j * Nrowrt - rt_offset);
                    BlasLike.dswap(ncolz + 1 - kmax, RT, Nrowrt, RT, Nrowrt, kmax + kmax * Nrowrt - rt_offset, j + kmax * Nrowrt - rt_offset);

                    RT[kmax + kmax * Nrowrt - rt_offset] = RT[j + j * Nrowrt - rt_offset];
                }
                /*set the diagonal element of R*/
                d = Math.Sqrt(dmax_);
                RT[j + j * Nrowrt - rt_offset] = d;
                if (j == ncolz) continue;
                /*
                set the above-diagonal elements of the k-th row of  r,
                and update the elements of all remaining rows
                */
                i = j + 1;
                for (k = i; k <= ncolz; ++k)
                {
                    t = RT[j + k * Nrowrt - rt_offset] / d;
                    RT[j + k * Nrowrt - rt_offset] = t;
                    /*R(I,K)  =  R(I,K)  - T * R(J,I),   I = i, k. */
                    if (t != 0) BlasLike.daxpy(k - j, -t, RT, Nrowrt, RT, 1, j + i * Nrowrt - rt_offset, i + k * Nrowrt - rt_offset);
                }
            }
            if (ncolr != ncolz && msg >= 80)
                lm_wmsg("\n//QPCRSH//  INDEFINITE PROJECTED HESSIAN.\n//QPCRSH// NCOLR=%6ld      NCOLZ=%6ld",
                    ncolr, ncolz);
            return ncolr;
        }
        short dqpgrad(short mode, int n, int nactiv, int nfree,
               ref int nhess, int Nq, int nrowh, int ncolh, int jadd,
                double alfa, ref double objqp, ref double gfixed, double gtp, double[] wrk2)
        {

            /*
                qpgrad	computes or updates..
                (1)	objqp, the value of the quadratic objective function, and
                (2)	the vectors  q(free)(t)g(free)  and  g(fixed),  where  q(free)
                    is the orthogonal factor of the	 a(free)  and  a  is the matrix
                    of constraints in the working set.  these vectors are stored in
                    elements  1,2,...,nfree	 and  nfree+1,...,n,  respectively,  of
                    the array  qtg
                (3)	the component of the gradient vector corresponding to a bound
                    constraint that has just been added to the working set
                */
            double deltaf;
            int jthcol;

            //   --wrk2;
            // --wrk1;
            // --x;
            //      --scale;
            //       --qtg;
            //        --p;
            //          --cvec;
            //       --kfree;
            //         --kactiv;

            jthcol = 0;
            switch (mode)
            {
                case 1:
                    /*
                        MODE = 1  ---	COMPUTE THE OBJECTIVE FUNCTION AND GRADIENT FROM
                                SCRATCH.  ALLOW FOR A DIAGONAL SCALING OF  X
                    */
                    BlasLike.dcopyvec(n, W, WRK);
                    if (scldqp) Factorise.ddmxmulv(n, LWRK, 1, WRK, 1);
                    qphess(n, nrowh, ncolh, jthcol, Q, WRK, QTG);
                    objqp = 0.5 * BlasLike.ddotvec(n, QTG, WRK) + BlasLike.ddotvec(n, c, WRK);
                    BlasLike.daxpy(n, 1, c, 1, QTG, 1);
                    if (scldqp) Factorise.ddmxmulv(n, LWRK, 1, QTG, 1);
                    /*COMPUTE Q(FREE)(T)(G(FREE) & G(FIXED).  discard ELEMENTS OF G(FREE)*/
                    dzyprod(6, n, nactiv, 0, nfree, Nq, KACTV, KFREE, QTG, WRK);
                    break;
                case 2:
                    /*
                            IF THE QP OBJECTIVE FUNCTION IS REDUCED BY A POSITIVE
                            STEP  ALFA,  OR  ALFA  IS NEGATIVE, UPDATE  OBJF
                            Q(FREE)(T)G(FREE)  AND  G(FIXED)  CORRESPONDING TO
                            THE CHANGE,  X = X + ALFA P
                    */
                    qphess(n, nrowh, ncolh, jthcol, Q, PX, WRK);
                    if (scldqp) Factorise.ddmxmulv(n, LWRK, 1, WRK, 1);
                    /*UPDATE  OBJQP. */
                    deltaf = alfa * gtp + 0.5 * alfa * alfa * BlasLike.ddotvec(n, PX, WRK);
                    if (deltaf > BlasLike.lm_eps && alfa > BlasLike.lm_eps)
                        /*THE STEP  ALFA  DOES NOT DECREASE THE OBJECTIVE FUNCTION. */
                        //			AddLog((char*)"Gone up by %-.10e baby\n",deltaf);
                        return -1;
                    objqp += deltaf;
                    /*UPDATE QTG.  USE P  AS TEMPORARY WORK SPACE*/
                    dzyprod(6, n, nactiv, 0, nfree, Nq, KACTV, KFREE,
                    WRK, wrk2);
                    BlasLike.daxpy(n, alfa, WRK, 1, QTG, 1);
                    break;
                case 3:
                    /*COMPUTE THE  JADD-TH COMPONENT OF THE GRADIENT VECTOR.*/
                    jthcol = jadd;
                    BlasLike.dzerovec(n, wrk2);
                    wrk2[jthcol] = 1;
                    qphess(n, nrowh, ncolh, jthcol, Q, wrk2, WRK);
                    if (scldqp)
                    {
                        Factorise.ddmxmulv(n, LWRK, 1, WRK, 1);
                        gfixed = LWRK[jadd - 1] * (BlasLike.ddotvec(n, WRK, W) + c[jadd - 1]);
                    }
                    else gfixed = BlasLike.ddotvec(n, WRK, W) + c[jadd - 1];
                    break;
            }
            ++nhess;
            return 0;
        }
        void dqpprt(int orthog, int isdel, int iter, int jadd, int jdel, int nactiv, int ncolz, int nfree, int n, int nclin, int nhess, double alfa, double condh, double condt, double obj, double gfnorm, double ztgnrm, double emax, double[] wrk2)
        {

            /*
                 dqpprt  prints various levels of output for  qpcore
                       msg    cumulative result
                       ---    -----------------
                    <=   0    no output
                    ==   1    nothing now (but full output later)
                    ==   5    one terse line of output
                    >=  10    same as 5 (but full output later)
                    >=  15    nothing more if  iter < istart
                              otherwise,  x,  istate  and  kfree
                    >=  20    multipliers (printed outside qpprt)
                              the array  ax
                    >=  30    diagonals of  t  and  r
                    >=  80    debug output
                    ==  99    cvec  and  hess  (called from qpdump)
            */
            char[] lstate =/*" LUETV";*/{ ' ', 'L', 'U', 'E', 'T', 'V' };
            string hdrfmt = "\n\n  ITN JDEL  JADD       STEP NHESS   OBJECTIVE NCOLZ %s NORM ZTG   COND T COND ZHZ  HESS MOD";
            string hdrfmt2 = "%5ld%5ld%c%5ld%c%10.2lg%6ld%12.4lg%6ld%11.2lg%10.2lg%9.1lg%9.1lg%10.2lg";

            var ladd = new char[1];
            var ldel = new char[1];
            int k, laddi, ldeli;
            //    --x;
            //--kfree;
            // --istate;

            if (msg < 5) return;
            ldeli = 0;
            laddi = 0;
            if (jdel > 0) ldeli = isdel;
            if (jdel < 0)
            {
                ldeli = 5;
                jdel = -jdel;
            }
            if (jadd > 0) laddi = ISTATE[jadd - 1];
            ldel[0] = lstate[ldeli];
            ladd[0] = lstate[laddi];
            if (msg < 15)
            {
                /*print heading (possibly) and terse line*/
                if (iter == 0 && jdel == 0)
                    lm_wmsg(hdrfmt, (orthog != 0 ? "NORM GFREE" : "  NORM QTG"));
                lm_wmsg(hdrfmt2, iter, jdel, ldel[0], jadd,
                    ladd[0], alfa, nhess, obj,
                    ncolz, gfnorm, ztgnrm, condt, condh, emax);
                return;
            }
            /*print terse line,  x,  istate,  kfree*/
            itrwri("QP", iter);
            lm_wmsg(hdrfmt, (orthog != 0 ? "NORM GFREE" : "  NORM QTG"));
            lm_wmsg(hdrfmt2, iter, jdel, ldel[0], jadd,
                ladd[0], alfa, nhess, obj, ncolz, gfnorm, ztgnrm, condt, condh, emax);
            lm_mdvwri("\nQP VARIABLES", n, W);
            lm_mdvwri("\nSTATUS OF THE QP BOUND CONSTRAINTS", n, ISTATE);
            if (nclin > 0) lm_mdvwri("\nSTATUS OF THE QP GENERAL CONSTRAINTS", nclin, ISTATE, n);
            if (nfree > 0) lm_mdvwri("\nLIST OF FREE QP VARIABLES", nfree, KFREE);
            /*compute and print  ax.  use  work  to avoid side effects*/
            if (msg < 20) return;
            if (nclin > 0)
            {
                for (k = 0; k < (int)nclin; ++k) wrk2[k] = BlasLike.ddot(n, A, NROWA, W, 1, k);
                lm_mdvwri("\nVALUES OF QP GENERAL LINEAR CONSTRAINTS", nclin, wrk2);
            }
            /*print all the diagonals of  t  and  r*/
            if (msg < 30) return;
            if (nactiv > 0)
            {
                BlasLike.dcopy(nactiv, RT, NROWRT - 1, WRK, 1, nactiv + ncolz * NROWRT - 1);
                lm_mdvwri("\nDIAGONALS OF QP WORKING SET FACTOR  T", nactiv, WRK);
            }
            if (ncolz > 0)
            {
                ColourConsole.WriteInfo("\nDIAGONALS OF QP PRJ. HESSIAN FACTOR  R");
                lm_gdvwri(ncolz, RT, NROWRT + 1);
            }
        }
        void dqpcolr(ref bool nocurv, ref bool posdef, ref bool renewr, int n, ref int ncolr, ref int nfree, ref int nrowh, ref int ncolh, ref int nhess, int[] kfree, ref double cslast, ref double snlast, ref double drmax, ref double emax, ref double hsize, ref double rdlast, double[] hz1)
        {
            /*
                dqpcolr  is used to compute elements of the  (ncolr)-th  column of
                r,  the cholesky factor of the projected hessian.  if  renewr  is
                true  on entry,  the complete column is to be computed
                otherwise, only the last diagonal element is required
                if the resulting projected hessian is singular or indefinite, its
                last diagonal element is increased by an amount  emax  that
                ensures positive definiteness.  this diagonal modification will
                alter the scale of the qp search vector  p, but not its
                direction

                on exit,  qpcolr will have stored the  ncolr  elements of the new
                column of  r  in the array  rt,  and set the variables  nocurv,
                posdef,  renewr,  drmax,  emax  and  hsize
            */
            int hess_dim1, hess_offset, rt_dim1, rt_offset, zy_dim1, zy_offset,
                i__1;


            double rdsq, zthz1;
            int j, k;
            short idiag;
            double s;
            double rnorm;
            int ncolr1;
            int jthcol;
            double rdsmin, rdsmax;
            //          --wrk;
            //  --hz1;
            zy_dim1 = nq;
            zy_offset = zy_dim1 + 1;
            //    zy -= zy_offset;
            // --scale;
            rt_dim1 = NROWRT;
            rt_offset = rt_dim1 + 1;
            //    rt -= rt_offset;
            hess_dim1 = nrowh;
            hess_offset = hess_dim1 + 1;
            //   hess -= hess_offset;
            //   --kfree;

            if (renewr)
            {
                goto L20;
            }
            /*

            only the last element of the new column of  r  need be computed
            this situation can only occur when a constraint is added to the
            working set with  zthz1  not positive definite

            the last diagonal element of  r  is that of  zthz1  plus a diagonal
            modification.  the square of the true diagonal is recovered from
            the rotations used to update  r  when the constraint was added to
            the working set
            */
            rdlast = RT[ncolr + ncolr * rt_dim1 - rt_offset];
            s = Math.Abs(snlast);
            rdsq = (cslast - s) * rdlast * ((cslast + s) * rdlast);
            goto L120;
        /*
        the projected hessian is expanded by a row and column.  compute
        the first  (ncolr - 1)  elements of the new column of the
        cholesky factor r.  also, compute  rdsq,  the square of the last
        diagonal element
        */
        L20:
            BlasLike.dzerovec(n, WRK);

            if (!UNITQ)
            {
                /*
                     expand the new column of z in to an n-vector
                */
                i__1 = nfree;
                for (k = 1; k <= i__1; ++k)
                {
                    j = kfree[k - 1];
                    WRK[j - 1] = ZY[k + ncolr * zy_dim1 - zy_offset];
                    /* L40: */
                }
                if (scldqp) Factorise.ddmxmulv(n, LWRK, 1, WRK, 1);
                jthcol = 0;
            }
            else
            {
                /*
                     only bounds are in the working set (nfree  is equal to  ncolz)
                     the (ncolr)-th  column of  z  is just a column of the identity
                     matrix
                */

                jthcol = kfree[ncolr - 1];
                WRK[jthcol - 1] = 1;
            }
            /*
                 compute the hessian times the last column of z
            */

            qphess(n, nrowh, ncolh, jthcol, Q, WRK, hz1);
            ++nhess;
            if (UNITQ && scldqp) BlasLike.dscalvec(n, LWRK[jthcol - 1], hz1);
            if (scldqp) Factorise.ddmxmulv(n, LWRK, 1, hz1, 1);
            /*
                 compute the  (ncolr)-th  column of  z(t)h z            */
            dzyprod(4, n, nfree, ncolr, nfree, nq, KFREE, KFREE, hz1, WRK);
            BlasLike.dcopyvec(ncolr, hz1, RT, 0, ncolr * rt_dim1 + 1 - rt_offset);
            /*
            compute the first  (ncolr - 1)  elements of the last column of  r.
            */
            ncolr1 = ncolr - 1;
            zthz1 = RT[ncolr + ncolr * rt_dim1 - rt_offset];
            rdsq = zthz1;
            if (ncolr1 != 0)
            {
                idiag = dtmxsolve((short)-1, ncolr1, RT, NROWRT, RT, (short)1, 0, ncolr * rt_dim1 + 1 - rt_offset);
                //  rnorm = dnrm2vec(ncolr1, &pRT[*ncolr * rt_dim1 + 1 - rt_offset]);
                rnorm = dnrm2vec(ncolr1, RT, ncolr * rt_dim1 + 1 - rt_offset);

                /*
                compute the square of the last diagonal element of  r
                */
                rdsq = zthz1 - rnorm * rnorm;
            }
            /*
            update the estimate of the norm of the hessian
            */

            hsize = Math.Max(hsize, zthz1);
        /*
        compute  rdlast, the last diagonal of  r.  the variables posdef
        and  nocurv  are set here.  they are used to indicate if the new
        projected hessian is positive definite or singular.  if  posdef
        is set to false,  rdlast  will be that of  zthz1  plus a diagonal
        modification. if the required diagonal modification is large,
        renewr  will be set to be  true,  indicating that the last row
        and column of  r  must be recomputed when a constraint is added
        to the working set during the next iteration
        */
        L120:
            nocurv = false;
            renewr = false;
            emax = 0;
            /*
            rdsmin  is the square of the smallest allowable diagonal element
            for a positive-definite cholesky factor.  note that the test for a
            singular matrix is scale dependent
            */
            rdsmin = (ncolr > 1 ? drmax * drmax : hsize) * BlasLike.lm_eps;
            posdef = rdsq > rdsmin;
            if (posdef)
            {
                goto L180;
            }
            if (rdsq < -rdsmin)
            {
                goto L140;
            }
            /*
                 the projected hessian is singular

            the quadratic has no curvature along at least one direction.  the
            perturbation  emax  is chosen to make the new eigenvalue of  zthz1
            small and positive
            */
            emax = rdsmin - rdsq;
            rdsq = rdsmin;
            nocurv = true;
            goto L180;
        /*
        the projected hessian is indefinite.  there are two cases

        case 1.  the modulus of the new last diagonal of  r  is not too
        large.  the modulus of  rdsq  is used for the square root
        */
        L140:
            rdsmax = 10 * hsize;
            if (rdsq < -rdsmax)
            {
                goto L160;
            }
            emax = -2 * rdsq;
            goto L180;
        /*
        case 2.  the modulus of the last diagonal of  r  is judged to be
        too large (some loss of precision may have occurred).  set
        renewr  so that the last column is recomputed later
        */
        L160:
            emax = rdsmax - rdsq;
            renewr = true;
            rdsq = rdsmax;
        /*
        compute the last diagonal element
        */
        L180:
            rdlast = Math.Sqrt((Math.Abs(rdsq)));
            RT[ncolr + ncolr * rt_dim1 - rt_offset] = rdlast;
            if (msg >= 80 && !posdef)
                lm_wmsg(
            "\n//QPCOLR//  POSDEF NOCURV          EMAX        RDLAST\n//QPCOLR//%8ld%7ld%14.4lg%14.4lg",
                posdef, nocurv, emax, rdlast);
            return;
        }
        void dqpchkp(int n, int nclin, int issave, int jdsave)
        {/*
	dqpchkp  is called when a constraint has just been deleted and the
	sign of the search direction  p  may be incorrect because of
	rounding errors in the computation of the projected gradient ztg.
	the sign of the search direction (and therefore the product  ap)
	is fixed by forcing p to satisfy the constraint (with index jdsave)
	that was just deleted.  variables that were held temporarily fixed
	(with istate = 4)  are not checked for feasibility
*/
            double atp = (jdsave <= n ? PX[jdsave - 1] : AP[jdsave - n - 1]);

            if (issave == 4) return;

            if (msg >= 80) lm_wmsg("\n//QPCHKP//  JDSAVE ISSAVE            ATP\n//QPCHKP//%8ld%7ld%15.5lg",
                    jdsave, issave, atp);

            if ((issave == 2 && atp <= 0) || (issave == 1 && atp >= 0)) return;

            /*REVERSE THE DIRECTION OF  P  AND  AP. */
            BlasLike.dnegvec(n, PX);
            if (nclin > 0) BlasLike.dnegvec(nclin, AP);
        }
        public static void printV<T>(string name, T[] a, int upto = -1, int astart = 0)
        {
            ColourConsole.WriteInfo(name);
            if (upto == -1) upto = a.Length;
            for (int i = astart; i < upto; ++i)
            {
                var p = a[i].GetType();
                if (p.FullName == "System.Double")
                    ColourConsole.Write($"{a[i],11:F8} ", ConsoleColor.Magenta);
                else
                    ColourConsole.Write($"{a[i]} ", ConsoleColor.DarkGreen);
                if ((i - astart) % 10 == 9) Console.Write("\n");
            }
            Console.Write("\n");
        }
        public static short LPopt(int n, int m, double[] ww, double[] LL, double[] UU, double[] AA, double[] cc, ref double objective, ref int iter, double[] LAMBDAback = null)
        {
            var opt = new Optimise();
            var lp = 1;
            var itmax = 20000;
            var orthog = 1;
            var nclin = m;
            var nctotl = n + m;
            var nrowa = m;
            double obj = 1e10;
            var featolv = small_round(1e-8);
            int cold = 1;
            short msglvl = -1000;
            opt.ISTATE = new int[n + m + n + n];
            var lwrk = 2 * (n * (n + 2) + m) + m;
            opt.LAMBDA = new double[n + m];
            opt.FEATOL = new double[n + m];
            opt.LWRK = new double[lwrk];
            opt.A = AA;
            opt.L = LL;
            opt.U = UU;
            opt.c = cc;
            opt.W = ww;
            BlasLike.dsetvec(n + m, 0, opt.LAMBDA);
            BlasLike.dsetvec(n + m, featolv, opt.FEATOL);
            short ifail = 89;
            short back;
            opt.clocker(true);
            back = opt.dqpsol(itmax, msglvl, n, m, n + m, m,
      n + n, 1, cold, lp, orthog, ref
       iter, ref obj, n + n, lwrk, ifail);
            var tt = opt.clocker();
            ColourConsole.WriteLine($"Time elapsed {tt} m secs", ConsoleColor.DarkCyan);
            objective = obj;
            if (LAMBDAback != null) BlasLike.dcopyvec(LAMBDAback.Length, opt.LAMBDA, LAMBDAback);
            return back;
        }
        public short QPopt(int n, int m, double[] ww, double[] LL, double[] UU, double[] AA, double[] cc, double[] QQ, ref double objective, ref int iter, int lp = 0, double[] LAMBDAback = null)
        {
            var opt = this;
            if (opt.h == null) opt.h = opt.qphess1;
            var itmax = (short)20000;
            var orthog = 1;
            var nclin = m;
            var nctotl = n + m;
            var nrowa = m;
            double obj = 1e10;
            var featolv = small_round(1e-8);
            int cold = 1;
            short msglvl = -1000;
            opt.ISTATE = new int[n + m + n + n];
            var lwrk = 2 * (n * (n + 2) + m) + m;
            opt.LAMBDA = new double[n + m];
            opt.FEATOL = new double[n + m];
            opt.LWRK = new double[lwrk];
            BlasLike.dsetvec(n + m, 0, opt.LAMBDA);
            BlasLike.dsetvec(n + m, featolv, opt.FEATOL);
            short ifail = 89;
            short back;
            opt.A = AA;
            opt.L = LL;
            opt.U = UU;
            opt.c = cc;
            opt.W = ww;
            opt.Q = QQ;
            opt.clocker(true);
            back = opt.dqpsol(itmax, msglvl, n, m, n + m, m,
      n + n, 1, cold, lp, orthog, ref
       iter, ref obj, n + n, lwrk, ifail);
            objective = obj;
            var tt = opt.clocker();
            ColourConsole.WriteLine($"Time elapsed {tt} m secs", ConsoleColor.DarkCyan);
            if (LAMBDAback != null) BlasLike.dcopyvec(LAMBDAback.Length, opt.LAMBDA, LAMBDAback);
            return back;
        }
        int badboundcount()
        {
            var back = 0;
            for (int i = 0; i < L.Length; ++i) back += ((L[i] <= U[i]) ? 0 : 1);
            return back;
        }
    }
}
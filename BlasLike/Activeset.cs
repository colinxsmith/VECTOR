using System;
using Blas;
using Solver;
using System.Diagnostics;
namespace ActiveSet
{
    public class Optimise
    {
        public static int DAS_Sol_nout = 0;
        public static int DAS_Sol_msg = 0;
        public static int DAS_Sol_istart = 0;
        public static int DAS_Sol_asize = 0;
        public static int DAS_Sol_dtmax = 0;
        public static int DAS_Sol_dtmin = 0;
        public static int DAS_Sol_nrowrt = 0;
        public static int DAS_Sol_ncolrt = 0;
        public static int DAS_Sol_nq = 0;
        public static int DAS_Sol_ncqp = 0;
        public static int DAS_Sol_nrowqp = 0;
        public static int DAS_Sol_scldqp = 0;
        public static double DAS_Sol_zgfacc = -1.0;
        public static unsafe void** DAS_Sol_ploc;
        public static int msg;
        public static bool scldqp;
        public unsafe static int nrowqp;
        public unsafe static int* pKACTV;
        public unsafe static int* pKFREE;
        public unsafe static int[] loclc;
        public unsafe static int[] locnp;
        public static int istart;
        public static int[] parm;
        public static int nq;
        public static double asize;
        public static int nrowrt;
        public static int ncolrt;
        public static double d__alfa;
        public static double pnorm1;
        public static double dtmax;
        public static double dtmin;
        public unsafe static double* pANORM;
        public unsafe static double* pWRK;
        public unsafe static double* pRT;
        public unsafe static double* pZY;
        public unsafe static double* pPX;
        public unsafe static double* pQTG;
        public unsafe static double* pRLAM;
        public unsafe static double* pAP;
        public static double CR(double x) { return x; }
        public static T CS<T>(T x) { return x; }
        public static int CN(int x) { return x; }
        public static int CI(int x) { return x; }
        public static byte CB(byte x) { return x; }
        public static bool CB(bool x) { return x; }
        public static T CL<T>(T x) { return x; }
        public static long timebase;
        public static double timeaquired;
        public static void wrexit(string name, int inform, int iter)
        {
            Console.WriteLine($"{name} Inform={inform} Iter={iter}");
        }
        public static void wdinky(string name, double ztgnrm, double dinky)
        {
            Console.WriteLine($"\n//{name}//         ZTGNRM         DINKY\n//{name}//{ztgnrm}{dinky}");
            //lm_wmsg("\n//%s//         ZTGNRM         DINKY\n//%s//%14.5lg%14.5lg",
            //name,name,ztgnrm,dinky);
        }
        public static double clocker(bool start = false)
        {
            long t = DateTimeOffset.UtcNow.ToUnixTimeSeconds();
            if (start) timeaquired = 0;
            else timeaquired += (double)(t - timebase);
            timebase = t;
            return timeaquired;
        }
        public unsafe static void dlpcore(bool lp, int minsum, bool orthog, int* unitq, int vertex, int* inform, int* iter,
                int itmax, byte lcrash, int n, int* nclin, int* nctotl, int* nrowa, int* nactiv,
                int* nfree, int* numinf, int* istate, int* kactiv, int* kfree, double* obj, double* xnorm,
                double* a, double* ax, double* bl, double* bu, double* clamda, double* cvec, double* featol, double* x, int* iw, double* w)
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
            int ifix;
            bool prnt;
            bool added;
            double palfa;
            int ifail;
            double bigdx;
            int isdel = 0;
            double objlp, condt;
            double anorm, dinky;
            int ncnln;
            bool stall;
            double wgfix = 1e9;
            int ncolz;
            double pnorm;
            bool nullr;
            //    int nrowj;
            int nclin0, kb;
            double bigalf, bigbnd, epspt9;
            int was_is;
            double feamax, feamin;
            bool delete_;
            int nfixed;
            int jbigst, kbigst;
            double tolact;
            bool modfyg;
            double condmx, atphit, cslast, rdlast, objsiz, snlast, suminf,
                 trulam;
            int idummy, msglvl;
            int jsmlst, ksmlst;
            double smllst;
            int mstall;
            double ztgnrm;
            int nstall;
            bool firstv;
            int hitlow;
            bool unitpg;
            double bnd;
            double gtp;

            --w;
            --iw;
            --x;
            --featol;
            --cvec;
            --clamda;
            --bu;
            --bl;
            --ax;
            a_dim1 = *nrowa;
            a_offset = a_dim1 + 1;
            a -= a_offset;
            --kfree;
            --kactiv;
            --istate;
            /*     SPECIFY MACHINE-DEPENDENT PARAMETERS. */
            /*INITIALIZE */
            ncnln = 0;
            nclin0 = Math.Max(*nclin, 1);
            //    nrowj = 1;
            *inform = 0;
            *iter = 0;
            jadd = 0;
            jdel = 0;
            ndel = 0;
            nstall = 0;
            *numinf = 1;
            msglvl = msg;
            msg = 0;
            if (*iter >= istart)
            {
                msg = msglvl;
            }
            bigbnd = parm[0];
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
            BlasLike.dxminmax(*nctotl, &featol[1], 1, &feamax, &feamin);
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
            var nq_c = nq;
            var nrowrt_c = nrowrt;
            var ncolrt_c = ncolrt;
            dlpcrsh(orthog, unitq, vertex, lcrash, n, nclin, nctotl, &nq_c,
                nrowa, &nrowrt_c, &ncolrt_c, nactiv, &ncolz, nfree, &istate[1], &
                kactiv[1], &kfree[1], &bigbnd, &tolact, xnorm, &a[a_offset],
                pANORM, &ax[1], &bl[1], &bu[1], &x[1], pQTG, pRT, pZY,
                pPX, pWRK, pRLAM);
            nq = nq_c;
            nrowrt = nrowrt_c;
            ncolrt = ncolrt_c;
            dlpgrad(lp, n, CN(*nctotl), CN(*nrowa), bigbnd, feamin, numinf, &suminf, &istate[
                1], &a[a_offset], &bl[1], &bu[1], &cvec[1], &featol[1], pQTG,
                &x[1]);
            dzyprod(6, n, *nactiv, ncolz, *nfree, nq, *unitq, &kactiv[1], &kfree[1]
                , pQTG, pZY, pWRK);
            *obj = suminf;
            if (lp) objlp = BlasLike.ddotvec(n, &cvec[1], &x[1]);
            if (lp && *numinf == 0) *obj = objlp;
            if (*numinf == 0 && !(lp)) goto L320;
            /* .......................START OF THE MAIN LOOP........................ 
            */
            /*     DEFINE SMALL QUANTITIES THAT REFLECT THE MAGNITUDE OF  C,  X, */
            /*     AND THE NORM OF THE CONSTRAINTS IN THE WORKING SET. */
            L20:
            clocker();
            objsiz = (1 + Math.Abs(*obj)) / (1 + *xnorm);
            if (*numinf == 0)
            {
                objsiz = (BlasLike.lm_eps + Math.Abs(*obj)) / (BlasLike.lm_eps + *xnorm);
            }
            anorm = 0;
            if (*nactiv > 0)
            {
                anorm = Math.Abs(dtmax);
            }
            dinky = epspt9 * Math.Max(anorm, objsiz);
            /*     COMPUTE THE NORMS OF THE PROJECTED GRADIENT AND THE GRADIENT WITH 
            */
            /*     RESPECT TO THE FREE VARIABLES. */
            ztgnrm = 0;
            if (ncolz > 0) ztgnrm = dnrm2vec(ncolz, pQTG);

            if (msg >= 80) wdinky("LPCORE", ztgnrm, dinky);
            delete_ = ztgnrm <= dinky;
            /*     PRINT THE DETAILS OF THIS ITERATION. */
            prnt = added || ndel > 1;
            if (!prnt) goto L40;
            var dtmax_c = dtmax;
            var dtmin_c = dtmin;

            condt = dprotdiv(&dtmax_c, &dtmin_c, &ifail);
            dtmax = dtmax_c;
            dtmin = dtmin_c;
            if (ifail != 0 && dtmax == 0) condt = BlasLike.lm_max;
            dlpprt(lp, CN(*nrowa), nrowrt, n, CN(*nclin), CN(*nfree), isdel, CN(*nactiv),
                CN(ncolz), CN(*iter), jadd, jdel, alfa, condt, CN(*numinf),
                suminf, objlp, &istate[1], &kfree[1], &a[a_offset], pRT, &x[1], pWRK, pAP);
            added = false;
            jadd = 0;
            jdel = 0;
        L40:
            if (*numinf == 0 && !lp)
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
            dgetlamd(lprob, n, *nactiv, ncolz, *nfree, *nrowa, nrowrt,
                &jsmlst, &ksmlst, &smllst, &istate[1], &kactiv[1], &a[
                a_offset], pANORM, pQTG, pRLAM, pRT);
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
            isdel = istate[jdel];
            istate[jdel] = 0;
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
            if (*numinf == 0 || minsum != 0)
            {
                goto L280;
            }
            /*     FIND THE BIGGEST MULTIPLIER LARGER THAN UNITY. */
            /*     FOR THE PURPOSES OF THE TEST,  THE  J-TH  MULTIPLIER IS SCALED */
            /*     BY  FEATOL(J)/FEAMIN.  THIS FORCES CONSTRAINTS WITH LARGER  FEATOL 
            */
            /*     VALUES TO BE DELETED FIRST. */
            dlpbgst(n, CN(*nactiv), CN(*nfree), &jbigst, &kbigst, &istate[1], &kactiv[1], dinky, feamin,
                &trulam, &featol[1], pRLAM);
            if (jbigst == 0)
            {
                goto L280;
            }
            jdel = jbigst;
            kdel = kbigst;
            isdel = istate[jbigst];

            was_is = trulam > 0 ? -2 : -1;
            istate[jbigst] = was_is;
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
            ddelcon(modfyg, orthog, *unitq, CN(jdel), CN(kdel), CN(*nactiv), CN(ncolz), CN(*nfree), n,
                CN(nq), CN(*nrowa), CN(nrowrt), &kactiv[1], &kfree[1], &a[a_offset],
                pQTG, pRT, pZY);
            ++ncolz;
            if (jdel <= (int)n)
            {
                ++(*nfree);
            }
            if (jdel > (int)n)
            {
                --(*nactiv);
            }
            goto L20;
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE SEARCH DIRECTION,  P = - Z*(PROJECTED GRADIENT). */
        /* --------------------------------------------------------------------- 
        */
        L100:
            if (*iter >= itmax)
            {
                goto L400;
            }
            ++(*iter);
            if (*iter >= istart)
            {
                msg = msglvl;
            }
            dfindp(&nullr, &unitpg, unitq, n, nclin, &nq_c, nrowa, &
                nrowrt_c, &ncolz, &ncolz, nfree, &istate[1], &kfree[1],
                delete_, &gtp, &pnorm, &rdlast, &a[a_offset], pAP,
                pPX, pQTG, pRT, pWRK, pZY, pWRK);
            nq = nq_c;
            nrowrt = nrowrt_c;
            /* --------------------------------------------------------------------- 
            */
            /*     FIND THE CONSTRAINT WE BUMP INTO ALONG  P. */
            /*     UPDATE  X  AND  AX  IF THE STEP  ALFA  IS NONZERO. */
            /* --------------------------------------------------------------------- 
            */
            /*     ALFA  IS INITIALIZED TO  BIGALF.  IF IT REMAINS THAT WAY AFTER */
            /*     THE CALL TO BNDALF, IT WILL BE REGARDED AS INFINITE. */
            bigalf = dprotdiv(&bigdx, &pnorm, &ifail);
            if (ifail != 0 && bigdx == 0) bigalf = BlasLike.lm_max;
            *inform = dbndalf(firstv, &hitlow, &istate[1], &jadd, n,
                CN(*nctotl), CN(*numinf), &alfa, &palfa, &atphit, &bigalf, &bigbnd,
                &pnorm, pANORM, pAP, &ax[1], &bl[1], &bu[1], &featol[1], pPX, &x[1]);
            if (*inform != 0 || jadd == 0)
            {
                goto L300;
            }
            /*     TEST IF  ALFA*PNORM  IS NEGLIGIBLE. */
            stall = (Math.Abs(d__alfa * pnorm1)) <= epspt9 * *xnorm;
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
            BlasLike.daxpyvec(n, alfa, pPX, &x[1]);
            if (*nclin > 0)
                BlasLike.daxpy(*nclin, alfa, pAP, 1, &ax[1], 1);
            *xnorm = dnrm2vec(n, &x[1]);
            if (lp) objlp = BlasLike.ddotvec(n, &cvec[1], &x[1]);
            /*     IF  X  IS NOT YET FEASIBLE,  COMPUTE  OBJ  AND  GRAD  AS THE VALUE 
            */
            /*     AND GRADIENT OF THE SUM OF INFEASIBILITIES (IF  X  IS FEASIBLE, */
            /*     THE VECTOR  QTG  IS UPDATED AND  GRAD  NEED NOT BE COMPUTED). */
            L140:
            if (*numinf == 0) goto L160;
            dlpgrad(lp, n, CN(*nctotl), CN(*nrowa), bigbnd, feamin, numinf, &suminf, &istate[
                1], &a[a_offset], &bl[1], &bu[1], &cvec[1], &featol[1], pQTG,
                &x[1]);
            if (!orthog && jadd <= (int)n) wgfix = pQTG[jadd - 1];
            dzyprod(6, n, *nactiv, ncolz, *nfree, nq, *unitq, &kactiv[1], &kfree[1]
                , pQTG, pZY, pWRK);
            *obj = suminf;
        /* --------------------------------------------------------------------- 
        */
        /*     ADD A CONSTRAINT TO THE WORKING SET. */
        /* --------------------------------------------------------------------- 
        */
        /*     UPDATE  ISTATE. */
        L160:
            if (lp && *numinf == 0)
            {
                *obj = objlp;
            }
            if (hitlow != 0)
            {
                istate[jadd] = 1;
            }
            if (hitlow == 0)
            {
                istate[jadd] = 2;
            }
            if (bl[jadd] == bu[jadd])
            {
                istate[jadd] = 3;
            }
            /*     IF A BOUND IS TO BE ADDED, MOVE  X  EXACTLY ONTO IT, EXCEPT WHEN */

            /*     A NEGATIVE STEP WAS TAKEN.  (BNDALF  MAY HAVE HAD TO MOVE TO SOME 
            */
            /*     OTHER CLOSER CONSTRAINT.) */
            iadd = jadd - n;
            if (jadd > (int)n) goto L200;
            bnd = (hitlow != 0 ? bl : bu)[jadd];
            if (alfa >= 0) x[jadd] = bnd;
            for (ifix = 1; ifix <= *nfree; ++ifix)
            {
                if (kfree[ifix] == jadd)
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
            nq_c = nq;
            nrowrt_c = nrowrt;
            *inform = daddcon(modfyg, false, orthog, unitq, &ifix, &iadd, &jadd,
                nactiv, &ncolz, &ncolz, nfree, n, &nq_c, nrowa, &nrowrt_c, &
                kfree[1], &condmx, &cslast, &snlast, &a[a_offset], pQTG,
                pRT, pZY, pWRK, pPX);
            nq = nq_c;
            nrowrt = nrowrt_c;
            --ncolz;
            nfixed = n - *nfree;
            if (nfixed == 0)
            {
                goto L240;
            }
            kb = *nactiv + nfixed;
            for (idummy = 1; idummy <= nfixed; ++idummy)
            {
                kactiv[kb + 1] = kactiv[kb];
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
                --(*nfree);
                kactiv[*nactiv + 1] = jadd;
                if (!orthog)
                {
                    if (*numinf > 0) pQTG[*nfree] = wgfix;
                    else if (lp) pQTG[*nfree] = cvec[jadd];
                }
            }
            else
            {
                /*ADD A GENERAL LINEAR CONSTRAINT*/
                ++(*nactiv);
                kactiv[*nactiv] = iadd;
            }
            goto L20;
        /* .........................END OF MAIN LOOP............................ 
        */
        /*     NO CONSTRAINTS TO DROP. */
        L280:
            if (*numinf > 0)
            {
                goto L340;
            }
            goto L320;
        /*     ERROR IN  BNDALF  --  PROBABLY UNBOUNDED LP. */
        L300:
            if (*numinf == 0)
            {
                goto L360;
            }
            goto L340;
        /*     FEASIBLE SOLUTION FOUND, OR OPTIMAL LP SOLUTION. */
        L320:
            *inform = 0;
            goto L420;
        /*     THE LINEAR CONSTRAINTS AND BOUNDS APPEAR TO BE INFEASIBLE. */
        L340:
            *inform = 1;
            goto L420;
        /*     UNBOUNDED LP. */
        L360:
            *inform = 2;
            goto L420;
        /*     TOO MANY ITERATIONS WITHOUT CHANGING  X. */
        L380:
            *inform = 3;
            goto L420;
        /*     TOO MANY ITERATIONS. */
        L400:
            *inform = 4;
        /* --------------------------------------------------------------------- 
        */
        /*     PRINT FULL SOLUTION.  IF NECESSARY, RECOMPUTE THE MULTIPLIERS. */
        /* --------------------------------------------------------------------- 
        */
        L420:
            msg = msglvl;
            if (msg >= 1)
            {
                wrexit("LP", *inform, *iter);
            }
            if (*inform > 0) dgetlamd(lprob, n, *nactiv, ncolz, *nfree, *nrowa,
                nrowrt, &jsmlst, &ksmlst, &smllst, &istate[1], &
                kactiv[1], &a[a_offset], pANORM, pQTG, pRLAM,
                pRT);
            if (!lp && *inform == 0) BlasLike.dzerovec(n, pRLAM);
            dprtsol(*nfree, *nrowa, n, *nclin, ncnln, *nctotl, bigbnd,
                *nactiv, &istate[1], &kactiv[1], &a[a_offset],
                &bl[1], &bu[1], &x[1], &clamda[1], pRLAM, &x[1]);
        }
        public unsafe static double dnrm2vec(int n, double* x)
        {
            if (n == 1) return (x[0] < 0.0 ? -x[0] : x[0]);
            else
            {
                double scale = 0.0, ssq = 1.0;
                dsssqvec(n, x, &scale, &ssq);
                return sc_norm(scale, ssq);
            }
        }
        public static double sc_norm(double scale, double ssq)
        {
            ssq = Math.Sqrt(ssq);
            return (scale < BlasLike.lm_max / ssq ? scale * ssq : BlasLike.lm_max);
        }
        public static void dsssqvec(int n, double[] x, double[] pscale, double[] psumsq, int px = 0)
        {
            if (n > 0)
            {
                double absxi, d, sumsq = psumsq[0], scale = pscale[0];
                Debug.Assert(scale >= 0);
                for (int i = 0; i < n; ++i)
                {
                    absxi = x[i + px];
                    if (absxi == 0) continue;
                    if (absxi < 0) absxi = -absxi;
                    if (scale < absxi)
                    {
                        d = scale / absxi;
                        sumsq = sumsq * (d * d) + 1;
                        scale = absxi;
                    }
                    else
                    {
                        d = absxi / scale;
                        sumsq += d * d;
                    }
                }
                pscale[0] = scale;
                psumsq[0] = sumsq;
            }
        }

        public unsafe static void dsssqvec(int n, double* x, double* pscale, double* psumsq)
        {
            if (n > 0)
            {
                double absxi, d, sumsq = psumsq[0], scale = pscale[0];
                Debug.Assert(scale >= 0);
                for (int i = 0; i < n; ++i)
                {
                    absxi = x[i];
                    if (absxi == 0) continue;
                    if (absxi < 0) absxi = -absxi;
                    if (scale < absxi)
                    {
                        d = scale / absxi;
                        sumsq = sumsq * (d * d) + 1;
                        scale = absxi;
                    }
                    else
                    {
                        d = absxi / scale;
                        sumsq += d * d;
                    }
                }
                pscale[0] = scale;
                psumsq[0] = sumsq;
            }
        }
        public unsafe static double dprotdiv(double* a, double* b, int* fail)
        {
            /*
                dprotdiv returns the value div given by
                    div =	( a/b                 if a/b does not overflow,
                            (
                        ( 0.0                 if a == 0.0,
                        (
                        ( sign( a/b )*flmax   if a != 0.0  and a/b would overflow,

                where  flmax  is a large value. in addition if
                a/b would overflow then  fail is returned as 1, otherwise  fail is
                returned as 0
                note that when  a and b  are both zero, fail is returned as 1, but
                div  is returned as  0.0. in all other cases of overflow  div is such
                that  abs( div ) = flmax

                when  b = 0  then  sign( a/b )  is taken as  sign( a )
            */
            double absb, div;
            int dfail;

            if (fail == null) fail = &dfail;

            if (*a == 0.0)
            {
                div = 0.0;
                *fail = *b == 0 ? 1 : 0;
            }
            else if (*b == 0.0)
            {
                div = BlasLike.dsign(BlasLike.lm_rsafe_range, *a);
                *fail = 1;
            }
            else
            {
                absb = Math.Abs(*b);
                if (absb >= 1.0)
                {
                    *fail = 0;
                    div = (Math.Abs(*a) >= absb * BlasLike.lm_safe_range ? *a / *b : 0.0);
                }
                else if (Math.Abs(*a) <= absb * BlasLike.lm_rsafe_range)
                {
                    *fail = 0;
                    div = *a / *b;
                }
                else
                {
                    *fail = 1;
                    div = BlasLike.lm_rsafe_range;
                    if ((*a < 0.0 && *b > 0.0) || (*a > 0.0 && *b < 0.0)) div = -div;
                }
            }
            return div;
        }
        public unsafe static void dlpprt(bool lp, int nrowa, int Nrowrt, int n, int nclin, int nfree, int isdel, int nactiv,
                int ncolz, int iter, int jadd, int jdel, double alfa, double condt,
                int numinf, double suminf, double objlp, int* istate, int* kfree, double* a,
                double* rt, double* x, double* wrk1, double* wrk2)
        {

            char[] lstate = { ' ', 'L', 'U', 'E', 'T' };

            char ladd, ldel;
            int inct, k, laddi, ldeli;

            --wrk2;
            --wrk1;
            --x;
            rt -= Nrowrt + 1;
            a -= nrowa + 1;
            --kfree;
            --istate;

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
                laddi = istate[jadd];
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
            lm_mdvwri("\nLP VARIABLES", n, &x[1]);
            lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", n, &istate[1]);
            if (nclin > 0)
                lm_mdvwri("\nSTATUS OF THE LP GENERAL CONSTRAINTS", nclin,
                    &istate[n + 1]);
            if (nfree > 0) lm_mdvwri("\nLIST OF FREE LP VARIABLES", nfree, &kfree[1]);
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
                    wrk2[k] = BlasLike.ddot(n, &a[k + nrowa], nrowa, &x[1], 1);
                }
                lm_mdvwri("\nVALUES OF LP GENERAL LINEAR CONSTRAINTS", nclin, &wrk2[1]);
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
                    BlasLike.dcopy(nactiv, &rt[nactiv + (ncolz + 1) * Nrowrt], inct, &wrk1[1], 1);
                    lm_mdvwri("\nDIAGONALS OF LP WORKING SET FACTOR  T", nactiv, &wrk1[1]);
                }
            }
        }
        public unsafe static void lm_mdvwri(string nn, int na, double* wrk)
        {
            Console.WriteLine(nn);
            for (int i = 0; i < na; ++i)
            {
                Console.Write($"{wrk[i]} ");
                if (i % 6 == 5) Console.Write("\n");
            }
            Console.Write("\n");
        }
        public unsafe static void lm_mdvwri(string nn, int na, int* wrk)
        {
            Console.WriteLine(nn);
            for (int i = 0; i < na; ++i)
            {
                Console.Write($"{wrk[i]} ");
                if (i % 6 == 5) Console.Write("\n");
            }
            Console.Write("\n");
        }
        public static void itrwri(string name, int n)
        {
            Console.WriteLine($"{name} iteration {n}");
        }
        public unsafe static void dlpcrsh(bool orthog, int* unitq, int vertex, byte lcrash, int n, int* nclin, int* nctotl,
        int* Nq, int* nrowa, int* Nrowrt, int* Ncolrt, int* nactiv, int* ncolz,
        int* nfree, int* istate, int* kactiv, int* kfree, double* bigbnd, double* tolact,
        double* xnorm, double* a, double* anorm, double* ax, double* bl, double* bu, double* x, double* qtg, double* rt,
        double* zy, double* p, double* wrk1, double* wrk2)
        {

            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1,
                i__2;
            int c__1 = 1;
            int n__ = n;


            int iadd, jadd;
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
            double condmx, cslast;
            int inform;
            double resmin, colsiz, snlast;
            int idummy;
            double rowmax;
            double bnd, res;


            --wrk2;
            --wrk1;
            --p;
            zy_dim1 = *Nq;
            zy_offset = zy_dim1 + 1;
            zy -= zy_offset;
            i__1 = *Nrowrt * *Ncolrt;
            BlasLike.dzerovec(i__1, rt);
            rt_dim1 = *Nrowrt;
            rt_offset = rt_dim1 + 1;
            rt -= rt_offset;
            --qtg;
            --x;
            --bu;
            --bl;
            --ax;
            --anorm;
            a_dim1 = *nrowa;
            a_offset = a_dim1 + 1;
            a -= a_offset;
            --kfree;
            --kactiv;
            --istate;

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
                Console.WriteLine($"{lcrash},{CL(*nclin)}, {CL(*nctotl)}");
                lm_mdvwri("\nLP VARIABLES BEFORE CRASH...", n, &x[1]);
                lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", CN(*nctotl), &istate[1]);
            }
            nfixed = 0;
            *nactiv = 0;
            nartif = 0;
            /*     IF A COLD START IS BEING MADE, INITIALIZE  ISTATE. */
            /*     IF  BL(J) = BU(J),  SET  ISTATE(J)=3  FOR ALL VARIABLES AND LINEAR 
            */
            /*     CONSTRAINTS. */
            if (lcrash > 0) goto L60;
            i__1 = *nctotl;
            for (j = 1; j <= i__1; ++j)
                istate[j] = bl[j] == bu[j] ? 3 : 0;
            /* L40: */
            /*     INITIALIZE  NFIXED,  NACTIV  AND  KACTIV. */
            /*     ENSURE THAT THE NUMBER OF BOUNDS AND GENERAL CONSTRAINTS IN THE */
            /*     WORKING SET DOES NOT EXCEED  N. */
            L60:
            i__1 = *nctotl;
            for (j = 1; j <= i__1; ++j)
            {
                if (nfixed + *nactiv == (int)n || istate[j] == 4) istate[j] = 0;
                if (istate[j] > 0)
                {
                    if (j <= (int)n)
                    {
                        ++nfixed;
                        if (istate[j] == 1) x[j] = bl[j];
                        if (istate[j] >= 2) x[j] = bu[j];
                    }
                    else
                    {
                        ++(*nactiv);
                        if (lcrash < 2) kactiv[*nactiv] = j - n;
                    }
                }
            }
            *nfree = n - nfixed;
            *ncolz = *nfree - *nactiv;
            if (msg >= 80)
            {
                lm_mdvwri("\nLP VARIABLES AFTER CRASH INIT...", n, &x[1]);
                lm_mdvwri("\nSTATUS OF THE LP BOUND   CONSTRAINTS", CN(*nctotl), &istate[1]);
            }
            /*if a hot start is required, the tq factorization is already known*/
            if (lcrash > 1)
            {
                goto L600;
            }
            dtmax = 1;
            dtmin = 1;
            *unitq = 1;
            /*     COMPUTE THE 2-NORMS OF THE CONSTRAINT ROWS. */
            asize = 1;
            if (*nclin == 0)
            {
                goto L140;
            }
            i__1 = *nclin;
            for (j = 1; j <= i__1; ++j)
            {
                anorm[j] = dnrm2(n, &a[j + a_dim1], *nrowa);
                /* L120: */
            }
            var asize_c = asize;
            BlasLike.dxminmax(*nclin, &anorm[1], 1, &asize_c, &amin);
            asize = asize_c;
        L140:
            if (lcrash > 0) goto L320;
            /*
            ---------------------------------------------------------------------
            if a cold start is required, an attempt is made to add as many
            constraints as possible to the working set
            ---------------------------------------------------------------------
            */
            if (nfixed + *nactiv == (int)n)
            {
                goto L460;
            }
            /* see if any variables are outside their bounds */
            for (j = 1; j <= (int)n; ++j)
            {
                if (istate[j] != 0) goto L200;
                b1 = bl[j];
                b2 = bu[j];
                nolow = b1 <= -(*bigbnd);
                noupp = b2 >= *bigbnd;

                was_is = 0;
                if (nolow)
                {
                    goto L160;
                }
                if (x[j] - b1 <= (1 + Math.Abs(b1)) * *tolact)
                {
                    was_is = 1;
                }
            L160:
                if (noupp)
                {
                    goto L180;
                }
                if (b2 - x[j] <= (1 + Math.Abs(b2)) * *tolact)
                {

                    was_is = 2;
                }
            L180:
                if (was_is == 0)
                {
                    goto L200;
                }
                /*        SET VARIABLE EQUAL TO ITS BOUND. */
                istate[j] = was_is;
                if (was_is == 1) x[j] = b1;
                if (was_is == 2) x[j] = b2;
                ++nfixed;
                if (nfixed + *nactiv == (int)n) goto L460;
                L200:
                ;
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
            if (*nclin == 0) goto L320;
            i__1 = *nclin;
            for (i = 1; i <= i__1; ++i)
            {
                j = n + i;
                if (istate[j] > 0)
                {
                    goto L220;
                }
                ax[i] = BlasLike.ddot(n, &a[i + a_dim1], CI(*nrowa), &x[1], 1);
            L220:
                ;
            }
            toobig = *tolact + *tolact;
            i__1 = n;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                resmin = toobig;
                was_is = 0;
                i__2 = *nclin;
                for (i = 1; i <= i__2; ++i)
                {
                    j = n + i;
                    if (istate[j] > 0)
                    {
                        goto L280;
                    }
                    b1 = bl[j];
                    b2 = bu[j];
                    nolow = b1 <= -(*bigbnd);
                    noupp = b2 >= *bigbnd;
                    resl = toobig;
                    resu = toobig;
                    if (nolow)
                    {
                        goto L240;
                    }
                    resl = Math.Abs(ax[i] - b1) / (1 + Math.Abs(b1));
                L240:
                    if (noupp)
                    {
                        goto L260;
                    }
                    resu = Math.Abs(ax[i] - b2) / (1 + Math.Abs(b2));
                L260:
                    res = Math.Min(resl, resu);
                    if (res >= *tolact)
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
                ++(*nactiv);
                kactiv[*nactiv] = imin;
                j = n + imin;
                istate[j] = was_is;
                if (nfixed + *nactiv == (int)n)
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
            *ncolz = n - nfixed - *nactiv;
            if (vertex == 0 || *ncolz == 0)
            {
                goto L460;
            }
            /*     COMPUTE LENGTHS OF COLUMNS OF SELECTED LINEAR CONSTRAINTS */
            /*     (JUST THE ONES CORRESPONDING TO FREE VARIABLES). */
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (istate[j] != 0)
                {
                    goto L380;
                }
                colsiz = 0;
                if (*nclin == 0)
                {
                    goto L360;
                }
                i__2 = *nclin;
                for (k = 1; k <= i__2; ++k)
                {
                    i = n + k;
                    if (istate[i] > 0)
                    {
                        colsiz += Math.Abs(a[k + j * a_dim1]);
                    }
                    /* L340: */
                }
            L360:
                wrk1[j] = colsiz;
            L380:
                ;
            }
            /*     FIND THE  NARTIF  SMALLEST SUCH COLUMNS. */
            /*     THIS IS AN EXPENSIVE LOOP.  LATER WE CAN REPLACE IT */
            /*     BY A 4-PASS PROCESS (SAY), ACCEPTING THE FIRST COL THAT */
            /*     IS WITHIN  T  OF  COLMIN, WHERE  T = 0.0, 0.001, 0.01, 0.1 (SAY). 
            */
            i__1 = *ncolz;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                colmin = BlasLike.lm_max;
                i__2 = n;
                for (j = 1; j <= i__2; ++j)
                {
                    if (istate[j] != 0)
                    {
                        goto L400;
                    }
                    if (*nclin == 0)
                    {
                        goto L420;
                    }
                    colsiz = wrk1[j];
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
                istate[j] = 4;
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
            *nfree = 0;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (istate[j] != 0)
                {
                    goto L480;
                }
                ++(*nfree);
                kfree[*nfree] = j;
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
            *ncolz = *nfree;
            if (*nactiv == 0)
            {
                goto L540;
            }
            nact1 = *nactiv;
            *nactiv = 0;
            dtqadd(orthog, unitq, &inform, &c__1, &nact1, nactiv, ncolz, nfree, &n__,
                Nq, nrowa, Nrowrt, Ncolrt, &istate[1], &kactiv[1], &kfree[
                1], &condmx, &a[a_offset], &qtg[1], &rt[rt_offset], &zy[zy_offset]
                , &wrk1[1], &wrk2[1]);
            /*     IF A VERTEX IS REQUIRED BUT  TQADD  WAS UNABLE TO ADD ALL OF THE */

            /*     SELECTED GENERAL CONSTRAINTS, ADD MORE TEMPORARY BOUNDS. */
            if (vertex != 0 || *ncolz == 0)
            {
                goto L540;
            }
            ncolz1 = *ncolz;
            i__1 = ncolz1;
            for (idummy = 1; idummy <= i__1; ++idummy)
            {
                rowmax = 0;
                i__2 = *nfree;
                for (i = 1; i <= i__2; ++i)
                {
                    rnorm = dnrm2(*ncolz, &zy[i + zy_dim1], *Nq);
                    if (rowmax >= rnorm)
                    {
                        goto L500;
                    }
                    rowmax = rnorm;
                    ifix = i;
                L500:
                    ;
                }
                jadd = kfree[ifix];
                inform = daddcon(false, false, orthog, unitq, &ifix, &iadd, &jadd,
                    nactiv, ncolz, ncolz, nfree, n, Nq, nrowa, Nrowrt, &
                    kfree[1], &condmx, &cslast, &snlast, &a[a_offset], &qtg[1], &
                    rt[rt_offset], &zy[zy_offset], &wrk1[1], &wrk2[1]);
                --(*nfree);
                --(*ncolz);
                ++nartif;
                istate[jadd] = 4;
                /* L520: */
            }
        /*     SET ELEMENTS  NACTIV + 1, ......, NACTIV + NFIXED  OF  KACTIV TO */

        /*     POINT TO THE FIXED VARIABLES. */
        L540:
            kb = *nactiv;
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (istate[j] == 0) continue;
                ++kb;
                kactiv[kb] = j;
            }
            /* --------------------------------------------------------------------- 
            */
            /*     THE TQ FACTORIZATION HAS BEEN COMPUTED.  FIND THE POINT CLOSEST TO 
            */
            /*     THE USER-SUPPLIED  X  THAT LIES ON THE INITIAL WORKING SET. */
            /* --------------------------------------------------------------------- 
            */
            /*     SET WRK1 = RESIDUALS FOR CONSTRAINTS IN THE WORKING SET. */
            if (*nactiv == 0)
            {
                goto L600;
            }
            i__1 = *nactiv;
            for (i = 1; i <= i__1; ++i)
            {
                k = kactiv[i];
                j = n + k;
                bnd = bl[j];
                if (istate[j] > 1)
                {
                    bnd = bu[j];
                }
                wrk1[i] = bnd - BlasLike.ddot(n, &a[k + a_dim1], CI(*nrowa), &x[1], 1);
                /* L580: */
            }
            /*     SOLVE FOR P, THE SMALLEST CORRECTION TO X THAT GIVES A POINT */
            /*     ON THE CONSTRAINTS IN THE WORKING SET. */
            /*     FIRST SOLVE  T*WRK1 = RESIDUALS, THEN GET  P = Y*WRK1. */
            idiag = 1;
            drtmxsolve(2, CN(*nactiv), &rt[(*ncolz + 1) * rt_dim1 + 1], CN(*Nrowrt), &wrk1[1],
                &idiag);
            BlasLike.dzerovec(n, &p[1]);
            BlasLike.dcopyvec(*nactiv, &wrk1[1], &p[*ncolz + 1]);
            dzyprod(2, n, *nactiv, *ncolz, *nfree, *Nq, *unitq, &kactiv[1], &kfree[1],
                &p[1], &zy[zy_offset], &wrk1[1]);
            BlasLike.daxpy(n, CR(1), &p[1], 1, &x[1], 1);
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE 2-NORM OF  X. */
        /*     INITIALIZE  AX  FOR ALL GENERAL CONSTRAINTS. */
        /* --------------------------------------------------------------------- 
        */
        L600:
            *xnorm = dnrm2vec(n, &x[1]);
            if (*nclin == 0)
            {
                goto L640;
            }
            BlasLike.dzerovec(*nclin, &ax[1]);
            i__1 = n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[j] != 0)
                {
                    BlasLike.daxpy(*nclin, x[j], &a[j * a_dim1 + 1], 1, &ax[1], 1);
                }
                /* L620: */
            }
        /*     A POINT THAT SATISFIES THE INITIAL WORKING SET HAS BEEN FOUND. */
        L640:
            *ncolz = *nfree - *nactiv;
            nfixed = n - *nfree;
            if (msg >= 80)
            {
                //  lm_wmsg(
                //"\nLPCRSH. WORKING SET SELECTED ...\nBOUNDS = %ld TEMPORARY BOUNDS = %ld GENERAL LINEAR = %ld",
                //  CL(nfixed), CL(nartif), CL(*nactiv));
                Console.WriteLine($"{nfixed},{nartif}, {CL(*nactiv)}");
                lm_mdvwri("\nLP VARIABLES AFTER  CRASH...", n, &x[1]);
            }
        }
        public unsafe static void dlpgrad(bool lp, int n, int nctotl, int nrowa, double bigbnd, double feamin,
        int* numinf, double* suminf, int* istate, double* a, double* bl,
        double* bu, double* cvec, double* featol, double* grad, double* x)
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

            --x;
            --grad;
            --featol;
            --cvec;
            --bu;
            --bl;
            a -= nrowa + 1;
            --istate;

            if (*numinf != 0)
            {
                *numinf = 0;
                *suminf = 0;
                BlasLike.dzerovec(n, &grad[1]);
                for (j = 1; j <= (int)nctotl; ++j)
                {
                    /* do nothing if the variable or constraint is at a bound */
                    if (istate[j] > 0) goto L60;
                    feasj = featol[j];
                    nolow = bl[j] <= -bigbnd;
                    noupp = bu[j] >= bigbnd;
                    k = j - n;
                    if (j <= (int)n)
                    {
                        atx = x[j];
                    }
                    if (j > (int)n)
                    {
                        atx = BlasLike.ddot(n, &a[k + nrowa], nrowa, &x[1], 1);
                    }
                    istate[j] = 0;
                    /* see if the lower bound is violated */
                    if (nolow)
                    {
                        goto L20;
                    }
                    s = bl[j] - atx;
                    if (s <= feasj)
                    {
                        goto L20;
                    }
                    istate[j] = -2;
                    weight = -feamin / feasj;
                    goto L40;
                /*see if the upper bound is violated*/
                L20:
                    if (noupp)
                    {
                        goto L60;
                    }
                    s = atx - bu[j];
                    if (s <= feasj)
                    {
                        goto L60;
                    }
                    istate[j] = -1;
                    weight = feamin / feasj;
                /* add the infeasibility */
                L40:
                    ++(*numinf);
                    *suminf += Math.Abs(weight) * s;
                    if (j <= (int)n)
                    {
                        grad[j] = weight;
                    }
                    if (j > (int)n) BlasLike.daxpy(n, weight, &a[k + nrowa], nrowa, &grad[1], 1);
                    L60:
                    ;
                }
            }
            /* if feasible, install true objective */
            if (lp && *numinf == 0) BlasLike.dcopyvec(n, &cvec[1], &grad[1]);
        }
        public unsafe static void dzyprod(short mode, int n, int nactiv, int ncolz, int nfree, int nq, int unitq, int* kactiv, int* kfree, double* v, double* zy, double* wrk)
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
            --wrk;
            zy -= nq + 1;
            --v;
            --kfree;
            --kactiv;

            nfixed = n - nfree;
            j1 = 1;
            j2 = nfree;
            if (mode == 1 || mode == 4) j2 = ncolz;
            else if (mode == 2 || mode == 5 || mode == 7) j1 = ncolz + 1;
            lenv = j2 - j1 + 1;
            if (mode < 4)
            {
                /*MODE = 1, 2  OR  3*/
                if (nfree > 0) BlasLike.dzerovec(nfree, &wrk[1]);
                /*COPY  V(FIXED)  INTO THE END OF  WRK. */
                if (mode != 1 && nfixed != 0)
                    BlasLike.dcopyvec(nfixed, &v[nfree + 1], &wrk[nfree + 1]);
                /*SET  WRK  =  RELEVANT PART OF  ZY * V. */
                if (lenv > 0)
                {
                    if (unitq==1) BlasLike.dcopyvec(lenv, &v[j1], &wrk[j1]);
                    else for (j = j1; j <= j2; ++j)
                            if (v[j] != 0)
                                BlasLike.daxpy(nfree, v[j], &zy[j * nq + 1], 1, &wrk[1], 1);
                }
                /*EXPAND  WRK  INTO  V  AS A FULL N-VECTOR. */
                BlasLike.dzerovec(n, &v[1]);
                if (nfree > 0)
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k];
                        v[j] = wrk[k];
                    }
                /*COPY  WRK(FIXED)  INTO THE APPROPRIATE PARTS OF  V*/
                if (mode == 1 || nfixed == 0) return;
                for (l = 1; l <= nfixed; ++l)
                {
                    kw = nfree + l;
                    ka = nactiv + l;
                    j = kactiv[ka];
                    v[j] = wrk[kw];
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
                        j = kactiv[ka];
                        wrk[kw] = v[j];
                    }
                /*PUT THE FREE  COMPONENTS OF  V  INTO THE BEGINNING OF  WRK. */
                if (nfree != 0)
                {
                    for (k = 1; k <= (int)nfree; ++k)
                    {
                        j = kfree[k];
                        wrk[k] = v[j];
                    }
                    /*SET  V  =  RELEVANT PART OF  ZY(T) * WRK*/
                    if (lenv > 0)
                    {
                        if (unitq==1) BlasLike.dcopyvec(lenv, &wrk[j1], &v[j1]);
                        else for (j = j1; j <= j2; ++j)
                                v[j] = BlasLike.ddotvec(nfree, zy + j * nq + 1, wrk + 1);
                    }
                }
                /*COPY THE FIXED COMPONENTS OF WRK INTO THE END OF  V*/
                if (mode != 4 && mode <= 6 && nfixed != 0)
                    BlasLike.dcopyvec(nfixed, &wrk[nfree + 1], &v[nfree + 1]);
            }
        }
        public unsafe static void dgetlamd(string lprob, int n, int nactiv, int ncolz, int nfree, int nrowa, int Nrowrt, int* jsmlst, int* ksmlst, double* smllst, int* istate, int* kactiv, double* a, double* anorm, double* qtg, double* rlamda, double* rt)
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

            rt -= Nrowrt + 1;
            --rlamda;
            --qtg;
            --anorm;
            a -= nrowa + 1;
            --kactiv;
            --istate;

            /*
            first, compute the lagrange multipliers for the general
            constraints in the working set, by solving
            t(transpose)*rlamda = y(t)*grad
            */
            nfixed = n - nfree;
            nlam = nfixed + nactiv;
            if (nactiv > 0)
            {
                BlasLike.dcopyvec(nactiv, &qtg[ncolz + 1], &rlamda[1]);
                idiag = 1;
                drtmxsolve(-2, nactiv, &rt[(ncolz + 1) * Nrowrt + 1], Nrowrt,
                    &rlamda[1], &idiag);
            }
            /*
            now set elements nactiv, nactiv+1,... of rlamda equal to the
            multipliers for the bound constraints in the working set
            */
            for (l = 1; l <= nfixed; ++l)
            {
                kb = nactiv + l;
                j = kactiv[kb];
                jgfxd = nfree + l;
                blam = qtg[jgfxd];
                for (ka = 1; ka <= nactiv; ++ka)
                {
                    i = kactiv[ka];
                    blam -= a[i + j * nrowa] * rlamda[ka];
                }
                rlamda[kb] = blam;
            }

            /*find  allmax and smllst*/
            *smllst = BlasLike.lm_max;
            *jsmlst = 0;
            *ksmlst = 0;
            for (k = 1; k <= nlam; ++k)
            {
                j = kactiv[k];
                if (k > nactiv) anormj = 1;
                else
                {
                    anormj = anorm[j];
                    j += n;
                }
                was_is = istate[j];
                /*
                change the sign of the estimate if the constraint is in the
                working set (or violated) at its upper bound
                if (was_is == 3) rlam = Math.Abs(rlam);
                */
                if (was_is != 3)
                {
                    rlam = rlamda[k] * anormj;
                    /*not a fixed variable or an equality constraint*/
                    if (was_is == 2) rlam = -rlam;
                    else if (was_is == 4) rlam = -Math.Abs(rlam);
                    /*find the smallest multiplier for the inequalities*/
                    if (rlam < *smllst)
                    {
                        *smllst = rlam;
                        *jsmlst = j;
                        *ksmlst = k;
                    }
                }
            }

            /*if required, print the multipliers*/
            if (msg < 20) return;
            if (nactiv > 0) w_lam(lprob, "CONSTRAINTS...", nactiv, &kactiv[1], &rlamda[1]);
            l = nactiv + 1;
            if (l <= nlam) w_lam(lprob, "BOUND CONSTRAINTS...", nlam - nactiv, &kactiv[l], &rlamda[l]);
            if (msg >= 80)
                lm_wmsg("\n//dgetlam//  JSMLST     SMLLST     KSMLST\n//dgetlam//%8ld%11.2lg%11ld",
                    CL(*jsmlst), *smllst, CL(*ksmlst));
        }
        public unsafe static void w_lam(string msg1, string msg2, int n, int* a, double* r)
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
        public static void lm_wmsg<T>(string mess, T n1)
        {
            Console.WriteLine($"{mess} {n1}");
        }
        public static void lm_wmsg(string mess, int n1, double n2, int n3)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3}");
        }
        public static void lm_wmsg(string mess, double n1, double n2, double n3)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3}");
        }
        public static void lm_wmsg(string mess, double n1, double n2, double n3, double n4)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4}");
        }
        public static void lm_wmsg<T>(string mess, int n1, int n2, int n3, int n4, int n5, T n6)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4} {n5} {n6}");
        }
        public static void lm_wmsg(string mess, string n1, string n2, string n3, int n4, double n5, double n6)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4} {n5} {n6}");
        }
        public static void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, bool n6, int n7, int n8, int n9, int n10)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9} {n10}");
        }
        public static void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, bool n6, double n7, int n8, double n9)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9}");
        }
        public static void lm_wmsg(string mess, int n1, int n2, double n3, double n4, double n5, int n6, double n7, int n8, double n9)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} {n9}");
        }
        public static void lm_wmsg(string mess, string n1, string n2, int n3)
        {
            Console.WriteLine($"{mess} {n1} {n2} {n3}");
        }
        public static void lm_wmsg<T>(string mess, string n1, T n2)
        {
            Console.WriteLine($"{mess} {n1} {n2}");
        }
        public static void lm_wmsg(string mess, double n1, double n2)
        {
            Console.WriteLine($"{mess} {n1} {n2}");
        }
        public unsafe static void dlpbgst(int n, int nactiv, int nfree, int* jbigst, int* kbigst,
        int* istate, int* kactiv, double dinky, double feamin, double* trulam, double* featol, double* rlamda)
        {
            double rlam, biggst;
            int nlam, j, k, was_is, nfixed;

            --rlamda;
            --featol;
            --kactiv;
            --istate;

            *jbigst = 0;
            nfixed = n - nfree;
            nlam = nfixed + nactiv;
            if (nlam == 0) return;
            biggst = 1 + dinky;
            for (k = 1; k <= nlam; ++k)
            {
                j = kactiv[k];
                if (k <= nactiv) j += n;

                was_is = istate[j];
                if (was_is >= 1)
                {
                    rlam = rlamda[k];
                    if (was_is == 2) rlam = -rlam;
                    if (was_is == 3) rlam = Math.Abs(rlam);
                    rlam = featol[j] / feamin * rlam;
                    if (biggst < rlam)
                    {
                        biggst = rlam;
                        *trulam = rlamda[k];
                        *jbigst = j;
                        *kbigst = k;
                    }
                }
            }
            if (msg >= 80) lm_wmsg("\n//LPBGST// JBIGST         BIGGST\n//LPBGST//%7ld%15.4lg", CL(*jbigst), biggst);
        }
        public unsafe static void ddelcon(bool modfyg, bool orthog, int unitq, int jdel, int kdel, int nactiv, int ncolz, int nfree, int n, int Nq, int nrowa, int Nrowrt, int* kactiv, int* kfree, double* a, double* qtg, double* rt, double* zy)
        {/*
	ddelcon updates the factorization of the matrix of
	constraints in the working set,  A(free)*(Z Y) = (0 T)
	if there are no general constraints in the working set and the
	matrix  Q = (Z Y)  is the identity, Q will not be touched
*/
            int i, j, k, l, ldiag;
            int nfree1, nactp1, nactv1, ka;
            double store;
            double cs, sn;
            int ibegin, ifreed;
            int nfreei, nactpi, istore;


            zy -= Nq + 1;
            rt -= Nrowrt + 1;
            --qtg;
            a -= nrowa + 1;
            --kfree;
            --kactiv;

            if (jdel <= n)
            {
                /*A SIMPLE BOUND IS BEING DELETED FROM THE WORKING SET. */
                ifreed = kdel - nactiv;
                if (msg >= 80)
                    lm_wmsg("BOUND DELETED", CL(nactiv), CL(ncolz), CL(nfree), CL(ifreed),
                        CL(jdel), unitq);
                nactv1 = nactiv;
                nfree1 = nfree + 1;
                ibegin = 1;
                kfree[nfree1] = jdel;

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
                    store = qtg[nfree1];
                    qtg[nfree1] = qtg[nfreei];
                    qtg[nfreei] = store;
                    istore = kactiv[nactp1];
                    kactiv[nactp1] = kactiv[nactpi];
                    kactiv[nactpi] = istore;
                }

                /*COPY THE INCOMING COLUMN OF A INTO THE END OF  T. */
                if (unitq!=0)
                {
                    for (ka = 1; ka <= nactiv; ++ka)
                    {
                        i = kactiv[ka];
                        rt[ka + nfree1 * Nrowrt] = a[i + jdel * nrowa];
                    }
                    /*EXPAND  Q  BY ADDING A UNIT ROW AND COLUMN. */
                    BlasLike.dzero(nfree, &zy[nfree1 + Nq], Nq);
                    BlasLike.dzerovec(nfree, &zy[nfree1 * Nq + 1]);
                    zy[nfree1 + nfree1 * Nq] = 1;
                }
            }
            else
            {
                /*A GENERAL CONSTRAINT IS BEING DELETED FROM THE WORKING SET*/
                if (msg >= 80)
                    lm_wmsg("CONSTRAINT DELETED", CL(nactiv), CL(ncolz), CL(nfree), CL(kdel),
                        CL(jdel), unitq);
                nactv1 = nactiv - 1;
                nfree1 = nfree;
                ibegin = kdel;
                /*DELETE A ROW OF  T  AND MOVE THE ONES BELOW IT UP. */
                for (i = kdel; i <= nactv1; ++i)
                {
                    j = i + 1;
                    kactiv[i] = kactiv[j];
                    ldiag = nfree - i;
                    BlasLike.dcopy(j, rt + j + ldiag * Nrowrt, Nrowrt, rt + i + ldiag * Nrowrt, Nrowrt);
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
                    delmgen(orthog, rt + i + (k + 1) * Nrowrt, rt + i + k * Nrowrt, &cs, &sn);
                    if (l > 0) delm(orthog, l, rt + i + 1 + (k + 1) * Nrowrt, 1, rt + i + 1 + k * Nrowrt, 1, cs, sn);
                    if (nactv1 > 0) delm(orthog, nfree1, zy + (k + 1) * Nq + 1, 1, zy + k * Nq + 1, 1, cs, sn);
                    if (modfyg) delm(orthog, 1, qtg + k + 1, 1, qtg + k, 1, cs, sn);
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
                    kactiv[j] = kactiv[j + 1];
                    ++j;
                }
            }
            /*ESTIMATE THE CONDITION NUMBER OF  T. */
            if (nactv1 > 0)
            {
                var dtmax_c = dtmax;
                var dtmin_c = dtmin;
                BlasLike.dxminmax(nactv1, &rt[nactv1 + (ncolz + 2) * Nrowrt], Nrowrt - 1, &dtmax_c, &dtmin_c);
                dtmax = dtmax_c;
                dtmin = dtmin_c;
            }
        }
        public unsafe static void dfindp(bool* nullr, bool* unitpg, int* unitq, int n, int* nclin, int* Nq, int* nrowa, int* Nrowrt, int* ncolr, int* ncolz, int* nfree, int* istate, int* kfree, bool negligible, double* gtp, double* pnorm, double* rdlast, double* a, double* ap, double* p, double* qtg, double* rt, double* v, double* zy, double* work)
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

            --work;
            zy_dim1 = *Nq;
            zy_offset = zy_dim1 + 1;
            zy -= zy_offset;
            --v;
            rt_dim1 = *Nrowrt;
            rt_offset = rt_dim1 + 1;
            rt -= rt_offset;
            --qtg;
            --p;
            --ap;
            a_dim1 = *nrowa;
            a_offset = a_dim1 + 1;
            a -= a_offset;
            --kfree;
            --istate;

            BlasLike.dcopyvec(*ncolr, &qtg[1], &p[1]);
            BlasLike.dscalvec(*ncolr, CR(-1), &p[1]);
            if (*nullr) goto L60;
            *rdlast = rt[*ncolr + *ncolr * rt_dim1];
            /*     *** */
            /*     CORRECTION INSERTED BY MHW, 22 OCT 1985. */
            /*     THIS ENSURES A NON-ZERO SEARCH DIRECTION. */
            /*     *** */
            if (*ncolr < *ncolz && negligible) p[*ncolr] = *rdlast;
            /* --------------------------------------------------------------------- 
            */
            /*     SOLVE THE SYSTEM   R(T)R (PZ) = - Z(T)G(FREE). */
            /* --------------------------------------------------------------------- 
            */
            if (*unitpg)
            {
                goto L20;
            }
            /*     PERFORM THE FORWARD SUBSTITUTION  R(T)V = - Z(T)G(FREE). */
            idiag = dtmxsolve(CS((short)-1), *ncolr, &rt[rt_offset], *Nrowrt, &p[1], CS((short)1));
            goto L40;
        /*     THE PROJECTED GRADIENT IS A MULTIPLE OF THE UNIT VECTOR, THE */
        /*     FORWARD SUBSTITUTION MAY BE AVOIDED. */
        L20:
            if (negligible) p[*ncolr] = -1;
            else p[*ncolr] /= *rdlast;
            /*     PERFORM THE BACKWARD SUBSTITUTION   R(PZ) = P. */
            L40:
            BlasLike.dcopyvec(*ncolr, &p[1], &v[1]);
            idiag = dtmxsolve(CS((short)1), *ncolr, &rt[rt_offset], *Nrowrt, &p[1], CS((short)1));
        /* --------------------------------------------------------------------- 
        */
        /*     THE VECTOR  (PZ)  HAS BEEN COMPUTED. */
        /* --------------------------------------------------------------------- 
        */
        /*     COMPUTE THE DIRECTIONAL DERIVATIVE  G(T)P = (GZ)(T)(PZ). */
        L60:
            *gtp = BlasLike.ddotvec(*ncolr, &qtg[1], &p[1]);
            /* --------------------------------------------------------------------- 
            */
            /*     COMPUTE  P = Z * PZ. */
            /* --------------------------------------------------------------------- 
            */
            /*     NACTIV  AND  KACTIV  ARE NOT USED IN  ZYPROD.  N  AND  KFREE */
            /*     SERVE AS ARGUMENTS FOR  NACTIV  AND  KACTIV. */
            dzyprod(1, n, n, *ncolr, *nfree, *Nq, *unitq, &kfree[1], &kfree[1], &p[1],
                 &zy[zy_offset], &work[1]);
            *pnorm = dnrm2vec(*nfree, &work[1]);
            if (msg >= 80) lm_mdvwri("\n//FINDP//   P ... ", n, &p[1]);
            /* --------------------------------------------------------------------- 
            */
            /*     COMPUTE  AP. */
            /* --------------------------------------------------------------------- 
            */
            if (*nclin > 0)
            {
                BlasLike.dzerovec(*nclin, &ap[1]);
                for (j = 1; j <= (int)n; ++j)
                {
                    if (istate[j] <= 0)
                        BlasLike.daxpy(*nclin, p[j], &a[j * a_dim1 + 1], 1, &ap[1], 1);
                    /* L100: */
                }
                if (msg >= 80) lm_mdvwri("\n//FINDP//  AP ... ", CN(*nclin), &ap[1]);
            }
            return;
        }
        public unsafe static short dbndalf(bool firstv, int* hitlow, int* istate, int* jadd, int n, int nctotl, int numinf, double* alfa, double* palfa, double* atphit, double* bigalf, double* bigbnd, double* pnorm, double* anorm, double* ap, double* ax, double* bl, double* bu, double* featol, double* p, double* x)
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
            int jadd1;
            double alfa1, alfa2;
            int jadd2;
            bool hlow1, hlow2, lastv;
            int i, j;
            double palfa1, palfa2, apmax1, apmax2;
            int jsave1, jsave2;
            double epspt9;
            int js;
            double absatp;
            double rownrm, atp, res, atx, atp1, atp2;
            --x;
            --p;
            --featol;
            --bu;
            --bl;
            --ax;
            --ap;
            --anorm;
            --istate;

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
            dbdpert(firstv, false, bigalf, bigbnd, pnorm, &jadd1, &jadd2, &palfa1, &
                palfa2, &istate[1], n, nctotl, &anorm[1], &ap[1], &
                ax[1], &bl[1], &bu[1], &featol[1], &p[1], &x[1]);
            jsave1 = jadd1;
            jsave2 = jadd2;

            /*
            SECOND PASS -- RECOMPUTE STEP-LENGTHS WITHOUT PERTURBATION
            AMONGST CONSTRAINTS THAT ARE CLOSE TO THE PERTURBED STEPS
            CHOOSE THE ONE (OF EACH TYPE) THAT MAKES THE LARGEST ANGLE
            WITH THE SEARCH DIRECTION
            */
            if (msg == 99) Console.WriteLine(
        "BNDALF ENTERED\n    J  JS         FEATOL         AX             AP     JADD1        ALFA1     JADD2        ALFA2");
            alfa1 = *bigalf;
            alfa2 = 0.0;
            if (firstv) alfa2 = *bigalf;
            apmax1 = 0.0;
            apmax2 = 0.0;
            atp1 = 0.0;
            atp2 = 0.0;
            hlow1 = false;
            hlow2 = false;
            lastv = !firstv;
            for (j = 1; j <= nctotl; ++j)
            {
                js = istate[j];
                if (js > 0) continue;
                if (j > n)
                {
                    /*GENERAL LINEAR CONSTRAINT. */
                    i = j - n;
                    atx = ax[i];
                    atp = ap[i];
                    /*			lm_wmsg((char*)"atx %e  atp %e  istate[%d]=%d",atx,atp,j,js); */
                    rownrm = anorm[i] + 1.0;
                }
                else
                {
                    /*BOUND CONSTRAINT. */
                    atx = x[j];
                    atp = p[j];
                    rownrm = 1.0;
                }
                if (Math.Abs(atp) <= epspt9 * rownrm * *pnorm) res = -1.0;
                else if (atp > BlasLike.lm_eps)
                {
                    /*ATX IS INCREASING. */
                    /*TEST FOR SMALLER ALFA1 IF UPPER BOUND IS SATISFIED. */
                    if (js != -1)
                    {
                        if (bu[j] < *bigbnd)
                        {
                            res = bu[j] - atx;
                            if (((int)j == jsave1 || palfa1 * atp >= res)
                                && apmax1 * rownrm * *pnorm < atp)
                            {
                                apmax1 = atp / (rownrm * *pnorm);
                                alfa1 = res / atp;
                                jadd1 = j;
                                atp1 = atp;
                                hlow1 = false;
                            }
                        }
                        /*TEST FOR BIGGER ALFA2 IF LOWER BOUND IS VIOLATED. */
                        if (js == -2)
                        {
                            res = bl[j] - atx;
                            if ((firstv || (int)j == jsave2 || palfa2 * atp <= res)
                                && (lastv || (int)j == jsave2 || palfa2 * atp >= res)
                                && (apmax2 * rownrm * *pnorm < atp))
                            {
                                apmax2 = atp / (rownrm * *pnorm);
                                if (atp >= 1.0) alfa2 = res / atp;
                                else if (res < *bigalf * atp) alfa2 = res / atp;
                                else alfa2 = *bigalf;
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
                    if (bl[j] > -(*bigbnd))
                    {
                        res = atx - bl[j];
                        if (((int)j == jsave1 || palfa1 * absatp >= res)
                            && (apmax1 * rownrm * *pnorm < absatp))
                        {
                            apmax1 = absatp / (rownrm * *pnorm);
                            alfa1 = res / absatp;
                            jadd1 = j;
                            atp1 = atp;
                            hlow1 = true;
                        }
                    }
                    /*TEST FOR BIGGER ALFA2 IF UPPER BOUND IS VIOLATED*/
                    if (js == -1)
                    {
                        res = atx - bu[j];
                        if ((firstv || (int)j == jsave2 || palfa2 * absatp <= res)
                            && (lastv || (int)j == jsave2 || palfa2 * absatp >= res)
                            && (apmax2 * rownrm * *pnorm < absatp))
                        {
                            apmax2 = absatp / (rownrm * *pnorm);
                            if (absatp >= 1.0) alfa2 = res / absatp;
                            else if (res < *bigalf * absatp) alfa2 = res / absatp;
                            else alfa2 = *bigalf;
                            jadd2 = j;
                            atp2 = atp;
                            hlow2 = false;
                        }
                    }
                }
                if (msg == 99)
                    lm_wmsg("%5ld%4ld%15.5lg%15.5lg%15.5lg%6ld%17.7lg%6ld%17.7lg",
                        CL(j), CL(js), featol[j], atx, atp, CL(jadd1), alfa1,
                        CL(jadd2), alfa2);
            }

            /*IF FEASIBLE, ONLY ALFA1 WILL HAVE BEEN SET. */
            *alfa = alfa1;
            *palfa = palfa1;
            *jadd = jadd1;
            *atphit = atp1;
            *hitlow = hlow1 ? 1 : 0;
            if (numinf != 0 && jadd2 != 0 && (alfa2 < alfa1 || (alfa2 <= palfa1 && apmax2 >= apmax1)))
            {
                /*
                INFEASIBLE -- SEE IF WE STEP TO THE FURTHEST VIOLATED CONSTRAINT. 
                BE PREPARED TO STEP IN THE RANGE  (ALFA1, PALFA1)  IF THE VIOLATED 
                CONSTRAINT HAS A LARGER VALUE OF  AP
                */
                *alfa = alfa2;
                *jadd = jadd2;
                *atphit = atp2;
                *hitlow = hlow2 ? 1 : 0;
            }
            else if (*alfa < -BlasLike.lm_eps)
            {
                /*
                NEGATIVE STEP
                JADD  WILL RETAIN ITS CURRENT VALUE, BUT WE MAY SHORTEN  ALFA
                TO BE  - PALFA1,  THE STEP TO THE NEAREST PERTURBED SATISFIED
                CONSTRAINT ALONG THE DIRECTION -P
                */
                dbdpert(firstv, true, bigalf, bigbnd, pnorm, &jadd1, &jadd2, &palfa1, &
                    palfa2, &istate[1], n, nctotl, &anorm[1], &ap[1], &
                    ax[1], &bl[1], &bu[1], &featol[1], &p[1], &x[1]);
                if (msg >= 80) lm_wmsg("NEGATIVE STEP", *alfa, palfa1);
                d__1 = Math.Abs(*alfa);
                *alfa = -Math.Min(d__1, palfa1);
            }

            /*
            TEST FOR UNDEFINED OR INFINITE STEP.  THIS SHOULD MEAN THAT THE
            SOLUTION IS UNBOUNDED
            */
            if (*jadd == 0)
            {
                *alfa = *bigalf;
                *palfa = *bigalf;
                inform = 2;
            }
            if (*alfa >= *bigalf) inform = 3;
            if (msg >= 80 && inform > 0) lm_wmsg(
        "\n//BNDALF//  UNBOUNDED STEP.\n//BNDALF//  JADD          ALFA\n//BNDALF//  %4ld%15.5lg",
                CL(*jadd), *alfa);
            return inform;
        }
        public unsafe static
short daddcon(bool modfyg, bool modfyr, bool orthog, int* unitq, int* ifix, int* iadd, int* jadd, int* nactiv, int* ncolr, int* ncolz, int* nfree, int n, int* Nq, int* nrowa, int* Nrowrt, int* kfree, double* condmx, double* cslast, double* snlast, double* a, double* qtg, double* rt, double* zy, double* wrk1, double* wrk2)
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

            double beta, cond;
            int nelm, inct, nact1;
            double d;
            int i, j, k;
            int ldiag;
            int ifail;
            double delta;
            double dtnew = -9999999999;
            int iswap;
            int lrowr, nfree1 = 0, ncolz1;
            double cs, condbd, sn;
            double tdtmin, tdtmax;
            int itrans, kp1;

            --wrk2;
            --wrk1;
            zy_dim1 = *Nq;
            zy_offset = zy_dim1 + 1;
            zy -= zy_offset;
            rt_dim1 = *Nrowrt;
            rt_offset = rt_dim1 + 1;
            rt -= rt_offset;
            --qtg;
            a_dim1 = *nrowa;
            a_offset = a_dim1 + 1;
            a -= a_offset;
            --kfree;

            /*
            if the condition estimator of the updated factors is greater than
            condbd,  a warning message is printed
            */
            condbd = Math.Pow(BlasLike.lm_eps, -0.9);
            ncolz1 = *ncolz - 1;
            if (*jadd > (int)n)
            {
                goto L60;
            }
            /*
                 a simple bound has entered the working set.  iadd  is not used
            */
            if (msg >= 80)
                lm_wmsg(
            "\n//ADDCON//  SIMPLE BOUND ADDED.\n//ADDCON// NACTIV NCOLZ NFREE  IFIX  JADD UNITQ\n//ADDCON//%7ld%6ld%6ld%6ld%6ld%6ld",
                    CL(*nactiv), CL(*ncolz), CL(*nfree), CL(*ifix),
                    CL(*jadd), CL(*unitq));
            /*
                 set  wrk1 = appropriate row of  q
                 reorder the elements of  kfree (this requires reordering the
                 corresponding rows of  q)
            */
            nfree1 = *nfree - 1;
            nact1 = *nactiv;
            if (*unitq!=0)
            {
                goto L20;
            }
            /*
                 q  is stored explicitly.  interchange components  ifix  and  nfree
                 of  kfree  and swap the corresponding rows of  q
            */
            BlasLike.dcopy(*nfree, &zy[*ifix + zy_dim1], *Nq, &wrk1[1], 1);
            if (*ifix == *nfree) goto L180;
            kfree[*ifix] = kfree[*nfree];
            BlasLike.dcopy(*nfree, &zy[*nfree + zy_dim1], *Nq, &zy[*ifix + zy_dim1], *Nq);
            goto L180;
        /*
             q  is not stored, but  kfree  defines an ordering of the columns
             of the identity matrix that implicitly define  z
             reorder  kfree  so that variables  ifix+1,...,nfree  are moved one
             position to the left
        */
        L20:
            BlasLike.dzerovec(*nfree, &wrk1[1]);
            wrk1[*ifix] = 1.0;
            if (*ifix == *nfree)
            {
                goto L180;
            }
            i__1 = nfree1;
            for (i = *ifix; i <= i__1; ++i)
            {
                kfree[i] = kfree[i + 1];
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
                    CL(*nactiv), CL(*ncolz), CL(*nfree), CL(*iadd), CL(*jadd), CL(*unitq));
            nact1 = *nactiv + 1;
            /*
                 transform the incoming row of  a  by  q(t)
            */
            BlasLike.dcopy(n, &a[*iadd + a_dim1], *nrowa, &wrk1[1], 1);
            dzyprod(8, n, *nactiv, *ncolz, *nfree, *Nq, *unitq, &kfree[1], &kfree[1], &
                wrk1[1], &zy[zy_offset], &wrk2[1]);
            if ((*unitq)==0)
            {
                goto L100;
            }
            /*
                 this is the first general constraint to be added  --  set  q = i.
            */
            i__1 = *nfree;
            for (j = 1; j <= i__1; ++j)
            {
                BlasLike.dzerovec(*nfree, &zy[j * zy_dim1 + 1]);
                zy[j + j * zy_dim1] = 1.0;
            }
            *unitq = 0;
        /*
             check that the incoming row is not dependent upon those
             already in the working set
        */
        L100:
            dtnew = dnrm2vec(*ncolz, &wrk1[1]);
            if (nact1 > 1)
            {
                goto L140;
            }
            /*
                 this is the only general constraint in the working set
            */
            var asize_c = asize;
            cond = dprotdiv(&asize_c, &dtnew, &ifail);
            asize = asize_c;
            if (ifail != 0 && asize == 0)
            {
                cond = BlasLike.lm_max;
            }
            if (cond >= *condmx)
            {
                goto L480;
            }
            if (cond >= condbd && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n", CL(*jadd));
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
            cond = dprotdiv(&tdtmax, &tdtmin, &ifail);
            if (ifail != 0 && tdtmax == 0)
            {
                cond = BlasLike.lm_max;
            }
            if (cond >= *condmx)
            {
                goto L480;
            }
            if (cond >= condbd && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n", CL(*jadd));
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
            if (modfyr || *unitq!=0)
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
            detagen(CN(ncolz1), &wrk1[*ncolz], &wrk1[1], CI(1), &iswap, &itrans)
                ;
            if (iswap > 0)
                delm(orthog, CN(*nfree), &zy[*ncolz * zy_dim1 + 1], 1, &zy[
                    iswap * zy_dim1 + 1], 1, zero, zero);
            if (itrans == 0)
            {
                goto L220;
            }
            i__1 = ncolz1;
            for (j = 1; j <= i__1; ++j)
            {
                d = wrk1[j];
                if (d == 0)
                {
                    goto L200;
                }
                BlasLike.daxpy(CN(*nfree), d, &zy[*ncolz * zy_dim1 + 1], 1, &zy[j * zy_dim1 + 1], 1);
            L200:
                ;
            }
        L220:
            if (!modfyg)
            {
                goto L360;
            }
            if (iswap > 0)
                delm(orthog, 1, &qtg[*ncolz], 1, &qtg[iswap], 1, zero, zero);
            if (itrans > 0)
            {
                BlasLike.daxpy(ncolz1, qtg[*ncolz], &wrk1[1], 1, &qtg[1], 1);
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
            delta = wrk1[*ncolz];
            var c__1 = 1;
            var lm_eps = BlasLike.lm_eps;
            dhhrflctgen(&ncolz1, &delta, &wrk1[1], &c__1, &lm_eps, &beta);
            if (beta != 0)
            {
                wrk1[*ncolz] = beta;
            }
            if (beta <= 0)
            {
                goto L360;
            }
            BlasLike.dzerovec(*nfree, &wrk2[1]);
            i__1 = *ncolz;
            for (j = 1; j <= i__1; ++j)
            {
                d = wrk1[j];
                if (d == 0)
                {
                    goto L260;
                }
                BlasLike.daxpy(*nfree, d, &zy[j * zy_dim1 + 1], 1, &wrk2[1], 1);
            L260:
                ;
            }
            i__1 = *ncolz;
            for (j = 1; j <= i__1; ++j)
            {
                d = wrk1[j];
                if (d == 0)
                {
                    goto L280;
                }
                d = -d / beta;
                BlasLike.daxpy(*nfree, d, &wrk2[1], 1, &zy[j * zy_dim1 + 1], 1);
            L280:
                ;
            }
            if (!modfyg)
            {
                goto L300;
            }
            d = BlasLike.ddotvec(*ncolz, &wrk1[1], &qtg[1]);
            d = -d / beta;
            BlasLike.daxpy(*ncolz, d, &wrk1[1], 1, &qtg[1], 1);
        L300:
            wrk1[*ncolz] = delta;
            goto L360;

        /*r  has to be modified.  use a sequence of 2*2 transformations*/
        L320:
            lrowr = *ncolr;
            i__1 = ncolz1;
            for (k = 1; k <= i__1; ++k)
            {
                /*
                compute the transformation that reduces wrk1(k) to zero,
                then apply it to the relevant columns of  z  and  grad(t)q.
                */
                kp1 = k + 1;
                delmgen(orthog, &wrk1[kp1], &wrk1[k], &cs, &sn);
                if ((*unitq)==0)
                {
                    delm(orthog, CN(*nfree), &zy[kp1 * zy_dim1 + 1], 1, &zy[k *
                        zy_dim1 + 1], 1, cs, sn);
                }
                if (modfyg)
                    delm(orthog, 1, &qtg[kp1], 1, &qtg[k], 1, cs, sn);
                /*
                apply the same transformation to the cols of  r  if relevant
                this generates a subdiagonal element in  r  which must be
                eliminated by a row rotation.  the last such row rotation
                is needed by  qpcore
                */
                if (!(modfyr && k < *ncolr))
                {
                    goto L340;
                }
                rt[kp1 + k * rt_dim1] = 0;
                delm(orthog, CN(kp1), &rt[kp1 * rt_dim1 + 1], 1, &rt[k *
                    rt_dim1 + 1], 1, cs, sn);
                BlasLike.drotg(&rt[k + k * rt_dim1], &rt[kp1 + k * rt_dim1], cslast, snlast);

                rt[kp1 + k * rt_dim1] = 0;
                --lrowr;
                dsymplanerotate(CN(lrowr), &rt[k + kp1 * rt_dim1],
                        CI(*Nrowrt), &rt[kp1 + kp1 * rt_dim1],
                        CI(*Nrowrt), *cslast, *snlast);
            L340:
                ;
            }
        /* if adding a general constraint, insert the new row of  t  and exit */
        L360:
            if (*jadd > (int)n)
            {
                BlasLike.dcopy(nact1, &wrk1[*ncolz], 1, &rt[nact1 + *ncolz * rt_dim1], *Nrowrt);
                return 0;
            }
            /*
                 we are adding a bound.  continue reducing the elements of  wrk1
                 to zero.  this affects  y,  t  and  qtg
                 first, set the super-diagonal elements of t to zero
            */
            if (*nactiv == 0)
            {
                goto L440;
            }
            BlasLike.dzero(*nactiv, &rt[*nactiv + *ncolz * rt_dim1], *Nrowrt - 1);
            nelm = 1;
            ldiag = *nactiv;
            i__1 = nfree1;
            for (k = *ncolz; k <= i__1; ++k)
            {
                delmgen(orthog, &wrk1[k + 1], &wrk1[k], &cs, &sn);
                delm(orthog, CN(*nfree), &zy[(k + 1) * zy_dim1 + 1], 1, &zy[k *
                    zy_dim1 + 1], 1, cs, sn);
                delm(orthog, CN(nelm), &rt[ldiag + (k + 1) * rt_dim1], 1,
                    &rt[ldiag + k * rt_dim1], 1, cs, sn);
                if (modfyg)
                    delm(orthog, 1, &qtg[k + 1], 1, &qtg[k], 1, cs, sn);
                ++nelm;
                --ldiag;
            }
            /*
            the diagonals of  t  have been altered.  recompute the largest and
            smallest values
            */
            inct = *Nrowrt - 1;
            var dtmax_c = dtmax;
            var dtmin_c = dtmin;
            BlasLike.dxminmax(*nactiv, &rt[*nactiv + (ncolz1 + 1) * rt_dim1], inct, &dtmax_c, &
                dtmin_c);
            dtmax = dtmax_c;
            dtmin = dtmin_c;
            if (dtmin / dtmax * *condmx < 1.0) goto L480;
            if (dtmin / dtmax * condbd < 1.0 && msg >= 0)
                lm_wmsg("\n*** WARNING\n *** SERIOUS ILL-CONDITIONING IN THE WORKING SET AFTER ADDING CONSTRAINT %5ld\n *** OVERFLOW MAY OCCUR IN SUBSEQUENT ITERATIONS\n\n",
                CL(*jadd));
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
                qtg[*nfree] /= wrk1[*nfree];
            }
            /*the factorization has been successfully updated*/
            return 0;

        /*THE PROPOSED WORKING SET APPEARS TO BE LINEARLY DEPENDENT*/
        L480:
            if (msg >= 80)
            {
                Console.WriteLine("\n//ADDCON//  DEPENDENT CONSTRAINT REJECTED");
                if (*jadd <= (int)n)
                    lm_wmsg(
                    "\n//ADDCON//     ASIZE     DTMAX     DTMIN\n//ADDCON//%10.2le%10.2le%10.2le",
                        asize, dtmax, dtmin);
                else lm_wmsg(
             "\n//ADDCON//     ASIZE     DTMAX     DTMIN     DTNEW\n//ADDCON//%10.2le%10.2le%10.2le",
                 asize, dtmax, dtmin, dtnew);
            }
            return 1;
        }
        public unsafe static
        void dprtsol(int nfree, int nrowa, int n, int nclin, int ncnln, int nctotl,
                    double bigbnd, int nactiv, int* istate, int* kactiv,
                    double* a, double* bl, double* bu, double* c, double* clamda, double* rlamda, double* x)
        {

            /*
                dprtsol	expands the lagrange multipliers into  clamda
                if  msg >= 10  or  msg == 1,  prtsol  then prints  x, a*x,
                c(x), their bounds,  the multipliers, and the residuals
                (distance to the nearest bound)

                dprtsol	is called by  lpcore, qpcore, lccore and npcore	 just
                before they exit
            */
            string id = "VLN";
            string lstate = "--++FRLLULEQTB";


            int nlam, nplin, nfixed, j, k, ip, was_is;
            double v, b1, b2, res, res2, wlam;
            string ls;
            char[] id3 = new char[1];

            --x;
            --rlamda;
            --clamda;
            --c;
            --bu;
            --bl;
            a -= nrowa + 1;
            --kactiv;
            --istate;

            nplin = n + nclin;
            /*EXPAND BOUND, LINEAR AND NONLINEAR MULTIPLIERS INTO CLAMDA*/
            BlasLike.dzerovec(nctotl, &clamda[1]);
            nfixed = n - nfree;
            nlam = nactiv + nfixed;
            for (k = 1; k <= nlam; ++k)
            {
                j = kactiv[k];
                if (k <= nactiv) j += n;
                clamda[j] = rlamda[k];
            }
            if (msg < 10 && msg != 1) return;
            Console.WriteLine("\n\nVARBL STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
            id3[0] = id[0];
            for (j = 1; j <= nctotl; ++j)
            {
                b1 = bl[j];
                b2 = bu[j];
                wlam = clamda[j];
                was_is = istate[j];
                ls = lstate.Substring(((was_is + 2) << 1)); //IS THIS RIGHT
                                                            //		ls = lstate + ((was_is + 2) << 1);
                if (j <= n)
                {
                    /* SECTION 1 -- THE VARIABLES  X. */
                    /* ------------------------------ */
                    k = j;
                    v = x[j];
                }
                else if (j <= nplin)
                {
                    /* SECTION 2 -- THE LINEAR CONSTRAINTS  A*X. */
                    /* ----------------------------------------- */
                    if (j == n + 1)
                    {
                        Console.WriteLine("\n\nLNCON STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
                        id3[0] = id[1];
                    }
                    k = j - n;
                    v = BlasLike.ddot(n, &a[k + nrowa], nrowa, &x[1], 1);
                }
                else
                {
                    /*        SECTION 3 -- THE NONLINEAR CONSTRAINTS  C(X). */
                    /*        --------------------------------------------- */
                    if (ncnln <= 0) continue;
                    if (j == nplin + 1)
                    {
                        Console.WriteLine("\n\nNLCON STATE     VALUE      LOWER BOUND    UPPER BOUND    LAGR MULT   RESIDUAL");
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
                Console.WriteLine($"{id3[0]},{CL(k)},{ls[0]},{ls[1]},{v}");
                if (ip != 0) Console.WriteLine($"{b1}");
                else Console.WriteLine("     NONE      ");
                if (ip < 3) Console.WriteLine($"{b2}");
                else Console.WriteLine("     NONE      ");
                Console.WriteLine($"{wlam},{res}");
            }
        }
        public unsafe static double dnrm2(int n, double* x, int incx)
        {

            if (n == 1) return (*x < 0.0 ? -*x : *x);
            else
            {
                double scale = 0.0, ssq = 1.0;
                if (incx == 1) dsssqvec(n, x, &scale, &ssq);
                else BlasLike.dsssq(n, x, incx, &scale, &ssq);
                return sc_norm(scale, ssq);
            }
        }
        public unsafe static void dtqadd(bool orthog, int* unitq, int* inform, int* k1, int* k2, int* nactiv, int* ncolz, int* nfree, int* n, int* nq, int* nrowa, int* nrowrt, int* ncolrt, int* istate, int* kactiv, int* kfree, double* condmx, double* a, double* qtg, double* rt, double* zy, double* wrk1, double* wrk2)
        {
            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1;

            int iadd, jadd, ifix, i, k, l, iswap;
            double cslast, snlast;

            /*     TQADD INCLUDES GENERAL LINEAR CONSTRAINTS  K1  THRU  K2  AS NEW */
            /*     COLUMNS OF THE TQ FACTORIZATION STORED IN  RT, ZY. */
            --wrk2;
            --wrk1;
            zy_dim1 = *nq;
            zy_offset = zy_dim1 + 1;
            zy -= zy_offset;
            rt_dim1 = *nrowrt;
            rt_offset = rt_dim1 + 1;
            rt -= rt_offset;
            --qtg;
            a_dim1 = *nrowa;
            a_offset = a_dim1 + 1;
            a -= a_offset;
            --kfree;
            --kactiv;
            --istate;

            i__1 = *k2;
            for (k = *k1; k <= i__1; ++k)
            {
                iadd = kactiv[k];
                jadd = *n + iadd;
                if (*nactiv == *nfree)
                {
                    goto L20;
                }
                *inform = daddcon(false, false, orthog, unitq, &ifix, &iadd, &jadd,
                    nactiv, ncolz, ncolz, nfree, CN(*n), nq, nrowa, nrowrt, &
                    kfree[1], condmx, &cslast, &snlast, &a[a_offset], &qtg[1], &
                    rt[rt_offset], &zy[zy_offset], &wrk1[1], &wrk2[1]);
                if (*inform > 0)
                {
                    goto L20;
                }
                ++(*nactiv);
                --(*ncolz);
                goto L40;
            L20:
                istate[jadd] = 0;
                kactiv[k] = -kactiv[k];
            L40:
                ;
            }
            if (*nactiv == *k2)
            {
                return;
            }
            /*     SOME OF THE CONSTRAINTS WERE CLASSED AS DEPENDENT AND NOT INCLUDED 
            */
            /*     IN THE FACTORIZATION.  MOVE ACCEPTED INDICES TO THE FRONT OF */
            /*     KACTIV  AND SHIFT REJECTED INDICES (WITH NEGATIVE VALUES) TO */
            /*     THE END. */
            l = *k1 - 1;
            i__1 = *k2;
            for (k = *k1; k <= i__1; ++k)
            {
                i = kactiv[k];
                if (i < 0)
                {
                    goto L60;
                }
                ++l;
                if (l == k)
                {
                    goto L60;
                }
                iswap = kactiv[l];
                kactiv[l] = i;
                kactiv[k] = iswap;
            L60:
                ;
            }
        }
        public unsafe static void drtmxsolve(int job, int n, double* t, int nrt, double* b, int* idiag)
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

            int fail;
            double temp;
            int j, k;
            int jlast;

            --b;
            t_dim1 = nrt;
            t_offset = t_dim1 + 1;
            t -= t_offset;

            if (n >= 1 && nrt >= n && Math.Abs(job) <= 2)
            {
                goto L20;
            }
            *idiag = (int)lm_check_fail((short)(*idiag), CS((short)1), "drtmxsolve");
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
                b[j] = dprotdiv(&b[j], &t[j + k * t_dim1], &fail);
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
                b[j] = dprotdiv(&b[j], &t[j + k * t_dim1], &fail);
                if (fail != 0)
                {
                    goto L280;
                }
                if (j > 1)
                {
                    BlasLike.daxpy(j - 1, -b[j], &t[k * t_dim1 + 1], 1, &b[1], 1);
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
                b[j] = dprotdiv(&b[j], &t[j + k * t_dim1], &fail);
                if (fail != 0)
                {
                    goto L280;
                }
                if (j < (int)n)
                {
                    BlasLike.daxpy(n - j, -b[j], &t[j + 1 + k * t_dim1], 1, &b[j + 1], 1);
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
                temp = b[j];
                b[j] = b[k];
                b[k] = temp;
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
                if (j > 1) b[j] -= BlasLike.ddotvec((int)(j - 1), &t[k * t_dim1 + 1], &b[1]);
                b[j] = dprotdiv(&b[j], &t[j + k * t_dim1], &fail);
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
                if (j < (int)n) b[j] -= BlasLike.ddotvec((int)(n - j), &t[j + 1 + k * t_dim1], &b[j + 1]);
                b[j] = dprotdiv(&b[j], &t[j + k * t_dim1], &fail);
                if (fail != 0)
                {
                    goto L280;
                }
                --j;
                /* L240: */
            }
        L260:
            *idiag = 0;
            return;
        L280:
            *idiag = (int)lm_check_fail((short)(*idiag), (short)(j + 1), "drtmxsolve");
        }
        public static short lm_check_fail(short ifail, short ierror, string srname)
        {
            if (ierror != 0)
            {
                /*Abnormal exit from calling routine */
                if (ifail <= 0)
                {
                    /*Noisy exit */
                    Console.WriteLine($"Mathematics routine {srname}: exited with ifail={ierror}",
                    srname, ierror);
                    if (ifail != 0) { Console.WriteLine("Hard failure"); return -50; }
                }
                //#if 0
                else    /*Soft failure*/
                    Console.WriteLine(" ** Soft failure - control returned");
                //#endif
            }
            return ierror;
        }
        public unsafe static

void delmgen(bool orthog, double* x, double* y, double* cs, double* sn)
        {/*
	If orthog delmgen generates a plane rotation.  Otherwise,
	delmgen  generates an elimination transformation  E  such that
	(X Y)*E  =  (X  0)   OR   (Y  0),  depending on the relative
	sizes of  X  &  Y.
*/
            if (orthog) BlasLike.drotg(x, y, cs, sn);
            else
            {
                *cs = 1;
                *sn = 0;
                if (*y != 0)
                {
                    if (Math.Abs(*x) < Math.Abs(*y))
                    {
                        *cs = 0;
                        *sn = -(*x) / *y;
                        *x = *y;
                    }
                    else *sn = -(*y) / *x;
                }
            }
            *y = 0;
        }
        public unsafe static void delm(bool orthog, int n, double* x, int incx, double* y, int incy, double cs, double sn)
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
                if (cs <= 0) BlasLike.dswap(n, x, incx, y, incy);
                if (sn != 0) BlasLike.daxpy(n, sn, x, incx, y, incy);
            }
            else dsymplanerotate(n, x, incx, y, incy, cs, sn);
        }
        public unsafe static void dsymplanerotate(int n, double* x, int incx, double* y, int incy, double c, double s)
        {/*
	dsymplanerotate performs the symmetric plane rotation
	( x  y ) = ( x  y )*( c   s )   s != 0
                    	    ( s  -c )

	If s is supplied as zero then x and y are unaltered
*/
            int i__1, i__2;

            double temp1;
            int i, ix, iy;
            --y;
            --x;

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
                            temp1 = x[ix];
                            x[ix] = y[ix];
                            y[ix] = temp1;
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
                                temp1 = x[ix];
                                x[ix] = y[iy];
                                y[iy] = temp1;
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[ix];
                                x[ix] = y[iy];
                                y[iy] = temp1;
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
                            temp1 = -x[ix];
                            x[ix] = -y[ix];
                            y[ix] = temp1;
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
                                temp1 = -x[ix];
                                x[ix] = -y[iy];
                                y[iy] = temp1;
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = -x[ix];
                                x[ix] = -y[iy];
                                y[iy] = temp1;
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
                            temp1 = x[ix];
                            x[ix] = c * temp1 + s * y[ix];
                            y[ix] = s * temp1 - c * y[ix];
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
                                temp1 = x[ix];
                                x[ix] = c * temp1 + s * y[iy];
                                y[iy] = s * temp1 - c * y[iy];
                                iy += incy;
                            }
                        }
                        else
                        {
                            ix = 1 - (n - 1) * incx;
                            i__1 = n;
                            for (i = 1; i <= i__1; ++i)
                            {
                                temp1 = x[ix];
                                x[ix] = c * temp1 + s * y[iy];
                                y[iy] = s * temp1 - c * y[iy];
                                ix += incx;
                                iy += incy;
                            }
                        }
                    }
                }
            }
        }
        public unsafe static short dtmxsolve(short job, int n, double* t, int nrt, double* b, short idiag)
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

            --b;
            t -= nrt + 1;

            if (n < 1 || nrt < n || Math.Abs(job) > 2) return lm_check_fail(idiag, CS((short)1), "dtmxsolve");

            if (job == 1 || job == -2)
            {
                for (k = n; k >= 1; --k) if ((fail = (b[k] != 0.0 ? 1 : 0)) != 0) break;
            }
            else
            {
                for (k = 1; k <= n; ++k) if ((fail = (b[k] != 0.0 ? 1 : 0)) != 0) break;
            }

            if (fail != 0)
                switch (job)
                {
                    case 0:
                        for (; k <= n; ++k)
                        {
                            b[k] = dprotdiv(&b[k], &t[k + k * nrt], &fail);
                            if (fail != 0) break;
                        }
                        break;
                    case 1:
                        for (; k >= 1; --k)
                        {
                            b[k] = dprotdiv(&b[k], &t[k + k * nrt], &fail);
                            if (fail != 0) break;
                            if (k > 1) BlasLike.daxpy(k - 1, -b[k], &t[k * nrt + 1], 1, &b[1], 1);
                        }
                        break;
                    case 2:
                        for (; k <= n; ++k)
                        {
                            b[k] = dprotdiv(&b[k], &t[k + k * nrt], &fail);
                            if (fail != 0) break;
                            if (k < n) BlasLike.daxpy(n - k, -b[k], &t[k + 1 + k * nrt], 1, &b[k + 1], 1);
                        }
                        break;
                    case -1:
                        for (; k <= n; ++k)
                        {
                            if (k > 1) b[k] -= BlasLike.ddotvec((int)(k - 1), &t[k * nrt + 1], &b[1]);
                            b[k] = dprotdiv(&b[k], &t[k + k * nrt], &fail);
                            if (fail != 0) break;
                        }
                        break;
                    case -2:
                        for (; k >= 1; --k)
                        {
                            if (k < n) b[k] -= BlasLike.ddotvec((int)(n - k), &t[k + 1 + k * nrt], &b[k + 1]);
                            b[k] = dprotdiv(&b[k], &t[k + k * nrt], &fail);
                            if (fail != 0) break;
                        }
                        break;
                }

            return (short)(fail != 0 ? lm_check_fail(idiag, (short)(k + 1), "dtmxsolve") : 0);
        }
        public unsafe static void dbdpert(bool firstv, bool negstp, double* bigalf, double* bigbnd, double* pnorm, int* jadd1, int* jadd2, double* palfa1, double* palfa2, int* istate, int n, int nctotl, double* anorm, double* ap, double* ax, double* bl, double* bu, double* featol, double* p, double* x)
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
            int i, j;
            bool lastv;
            double epspt9;
            int js;
            double absatp, rownrm, atp, res, atx;
            --x;
            --p;
            --featol;
            --bu;
            --bl;
            --ax;
            --ap;
            --anorm;
            --istate;

            epspt9 = parm[3];
            if (msg == 99) Console.WriteLine(
        "\n   J  JS         FEATOL         AX             AP     JADD1       PALFA1     JADD2       PALFA2\n");
            lastv = !firstv;
            *jadd1 = 0;
            *jadd2 = 0;
            *palfa1 = *bigalf;
            *palfa2 = 0.0;
            if (firstv) *palfa2 = *bigalf;
            for (j = 1; j <= nctotl; ++j)
            {
                js = istate[j];
                if (js > 0) continue;
                if (j > n)
                {
                    /*GENERAL LINEAR CONSTRAINT. */
                    i = j - n;
                    atx = ax[i];
                    atp = ap[i];
                    rownrm = 1.0 + anorm[i];
                }
                else
                {
                    /*BOUND CONSTRAINT. */
                    atx = x[j];
                    atp = p[j];
                    rownrm = 1.0;
                }
                if (negstp) atp = -atp;
                if (Math.Abs(atp) <= epspt9 * rownrm * *pnorm) res = -1.0;
                else if (atp > BlasLike.lm_eps)
                {
                    /*AX IS INCREASING*/
                    /*TEST FOR SMALLER PALFA1 IF UPPER BOUND IS SATISFIED. */
                    if (js != -1)
                    {
                        if (bu[j] < *bigbnd)
                        {
                            res = bu[j] - atx + featol[j];
                            if (*bigalf * atp > Math.Abs(res) && *palfa1 * atp > res)
                            {
                                *palfa1 = res / atp;
                                *jadd1 = j;
                            }
                        }
                        /*TEST FOR DIFFERENT PALFA2 IF LOWER BOUND IS VIOLATED. */
                        if (js == -2)
                        {
                            res = bl[j] - atx - featol[j];
                            if (*bigalf * atp > Math.Abs(res)
                                && (firstv || *palfa2 * atp < res)
                                && (lastv || *palfa2 * atp > res))
                            {
                                *palfa2 = res / atp;
                                *jadd2 = j;
                            }
                        }
                    }
                }
                else if (js != -2)
                {
                    /*AX IS DECREASING TEST FOR SMALLER PALFA1 IF LOWER BOUND IS SATISFIED*/
                    absatp = -atp;
                    if (bl[j] > -(*bigbnd))
                    {
                        res = atx - bl[j] + featol[j];
                        if (*bigalf * absatp > Math.Abs(res) && *palfa1 * absatp > res)
                        {
                            *palfa1 = res / absatp;
                            *jadd1 = j;
                        }
                    }
                    /*TEST FOR DIFFERENT PALFA2 IF UPPER BOUND IS VIOLATED*/
                    if (js == -1)
                    {
                        res = atx - bu[j] - featol[j];
                        if (*bigalf * absatp > Math.Abs(res)
                            && (firstv || *palfa2 * absatp < res)
                            && (lastv || *palfa2 * absatp > res))
                        {
                            *palfa2 = res / absatp;
                            *jadd2 = j;
                        }
                    }
                }
                if (msg == 99)
                    lm_wmsg("%5ld%4ld%15.5lg%15.5lg%15.5lg%6ld%17.7lg%6ld%17.7lg",
                    CL(j), CL(js), featol[j], atx, atp,
                    CL(*jadd1), *palfa1, CL(*jadd2), *palfa2);
            }
        }
        public unsafe static void detagen(int n, double* alpha, double* x, int incx, int* iswap, int* itrans)
        {
            /*
                detagen  generates an elimination transformation  e  such that
                    e ( alpha )  =  ( delta ) ,
                      (   x   )     (   0   )

                where  e  has the form
                    e  =	( 1    ) p
                        ( z  i )

                for some n-vector  z  and permutation matrix  p  of order  n + 1.
                in certain circumstances ( x  very small in absolute terms or
                x very small compared to  alpha),  e  will be the identity matrix.
                detagen  will then leave  alpha  and  x  unaltered, and will return
                iswap = 0,  itrans = 0

                more generally,  iswap  and  itrans  indicate the various possible
                forms of  p  and  z  as follows
                    if  iswap  =  0,  p = i
                    if  iswap  gt 0,  p  interchanges  alpha  and  x(iswap)
                    if  itrans =  0,  z = 0  and the transformation is just  e = p
                    if  itrans gt 0,  z  is nonzero.  its elements are returned in  x.

                detagen  guards against overflow and underflow
                it is assumed that  flmin < epsmch**2 (i.e. rtmin < epsmch).
            */
            int imax = 1000000000;
            int nzero;
            double xmax, absalf, tol, axi;
            double* v, vlim;

            *iswap = 0;
            *itrans = 0;
            if (n < 1) return;
            absalf = Math.Abs(*alpha);
            xmax = 0;

            for (v = x, vlim = x + n * incx; v != vlim; v += incx)
                if (xmax < (axi = Math.Abs(*v)))
                {
                    xmax = axi;
                    imax = (int)(v - x);
                }

            /* exit if  x  is very small */
            if (xmax <= BlasLike.lm_rootmin) return;

            /* see if an interchange is needed for stability */
            if (absalf < xmax)
            {
                *iswap = imax + 1;
                xmax = x[imax];
                x[imax] = *alpha;
                *alpha = xmax;
            }

            /*
                 form the multipliers in  x.  they will be no greater than one
                 in magnitude.  change negligible multipliers to zero
            */
            tol = Math.Abs(*alpha) * BlasLike.lm_eps;
            nzero = 0;
            for (v = x; v != vlim; v += incx)
                if (Math.Abs(*v) > tol) *v = -*v / *alpha;
                else
                {
                    *v = 0;
                    ++nzero;
                }
            /*z is zero only if nzero=n*/
            if (nzero < n) *itrans = 1;
        }
        public unsafe static void dhhrflctgen(int* n, double* alpha, double* x, int* incx, double* tol, double* z1)
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
            double beta, work; //work was work[1] makes no sense!
            double scale, tl, ssq;

            --x;

            if (*n < 1)
            {
                *z1 = 0;
            }
            else
            {
                if (*tol <= 0 || *tol > 1)
                {
                    tl = 0;
                }
                else
                {
                    tl = Math.Abs(*alpha) * *tol;
                }
                ssq = 1;
                scale = 0;
                BlasLike.dsssq(*n, &x[1], *incx, &scale, &ssq);
                if (scale == 0 || scale < tl)
                {
                    *z1 = 0;
                }
                else
                {
                    if (*alpha != 0)
                    {
                        work = *alpha;//I MADE CHANGES  TO WORK; compare with libsafeqp
                        BlasLike.dsssqvec(1, &work, &scale, &ssq);
                        beta = -BlasLike.dsign(sc_norm(scale, ssq), *alpha);
                        *z1 = (beta - *alpha) / beta;
                    }
                    else
                    {
                        beta = -sc_norm(scale, ssq);
                        *z1 = 1;
                    }
                    BlasLike.dscal(*n, -1 / beta, &x[1], *incx);
                    *alpha = beta;
                }
            }
        }
        public unsafe static short dqpsol(short itmax, short msglvl, int n, int nclin, int nctotl, int nrowa, int nrowh, int ncolh, double* bigbnd, double* a, double* bl, double* bu, double* cvec, double* featol, double* hess, int cold, int lp, int orthog, double* x, int* istate, short* iter, double* obj, double* clamda, int* iw, int leniw, double* w, int lenw, short ifail)
        {

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
            int nfree, ncnln;
            int unitq;
            double xnorm;// epspt9;
            int lscale, maxact, minact;
            byte lcrash;
            //double tolact;
            int minfxd, inform, mxfree, nactiv, numinf;
            int litotl;
            short nerror;
            int minsum;
            int mxcolz;
            int vertex;
            int lwtotl;
            int lax;

            //#define NCLIN &nclin_
            int nclin_ = nclin;
            //#define NCTOTL &nctotl_
            int nctotl_ = nctotl;
            //#define NROWA &nrowa_
            int nrowa_ = nrowa;
            int iter_;
            //#define NROWH &nrowh_
            int nrowh_ = nrowh;
            //#define NCOLH &ncolh_
            int ncolh_ = ncolh;
            --w;
            --iw;
            --clamda;
            --istate;
            --x;
            --featol;
            --cvec;
            --bu;
            --bl;

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
            nrowrt = Math.Max(mxcolz, maxact);
            ncolrt = Math.Max(1, mxfree);
            ncnln = 0;

            /*allocate certain arrays that are not done in  alloc*/
            litotl = 0;
            lax = 1;
            lwtotl = lax + nrowa - 1;
            /*allocate remaining work arrays*/
            dalloc(2, n, nclin, ncnln, nctotl, iw, w, &litotl, &lwtotl);
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
                vertex, &istate[1], a, &w[lax], &bl[1], &bu[1],
                &cvec[1], &x[1]);
            if (msglvl == 99)
                dqpdump(n, nrowh, ncolh, &cvec[1], hess, pWRK, pPX);

            nerror = dchkdat(leniw, lenw, litotl, lwtotl, nrowa, n, nclin,
                    nctotl, &istate[1], pKACTV, lcrash,
                    bigbnd, a, &bl[1], &bu[1], &featol[1], &x[1]);
            *iter = 0;
            if (nerror != 0)
            {
                if (msglvl > 0)
                    lm_wmsg(
            " EXIT QPSOL-%6ld ERRORS FOUND IN THE INPUT PARAMETERS.  PROBLEM ABANDONED.",
                    CL(nerror));
                return lm_check_fail((short)CS(ifail), CS((short)9), "QPSOL");
            }

            /*
            no scaling is provided by this version of  dqpsol
            give a fake value for the start of the scale array
            */
            scldqp = false;
            lscale = 1;

            /*
            ---------------------------------------------------------------------
            call  lpcore  to obtain a feasible point, or solve a linear
            problem
            ---------------------------------------------------------------------
            */
            dlpcore((lp & 1) != 0, minsum, orthog != 0, &unitq, vertex, &inform, &iter_,
            itmx, lcrash, n, &nclin, &nctotl, &nrowa, &nactiv, &nfree, &numinf,
                &istate[1], pKACTV, pKFREE, obj, &xnorm, a, &
                w[lax], &bl[1], &bu[1], &clamda[1], &cvec[1], &featol[1], &x[1], &
                iw[1], &w[1]);
            *iter = (short)iter_;
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
                var nrowrt_c = nrowrt;
                dqpcore(orthog, &unitq, &inform, &iter_, &itmx, n, &nclin, &nctotl,
                        &nrowrt_c, &nrowh, &ncolh, &nactiv, &nfree, &istate[1],
                        pKACTV, pKFREE, obj, &xnorm, a, &w[lax], &bl[1], &bu[1], &clamda[1], &cvec[1], &featol[1], hess,
                        &w[lscale], &x[1], &iw[1], &w[1]);
                nrowrt = nrowrt_c;
                *iter = (short)iter_;
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
                        Console.WriteLine("WEAK LOCAL MINIMUM.");
                        break;
                    case 2:
                        lm_wmsg("\nEXIT dqpsol- %.2s SOLUTION IS UNBOUNDED.", l);
                        break;
                    case 3:
                        Console.WriteLine("ZERO MULTIPLIERS.");
                        break;
                    case 4:
                        Console.WriteLine("TOO MANY ITERATIONS WITHOUT CHANGING X.");
                        break;
                    case 5:
                        Console.WriteLine("TOO MANY ITERATIONS.");
                        break;
                    case 6:
                        Console.WriteLine("CANNOT SATISFY THE LINEAR CONSTRAINTS.");
                        break;
                    case 7:
                        Console.WriteLine("TOO MANY ITERATIONS WITHOUT CHANGING X IN THE LP PHASE.");
                        break;
                    case 8:
                        Console.WriteLine("TOO MANY ITERATIONS DURING THE LP PHASE.");
                        break;
                }
                if (numinf == 0) lm_wmsg("\n FINAL %.2s OBJECTIVE VALUE =%20.9lg", l, *obj);
                else lm_wmsg("\n FINAL SUM OF INFEASIBILITIES =%20.9lg", *obj);
            }

            return (short)(inform == 0 ? 0 : lm_check_fail((short)CS(ifail), (short)CS(inform), "QPSOL"));
        }
        public unsafe static void dlpdump(int n, int nclin, int nctotl, int nrowa, int lcrash, int lp, int minsum, int vertex, int* istate, double* a, double* ax, double* bl, double* bu, double* cvec, double* x)
        {
            int j, k;
            double atx;

            Console.WriteLine("\n\n\n\n\n\nOUTPUT FROM LPDUMP\n ******************");
            lm_wmsg("\nLCRASH =%d LP=%d MINSUM=%d VERTEX=%d", lcrash, lp, minsum, vertex);

            /*PRINT  A  BY ROWS AND COMPUTE  AX = A*X. */
            for (k = 0; k < nclin; ++k)
            {
                lm_wmsg("\nROW%6ld OF A ...", CL(k + 1));
                lm_gdvwri(n, a + k, nrowa);
                ax[k] = BlasLike.ddot(n, a + k, nrowa, x, 1);
            }

            /*PRINT  BL, BU  AND  X OR AX. */
            Console.WriteLine("\n              J      BL(J)          BU(J)           X(J)");
            for (j = 0; j < nctotl; ++j)
            {
                if (j < n)
                {
                    k = j;
                    atx = x[j];
                }
                else
                {
                    k = j - n;
                    atx = ax[k];
                    if (k != 0)
                        Console.WriteLine("\n              I    BL(N+I)        BU(N+I)         A(I)*X");
                }
                lm_wmsg("", CL(k + 1), bl[j], bu[j], atx);
            }
            if ((lp & 1) != 0) lm_mdvwri("\nCVEC ...", n, cvec);
            if (lcrash != 0) lm_mdvwri("\nISTATE ...", nctotl, &istate[1]);
        }
        public unsafe static void lm_gdvwri(int n, double* x, int inc)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                Console.Write($"{*x} ");
                x += inc;
                if (i % 6 == 5) Console.Write("\n");
            }
            Console.Write("\n");
        }
        public unsafe static void set_addr(int i, int l, int* iloc, void* a, int sa, void** ploc)
        {
            while (i <= l)
            {
                ploc[i] = ((byte*)a) + iloc[i] * sa;
                i++;
            }
        }
        public unsafe static void dalloc(byte nalg, int n, int nclin, int ncnln, int nctotl, int* iw, double* w, int* litotl, int* lwtotl)
        {
            int ladx, laqp, lrho, lslk, lqtg, lwrk, lcsl1, lmax1, lmax2, lslk1,
                lztg2, ldlam, lcjdx, lrlam, ldslk, lxbwd, lxfwd, lqpdx, lgrad2,
                lkfree, lcslam, lsigma, lshare, lenaqp = 1000000000, lkactv, lanorm, lqpadx,
                lg1, lg2, lqptol, liqpst, lnpwrk, lqpwrk, lx1, lx2, lbl, lap, lbu,
                ldx, lrt, lpx, lzy, lcs1, lcs2;
            int* loclp = null;

            switch (nalg)
            {
                case 1:
                case 2:
                    /*
                    allocate the addresses for  lpcore  and  qpcore
                    */
                    lkactv = *litotl + 1;
                    lkfree = lkactv + n;
                    *litotl = lkfree + n - 1;
                    lanorm = *lwtotl + 1;
                    lap = lanorm + nclin;
                    lpx = lap + nclin;
                    lqtg = lpx + n;
                    lrlam = lqtg + n;
                    lrt = lrlam + n;
                    lzy = lrt + nrowrt * ncolrt;
                    lwrk = lzy + nq * nq;
                    *lwtotl = lwrk + n - 1;
                    loclp[0] = lkactv;
                    loclp[1] = lkfree;
                    loclp[2] = lanorm;
                    loclp[3] = lap;
                    loclp[4] = lpx;
                    loclp[5] = lqtg;
                    loclp[6] = lrlam;
                    loclp[7] = lrt;
                    loclp[8] = lzy;
                    loclp[9] = lwrk;
                    set_addr(0, 1, loclp, iw, sizeof(int), DAS_Sol_ploc);
                    set_addr(2, 9, loclp, w, sizeof(double), DAS_Sol_ploc);
                    break;
                case 3:
                    /*
                    allocate the addresses for npcore
                    */
                    lkactv = *litotl + 1;
                    lkfree = lkactv + n;
                    liqpst = lkfree + n;
                    *litotl = liqpst + nctotl - 1;
                    /*
                    variables used not only by  dnpcore,  but also dlpcore and  dqpcore
                    */
                    lanorm = *lwtotl + 1;
                    lqtg = lanorm + nrowqp;
                    lrlam = lqtg + n;
                    lrt = lrlam + n;
                    lzy = lrt + nrowrt * ncolrt;
                    loclp[1] = lkactv;
                    loclp[2] = lkfree;
                    loclp[3] = lanorm;
                    loclp[6] = lqtg;
                    loclp[7] = lrlam;
                    loclp[8] = lrt;
                    loclp[9] = lzy;
                    /*
                    assign the addresses for the workspace arrays used by  dnpiqp
                    */
                    lqpadx = lzy + nq * nq;
                    lqpdx = lqpadx + nrowqp;
                    lqpwrk = lqpdx + n;
                    loclp[4] = lqpadx;
                    loclp[5] = lqpdx;
                    loclp[10] = lqpwrk;
                    /*
                    assign the addresses for arrays used in  npcore
                    */
                    if (ncnln == 0) lenaqp = 0;
                    if (ncnln > 0) lenaqp = nrowqp * n;
                    laqp = lqpwrk + n;
                    ladx = laqp + lenaqp;
                    lbl = ladx + nrowqp;
                    lbu = lbl + nctotl;
                    ldx = lbu + nctotl;
                    lg1 = ldx + n;
                    lg2 = lg1 + n;
                    lqptol = lg2 + n;
                    lx1 = lqptol + nctotl;
                    lnpwrk = lx1 + n;
                    locnp[0] = liqpst;
                    locnp[1] = laqp;
                    locnp[2] = ladx;
                    locnp[3] = lbl;
                    locnp[4] = lbu;
                    locnp[5] = ldx;
                    locnp[6] = lg1;
                    locnp[7] = lg2;
                    locnp[8] = lqptol;
                    locnp[9] = lx1;
                    locnp[10] = lnpwrk;
                    lcs1 = lnpwrk + nctotl;
                    lcs2 = lcs1 + ncnln;
                    lcsl1 = lcs2 + ncnln;
                    lcslam = lcsl1 + ncnln;
                    lcjdx = lcslam + ncnln;
                    ldlam = lcjdx + ncnln;
                    ldslk = ldlam + ncnln;
                    lrho = ldslk + ncnln;
                    lsigma = lrho + ncnln;
                    lslk1 = lsigma + ncnln;
                    lslk = lslk1 + ncnln;
                    locnp[11] = lcs1;
                    locnp[12] = lcs2;
                    locnp[13] = lcsl1;
                    locnp[14] = lcslam;
                    locnp[15] = lcjdx;
                    locnp[16] = ldlam;
                    locnp[17] = ldslk;
                    locnp[18] = lrho;
                    locnp[19] = lsigma;
                    locnp[20] = lslk1;
                    locnp[21] = lslk;
                    *lwtotl = lslk + ncnln - 1;
                    break;
                case 4:
                    /*
                    allocate the addresses for  lccore
                    */
                    lkactv = *litotl + 1;
                    lkfree = lkactv + n;
                    *litotl = lkfree + n - 1;
                    lztg2 = *lwtotl + 1;
                    loclc[0] = lztg2;
                    /*
                    arrays used not only by  dlccore,  but also  dlpcore
                    */
                    lanorm = lztg2 + n;
                    lap = lanorm + nclin;
                    lpx = lap + nclin;
                    lqtg = lpx + n;
                    lrlam = lqtg + n;
                    lrt = lrlam + n;
                    lzy = lrt + nrowrt * ncolrt;
                    lwrk = lzy + nq * nq;
                    loclp[1] = lkactv;
                    loclp[2] = lkfree;
                    loclp[3] = lanorm;
                    loclp[4] = lap;
                    loclp[5] = lpx;
                    loclp[6] = lqtg;
                    loclp[7] = lrlam;
                    loclp[8] = lrt;
                    loclp[9] = lzy;
                    loclp[10] = lwrk;
                    lshare = lwrk + n;
                    /*
                    assign the addresses of the workspace used by  dlcsrch
                    this workspace is shared by  dlcappg
                    */
                    lx2 = lshare;
                    lgrad2 = lx2 + n;
                    lmax1 = lgrad2 + n - 1;
                    /*
                    assign the addresses of the workspace used by  dlcappg
                    this workspace is shared by  dlcsrch
                    */
                    lxfwd = lshare;
                    lxbwd = lxfwd + n;
                    lmax2 = lxbwd + n - 1;
                    *lwtotl = Math.Max(lmax1, lmax2);
                    loclc[1] = lx2;
                    loclc[2] = lgrad2;
                    loclc[3] = lxfwd;
                    loclc[4] = lxbwd;
                    break;
            }
        }
        public unsafe static void dqpdump(int n, int nrowh, int ncolh, double* cvec, double* hess, double* wrk, double* hx)
        {
            int j, i;

            Console.WriteLine("\n\n\n\n\n\nOUTPUT FROM QPDUMP\n******************");
            lm_mdvwri("\nCVEC ...", n, cvec);

            /*PRINT  HESS  UNLESS IT APPEARS TO BE IMPLICIT. */
            lm_wmsg("\nNROWH =%6ld NCOLH =%6ld", CL(nrowh), CL(ncolh));
            if (nrowh > 1 || ncolh > 1)
            {
                if (ncolh == 1) lm_mdvwri("\nHESS ...", nrowh, hess);
                else
                {
                    i = Math.Min(ncolh, n);
                    for (j = 1; j <= i; ++j)
                    {
                        lm_wmsg("\nCOLUMN%6ld OF  HESS ...", CL(j));
                        lm_gdvwri(i, hess + j - 1, CI(nrowh));
                    }
                }
            }
            /*CALL  QPHESS  TO COMPUTE EACH COLUMN OF THE HESSIAN. */
            Console.WriteLine("\n\n THE FOLLOWING IS RETURNED BY  QPHESS.");
            BlasLike.dzerovec(n, wrk);
            for (j = 1; j <= n; ++j)
            {
                wrk[i = j - 1] = 1.0;
                qphess(n, nrowh, ncolh, j, hess, wrk, hx);
                lm_wmsg("\nCOLUMN%6ld FROM  QPHESS ...", CL(j));
                lm_gdvwri(n, hx, 1);
                wrk[i] = 0.0;
            }
        }
        public unsafe static void qphess(int n, int nrowh, int ncolh, int j, double* hess, double* wrk, double* hx)
        {
            Solver.Factorise.dsmxmulv(n, hess, wrk, hx);
        }
        public unsafe static short dchkdat(int liwork, int lwork, int litotl, int lwtotl, int nrowa, int n, int nclin, int nctotl, int* istate, int* kactiv, int lcrash, double* bigbnd, double* a, double* bl, double* bu, double* featol, double* x)
        {

            string id = "VARBL LNCON NLCON ";

            short nerror;
            int nplin = n + nclin;

            double ftol;
            double test;
            int j, k;
            double b1, b2;
            //    int l1;
            bool ok;
            int was_is;

            /*     CHKDAT  CHECKS THE DATA INPUT TO THE VARIOUS OPTIMIZERS. */
            /*     THE FOLLOWING QUANTITIES ARE NOT CHECKED... */
            /*     NROWA, N, NCLIN, NCTOTL */
            /*     KACTIV */
            /*     A, X */
            --x;
            --featol;
            --bu;
            --bl;
            --kactiv;
            --istate;

            nerror = 0;
            /* --------------------------------------------------------------------- 
            */
            /*     CHECK THAT THERE IS ENOUGH WORKSPACE TO SOLVE THE PROBLEM. */
            /* --------------------------------------------------------------------- 
            */
            ok = litotl <= liwork && lwtotl <= lwork;
            if (ok && msg <= 0)
            {
                goto L20;
            }
            if (msg >= 0)
                lm_wmsg(
            "\nWORKSPACE PROVIDED IS     IW(%6ld),  W(%6ld).\nTO SOLVE PROBLEM WE NEED  IW(%6ld),  W(%6ld).",
                    CL(liwork), CL(lwork), CL(litotl), CL(lwtotl));
            if (ok)
                goto L20;
            ++nerror;
            if (msg >= 0)
                Console.WriteLine("\nXXX  NOT ENOUGH WORKSPACE TO SOLVE PROBLEM.");
            /* --------------------------------------------------------------------- 
            */
            /*     CHECK THE BOUNDS ON ALL VARIABLES AND CONSTRAINTS. */
            /* --------------------------------------------------------------------- 
            */
            L20:
            for (j = 1; j <= (int)nctotl; ++j)
            {
                b1 = bl[j];
                b2 = bu[j];
                ok = b1 <= b2;
                if (ok)
                {
                    goto L40;
                }
                ++nerror;
                k = j;
                //      l1 = 1;
                if (j > (int)n)
                {
                    k = j - n;
                }
                //   if (j > (int)n)
                //   {
                //      l1 = 4;
                //  }
                if (j > (int)nplin)
                {
                    k -= nclin;
                }
                //   if (j > (int)nplin)
                //   {
                //      l1 = 7;
                //  }
                if (msg >= 0)
                    lm_wmsg(
            "\nXXX THE BOUNDS ON %.2s%.2s%.1s%3ld ARE INCONSISTENT BL =%16.7lg BU =%16.7lg",
                    id, id + 2, id + 4, CL(k), b1,
                        b2);
                L40:
                ;
            }
            /* --------------------------------------------------------------------- 
            */
            /*     CHECK  BIGBND  AND  FEATOL. */
            /* --------------------------------------------------------------------- 
            */
            ok = *bigbnd > 0;
            if (ok) goto L60;
            ++nerror;
            if (msg >= 0)
                lm_wmsg("\nXXX BIGBND IS NOT POSITIVE...%20.9lg", *bigbnd);
            L60:
            for (j = 1; j <= (int)nctotl; ++j)
            {
                ftol = featol[j];
                test = 1 + ftol;
                ok = test > 1;
                if (!ok && msg >= 0)
                    lm_wmsg(
            "\n*** WARNING -- FEATOL(%4ld) IS LESS THAN MACHINE PRECISION...%16.6lg",
                    CL(j), ftol);
            }
            /* --------------------------------------------------------------------- 
            */
            /*     IF WARM START, CHECK  ISTATE. */
            /* --------------------------------------------------------------------- 
            */
            /* L100: */
            if (lcrash != 0)
            {
                for (j = 1; j <= (int)nplin; ++j)
                {
                    was_is = istate[j];
                    ok = was_is >= -2 && was_is <= 4;
                    if (!ok)
                    {
                        ++nerror;
                        if (msg >= 0)
                            lm_wmsg(
                        "\nXXX COMPONENT%5ld OF ISTATE IS OUT OF RANGE...%10ld",
                        CL(j), CL(was_is));
                    }
                }
            }
            return nerror;
        }
        public unsafe static void dqpcore(int orthog, int* unitq, int* inform, int* iter, int* itmax, int n, int* nclin, int* nctotl, int* nrowa, int* nrowh, int* ncolh, int* nactiv, int* nfree, int* istate, int* kactiv, int* kfree, double* objqp, double* xnorm, double* a, double* ax, double* bl, double* bu, double* clamda, double* cvec, double* featol, double* hess, double* scale, double* x, int* iw, double* w)
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
            int ifix;
            double palfa;
            int ifail;
            double condh, bigdx;
            int isdel = 0;
            double condt, drmin, anorm, dinky, drmax, hsize;
            int ncnln, ncolr, ncolz;
            double pnorm;
            //	int nrowj;
            int nhess;
            bool nullr, stall, uncon;
            int nclin0, kb;
            double bigalf, bigbnd, epspt9;
            double gfixed, alfhit;
            bool refine;
            int jdsave, nfixed;
            bool posdef;
            int modfyg, negligible;
            double condmx, atphit, cslast, gfnorm, rdlast, objsiz, snlast;
            int idummy, issave = 0, msglvl;
            int jsmlst, ksmlst;
            double smllst;
            int nstall, numinf;
            double ztgnrm;
            int firstv, hitlow;
            bool nocurv, renewr, unitpg, zerolm;
            int modfyr;
            double bnd;
            double gtp = 0;
            --w;
            --iw;
            --x;
            --scale;
            --featol;
            --cvec;
            --clamda;
            --bu;
            --bl;
            --ax;
            --kfree;
            --kactiv;
            --istate;

            /*INITIALIZE */
            *iter = 0;
            jadd = 0;
            jdel = 0;
            jdsave = 0;
            nclin0 = Math.Max(*nclin, 1);
            ncnln = 0;
            ncolz = *nfree - *nactiv;
            //	nrowj = 1;
            nstall = 0;
            nhess = 0;
            numinf = 0;
            msglvl = msg;
            msg = 0;
            if (istart == 0) msg = msglvl;
            bigbnd = parm[0];
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
            ncolr = dqpcrsh(*unitq, n, ncolz, *nfree, &nhess_conv,
                nq, *nrowh, *ncolh, nrowrt, &kfree[1], &hsize, hess,
                pRT, &scale[1], pZY, pRLAM, pWRK);
            dqpgrad(1, *unitq, n, CN(*nactiv), CN(*nfree), &nhess_conv, nq,
                CN(*nrowh), CN(*ncolh), jadd, &kactiv[1], &kfree[1], alfa, objqp, &gfixed,
                gtp, &cvec[1], hess, pPX, pQTG, &scale[1], &x[1], pZY, pWRK, pRLAM);
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
                clocker();
                /*
                ******* section i.  test for convergence *******
                compute the norms of the projected gradient and the gradient with
                respect to the free variables.
                */
                ztgnrm = 0;
                if (ncolr > 0) ztgnrm = dnrm2vec(ncolr, pQTG);
                gfnorm = ztgnrm;
                if (*nfree > 0 && *nactiv > 0) gfnorm = dnrm2vec(*nfree, pQTG);

                /*
                define small quantities that reflect the magnitude of  c,  x,  h
                and the matrix of constraints in the working set
                */
                objsiz = (BlasLike.lm_eps + Math.Abs(*objqp)) / (BlasLike.lm_eps + *xnorm);
                anorm = 0;
                if (*nactiv > 0) anorm = Math.Abs(dtmax);
                /*Computing MAX */
                dinky = Math.Max(anorm, objsiz);
                dinky = Math.Max(epspt9, DAS_Sol_zgfacc) * Math.Max(dinky, gfnorm);

                if (msg >= 80) wdinky("QPCORE", ztgnrm, dinky);

                /*
                print the details of this iteration
                use the largest and smallest diagonals of r to estimate the
                condition number of the projected hessian matrix
                */
                var dtmax_c = dtmax;
                var dtmin_c = dtmin;
                condt = dprotdiv(&dtmax_c, &dtmin_c, &ifail);
                dtmax = dtmax_c;
                dtmin = dtmin_c;
                if (ifail != 0 && dtmax == 0) condt = BlasLike.lm_max;
                if (ncolr > 0) BlasLike.dxminmax(ncolr, pRT, nrowrt + 1, &drmax, &drmin);
                condh = dprotdiv(&drmax, &drmin, &ifail);
                if (ifail != 0 && drmax == 0) condh = BlasLike.lm_max;
                if (condh >= BlasLike.lm_rootmax) condh = BlasLike.lm_max;
                if (condh < BlasLike.lm_rootmax) condh *= condh;
                var nrowrt_c = nrowrt;
                dqpprt(orthog, isdel, *iter, jadd, jdel, *nactiv, ncolz, *nfree,
                    n, *nclin, *nrowa, *&nrowrt_c, nhess, &
                    istate[1], &kfree[1], alfa, condh, condt, *objqp, gfnorm,
                    ztgnrm, emax, a, pRT, &x[1], pWRK, pAP);
                nrowrt = nrowrt_c;
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
                            dgetlamd(lprob, n, *nactiv, ncolz, *nfree, *nrowa, nrowrt,
                                &jsmlst, &ksmlst, &smllst, &istate[1], &kactiv[1], a,
                                pANORM, pQTG, pRLAM, pRT);

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
                            isdel = istate[jdel];
                            issave = isdel;
                            istate[jdel] = 0;

                            /*update the TQ factorization of the matrix of constraints in the working set*/

                            ddelcon(modfyg != 0, orthog != 0, *unitq, CN(jdel), CN(kdel),
                                CN(*nactiv), CN(ncolz), CN(*nfree), n,
                                CN(nq), CN(*nrowa), CN(nrowrt), &kactiv[1], &kfree[1], a,
                                pQTG, pRT, pZY);
                            ++ncolz;
                            if (jdel <= (int)n) ++(*nfree);
                            else --(*nactiv);
                        }

                        /*
                        the projected hessian is expanded by a row and column. compute
                        the elements of the new column of the cholesky factor r
                        use the vector p as temporary work space
                        */
                        renewr = true;
                        ++ncolr;
                        var nq_ccc = nq;
                        var nrowrt_ccc = nrowrt;
                        dqpcolr(&nocurv, &posdef, &renewr, unitq, n, &ncolr,
                            nfree, &nq_ccc, nrowh, ncolh, &nrowrt_ccc, &nhess, &kfree[1], &
                            cslast, &snlast, &drmax, &emax, &hsize, &rdlast, hess,
                            pRT, &scale[1], pZY, pPX, pWRK);
                        nq = nq_ccc;
                        nrowrt = nrowrt_ccc;
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
                if (*iter >= *itmax)
                {
                    /*too many iterations*/
                    i = 5;
                    break;
                }

                ++(*iter);
                if (*iter >= istart) msg = msglvl;
                var nq_c = nq;
                nrowrt_c = nrowrt;
                dfindp(&nullr, &unitpg, unitq, n, nclin, &nq_c, nrowa, &
                    nrowrt_c, &ncolr, &ncolz, nfree, &istate[1], &kfree[1],
                    negligible != 0, &gtp, &pnorm, &rdlast, a, pAP,
                    pPX, pQTG, pRT, pWRK, pZY, pWRK);
                nq = nq_c;
                nrowrt = nrowrt_c;

                /*
                if a constraint has just been deleted and the projected gradient
                is small (this can only occur here when the projected hessian is
                indefinite), the sign of  p  may be incorrect because of rounding
                errors in the computation of  ztg.  fix the sign of  p	by
                forcing it to satisfy the constraint that was just deleted
                */
                if ((jdsave > 0 && negligible != 0) || zerolm) dqpchkp(n, CN(*nclin), issave, jdsave, pAP, pPX);

                /*
                find the constraint we bump into along p
                update x and a*x if the step alfa is nonzero
                alfhit	is initialized to  bigalf.  if it remains that way after
                the call to bndalf, it will be regarded as infinite
                */
                bigalf = dprotdiv(&bigdx, &pnorm, &ifail);
                if (ifail != 0 && bigdx == 0) bigalf = BlasLike.lm_max;
                dbndalf(firstv != 0, &hitlow, &istate[1], &jadd, n,
                    CN(*nctotl), numinf, &alfhit, &palfa, &atphit, &bigalf, &
                    bigbnd, &pnorm, pANORM, pAP, &ax[1], &bl[1], &bu[1], &
                    featol[1], pPX, &x[1]);

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
                stall = Math.Abs(alfa * pnorm) <= epspt9 * *xnorm;
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
                    if (dqpgrad(2, *unitq, n, CN(*nactiv), CN(*nfree), &nhess_conv, nq,
                        CN(*nrowh), CN(*ncolh), jadd, &kactiv[1], &kfree[1], alfa, objqp, &
                        gfixed, gtp, &cvec[1], hess, pPX, pQTG, &
                        scale[1], &x[1], pZY, pWRK, pRLAM) < -BlasLike.lm_eps)
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
                    BlasLike.daxpy(n, alfa, pPX, 1, &x[1], 1);
                    if (*nclin > 0) BlasLike.daxpy(*nclin, alfa, pAP, 1, &ax[1], 1);
                    *xnorm = dnrm2vec(n, &x[1]);
                }

                /*if an unconstrained step was taken, repeat the main loop*/
                if (uncon) continue;

                /*add a constraint to the working set. update  istate*/
                if (hitlow != 0) istate[jadd] = 1;
                else istate[jadd] = 2;
                if (bl[jadd] == bu[jadd]) istate[jadd] = 3;

                /*
                if a bound is to be added, move x exactly onto it, except when
                a negative step was taken.  (bndalf  may have had to move to some
                other closer constraint.)
                */
                iadd = jadd - n;
                if (jadd <= (int)n)
                {
                    if (hitlow != 0) bnd = bl[jadd];
                    else bnd = bu[jadd];
                    if (alfa >= 0) x[jadd] = bnd;
                    i = *nfree;
                    for (ifix = 1; ifix <= (int)i; ++ifix) if (kfree[ifix] == jadd) break;
                }
                /*
                update the TQ factors of the matrix of constraints in the
                working set.  use the array  p	as temporary work space
                */
                var nq_cc = nq;
                var nrowrt_cc = nrowrt;
                daddcon(modfyg != 0, modfyr != 0, orthog != 0, unitq, &ifix, &iadd, &jadd,
                    nactiv, &ncolr, &ncolz, nfree, n, &nq_cc, nrowa, &nrowrt_cc,
                    &kfree[1], &condmx, &cslast, &snlast, a, pQTG,
                    pRT, pZY, pWRK, pPX);
                nq = nq_cc;
                nrowrt_c = nrowrt_cc;
                --ncolr;
                --ncolz;
                nfixed = n - *nfree;
                if (nfixed != 0)
                {
                    kb = *nactiv + nfixed;
                    i = nfixed;
                    for (idummy = 1; idummy <= nfixed; ++idummy)
                    {
                        kactiv[kb + 1] = kactiv[kb];
                        --kb;
                    }
                }
                if (jadd <= (int)n)
                {
                    /*
                    add a bound.  if stabilized eliminations are being used to update
                    the  TQ	 factorization,	 recompute the component of the gradient
                    corresponding to the newly fixed variable
                    use the array  p  as temporary work space
                    */
                    --(*nfree);
                    kactiv[*nactiv + 1] = jadd;
                    int nhess_conv_c = nhess;
                    if (orthog == 0)
                    {
                        dqpgrad(3, *unitq, n, CN(*nactiv), CN(*nfree), &nhess_conv, nq,
                            CN(*nrowh), CN(*ncolh), jadd, &kactiv[1], &kfree[1], alfa, objqp,
                            &pQTG[*nfree], gtp, &cvec[1], hess, pPX, pQTG, &
                            scale[1], &x[1], pZY, pWRK, pPX);
                        nhess = nhess_conv_c;
                    }
                }
                else
                {
                    /*add a general linear constraint*/
                    ++(*nactiv);
                    kactiv[*nactiv] = iadd;
                }

                /*
                repeat the main loop if the projected hessian that was used to
                compute this search direction was positive definite
                */
                if (ncolr == 0) posdef = true;
                if (ncolr == 0) emax = 0;

                if (!posdef)
                    /*
                    the projected hessian was not sufficiently positive definite
                    before the constraint was added.  either compute the true value
                    of the last diagonal of	 r  or	recompute the whole of its last
                    column. use the array  rlam  as temporary work space
                    */
                    nq_cc = nq;
                nrowrt_cc = nrowrt;
                dqpcolr(&nocurv, &posdef, &renewr, unitq, n, &ncolr,
                    nfree, &nq_cc, nrowh, ncolh, &nrowrt_cc, &nhess, &
                    kfree[1], &cslast, &snlast, &drmax, &emax, &hsize, &rdlast,
                    hess, pRT, &scale[1], pZY, pRLAM,
                    pWRK);
                nq = nq_cc;
                nrowrt = nrowrt_cc;
                /*.........................END OF MAIN LOOP............................*/
            }


            *inform = i;
            /*PRINT FULL SOLUTION*/
            msg = msglvl;
            if (msg >= 1) wrexit(lprob, i, *iter);
            if (i > 0) dgetlamd(lprob, n, *nactiv, ncolz, *nfree, *nrowa,
                     nrowrt, &jsmlst, &ksmlst, &smllst, &istate[1], &
                     kactiv[1], a, pANORM, pQTG, pRLAM, pRT);
            dprtsol(*nfree, *nrowa, n, *nclin, ncnln, *nctotl, bigbnd,
                *nactiv, &istate[1], &kactiv[1], a,
                &bl[1], &bu[1], &x[1], &clamda[1], pRLAM, &x[1]);
        }
        public unsafe static int dqpcrsh(int unitq, int n, int ncolz, int nfree, int* nhess, int Nq, int nrowh, int ncolh, int Nrowrt, int* kfree, double* hsize, double* hess, double* rt, double* scale, double* zy, double* hz1, double* wrk)
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
            int i, j, k, kmax, ksave, jthcol, ncolr;


            if (ncolz == 0) return 0;
            --wrk;
            zy_offset = Nq + 1;
            zy -= zy_offset;
            --scale;
            rt -= Nrowrt + 1;
            --kfree;
            ncolr = 0;
            /*
            compute  z(t) h z  and store the upper-triangular symmetric part
            in the first  ncolz  columns of  rt
            */
            for (k = 1; k <= ncolz; ++k)
            {

                BlasLike.dzerovec(n, &wrk[1]);
                if (unitq!=0)
                {
                    /* expand the column of  z  into an  n-vector */
                    for (i = 1; i <= nfree; ++i)
                    {
                        j = kfree[i];
                        wrk[j] = zy[i + k * Nq];
                    }
                    if (scldqp) Factorise.ddmxmulv(n, &scale[1], 1, &wrk[1], 1);
                    jthcol = 0;
                }
                else
                {
                    /*
                        only bounds are in the working set.  the  k-th column of  z is
                        just a column of the identity matrix
                    */
                    jthcol = kfree[k];
                    wrk[jthcol] = 1;
                }
                /*set  rt(*,k)  =  top of   h * (column of  z)*/
                qphess(n, nrowh, ncolh, jthcol, hess, &wrk[1], hz1);
                ++(*nhess);
                if (unitq!=0 && scldqp) BlasLike.dscalvec(n, scale[jthcol], hz1);
                if (scldqp) Factorise.ddmxmulv(n, &scale[1], 1, hz1, 1);
                dzyprod(4, n, nfree, ncolz, nfree, Nq, unitq, &kfree[1], &kfree[1],
                    hz1, &zy[zy_offset], &wrk[1]);
                BlasLike.dcopyvec(ncolz, hz1, &rt[k * Nrowrt + 1]);
                /*update an estimate of the size of the projected hessian*/
                t = Math.Abs(rt[k + k * Nrowrt]);
                if (t > *hsize) *hsize = t;
            }
            /*
                 form the cholesky factorization  r(t) r  =  z(t) h z  as far as
                 possible, using symmetric row and column interchanges
            */
            dmin_ = BlasLike.lm_eps * *hsize;
            for (j = 1; j <= ncolz; ++j)
            {
                /*FIND THE MAXIMUM REMAINING DIAGONAL*/
                kmax = j;
                dmax_ = rt[j + j * Nrowrt];
                for (k = j; k <= ncolz; ++k)
                {
                    d = rt[k + k * Nrowrt];
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
                    if (unitq!=0)
                    {
                        BlasLike.dcopyvec(nfree, &zy[kmax * Nq + 1], &wrk[1]);
                        BlasLike.dcopyvec(nfree, &zy[j * Nq + 1], &zy[kmax * Nq + 1]);
                        BlasLike.dcopyvec(nfree, &wrk[1], &zy[j * Nq + 1]);
                    }
                    else
                    {
                        /*Z is not stored explicitly*/
                        ksave = kfree[kmax];
                        kfree[kmax] = kfree[j];
                        kfree[j] = ksave;
                    }
                    /*interchange rows and columns of the projected hessian*/
//#if 1
			BlasLike.dswapvec( j, &rt[1+kmax*Nrowrt], &rt[1+j*Nrowrt]);
			BlasLike.dswap(kmax-j+1,&rt[j+kmax*Nrowrt], 1, &rt[j+j*Nrowrt], Nrowrt );
			BlasLike.dswap(ncolz+1-kmax,&rt[kmax+kmax*Nrowrt], Nrowrt, &rt[j+kmax*Nrowrt], Nrowrt );
/*#else
                    for (i = 1; i <= j; ++i)
                    {
                        t = rt[i + kmax * Nrowrt];
                        rt[i + kmax * Nrowrt] = rt[i + j * Nrowrt];
                        rt[i + j * Nrowrt] = t;
                    }
                    for (k = j; k <= kmax; ++k)
                    {
                        t = rt[k + kmax * Nrowrt];
                        rt[k + kmax * Nrowrt] = rt[j + k * Nrowrt];
                        rt[j + k * Nrowrt] = t;
                    }
                    for (k = kmax; k <= ncolz; ++k)
                    {
                        t = rt[kmax + k * Nrowrt];
                        rt[kmax + k * Nrowrt] = rt[j + k * Nrowrt];
                        rt[j + k * Nrowrt] = t;
                    }
#endif*/
                    rt[kmax + kmax * Nrowrt] = rt[j + j * Nrowrt];
                }
                /*set the diagonal element of R*/
                d = Math.Sqrt(dmax_);
                rt[j + j * Nrowrt] = d;
                if (j == ncolz) continue;
                /*
                set the above-diagonal elements of the k-th row of  r,
                and update the elements of all remaining rows
                */
                i = j + 1;
                for (k = i; k <= ncolz; ++k)
                {
                    t = rt[j + k * Nrowrt] / d;
                    rt[j + k * Nrowrt] = t;
                    /*R(I,K)  =  R(I,K)  - T * R(J,I),   I = i, k. */
                    if (t != 0) BlasLike.daxpy(k - j, -t, &rt[j + i * Nrowrt], Nrowrt, &rt[i + k * Nrowrt], 1);
                }
            }
            if (ncolr != ncolz && msg >= 80)
                lm_wmsg("\n//QPCRSH//  INDEFINITE PROJECTED HESSIAN.\n//QPCRSH// NCOLR=%6ld      NCOLZ=%6ld",
                    CL(ncolr), CL(ncolz));
            return ncolr;
        }

    }
}
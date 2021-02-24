using System;
using Blas;
using System.Diagnostics;
namespace ActiveSet
{
    public class Linear
    {
        public static int msg;
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
        public static int CN(int x) { return x; }
        public static int CI(int x) { return x; }
        public static byte CB(byte x) { return x; }
        public static bool CB(bool x) { return x; }
        public static int CL(int x) { return x; }
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
        public unsafe static void dlpcore(bool lp, int minsum, bool orthog, bool* unitq, int vertex, int* inform, int* iter,
                int* itmax, byte lcrash, int n, int* nclin, int* nctotl, int* nrowa, int* nactiv,
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
            int nrowj;
            int nclin0, kb;
            double bigalf, bigbnd, epspt9;
            int was_is;
            double feamax, feamin;
            bool delete_;
            int nfixed;
            int jbigst, kbigst;
            double tolact;
            int modfyg;
            double condmx, atphit, cslast, rdlast, objsiz, snlast, suminf,
                 trulam;
            int idummy, msglvl;
            int jsmlst, ksmlst;
            double smllst;
            int mstall;
            double ztgnrm;
            int nstall;
            bool firstv, hitlow;
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
            nrowj = 1;
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
            modfyg = 1;
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
            ddelcon(modfyg, orthog, CB(*unitq), CN(jdel), CN(kdel), CN(*nactiv), CN(ncolz), CN(*nfree), n,
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
            if (*iter >= *itmax)
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
            if (hitlow)
            {
                istate[jadd] = 1;
            }
            if (!hitlow)
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
            bnd = (hitlow ? bl : bu)[jadd];
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
            *inform = daddcon(modfyg, 0, orthog, unitq, &ifix, &iadd, &jadd,
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
        public unsafe static void dlpcrsh(bool orthog, bool* unitq, int vertex, byte lcrash, int n, int* nclin, int* nctotl,
        int* Nq, int* nrowa, int* Nrowrt, int* Ncolrt, int* nactiv, int* ncolz,
        int* nfree, int* istate, int* kactiv, int* kfree, double* bigbnd, double* tolact,
        double* xnorm, double* a, double* anorm, double* ax, double* bl, double* bu, double* x, double* qtg, double* rt,
        double* zy, double* p, double* wrk1, double* wrk2)
        {

            int a_dim1, a_offset, rt_dim1, rt_offset, zy_dim1, zy_offset, i__1,
                i__2;
            int c__1=1;
            int n__ = n;


            int iadd, jadd;
            double amin;
            int imin = 1000000000, jmin = 1000000000, ifix=-111111111;
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
            *unitq = true;
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
                    resl = Math.Abs(ax[i] - b1)/ (1 + Math.Abs(b1));
                L240:
                    if (noupp)
                    {
                        goto L260;
                    }
                    resu =  Math.Abs(ax[i] - b2) / (1 + Math.Abs(b2));
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
                        colsiz +=  Math.Abs(a[k + j * a_dim1]);
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
            if (!vertex || *ncolz == 0)
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
                inform = daddcon(0, 0, orthog, unitq, &ifix, &iadd, &jadd,
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
    }
}
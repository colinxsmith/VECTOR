using System;
using Blas;
using Solver;
using Ordering;
using System.IO;
using DataFile;
using System.Diagnostics;
namespace Portfolio
{
    public class Portfolio
    {
        int BACK;
        public void set_repeat<T>(int n, T p, T[] a)
        {
            //while(n--){*a++ = p;}
            for (var i = 0; i < n; ++i)
            {
                a[i] = p;
            }
        }
        public class roundstep
        {
            public double util;
            public double[] L;
            public double[] U;
            public double[] kL;
            public double[] kU;
            public double[] w;
            public OptParamRound info;
            public roundstep next;
            public roundstep prev;
            public int nround;
            public int back;
            public int[] bound_order;
            public int[] can_repeat;
            public int count;
            public bool success;
        }
        public double RoundInnerUtil(object info)
        {
            OptParamRound OP = (OptParamRound)info;
            INFO sendInput = (INFO)OP.MoreInfo;
            var w = new double[sendInput.n];
            var gradient = new double[sendInput.n];
            BlasLike.dcopyvec(sendInput.n, wback, w);
            int baskethere = -1, tradeshere = -1;
            return PortfolioUtility(sendInput.n, gamma, kappa, sendInput.buy, sendInput.sell, sendInput.alpha, w, gradient, ref baskethere, ref tradeshere, false);
        }
        public int RoundInnerOpt(object info)
        {
            OptParamRound OP = (OptParamRound)info;
            INFO vars = (INFO)OP.MoreInfo;
            double[] minholdlot = OP.minholdlot;
            double[] mintradelot = OP.mintradelot;
            OP.minholdlot = null;
            OP.mintradelot = null;
            OP.back = BACK = BasicOptimisation(vars.n, vars.m, vars.nfac, vars.A, OP.lower, OP.upper, gamma, kappa, vars.delta, vars.value, vars.valuel, vars.rmin, vars.rmax, vars.
                             alpha, vars.initial, vars.buy, vars.sell, vars.names, vars.useIP, vars.nabs, vars.A_abs, vars.L_abs, vars.U_abs, vars.mabs, vars.I_a);
            BlasLike.dcopyvec(vars.n, wback, OP.x);
            return BACK;
        }
        public delegate double UFUNC(int n, double gamma, double kappa, double[] buy, double[] sell, double[] alpha, double[] w, double[] gradient, ref int basket, ref int trades, bool print = true, double thresh = 1E-14);
        public delegate void GFUNC(int n, double[] w, double[] g);
        public delegate int FUNCGEN(object i);
        public delegate double UTILGEN(object i);
        public delegate int OFUNC(int n, int m, int nfac, double[] A, double[] L, double[] U, double gamma, double kappa, double delta, double value, double valuel, double rmin, double rmax, double[] alpha, double[] initial, double[] buy, double[] sell, string[] names, bool useIP = true, int nabs = 0, double[] A_abs = null, double[] L_abs = null, double[] U_abs = null, int mabs = 0, int[] I_a = null);
        public class OptParamRound
        {
            public int n;
            public int m;
            public double[] x;
            public double[] grad;
            public UTILGEN UtilityFunc;
            public FUNCGEN OptFunc;//This must set back member
            public GFUNC GradFunc;
            public object MoreInfo;
            public double[] lower;
            public double[] upper;
            public double[] c;
            public int back;
            public int mabs;
            public int lp;
            public int basket;
            public int trades;
            public int longbasket;
            public int shortbasket;
            public int tradebuy;
            public int tradesell;
            public double[] initial;
            public double equalbounds;
            public double dropfac;
            public double[] minholdlot;
            public double[] mintradelot;
            public object TimeOptData;
            public string dump;
        }
        double unround = 1e60;
        double round_eps = BlasLike.lm_eps8;
        double digitisei(double w, double initial, double minl, double sizl, double minlb = 0.0)
        {
            double ww = w - initial;
            double wa = Math.Abs(ww);
            double digit = 0, one = 1.0, p5 = 0.5;
            if (minlb == 0)
            {
                if (wa < BlasLike.lm_eps) { digit = 0; }
                else if (wa < minl) { digit = 0/*wa/minl*/; }
                else if (sizl < BlasLike.lm_eps) { digit = unround; }
                else if (wa >= minl) { digit = one + (wa - minl) / sizl; }
                if (ww < 0) { digit = -digit; }
            }
            else
            {//This isn't right or needed!
                if (ww < minl && ww > minlb)
                {
                    double init_check = (minl + minlb) * p5;
                    if (Math.Abs(init_check - initial) < 1e-15)
                    {
                        digit = initial;
                    }
                    else if (Math.Abs(init_check) < 1e-15)
                    {
                        digit = 0;
                    }
                }
                else { digit = unround; }
            }
            if (Math.Abs(digit) != unround)
            {

                digit = check_digit(digit);
            }

            return digit;
        }
        double check_digit(double digit)
        {
            double delta = Math.Abs(Math.Abs(digit - (long)(digit)) - 1);
            if (delta <= (BlasLike.lm_rooteps))
            {
                long ndelta = (long)(delta / BlasLike.lm_eps);
                if (digit > 0)
                {
                    digit += ndelta * BlasLike.lm_eps;
                }
                else if (digit < 0)
                {
                    digit -= ndelta * BlasLike.lm_eps;
                }
            }
            return digit;
        }
        double digit2w(double w, double initial, double d, double minl, double sizl, double minlb = 0.0)
        {
            d = check_digit(d);
            //Find nearest integer to d[i] (from digitise) and work out corresponding weight
            double roundw = w, p5 = 0.5;
            long p = 0, lone = 1;
            if (Math.Abs(d) == unround) { p = 1000000000; }
            else if (Math.Abs(d) <= BlasLike.lm_eps8) { p = 0; }
            /*	else if(d>0)p=(long)floor(d+p5);
                else if(d<0)p=(long)floor(d-p5);*/
            else p = (long)Math.Floor(Math.Abs(d) + p5);//This is the right way

            if (d < 0) p = -p;
            if (Math.Abs(d) == unround) { roundw = w - initial; }
            else if (Math.Abs(p) <= lone)
            {
                if (minlb == 0)
                    roundw = minl * p;
                else if (roundw >= 0)
                    roundw = minl * p;
                else if (roundw <= 0)
                    roundw = minlb * p;
            }
            else if (d > 0 && minl == sizl) { roundw = sizl * p; }
            else if (d < 0 && minl == sizl) { roundw = sizl * p; }
            else if (d > 0) { roundw = sizl * (p - lone) + minl; }
            else if (d < 0) { roundw = sizl * (p + lone) - minl; }
            roundw += initial;

            return roundw;
        }

        double round_weight(double x, double initial, double minl, double sizel, double minlb = 0.0)
        {
            return digit2w(x, initial, digitisei(x, initial, minl, sizel, minlb), minl, sizel, minlb);
        }
        public int round_check(int n, double[] w, double[] initial, double[] L, double[] U, double[] minlot, double[] sizelot, double eps = 0.0)
        {
            if (eps == 0.0) eps = BlasLike.lm_rooteps;
            int i, bad;
            bool badi;
            double init;
            long kk;
            for (i = 0, bad = 0; i < n; ++i)
            {
                init = initial != null ? initial[i] : 0;
                badi = false;
                if (Math.Abs(Math.Abs(w[i]) - Math.Abs(round_weight(w[i], init, minlot[i], sizelot != null ? sizelot[i] : 0))) > eps)
                {
                    if (Math.Abs(Math.Abs(w[i] - init) - minlot[i]) > eps)
                    {
                        if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                        {
                            kk = (long)((Math.Abs(w[i] - init) - minlot[i]) / sizelot[i]);
                            if (Math.Abs(kk * sizelot[i] + minlot[i] - Math.Abs(w[i] - init)) > eps)
                            {
                                badi = true;ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i]}[/yellow][red] BAD[/red]{L[i]} {w[i]} {U[i]}");
                            }
                        }
                        else if (Math.Abs(w[i] - init) - Math.Abs(minlot[i]) < -eps && Math.Abs(w[i] - init) > eps)
                            {badi = true;ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i]}[/yellow][red] BAD[/red]{L[i]} {w[i]} {U[i]}");}
                    }
                }
                if (L[i] == U[i]) continue;
                if (w[i] < L[i] - BlasLike.lm_eps8 || w[i] > U[i] + BlasLike.lm_eps8)
                {
                    badi = true;ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i]}[/yellow][red] BAD[/red]{L[i]} {w[i]} {U[i]}");
                }
                if (badi) bad++;
            }
            return n - bad;
        }///<summary>Check how well x is rounded</summary>
         ///<param name="n">Number of assets</param>
         ///<param name="x">Rounded weights</param>
         ///<param name="initial">Initial weights</param>
         ///<param name="minlot">Bottom "ladder rung"</param>
         ///<param name="sizelot">Other "ladder rungs"</param>
         ///<param name="shake">Output shake[i] = -1 if asset i is rounded, otherwise i"</param>
        public int roundcheck(int n, double[] x, double[] initial, double[] minlot, double[] sizelot, int[] shake)
        {
            var nround = 0;
            for (var i = 0; i < n; ++i)
            {
                var init = initial != null ? initial[i] : 0;
                var dd = 0L;
                double nwL, nwU;
                if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                {
                    dd = (long)digitisei(x[i], init, minlot[i], sizelot[i]);
                    if (x[i] > init)
                    {
                        nwL = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                        nwU = digit2w(x[i], init, dd + 1, minlot[i], sizelot[i]);
                    }
                    else
                    {
                        nwL = digit2w(x[i], init, dd - 1, minlot[i], sizelot[i]);
                        nwU = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                    }
                    if (Math.Abs(x[i] - nwL) < round_eps || Math.Abs(x[i] - nwU) < round_eps || Math.Abs(x[i] - init) < BlasLike.lm_eps)
                    {
                        shake[i] = -1;
                        nround++;
                    }
                    else shake[i] = i;
                }
                else if (minlot != null && minlot[i] > BlasLike.lm_eps)
                {
                    if (Math.Abs(x[i] - init) >= minlot[i] || Math.Abs(x[i] - init) < BlasLike.lm_eps)
                    {
                        shake[i] = -1;
                        nround++;
                    }
                    else shake[i] = i;
                }
                else
                {
                    shake[i] = -1;
                    nround++;
                }
            }
            return nround;
        }
        public void treenext(roundstep rstep, double[] initial, double[] minlot,
                                double[] sizelot, bool passedfromthresh = false, double[] thresh = null)
        {
            OptParamRound info = rstep.info;
            int n = info.n;
            int m = info.m;
            int firstlim = (n < 100) ? n : n, roundy = n;
            int stuck;
            double[] x = info.x;//,c=info.c,H=info.H;
            double[] bound_error = new double[n];

            if (rstep.prev != null)
            {
                BlasLike.dcopyvec(m + n, rstep.L, info.lower);
                BlasLike.dcopyvec(m + n, rstep.U, info.upper);
                //	info.OptSetup(basket,trades);
                info.OptFunc(info);
                //	rstep.util=info.utility_base(n,x,c,H);
                rstep.util = info.UtilityFunc(info);
                rstep.back = info.back;
                BlasLike.dcopyvec(n, x, rstep.w);
                BlasLike.dcopyvec(m + n, rstep.kL, info.lower);
                BlasLike.dcopyvec(m + n, rstep.kU, info.upper);
            }


            rstep.nround = 0;
            bool fixup = false;
            int i, j, i6 = 0;
            double init = 0, nwL = 0, nwU = 0;
            long dd;
            roundstep next = rstep.next = new roundstep(), roundstuck;
            if (next == null) return;
            next.can_repeat = rstep.can_repeat;
            next.success = rstep.success;

            next.kL = rstep.kL;
            next.kU = rstep.kU;
            next.L = new double[n + m];
            next.U = new double[n + m];
            next.w = new double[n];
            next.bound_order = new int[n];
            if (!(next.L != null && next.U != null && next.w != null && next.bound_order != null)) return;
            next.prev = rstep;
            next.info = info;
            next.next = null;
            next.nround = 0;
            next.util = BlasLike.lm_max;
            next.count = next.prev.count + 1;
            next.back = info.back;
            if (rstep.prev != null && rstep.prev.nround == n && rstep.back < 2)
            {
                rstep.back = 6; fixup = true;
            }
            else
                fixup = false;

            while (rstep.back == 6 && i6 < n)
            {
                for (j = i6; j < n; ++j)
                {
                    dd = 0;
                    if (rstep.prev != null && rstep.prev.nround == n)
                        i = rstep.prev.bound_order[j];
                    else
                        i = rstep.bound_order[j];
                    init = initial != null ? initial[i] : 0;
                    if (rstep.L[i] == rstep.kL[i] && rstep.U[i] == rstep.kU[i]) { i6++; continue; }
                    else
                    {
                        if (rstep.prev != null && rstep.prev.L[i] == rstep.kL[i] && rstep.prev.U[i] == rstep.kU[i])
                        {
                            if (rstep.U[i] == rstep.kU[i])
                            {
                                if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                                {
                                    dd = (long)digitisei(rstep.prev.w[i], init, minlot[i], sizelot[i]);
                                    if (rstep.prev.w[i] > init)
                                    {
                                        nwL = digit2w(rstep.prev.w[i], init, dd, minlot[i], sizelot[i]);
                                        nwU = digit2w(rstep.prev.w[i], init, dd + 1, minlot[i], sizelot[i]);
                                    }
                                    else
                                    {
                                        nwL = digit2w(rstep.prev.w[i], init, dd - 1, minlot[i], sizelot[i]);
                                        nwU = digit2w(rstep.prev.w[i], init, dd, minlot[i], sizelot[i]);
                                    }
                                }
                                else if (rstep.prev.w[i] >= init)
                                {
                                    nwL = init;
                                    nwU = minlot[i] + init;
                                }
                                else if (rstep.prev.w[i] < init)
                                {
                                    nwU = init;
                                    nwL = -minlot[i] + init;
                                }
                                rstep.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.kU[i]));
                                rstep.L[i] = rstep.kL[i];
                            }
                            else if (rstep.L[i] == rstep.kL[i])
                            {
                                if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                                {
                                    dd = (long)digitisei(rstep.prev.w[i], init, minlot[i], sizelot[i]);
                                    if (rstep.prev.w[i] > init)
                                    {
                                        nwL = digit2w(rstep.prev.w[i], init, dd, minlot[i], sizelot[i]);
                                        nwU = digit2w(rstep.prev.w[i], init, dd + 1, minlot[i], sizelot[i]);
                                    }
                                    else
                                    {
                                        nwL = digit2w(rstep.prev.w[i], init, dd - 1, minlot[i], sizelot[i]);
                                        nwU = digit2w(rstep.prev.w[i], init, dd, minlot[i], sizelot[i]);
                                    }
                                }
                                else if (rstep.prev.w[i] >= init)
                                {
                                    nwL = init;
                                    nwU = minlot[i] + init;
                                }
                                else if (rstep.prev.w[i] < init)
                                {
                                    nwU = init;
                                    nwL = -minlot[i] + init;
                                }
                                if (Math.Abs((nwU - (rstep.prev.w[i])) - (rstep.prev.w[i] - nwL)) < BlasLike.lm_rooteps)
                                {
                                    if (rstep.L[i] == rstep.kL[i])
                                    {
                                        rstep.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.L[i]));
                                        rstep.U[i] = rstep.kU[i];
                                    }
                                    else if (rstep.U[i] == rstep.kU[i])
                                    {
                                        rstep.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.U[i]));
                                        rstep.L[i] = rstep.kL[i];
                                    }
                                }
                                else
                                {
                                    rstep.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.kL[i]));
                                    rstep.U[i] = rstep.kU[i];
                                }
                            }
                        }
                        else
                        {
                            double wuse = fixup ? rstep.w[i] : rstep.prev.w[i];
                            if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                            {
                                dd = (long)digitisei(wuse, init, minlot[i], sizelot[i]);
                                if (wuse > init)
                                {
                                    nwL = Math.Max((digit2w(wuse, init, dd - 1, minlot[i], sizelot[i])), rstep.kL[i]);
                                    nwU = Math.Min((digit2w(wuse, init, dd + 1, minlot[i], sizelot[i])), rstep.kU[i]);
                                }
                                else
                                {
                                    nwL = Math.Max((digit2w(wuse, init, dd - 1, minlot[i], sizelot[i])), rstep.kL[i]);
                                    nwU = Math.Min((digit2w(wuse, init, dd + 1, minlot[i], sizelot[i])), rstep.kU[i]);
                                }
                            }
                            else if (wuse >= init)
                            {
                                nwL = init;
                                nwU = minlot[i] + init;
                            }
                            else if (wuse < init)
                            {
                                nwU = init;
                                nwL = -minlot[i] + init;
                            }
                            if (wuse > rstep.U[i])
                            {
                                rstep.U[i] = rstep.kU[i];
                                rstep.L[i] = Math.Max(nwU, rstep.kL[i]);
                            }
                            else if (wuse < rstep.L[i])
                            {
                                rstep.L[i] = rstep.kL[i];
                                rstep.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.kU[i]));
                            }
                            else if (rstep.prev != null && rstep.prev.nround == n)
                            {
                                if (rstep.L[i] == rstep.kL[i])
                                {
                                    rstep.U[i] = rstep.kU[i];
                                    rstep.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.kL[i]));
                                }
                                else if (rstep.U[i] == rstep.kU[i])
                                {
                                    rstep.L[i] = rstep.kL[i];
                                    rstep.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.kU[i]));
                                }
                                else if (wuse == rstep.L[i])
                                {
                                    rstep.L[i] = rstep.kL[i];
                                    rstep.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.kU[i]));
                                }
                                else if (wuse == rstep.U[i])
                                {
                                    rstep.U[i] = rstep.kU[i];
                                    rstep.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.kL[i]));
                                }
                            }
                        }
                        i6 = j + 1; break;
                    }
                }
                BlasLike.dcopyvec(m + n, rstep.L, info.lower);
                BlasLike.dcopyvec(m + n, rstep.U, info.upper);
                //		info.OptSetup(basket,trades);
                info.OptFunc(info);
                //		rstep.util=info.utility_base(n,x,c,H);
                rstep.util = info.UtilityFunc(info);
                rstep.back = info.back;
                if (rstep.back == 6)
                {
                    j = i6 - 1;
                    i = rstep.bound_order[j];
                    if (rstep.U[i] == rstep.kU[i])//&&rstep.prev.bound_order[j]!=rstep.bound_order[j])
                    {
                        rstep.U[i] = rstep.L[i];
                        rstep.L[i] = rstep.kL[i];
                    }
                    else if (rstep.L[i] == rstep.kL[i])//&&rstep.prev.bound_order[j]!=rstep.bound_order[j])
                    {
                        rstep.L[i] = rstep.U[i];
                        rstep.U[i] = rstep.kU[i];
                    }
                    BlasLike.dcopyvec(n, rstep.w, x);
                }
                BlasLike.dcopyvec(n, x, rstep.w);
                BlasLike.dcopyvec(m + n, rstep.kL, info.lower);
                BlasLike.dcopyvec(m + n, rstep.kU, info.upper);
            }
            BlasLike.dcopyvec(m + n, rstep.L, next.L);
            BlasLike.dcopyvec(m + n, rstep.U, next.U);
            double switch1 = 1;
            for (i = 0; i < n; ++i)
            {
                init = initial != null ? initial[i] : 0; dd = 0;
                if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                {
                    dd = (long)digitisei(x[i], init, minlot[i], sizelot[i]);
                    if (x[i] > init)
                    {
                        nwL = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                        nwU = digit2w(x[i], init, dd + 1, minlot[i], sizelot[i]);
                    }
                    else
                    {
                        nwL = digit2w(x[i], init, dd - 1, minlot[i], sizelot[i]);
                        nwU = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                    }
                    if (Math.Abs(x[i] - nwL) < round_eps || Math.Abs(x[i] - nwU) < round_eps || Math.Abs(x[i] - init) < BlasLike.lm_eps)
                    {
                        bound_error[i] = i;
                        rstep.nround++;
                        continue;
                    }
                    if ((nwU - (x[i])) / (x[i] - nwL) < switch1)
                    {
                        if (nwU >= rstep.L[i])
                            bound_error[i] = n + nwU - (x[i]);
                        else
                            bound_error[i] = n + x[i] - nwL;
                    }
                    else
                    {
                        if (nwL <= rstep.U[i])
                            bound_error[i] = n + x[i] - nwL;
                        else
                            bound_error[i] = n + nwU - (x[i]);
                    }
                }
                else
                {
                    if ((thresh != null && Math.Abs(x[i]) >= thresh[i]) || Math.Abs(x[i]) < BlasLike.lm_eps)
                    {
                        if (Math.Abs(x[i] - init) >= minlot[i] || Math.Abs(x[i] - init) < BlasLike.lm_eps)
                        {
                            bound_error[i] = i;
                            rstep.nround++;
                            continue;
                        }
                        else
                            bound_error[i] = n + Math.Abs(x[i] - init);
                    }
                    else
                    {
                        if (rstep.kL[i] == rstep.kU[i]) { rstep.nround++; bound_error[i] = i; }
                        else bound_error[i] = n + Math.Max(Math.Abs(x[i]), (Math.Abs(x[i] - init)));
                    }
                }
            }
            ColourConsole.WriteEmbeddedColourLine($"[green]first nround=[/green][cyan]{rstep.nround}[/cyan]");
            /*if(!sizelot)rstep.nround=thresh_check(n,x,initial,rstep.kL,rstep.kU,minlot,0,round_eps);
            else*/
            rstep.nround = round_check(n, x, initial, rstep.kL, rstep.kU, minlot, sizelot, round_eps);
            ColourConsole.WriteEmbeddedColourLine($"[green]then  nround=[/green][cyan]{rstep.nround}[/cyan]");
            Ordering.Order.getorder(n, bound_error, next.bound_order, null);//printorder(n,next.bound_order);
                                                                            //	for(j=rstep.nround;j<min(n/4+rstep.nround,n);++j)
            roundy = Math.Max(((int)(rstep.nround * .5 + n * .5)), (rstep.nround + 1));
            //	stuck=(rstep.prev&&(rstep.prev.nround==rstep.nround))?true:false;
            stuck = 0; roundstuck = rstep;
            while (roundstuck.prev != null && (roundstuck.prev.nround == roundstuck.nround))
            {
                stuck++;
                roundstuck = roundstuck.prev;
            }
            for (j = 0; j < Math.Min(Math.Max(1, roundy), n); ++j)
            {
                //		i=next.bound_order[n-j-1];
                i = next.bound_order[j]; dd = 0;
                init = initial != null ? initial[i] : 0;
                if (sizelot != null && sizelot[i] > BlasLike.lm_eps)
                {
                    dd = (long)digitisei(x[i], init, minlot[i], sizelot[i]);
                    if (x[i] > init)
                    {
                        if (!(j % 2 != 0 && next.count % 2 != 0))
                        {
                            nwL = digit2w(x[i], init, dd - 1, minlot[i], sizelot[i]);
                            nwU = digit2w(x[i], init, dd + (long)stuck, minlot[i], sizelot[i]);
                            if (nwL < rstep.kL[i])
                            {
                                nwL = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                                nwU = digit2w(x[i], init, dd + 1 + (long)stuck, minlot[i], sizelot[i]);
                            }
                        }
                        else
                        {
                            nwL = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                            nwU = digit2w(x[i], init, dd + 1 + (long)stuck, minlot[i], sizelot[i]);
                        }
                    }
                    else
                    {
                        if (!(j % 2 != 0 && next.count % 2 != 0))
                        {
                            nwL = digit2w(x[i], init, dd - (long)stuck, minlot[i], sizelot[i]);
                            nwU = digit2w(x[i], init, dd + 1, minlot[i], sizelot[i]);
                            if (nwU > rstep.kU[i])
                            {
                                nwL = digit2w(x[i], init, dd - 1 - (long)stuck, minlot[i], sizelot[i]);
                                nwU = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                            }
                        }
                        else
                        {
                            nwL = digit2w(x[i], init, dd - 1 - (long)stuck, minlot[i], sizelot[i]);
                            nwU = digit2w(x[i], init, dd, minlot[i], sizelot[i]);
                        }
                    }
                    //info.AddLog((char*)"closeness %d %e %e %e\n",i,x[i]-nwL,x[i]-nwU,x[i]);
                    if (Math.Abs(x[i] - nwL) < round_eps || Math.Abs(x[i] - nwU) < round_eps || Math.Abs(x[i] - init) < BlasLike.lm_eps)
                    {
                        continue;
                        //				break;
                    }
                    if (Math.Abs((nwU - (x[i])) - (x[i] - nwL)) < BlasLike.lm_rooteps)
                    {
                        if (false && !(j % 2 != 0 && next.count % 2 != 0))
                        {
                            next.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.L[i]));
                        }
                        /*    else if (false)
                            {
                                next.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.U[i]));
                            }*/
                    }
                    else if ((nwU - (x[i])) / (x[i] - nwL) < switch1)
                    {
                        if (nwU >= rstep.L[i] && nwU <= rstep.kU[i])
                            next.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.L[i]));
                        else
                            next.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.U[i]));
                    }
                    else
                    {
                        if (nwL <= rstep.U[i] && nwL >= rstep.kL[i])
                            next.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.U[i]));
                        else
                            next.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.L[i]));
                    }
                    if (next.U[i] < next.L[i])
                    {
                        if (x[i] > init)
                        {
                            next.L[i] = Math.Min(rstep.kU[i], Math.Max(nwU, rstep.kL[i]));
                            next.U[i] = rstep.kU[i];
                        }
                        else if (x[i] < init)
                        {
                            next.U[i] = Math.Max(rstep.kL[i], Math.Min(nwL, rstep.kU[i]));
                            next.L[i] = rstep.kL[i];
                        }
                        else
                        {
                            next.U[i] = Math.Max(rstep.kL[i], Math.Min(init, rstep.kU[i]));
                            next.L[i] = rstep.kL[i];
                        }
                        //printf((char*)"bound problem\n");
                    }
                    //			printf((char*)"%3d %20.8e %20.8e %20.8e trade %20.8e\n",i,next.L[i],next.U[i],x[i],x[i]-init);
                }
                else
                {
                    i = next.bound_order[j];
                    if ((thresh != null && Math.Abs(x[i]) >= thresh[i]) || Math.Abs(x[i]) < BlasLike.lm_eps)
                    {
                        if (!(rstep.nround == n && Math.Abs(Math.Abs(x[i] - init) - minlot[i]) < BlasLike.lm_eps8) && ((Math.Abs(x[i] - init) >= Math.Abs(minlot[i]) || Math.Abs(x[i] - init) < BlasLike.lm_eps)))
                        {
                            continue;
                        }
                        if (Math.Abs(minlot[i]) - Math.Abs(x[i] - init) < Math.Abs(x[i] - init) * switch1)
                        {
                            if (x[i] - init > 0)
                            {
                                next.L[i] = Math.Max(Math.Min(rstep.kU[i], minlot[i] + init), rstep.kL[i]);
                                next.L[i] = Math.Min(next.U[i], next.L[i]);
                            }
                            else
                            {
                                next.U[i] = Math.Min(Math.Max(rstep.kL[i], -minlot[i] + init), rstep.kU[i]);
                                next.U[i] = Math.Max(next.U[i], next.L[i]);
                            }
                        }
                        else
                        {
                            next.U[i] = Math.Min(Math.Max(rstep.kL[i], init), rstep.kU[i]);
                            next.L[i] = Math.Max(Math.Min(rstep.kU[i], init), rstep.kL[i]);
                        }
                        if (rstep.nround == n && Math.Abs(Math.Abs(x[i] - init) - minlot[i]) < BlasLike.lm_eps8 && rstep.can_repeat[i] != 0)
                        {
                            if (rstep.can_repeat[i] == 3)
                            {
                                rstep.can_repeat[i]--;
                                if (rstep.kL[i] < init + BlasLike.lm_eps8 && x[i] > init)
                                {
                                    next.L[i] = init + minlot[i];//rstep.kL[i] next time
                                    next.U[i] = Math.Min(init + minlot[i], rstep.kU[i]);
                                }
                                else if (rstep.kU[i] > init - BlasLike.lm_eps8 && x[i] < init)
                                {
                                    next.U[i] = init - minlot[i];//rstep.kU[i] next time
                                    next.L[i] = Math.Max(rstep.kL[i], init - minlot[i]);
                                }
                            }
                            else if (rstep.can_repeat[i] == 2)
                            {
                                rstep.can_repeat[i]--;
                                if (rstep.kL[i] < init + BlasLike.lm_eps8 && x[i] > init)
                                {
                                    next.L[i] = init + minlot[i];//rstep.kL[i] next time
                                    next.U[i] = rstep.kU[i];
                                }
                                else if (rstep.kU[i] > init - BlasLike.lm_eps8 && x[i] < init)
                                {
                                    next.U[i] = init - minlot[i];//rstep.kU[i] next time
                                    next.L[i] = rstep.kL[i];
                                }
                            }
                            else if (rstep.can_repeat[i] == 1)
                            {
                                rstep.can_repeat[i]--;
                                next.L[i] = rstep.kL[i];
                                next.U[i] = rstep.kU[i];
                            }
                        }
                    }
                    else if (Math.Abs(minlot[i]) - Math.Abs(x[i] - init) < Math.Abs(x[i] - init) * switch1)
                    {
                        if (Math.Abs(minlot[i]) - Math.Abs(x[i] - init) < Math.Abs(x[i] - init) * switch1 && (thresh != null && Math.Abs(thresh[i]) - Math.Abs(x[i]) < Math.Abs(x[i]) * switch1))
                        {
                            if (x[i] - init > 0)
                                next.L[i] = Math.Max(Math.Min(rstep.kU[i], Math.Max(thresh != null ? thresh[i] : 0, minlot[i] + init)), rstep.kL[i]);
                            else
                                next.U[i] = Math.Min(Math.Max(rstep.kL[i], Math.Min(thresh != null ? -thresh[i] : 0, -minlot[i] + init)), rstep.kU[i]);
                        }
                        else
                        {
                            if (x[i] - init > 0)
                                next.L[i] = Math.Max(Math.Min(rstep.kU[i], Math.Max(thresh != null ? thresh[i] : 0, init)), rstep.kL[i]);
                            else
                                next.U[i] = Math.Min(Math.Max(rstep.kL[i], Math.Min(thresh != null ? -thresh[i] : 0, init)), rstep.kU[i]);
                        }
                    }
                    else
                    {
                        if (x[i] - init > 0)
                            next.L[i] = Math.Max(Math.Min(rstep.kU[i], 0), rstep.kL[i]);
                        else
                            next.U[i] = Math.Min(Math.Max(rstep.kL[i], 0), rstep.kU[i]);
                    }
                }
            }
            if (rstep.nround == n && next.count == 2 && rstep.back <= 1) { rstep.util = info.UtilityFunc(info); return; }
            if (!next.success && rstep.nround == n && rstep.back <= 1)
                next.success = true;
            if (rstep.nround == n && passedfromthresh && next.count > 20) { rstep.util = info.UtilityFunc(info); return; }
            if ((rstep.nround < n && next.count < (firstlim * 2) && !next.success) || (next.count < firstlim/*&&info.TimeOptData==0*/))
            {
                ColourConsole.WriteEmbeddedColourLine($"[yellow]stage {next.count}[/yellow][green] {rstep.nround} rounded[/green]");
                //		_ASSERT(next.count<7);
                treenext(next, initial, minlot, sizelot, passedfromthresh, thresh);
            }
        }
        public bool treestart(OptParamRound info, bool passedfromthresh, double[] initial, double[] minlot,
                        double[] sizelot, double[] roundw, double[] thresh = null)
        {
            int i;
            int n = info.n;
            int m = info.m;
            double[] x = info.x;//At the start this is the unrounded optimal solution
            roundstep next = new roundstep();
            roundstep start, prev;
            next.can_repeat = new int[n];
            set_repeat(n, 3, next.can_repeat);

            next.success = false;
            next.count = 1;
            next.kL = new double[n + m];
            next.kU = new double[n + m];
            BlasLike.dcopyvec(n + m, info.lower, next.kL);
            BlasLike.dcopyvec(n + m, info.upper, next.kU);
            next.L = new double[n + m];
            next.U = new double[n + m];
            next.w = new double[n];
            next.bound_order = new int[n]; for (i = 0; i < n; ++i) next.bound_order[i] = i;
            next.prev = null;
            next.info = info;
            next.back = info.back;
            BlasLike.dcopyvec(n, x, next.w);
            BlasLike.dcopyvec(m + n, next.kL, next.L);
            BlasLike.dcopyvec(m + n, next.kU, next.U);
            treenext(next, initial, minlot, sizelot, passedfromthresh, thresh);//すぎの木
            start = next;
            prev = null;
            while (next != null)
            {
                prev = next;
                next = next.next;
            }
            double ulow = BlasLike.lm_max;
            bool back = false;
            int nround = 0;
            if (prev.prev != null && prev.prev.back < 2)
            {
                BlasLike.dcopyvec(n, prev.prev.w, roundw);
                ColourConsole.WriteEmbeddedColourLine($"[green]{prev.prev.nround} Rounded at level[/green]\t[yellow]{prev.prev.count}[/yellow]");
            }
            if (prev.nround == n && prev.back < 2)
            {
                BlasLike.dcopyvec(n, prev.w, roundw);
                ulow = prev.util;
                info.back = prev.back;
                back = true;
            }
            while (prev != null)
            {
                if (prev.nround > nround)
                {
                    BlasLike.dcopyvec(n, prev.w, roundw); nround = prev.nround;
                }
                if (prev.nround == n && ulow > prev.util && prev.back < 2)
                {
                    BlasLike.dcopyvec(n, prev.w, roundw);
                    ulow = prev.util;
                    info.back = prev.back;
                    back = true;
                }
                prev = prev.prev;
            }
            ColourConsole.WriteEmbeddedColourLine($"[green]{nround}[/green]\t[yellow]stocks were rounded properly[/yellow]");
            return back;
        }

        public void Rounding(int basket, int trades, double[] initial, double[] minlot,
                                double[] sizelot, double[] roundw, double[] minholdlot, double[] mintradelot, object info)
        {
            OptParamRound Op = new OptParamRound();
            Op.MoreInfo = info;
            Op.basket = basket;
            Op.trades = trades;
            Op.lower = (double[])((INFO)info).L.Clone();
            Op.m = ((INFO)info).m;
            Op.n = ((INFO)info).n;
            Op.OptFunc = RoundInnerOpt;
            Op.upper = (double[])((INFO)info).U.Clone();
            Op.UtilityFunc = RoundInnerUtil;
            Op.x = (double[])wback.Clone();
            Op.minholdlot = minholdlot;
            Op.mintradelot = mintradelot;
            var ffi = 0;
            for (var i = 0; i < Op.n; ++i)
            {
                if (Op.lower[i] == Op.upper[i])
                {
                    minlot[i] = 0;
                    sizelot[i] = 0;
                    ffi++;
                }
            }
            if (ffi > 0) ColourConsole.WriteEmbeddedColourLine($"[cyan]Changed lot for[/cyan] [red]{ffi}[/red][cyan] lots due to fixed bounds[/cyan]");
            for (var i = 0; i < Op.n; ++i)
            {
                if (Op.lower[i] >= 0 && minlot[i] >= Op.lower[i]) { ColourConsole.WriteEmbeddedColourLine($"[cyan]Increase Lower bound for {names[i]}[/cyan][green] {Op.lower[i]}[/green][red] due to minlot of {minlot[i]}[/red]"); Op.lower[i] = minlot[i]; }
                if (Op.upper[i] <= 0 && -minlot[i] <= Op.upper[i]) { ColourConsole.WriteEmbeddedColourLine($"[cyan]Decrease Upper bound for {names[i]}[/cyan][green] {Op.upper[i]}[/green][red] due to minlot of {-minlot[i]}[/red]"); Op.upper[i] = -minlot[i]; }
            }
            if (treestart(Op, false, initial, minlot, sizelot, roundw)) { BACK = 0; ((INFO)info).back = BACK; }
        }

        public class INFO
        {
            public int n; public int m; public int nfac; public double[] A; public double[] L; public double[] U; public
               double delta; public double value; public double valuel; public
             double rmin; public double rmax; public double[] alpha; public double[] initial; public double[] buy; public double[] sell; public
             string[] names; public bool useIP; public int nabs; public double[] A_abs; public double[] L_abs; public double[] U_abs; public
             int mabs; public int[] I_a;
            public double kappa = 0.5;
            public double[] bench;
            public double target;
            public int back = -1;
        }
        ///<summary> A function that returns risk - target risk. Use with Solve1D to do risk constraint</summary>
        ///<param name="gam">guess for gamma (for return risk utility) or (gamma and kappa for cost risk utility)</param>
        ///<param name="info">object for passing other information</param>
        public double CalcRisk(double gam, object info)
        {
            INFO vars = (INFO)info;
            vars.back = BasicOptimisation(vars.n, vars.m, vars.nfac, vars.A, vars.L, vars.U, gam, gam, vars.delta, vars.value, vars.valuel, vars.rmin, vars.rmax, vars.
                     alpha, vars.initial, vars.buy, vars.sell, vars.names, vars.useIP, vars.nabs, vars.A_abs, vars.L_abs, vars.U_abs, vars.mabs, vars.I_a);
            double[] www = (double[])wback.Clone();
            if (vars.bench != null) BlasLike.dsubvec(vars.n, www, vars.bench, www);
            var fix = nfixed;
            nfixed = 0;
            var backr = Math.Sqrt(Variance(www)) - vars.target;
            nfixed = fix;
            return backr;
        }

        ///<summary>Utility function that tries to drop assets to give desired basket and/or number of trades and also apply a risk constraint</summary>
        ///<param name="basket">Desired number of non zero assets</param>
        ///<param name="trades">Desired number of non zero trades</param>
        ///<param name="info">object for passing other information</param>
        ///<param name="targetRisk">Desired risk</param>
        public void DropRisk(int basket, int trades, double targetRisk, object info)
        {
            var baskethere = 0;
            var tradeshere = 0;
            var sendInput = (Portfolio.INFO)info;
            var newgamma = ActiveSet.Optimise.Solve1D(CalcRisk, 0, 1, 0, sendInput);
            var back = sendInput.back;
            if (newgamma > 10 || sendInput.back == 6) ColourConsole.WriteError("Infeasible target risk");
            else
            {
                kappa = gamma = newgamma;
                var riskh = CalcRisk(gamma, sendInput) + targetRisk;
                back = sendInput.back;
                if (back != 6) ColourConsole.WriteEmbeddedColourLine($"[green]risk for multiplier {gamma,20:e12} is[/green]\t[yellow]{riskh,20:e12}[/yellow]\t[cyan]target risk {targetRisk}[/cyan]");
                else ColourConsole.WriteError("INFEASIBLE");
                sendInput.useIP = false;
                back = BasicOptimisation(sendInput.n, sendInput.m, sendInput.nfac, sendInput.A, sendInput.L, sendInput.U, gamma, kappa, sendInput.delta, sendInput.value, sendInput.valuel, sendInput.rmin, sendInput.rmax, sendInput.
                 alpha, sendInput.initial, sendInput.buy, sendInput.sell, sendInput.names, sendInput.useIP, sendInput.nabs, sendInput.A_abs, sendInput.L_abs, sendInput.U_abs, sendInput.mabs, sendInput.I_a);
                var w = new double[sendInput.n];
                var gradient = new double[sendInput.n];
                BlasLike.dcopyvec(sendInput.n, wback, w);
                var utility = PortfolioUtility(sendInput.n, gamma, kappa, sendInput.buy, sendInput.sell, sendInput.alpha, w, gradient, ref baskethere, ref tradeshere);
                ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio Utility (standard form):\t[/magenta][green]{utility,20:e12}[/green]");
                if (back != 6) back = Dropper(sendInput.n, sendInput.m, sendInput.nfac, sendInput.A, sendInput.L, sendInput.U, gamma, kappa, sendInput.delta, sendInput.value, sendInput.valuel, sendInput.rmin, sendInput.rmax, sendInput.
                     alpha, sendInput.initial, sendInput.buy, sendInput.sell, sendInput.names, sendInput.useIP, sendInput.nabs, sendInput.A_abs, sendInput.L_abs, sendInput.U_abs, sendInput.mabs, sendInput.I_a, basket, baskethere, trades, tradeshere);
                if (back == 6) ColourConsole.WriteError("INFEASIBLE");
                BlasLike.dcopyvec(sendInput.n, wback, w);
                if (back != 6)
                {
                    BoundsSetToSign(sendInput.n, sendInput.L, sendInput.U, sendInput.initial, w, true);
                    //Try to get the risk constraint correct if possible
                    newgamma = ActiveSet.Optimise.Solve1D(CalcRisk, 0, 1, 0, sendInput);
                    if (newgamma > 10 || sendInput.back == 6) ColourConsole.WriteError("Infeasible target risk");
                    else
                    {
                        kappa = gamma = newgamma; riskh = CalcRisk(gamma, sendInput) + targetRisk;
                        ColourConsole.WriteEmbeddedColourLine($"[green]risk for multiplier {gamma,20:e12} is[/green]\t[yellow]{riskh,20:e12}[/yellow]\t[cyan]target risk {targetRisk}[/cyan]");
                    }
                    if (newgamma > 10) newgamma = gamma;
                    back = sendInput.back;
                    if (back == 6) ColourConsole.WriteError("INFEASIBLE");
                    else
                    {
                        kappa = gamma = newgamma;
                        sendInput.useIP = true;
                        BasicOptimisation(sendInput.n, sendInput.m, sendInput.nfac, sendInput.A, sendInput.L, sendInput.U, gamma, kappa, sendInput.delta, sendInput.value, sendInput.valuel, sendInput.rmin, sendInput.rmax, sendInput.
                         alpha, sendInput.initial, sendInput.buy, sendInput.sell, sendInput.names, sendInput.useIP, sendInput.nabs, sendInput.A_abs, sendInput.L_abs, sendInput.U_abs, sendInput.mabs, sendInput.I_a);
                        w = new double[sendInput.n];
                        gradient = new double[sendInput.n];
                        BlasLike.dcopyvec(sendInput.n, wback, w);
                        utility = PortfolioUtility(sendInput.n, gamma, kappa, sendInput.buy, sendInput.sell, sendInput.alpha, w, gradient, ref baskethere, ref tradeshere);
                        ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio Utility (standard form):\t[/magenta][green]{utility,20:e12}[/green]");
                    }
                }
            }
        }

        ///<summary>object for passing other information for risk constrained optimisation via goal seek method</summary>
        public class DropKeep
        {
            ///<summary>Utility for weights in w</summary>
            public double utility;
            ///<summary>Keep utility value for weights in w</summary>
            public double[] w;
        }
        public static void fixup_zero(double[] s, double thresh = 1e-14)
        {
            for (var i = 0; i < s.Length; ++i)
            {
                if (Math.Abs(s[i]) < thresh) s[i] = 0;
            }
        }
        ///<summary>
        /// Return the utility of a portfolio with weights w[i]
        /// Also return utility gradient
        /// and the number of non-zero w[i] (basket size) and non-zero w[i]-initial[i] (trades)
        ///</summary>
        ///<param name="n">Number of assets in portfolio</param>
        ///<param name="gamma">gamma/(1-gamma) is the multiplier for returns</param>
        ///<param name="kappa">kappa/(1-kappa) is the multiplier for costs</param>
        ///<param name="buy">buy[i] is the buy cost gradient for i'th asset</param>
        ///<param name="sell">sell[i] is the sell cost gradient for i'th asset</param>
        ///<param name="alpha">alpha[i] is the return for i'th asset</param>
        ///<param name="w">w[i] is the weight for i'th asset</param>
        ///<param name="gradient">gradient[i] is the gradient for i'th asset</param>
        ///<param name="print">print output if true</param>
        ///<param name="thresh">less than thresh means zero</param>
        public double PortfolioUtility(int n, double gamma, double kappa, double[] buy, double[] sell, double[] alpha, double[] w, double[] gradient, ref int basket, ref int trades, bool print = true, double thresh = 1e-14)
        {
            var nfixedo = nfixed;
            nfixed = 0;//Must set this here
            BlasLike.dzerovec(n, gradient);
            double back = 0;
            if (gamma != 0 && alpha != null)
            {
                BlasLike.daxpyvec(n, -gamma / (1 - gamma), alpha, gradient);
                back += BlasLike.ddotvec(n, w, gradient);
                if (bench != null) back -= BlasLike.ddotvec(n, bench, gradient);
            }
            if (kappa != 0 && buy != null && sell != null)
            {
                for (var i = 0; i < n; ++i)
                {
                    var costgrad = kappa / (1 - kappa) * (w[i] > initial[i] ? buy[i] : -sell[i]);
                    gradient[i] += costgrad;
                    back += costgrad * (w[i] - initial[i]);
                }
            }
            if (Q != null)
            {
                var implied = new double[n];
                hessmull(n, Q, w, implied);
                back += 0.5 * BlasLike.ddotvec(n, implied, w);
                BlasLike.daxpyvec(n, 1, implied, gradient);
                if (bench != null)
                {
                    hessmull(n, Q, bench, implied);
                    BlasLike.daxpyvec(n, -1, implied, gradient);
                    back -= 0.5 * BlasLike.ddotvec(n, implied, w);
                }
            }
            trades = 0; basket = 0;
            if (print) ColourConsole.WriteEmbeddedColourLine($"[green]{"Weight",14}[/green]\t[red]{"Gradient",14}[/red]\t[yellow]{"weight*gradient",15}[/yellow]");
            for (var i = 0; i < n; ++i)
            {
                if (Math.Abs(w[i]) >= thresh) basket++;
                if (Math.Abs(w[i] - initial[i]) >= thresh) trades++;
                if (print) ColourConsole.WriteEmbeddedColourLine($"[green]{w[i],14:E6}[/green]\t[red]{gradient[i],14:E6}[/red]\t[yellow]{w[i] * gradient[i],15:E6}[/yellow]");
            }
            ColourConsole.WriteEmbeddedColourLine($"[green]Basket size:{basket,4}[/green]\t[cyan]Trade size:{trades,4}[/cyan]");
            nfixed = nfixedo;
            return back;
        }
        public int Dropper(int n, int m, int nfac, double[] A, double[] L, double[] U,
        double gamma, double kappa, double delta, double value, double valuel,
        double rmin, double rmax, double[] alpha, double[] initial, double[] buy, double[] sell,
        string[] names, bool useIp = true, int nabs = 0, double[] A_abs = null, double[] Abs_L = null, double[] Abs_U = null,
        int mabs = 0, int[] I_A = null, int basket = -1, int baskethere = -1, int trades = -1, int tradeshere = -1)
        {
            if (basket < 0) { basket = n; baskethere = n; }
            if (trades < 0) { trades = n; tradeshere = n; }
            var twoside = (trades < n && basket < n);
            double[] gradient = new double[n];
            double[] w = (double[])wback.Clone();
            double utility = 0;
            if (basket < baskethere || trades < tradeshere)
            {
                int basketnow = 0, tradesnow = 0;
                var order = new int[w.Length];
                var orderT = new int[w.Length];
                var orderk = new int[w.Length];
                var orderkT = new int[w.Length];
                var dropbad = new byte[w.Length];
                for (var i = 0; i < w.Length; i++)
                {
                    if (L[i] > 0 || U[i] < 0) dropbad[i] = 1;
                    else if (L[i] == U[i] && (L[i] != initial[i] || L[i] != 0)) dropbad[i] = 1;
                }
                var oldL = (double[])L.Clone();
                var oldU = (double[])U.Clone();
                Ordering.Order.getorderabs(w.Length, w, order, dropbad);
                BlasLike.dsubvec(w.Length, w, initial, w);
                Ordering.Order.getorderabs(w.Length, w, orderT, dropbad);
                BlasLike.daddvec(w.Length, w, initial, w);
                var done = 0;
                if (basket < baskethere)
                    for (var i = basket; i < w.Length; ++i)
                    {
                        var ii = order[i];
                        L[ii] = U[ii] = 0;
                    }
                if (trades < tradeshere)
                {
                    for (int i = trades; i < w.Length; ++i)
                    {
                        var ii = orderT[i];
                        if (L[ii] == U[ii]) done++;
                        else L[ii] = U[ii] = initial[ii];
                    }
                }
                var back = BasicOptimisation(n, m, nfac, A, L, U, gamma, kappa, delta, value, valuel, rmin, rmax,
                    alpha, initial, buy, sell, names, useIp, nabs, A_abs, Abs_L, Abs_U, mabs, I_A);
                if (back != 6)
                {
                    BlasLike.dcopyvec(n, wback, w);
                    utility = PortfolioUtility(n, gamma, kappa, buy, sell, alpha, w, gradient, ref basketnow, ref tradesnow, false);
                    ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio Utility (standard form):\t[/magenta][green]{utility,20:e12}[/green]");
                }
                BlasLike.dcopyvec(L.Length, oldL, L);
                BlasLike.dcopyvec(U.Length, oldU, U);
                DropKeep[] overall = new DropKeep[0];
                if (back != 6 && basket >= basketnow && trades >= tradesnow)
                {
                    Array.Resize(ref overall, overall.Length + 1);
                    overall[overall.Length - 1] = new DropKeep();
                    overall[overall.Length - 1].w = (double[])w.Clone();
                    overall[overall.Length - 1].utility = utility;
                }
                else ColourConsole.WriteEmbeddedColourLine($"[green]NAIVE DROP was [/green][red]INFEASIBLE[/red]");
                int interimBasket, interimTrades;
                var scale = 0.5;
                basketnow = baskethere;
                tradesnow = tradeshere;
                interimBasket = (int)Math.Floor(basketnow * scale + (1.0 - scale) * basket);
                interimTrades = (int)Math.Floor(tradesnow * scale + (1.0 - scale) * trades);
                bool same = false;
                bool fast = true;
                var bad = 0;
                while ((basketnow > basket || tradesnow > trades) && bad < w.Length)
                {
                    ColourConsole.WriteEmbeddedColourLine($"[magenta]INTERIM BASKET[/magenta][green] {interimBasket}[/green]");
                    ColourConsole.WriteEmbeddedColourLine($"[darkcyan]INTERIM TRADES[/darkcyan][green] {interimTrades}[/green]");
                    if ((twoside && !fast) || basketnow > basket)
                        for (var i = interimBasket; i < w.Length; ++i)
                        {
                            var ii = order[i];
                            L[ii] = U[ii] = 0;
                        }
                    if ((twoside) || tradesnow > trades)
                    {
                        done = 0;
                        for (var i = w.Length - 1; i >= interimTrades - done; --i)
                        {
                            if (i == 0) break;
                            var ii = orderT[i];
                            if (L[ii] == U[ii] && L[ii] != initial[ii]) done++;
                            else L[ii] = U[ii] = initial[ii];
                        }
                    }
                    back = BasicOptimisation(n, m, nfac, A, L, U, gamma, kappa, delta, value, valuel, rmin, rmax,
                       alpha, initial, buy, sell, names, useIp, nabs, A_abs, Abs_L, Abs_U, mabs, I_A);
                    while (back == 6)
                    {
                        fast = false;
                        BlasLike.dcopyvec(L.Length, oldL, L);
                        BlasLike.dcopyvec(U.Length, oldU, U);
                        interimBasket = basketnow - 1;
                        interimTrades = tradesnow - 1;
                        if (interimBasket < basket) interimBasket = basket;
                        if (interimTrades < trades) interimTrades = trades;
                        ColourConsole.WriteEmbeddedColourLine($"[magenta]INTERIM BASKET[/magenta][green] {interimBasket}[/green]");
                        ColourConsole.WriteEmbeddedColourLine($"[magenta]INTERIM TRADES[/magenta][green] {interimTrades}[/green]");
                        if (basketnow > basket)
                            for (var i = interimBasket; i < w.Length; ++i)
                            {
                                var ii = order[i];
                                L[ii] = U[ii] = 0;
                            }
                        if (tradesnow > trades)
                        {
                            done = 0;
                            for (var i = w.Length - 1; i >= interimTrades - done; --i)
                            {
                                if (i == 0) break;
                                var ii = orderT[i];
                                if (L[ii] == U[ii] && L[ii] != initial[ii]) done++;
                                else L[ii] = U[ii] = initial[ii];
                            }
                        }
                        back = BasicOptimisation(n, m, nfac, A, L, U, gamma, kappa, delta, value, valuel, rmin, rmax,
                           alpha, initial, buy, sell, names, useIp, nabs, A_abs, Abs_L, Abs_U, mabs, I_A);
                        if (back != 6)
                        {
                            BlasLike.dcopyvec(n, wback, w);
                            utility = PortfolioUtility(n, gamma, kappa, buy, sell, alpha, w, gradient, ref basketnow, ref tradesnow, false);
                            ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio Utility (standard form):\t[/magenta][green]{utility,20:e12}[/green]");
                        }
                        if (basketnow > basket)
                        {
                            if (back == 6) dropbad[order[interimBasket]] = 1;
                            else
                            {
                                BlasLike.dcopyvec(n, wback, w); break;
                            }
                            Ordering.Order.getorderabs(w.Length, w, order, dropbad);
                        }
                        if (tradesnow > trades)
                        {
                            if (back == 6) dropbad[orderT[interimTrades]] = 1;
                            else
                            {
                                BlasLike.dcopyvec(n, wback, w); break;
                            }
                            BlasLike.dsubvec(w.Length, w, initial, w);
                            Ordering.Order.getorderabs(w.Length, w, orderT, dropbad);
                            BlasLike.daddvec(w.Length, w, initial, w);
                        }
                        bad = 0;
                        foreach (var k in dropbad)
                        {
                            bad += k;
                        }
                        ColourConsole.WriteEmbeddedColourLine($"[green]Number of assets that cannot be dropped:\t[/green][red]{bad,20}[/red]");
                        if (bad > basket) break;
                        if (bad > trades) break;
                    }
                    if (back == 6) break;
                    BlasLike.dcopyvec(n, wback, w);
                    utility = PortfolioUtility(n, gamma, kappa, buy, sell, alpha, w, gradient, ref basketnow, ref tradesnow, false);
                    ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio Utility (standard form):\t[/magenta][green]{utility,20:e12}[/green]");
                    BlasLike.dcopyvec(L.Length, oldL, L);
                    BlasLike.dcopyvec(U.Length, oldU, U);
                    if (fast) interimBasket = (int)Math.Floor(basketnow * scale + (1.0 - scale) * basket);
                    else interimBasket = basketnow - 1;
                    if (fast) interimTrades = (int)Math.Floor(tradesnow * scale + (1.0 - scale) * trades);
                    else interimTrades = tradesnow - 1;
                    if (basket < n && trades < n)
                    {
                        if (trades < basket) interimBasket = basket;
                        else interimTrades = trades;
                    }
                    Ordering.Order.getorderabs(w.Length, w, orderk, dropbad);
                    BlasLike.dsubvec(w.Length, w, initial, w);
                    Ordering.Order.getorderabs(w.Length, w, orderkT, dropbad);
                    BlasLike.daddvec(w.Length, w, initial, w);
                    same = true;
                    if (basketnow > basket)
                        for (var i = 0; i < w.Length; ++i)
                        {
                            same = same && (orderk[i] == order[i]);
                            order[i] = orderk[i];
                        }
                    if (tradesnow > trades)
                        for (var i = 0; i < w.Length; ++i)
                        {
                            same = same && (orderkT[i] == orderT[i]);
                            orderT[i] = orderkT[i];
                        }
                    if (same) ColourConsole.WriteEmbeddedColourLine($"[cyan]DROPPING ORDER DID NOT CHANGE THIS TIME[/cyan]");
                }
                BlasLike.dcopyvec(L.Length, oldL, L);
                BlasLike.dcopyvec(U.Length, oldU, U);
                if (back != 6 && basketnow <= basket && tradesnow <= trades)
                {
                    Array.Resize(ref overall, overall.Length + 1);
                    overall[overall.Length - 1] = new DropKeep();
                    overall[overall.Length - 1].w = (double[])w.Clone();
                    overall[overall.Length - 1].utility = utility;
                }
                var minutil = BlasLike.lm_max;
                var minId = -1;
                for (var i = 0; i < overall.Length; ++i)
                {
                    if (minutil > overall[i].utility)
                    {
                        minId = i;
                        minutil = overall[i].utility;
                    }
                }
                if (minId != -1)
                {
                    wback = overall[minId].w;
                    ColourConsole.WriteEmbeddedColourLine($"[magenta]Best drop utility[/magenta] ([darkmagenta]{minId + 1}[/darkmagenta]):\t[green]{overall[minId].utility,20:e12}[/green][cyan] (out of {overall.Length})[/cyan]");
                    return 0;
                }
                else { ColourConsole.WriteError("NUMBER CONSTRAINT COULD NOT BE MET"); return 6; }
            }
            return 0;
        }
        ///<summary>
        /// Project out the extra variables added to handle non-linear constraints
        /// leaving an effective model. Show that primal effective utility = dual effective utility
        ///</summary>
        ///<param name="LAMBDA">The Lagrangian multipliers as defined in Active Set</param>
        ///<param name="cextra">The addition to c due to benchmark</param>
        public void UtilityAnalysis(double[] LAMBDA, double[] cextra = null)
        {/*
            for (var i = 0; i < n + m; ++i)
            {
                Console.WriteLine($"{i,10}\t{LAMBDA[i],16:F8}");
            }*/
            var ntrue = this.ntrue - nfixed;
            var Ceff = new double[ntrue];
            var chere = cextra == null ? c : cextra;
            for (var i = 0; i < ntrue; ++i)
            {
                Ceff[i] = -BlasLike.ddot(m - mtrue, A, 1, LAMBDA, 1, mtrue + i * m, n + mtrue);
            }
            var implied = new double[n];
            if (Q != null) hessmull(n, Q, w, implied);
            BlasLike.daddvec(ntrue, Ceff, chere, Ceff);
            if (Q != null) BlasLike.daddvec(ntrue, Ceff, implied, Ceff);
            var cvals = new double[mtrue];
            for (var i = 0; i < mtrue; ++i)
            {
                cvals[i] = BlasLike.ddot(ntrue, A, m, w, 1, i);
                ColourConsole.WriteEmbeddedColourLine($"[red]LAMBDA {LAMBDA[i + n],20:f16}[/red][green] L {((L[i + n] == -BlasLike.lm_max) ? -100.0 : L[i + n]),12:F8}[/green][cyan] value {cvals[i],12:F8}[/cyan][green] U {((U[i + n] == BlasLike.lm_max) ? 100.0 : U[i + n]),12:F8}[/green]");
            }
            var dual = BlasLike.ddotvec(ntrue, LAMBDA, w) + BlasLike.ddotvec(mtrue, LAMBDA, cvals, n) - BlasLike.ddotvec(n - ntrue, LAMBDA, w, ntrue, ntrue) - 0.5 * BlasLike.ddotvec(ntrue, w, implied);
            var primal = BlasLike.ddotvec(ntrue, Ceff, w) - 0.5 * BlasLike.ddotvec(ntrue, w, implied);
            var old = Console.ForegroundColor;
            ColourConsole.WriteEmbeddedColourLine($"[red]Effective model with non-linear extra part projected out[/red]");
            ColourConsole.WriteEmbeddedColourLine($"\t\t\t\t\t[green]Weight[/green]\t\t\t[cyan]Effective Utility Gradient[/cyan]");
            for (var i = 0; i < ntrue; ++i)
            {
                ColourConsole.WriteEmbeddedColourLine($"[darkred]{i + 1,5}[/darkred][red]{names[i],30}[/red][green]{w[i],12:F8}[/green]\t\t\t[cyan]{Ceff[i],20:e12}[/cyan]");
            }
            ColourConsole.WriteEmbeddedColourLine($"Primal:\t[magenta]{primal,12:F8}[/magenta]");
            ColourConsole.WriteEmbeddedColourLine($"Dual:\t[cyan]{dual,12:F8}[/cyan]");
        }
        public virtual void WriteInputs(string filename)
        {
            using (StreamWriter writer = new StreamWriter(filename))
            {
                writer.WriteLine("n");
                writer.WriteLine(n);
                writer.WriteLine("ntrue");
                writer.WriteLine(ntrue);
                writer.WriteLine("m");
                writer.WriteLine(m);
                writer.WriteLine("mtrue");
                writer.WriteLine(mtrue);
                printVector("initial", initial, writer);
                printVector("WW", w, writer);
                printVector("LL", L, writer);
                printVector("UU", U, writer);
                printVector("bench", bench, writer);
                printVector("QQ", Q, writer, true);
                printVector("AA", A, writer, m);
                printVector("CC", c, writer);
            }
        }
        public Portfolio(string file)
        {
            inFile = file;
            if (inFile != "")
                using (var OptData = new InputSomeData())
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
                    c = OptData.mapDouble["alpha"];
                    BlasLike.dnegvec(c.Length, c);
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

        public void BuySellSetup(int n, int m, int nfac, double[] A, double[] L, double[] U, double gamma, double kappa, double delta, double[] alpha, double[] initial, double[] buy, double[] sell, string[] names, bool useIP = true)
        {
            if (delta < 0) delta = 2;
            // kappa = -1;
            var useCosts = kappa > 0.0;
            if (!useCosts) kappa = 0;
            this.ntrue = n;
            this.mtrue = m;
            makeQ();
            var bothsellbuy = false; //bothsellbuy = false means treat sell side only
            var N = n + n;
            var M = m + n + (delta < 1.0 ? 1 : 0);
            if (bothsellbuy)
            {
                N += n;
                M += 2 * n;
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
            if (useIP) BlasLike.dsetvec(n, BlasLike.lm_max, UU, n);
            else
            {
                BlasLike.dcopyvec(n, U, UU, 0, n);
            }
            if (bothsellbuy)
            {
                BlasLike.dsetvec(n, 0, LL, n + n);
                if (useIP) BlasLike.dsetvec(n, BlasLike.lm_max, UU, n + n);
                else
                {
                    BlasLike.dcopyvec(n, U, UU, 0, n + n);
                }
            }
            //Constraints
            BlasLike.dcopyvec(m, L, LL, n, N);
            BlasLike.dcopyvec(m, U, UU, n, N);
            for (var i = 0; i < m; i++)
            {
                BlasLike.dcopy(n, A, m, AA, M, i, i);
            }
            BlasLike.dcopyvec(n, initial, LL, 0, N + m);
            if (useIP) BlasLike.dsetvec(n, BlasLike.lm_max, UU, N + m);
            else
            {
                BlasLike.dsetvec(n, 1, UU, N + m);
            }
            for (var i = m; i < m + n; ++i)
            {
                BlasLike.dset(1, 1.0, AA, M, i + M * (i - m));
                BlasLike.dset(1, 1.0, AA, M, i + M * (n + i - m));
            }
            if (bothsellbuy)
            {
                BlasLike.dsccopyvec(n, 1.0, initial, UU, 0, N + m + n);
                if (useIP) BlasLike.dsetvec(n, -BlasLike.lm_max, LL, N + m + n);
                else
                {
                    BlasLike.dsetvec(n, -1.0, LL, N + m + n);
                }
                for (var i = m + n; i < m + n + n; ++i)
                {
                    BlasLike.dset(1, 1.0, AA, M, i + M * (i - m - n));
                    BlasLike.dset(1, -1.0, AA, M, i + M * (n + n + i - m - n));
                }


                BlasLike.dcopyvec(n, initial, UU, 0, N + m + n * 2);
                BlasLike.dcopyvec(n, initial, LL, 0, N + m + n * 2);
                for (var i = m + n * 2; i < m + n * 2 + n; ++i)
                {
                    BlasLike.dset(1, 1.0, AA, M, i + M * (i - m - n * 2));
                    BlasLike.dset(1, 1.0, AA, M, i + M * (n + i - m - n * 2));
                    BlasLike.dset(1, -1.0, AA, M, i + M * (n + n + i - m - n * 2));
                }
            }
            BlasLike.dsccopyvec(n, -gamma / (1 - gamma), alpha, CC);
            if (useCosts)
            {
                var mult = kappa / (1.0 - kappa);
                if (!bothsellbuy)
                {
                    BlasLike.daxpyvec(n, mult, buy, CC);
                    BlasLike.dsccopyvec(n, mult, buy, CC, 0, n);
                    BlasLike.daxpyvec(n, mult, sell, CC, 0, n);
                }
                else
                {
                    BlasLike.dsccopyvec(n, mult, sell, CC, 0, n);
                    BlasLike.dsccopyvec(n, mult, buy, CC, 0, 2 * n);
                }
            }
            else
            {
                BlasLike.dsetvec(n, 0, CC, n);
                if (bothsellbuy) BlasLike.dsetvec(n, 0, CC, n + n);
            }
            if (delta < 1.0)
            {
                LL[N + M - 1] = -BlasLike.lm_max * 0;// Proper lower bound <=0 is redundant
                if (bothsellbuy)
                {
                    BlasLike.dset(n * 2, 1.0, AA, M, M - 1 + M * n);
                    /* LL[N + M - 1] =*/
                    UU[N + M - 1] = delta * 2;
                }
                else
                {
                    BlasLike.dset(n, 1.0, AA, M, M - 1);
                    BlasLike.dset(n, 2.0, AA, M, M - 1 + M * n);
                    /*   LL[N + M - 1] = */
                    UU[N + M - 1] = 2.0 * delta + BlasLike.dsumvec(n, initial);
                }
            }
            this.L = LL;
            this.U = UU;
            this.A = AA;
            this.n = N;
            this.m = M;
            this.gamma = gamma;
            this.c = CC;
            if (useIP)
            {
                var back = InteriorOpt(1e-10, WW);
            }
            else
            {
                Q = null;
                BlasLike.dsetvec(WW.Length, 1.0 / n, WW);
                for (var i = 0; i < n; ++i)
                {
                    WW[i + n] = Math.Abs(initial[i]) < BlasLike.lm_eps ? 0 : Math.Max(0, (initial[i] - 1.0 / n));
                    if (bothsellbuy) WW[i + 2 * n] = Math.Abs(initial[i]) < BlasLike.lm_eps ? 0 : Math.Max(0, -(initial[i] - 1.0 / n));
                }
                this.initial = initial;
                this.w = WW;
                if (bothsellbuy) BlasLike.dcopyvec(n, initial, WW, 0, n + n);
                WriteInputs("./optinput1");
                var LAMBDAS = new double[N + M];
                var back = ActiveOpt(1, WW, LAMBDAS);
                Console.WriteLine($"back = {back}");
                makeQ();
                BlasLike.dsetvec(WW.Length, 1.0 / n, WW);
                for (var i = 0; i < n; ++i)
                {
                    WW[i + n] = Math.Abs(initial[i]) < BlasLike.lm_eps ? 0 : Math.Max(0, (initial[i] - 1.0 / n));
                    if (bothsellbuy) WW[i + 2 * n] = Math.Abs(initial[i]) < BlasLike.lm_eps ? 0 : Math.Max(0, -(initial[i] - 1.0 / n));
                }
                this.w = WW;
                WriteInputs("./optinput2");
                back = ActiveOpt(0, WW, LAMBDAS);
                Console.WriteLine($"back = {back}");
            }
            var turnover = 0.0;
            var cost = 0.0;
            for (var i = 0; i < n; ++i)
            {
                turnover += Math.Abs(WW[i] - initial[i]);
                if ((buy != null) && (sell != null))
                {
                    var diff = (WW[i] - initial[i]);
                    cost += diff > 0 ? diff * buy[i] : -diff * sell[i];
                }
                var c1 = BlasLike.ddot(N, AA, M, WW, 1, i + m);
                if (bothsellbuy)
                {
                    var c2 = BlasLike.ddot(N, AA, M, WW, 1, n + i + m);
                    if (WW[i] <= initial[i]) Console.WriteLine($"{names[i]}\t{(WW[i] - initial[i]),12:F8}\t{WW[i + n],12:F8} {c1,12:F8} {WW[i + n + n],12:F8} {c2,12:F8} {initial[i],12:F8}");
                    else Console.WriteLine($"{names[i]} {(WW[i] - initial[i]),12:F8}\t{WW[i + n],12:F8} {c1,12:F8} {WW[i + n + n],12:F8} {c2,12:F8} {initial[i],12:F8}");
                }
                else
                {
                    if (WW[i] <= initial[i]) Console.WriteLine($"{names[i]}\t{(WW[i] - initial[i]),12:F8}\t{WW[i + n],12:F8} {(c1 - initial[i]),12:F8}  {initial[i],12:F8}\t\t{(!useIP ? (UU[i + N + m] - c1) : 10):f2}");
                    else Console.WriteLine($"{names[i]} {(WW[i] - initial[i]),12:F8}\t{WW[i + n],12:F8} {(c1 - initial[i]),12:F8}  {initial[i],12:F8}\t\t{(!useIP ? (UU[i + N + m] - c1) : 10):f2}");
                }
            }
            var eret = BlasLike.ddotvec(n, alpha, WW);
            var variance = Variance(WW);
            var eretA = -BlasLike.ddotvec(n, CC, WW) + (bothsellbuy ? 0 : (kappa / (1 - kappa) * BlasLike.ddotvec(n, buy, WW)));
            var turn2 = BlasLike.dsumvec(n, WW, n) + (BlasLike.dsumvec(n, WW) - BlasLike.dsumvec(n, initial)) * 0.5;
            var costA = 0.0;
            if (!bothsellbuy) costA = BlasLike.ddotvec(n, WW, sell, n) + BlasLike.ddotvec(n, WW, buy, n) + BlasLike.ddotvec(n, WW, buy) - BlasLike.ddotvec(n, initial, buy);
            else costA = BlasLike.ddotvec(n, WW, sell, n) + BlasLike.ddotvec(n, WW, buy, n + n);
            Console.WriteLine($"Variance: {variance}");
            Console.WriteLine($"Return: {eret}: {eretA}");
            Console.WriteLine($"Turnover: {turnover * 0.5}: {turn2}");
            Console.WriteLine($"Cost: {cost}:  {costA}");
            var extra = 0.0;
            if (bench != null)
            {
                var implied = new double[ntrue];
                hessmull(n, 1, 1, 1, Q, bench, implied);
                extra = -BlasLike.ddotvec(ntrue, WW, implied);
            }
            var utility = -gamma / (1 - gamma) * eret + kappa / (1 - kappa) * (cost + (bothsellbuy ? 0 : BlasLike.ddotvec(n, initial, buy))) + 0.5 * variance + extra;
            var utilityA = BlasLike.ddotvec(N, CC, WW) + 0.5 * variance + extra;
            Console.WriteLine($"Utility: {utility}: {utilityA}");
            for (var i = 0; i < m; ++i)
            {
                var ccval = BlasLike.ddot(n, A, m, WW, 1, i);
                Console.WriteLine($"Portfolio constraint {i}: {ccval}");
            }
            //            ActiveSet.Optimise.printV("optimal weights", WW, n);
        }
        ///<summary>Portfolio Optimisation with BUY/SELL utility and LONG/SHORT constraints
        ///If a variable's upper and lower bounds are equal, this variable is re-ordered out of the optimisaion
        ///</summary>
        ///<param name="n">Number of Assets</param>
        ///<param name="m">Number of Constraints</param>
        ///<param name="nfac">Number of Risk Factors (-1 for Full Covariance)</param>
        ///<param name="A">Constraint Array A[i*m+j] exposure of i'th asset to j'th constraint</param>
        ///<param name="L">Array of lower bounds L[i] is for asset i L[n+i] for constraint i</param>
        ///<param name="U">Array of upper bounds U[i] is for asset i U[n+i] for constraint i</param>
        ///<param name="gamma">gamma/(1-gamma) is the multiplier for returns</param>
        ///<param name="kappa">kappa/(1-kappa) is the multiplier for costs</param>
        ///<param name="delta">desired turnover from initial</param>
        ///<param name="value">upper bound of value of the long side for long/short problem</param>
        ///<param name="valuel">lower bound of value of the long side for long/short problem</param>
        ///<param name="rmin">lower bound of short/long</param>
        ///<param name="rmax">upper bound of short/long</param>
        ///<param name="alpha">array of expected returns</param>
        ///<param name="initial">array of initial weights</param>
        ///<param name="buy">array of buy costs</param>
        ///<param name="sell">array of sell costs</param>
        ///<param name="names">array of asset names</param>
        ///<param name="useIP">if true use interior point method</param>
        ///<param name="nabs">number of absolute constraints in A_abs</param>
        ///<param name="A_abs">array of absolute constraints in long/short problem A_abs[i*nabs+j] exposure of i'th asset to j'th constraint</param>
        ///<param name="L_abs">array of lower bounds for absolute constraints (nabs+mabs)</param>
        ///<param name="U_abs">array of upper bounds for absolute constraints (nabs+mabs)</param>
        ///<param name="mabs">number of absolute constraints in I_a</param>
        ///<param name="I_a">integer array, if I_a[i] = k, the i'th constraint has data from the k'th constraint in A</param>
        public int BasicOptimisation(int n, int m, int nfac, double[] A, double[] L, double[] U,
        double gamma, double kappa, double delta, double value, double valuel,
        double rmin, double rmax, double[] alpha, double[] initial, double[] buy, double[] sell,
        string[] names, bool useIP = true, int nabs = 0, double[] A_abs = null, double[] L_abs = null, double[] U_abs = null,
        int mabs = 0, int[] I_a = null)
        {
            int back;
            nfixed = 0;
            ntrue = n;
            mtrue = m;
            fixedW = new double[n];
            makeQ();
            var mainorder = new int[n];
            var mainorderInverse = new int[n];
            var fixedSecondOrder = new double[n];
            for (var i = 0; i < n; ++i)
            {
                if (L[i] == U[i]) nfixed++;
                mainorder[i] = i;
            }
            var boundLU = new double[m];
            if (nfixed > 0)
            {
                int i = 0, I = n - 1, ifixed = 0;
                for (; I >= 0; --I)
                {
                    if (L[mainorder[I]] != U[mainorder[I]])
                    {
                        for (; i < I; ++i)
                        {
                            if (L[mainorder[i]] == U[mainorder[i]])
                            {
                                Order.swap(ref mainorder[i], ref mainorder[I]); ifixed++; i++; break;
                            }
                        }
                    }
                    else ifixed++;
                    if (ifixed == nfixed) break;
                }
                for (i = 0; i < n; ++i)
                {
                    mainorderInverse[mainorder[i]] = i;
                }
                Order.Reorder(n, mainorder, L);
                Order.Reorder(n, mainorder, U);
                Order.Reorder(n, mainorder, alpha);
                Order.Reorder(n, mainorder, initial);
                if (bench != null) Order.Reorder(n, mainorder, bench);
                if (buy != null) Order.Reorder(n, mainorder, buy);
                if (sell != null) Order.Reorder(n, mainorder, sell);
                if (names != null) Order.Reorder(n, mainorder, names);
                Order.Reorder_gen(n, mainorder, A, m, 1, true);
                if (A_abs != null) Order.Reorder_gen(n, mainorder, A_abs, nabs, 1, true);
                if (Q != null && nfac == -1) Order.ReorderSymm(n, mainorder, Q);
                else if (Q != null && nfac >= 0)
                {
                    Order.Reorder(n, mainorder, Q);
                    Order.Reorder_gen(n, mainorder, Q, nfac, 1, true, n);
                }
                if (debugLevel == 2)
                {
                    ActiveSet.Optimise.printV("L end before convert", L, -1, n - nfixed);
                    ActiveSet.Optimise.printV("U end before convert", U, -1, n - nfixed);
                }
                for (i = 0; i < m; ++i)
                {
                    boundLU[i] = BlasLike.ddot(nfixed, A, m, L, 1, (n - nfixed) * m + i, n - nfixed);
                }
                BlasLike.dzerovec(n - nfixed, fixedW);
                BlasLike.dcopyvec(nfixed, L, fixedW, n - nfixed, n - nfixed);
                Order.bound_reorganise(1, n, n - nfixed, m, L);
                if (debugLevel == 2) ActiveSet.Optimise.printV("L end", L, -1, n - nfixed);
                Order.bound_reorganise(1, n, n - nfixed, m, U);
                if (debugLevel == 2) ActiveSet.Optimise.printV("U end", U, -1, n - nfixed);
                var nfixedo = nfixed;
                nfixed = 0;
                hessmull(n, Q, fixedW, fixedSecondOrder);
                fixedVariance = 0.5 * BlasLike.ddotvec(n, fixedW, fixedSecondOrder);
                nfixed = nfixedo;
                n -= nfixed;
                BlasLike.dsubvec(m, L, boundLU, L, n, 0, n);
                BlasLike.dsubvec(m, U, boundLU, U, n, 0, n);
            }
            if (delta < 0) delta = 2;
            // kappa = -1;
            var useLS = (value > 0 || valuel >= 0 || rmax > 0 || rmin > 0 || nabs > 0 || mabs > 0);
            var useCosts = kappa > 0.0 && buy != null && sell != null;
            if (!useCosts) kappa = 0;
            var buysellvars = delta < 2 || useCosts;
            var buysellI = 0;
            var longshortI = 0;
            var longshortbuysell = 0;
            for (var i = 0; i < n; ++i)
            {
                if (buysellvars && (initial[i] > L[i] && initial[i] < U[i])) buysellI++;
                if (useLS && 0 > L[i] && 0 < U[i])
                {
                    longshortI++;
                    if (buysellI > 0 && Math.Abs(initial[i]) < BlasLike.lm_eps) longshortbuysell++;
                }
            }
            var longshortbuysellIndex = new int[longshortbuysell];
            var longshortbuysellIndex_inverse = new int[n];
            for (var i = 0; i < n; ++i) longshortbuysellIndex_inverse[i] = -1;
            var buysellIndex = new int[buysellI];
            var buysellIndex_inverse = new int[n];
            for (var i = 0; i < n; ++i) buysellIndex_inverse[i] = -1;
            var longshortIndex = new int[longshortI - longshortbuysell];
            var longshortIndex_inverse = new int[n];
            for (var i = 0; i < n; ++i) longshortIndex_inverse[i] = -1;
            buysellI = 0;
            longshortI = 0;
            longshortbuysell = 0;
            for (var i = 0; i < n; ++i)
            {
                if (buysellvars && (initial[i] > L[i] && initial[i] < U[i])) buysellIndex[buysellI++] = i;
                if (useLS && 0 > L[i] && 0 < U[i])
                {
                    if (buysellI > 0 && Math.Abs(initial[i]) < BlasLike.lm_eps) longshortbuysellIndex[longshortbuysell++] = i;
                    else longshortIndex[longshortI++] = i;
                }
            }
            for (var i = 0; i < buysellI; ++i) buysellIndex_inverse[buysellIndex[i]] = i;
            for (var i = 0; i < longshortI; ++i) longshortIndex_inverse[longshortIndex[i]] = i;
            for (var i = 0; i < longshortbuysell; ++i) longshortbuysellIndex_inverse[longshortbuysellIndex[i]] = i;
            var N = n + buysellI + longshortI;
            var M = m + buysellI + longshortI + (delta < 2.0 ? 1 : 0);
            if (value > 0) M++;
            if (rmax > 0 && rmin == rmax) M++;
            else if (rmax > 0) M++;
            if (rmin > 0 && rmin != rmax) M++;
            M += nabs;
            M += mabs;
            var CC = new double[N];
            var AA = new double[N * M];
            var LL = new double[N + M];
            var UU = new double[N + M];
            var WW = new double[N];
            //Variables
            BlasLike.dcopyvec(n, L, LL);
            BlasLike.dcopyvec(n, U, UU);
            BlasLike.dsetvec(buysellI, 0, LL, n);
            BlasLike.dsetvec(buysellI, useIP ? BlasLike.lm_max : 1, UU, n);
            BlasLike.dsetvec(longshortI, 0, LL, n + buysellI);
            BlasLike.dsetvec(longshortI, useIP ? BlasLike.lm_max : 1, UU, n + buysellI);
            //Constraints
            BlasLike.dcopyvec(m, L, LL, n, N);
            BlasLike.dcopyvec(m, U, UU, n, N);
            for (var i = 0; i < m; i++)
            {
                BlasLike.dcopy(n, A, m, AA, M, i, i);
            }
            for (var i = 0; i < buysellI; ++i)
            {
                LL[N + m + i] = initial[buysellIndex[i]];
            }
            for (var i = 0; i < longshortI; ++i)
            {
                LL[N + m + i + buysellI] = 0;
            }
            if (useIP) BlasLike.dsetvec(buysellI, BlasLike.lm_max, UU, N + m);
            else BlasLike.dsetvec(buysellI, 1, UU, N + m);
            if (useIP) BlasLike.dsetvec(longshortI, BlasLike.lm_max, UU, N + m + buysellI);
            else BlasLike.dsetvec(longshortI, 1, UU, N + m + buysellI);
            for (var i = m; i < m + buysellI; ++i)
            {
                var ind = buysellIndex[i - m];
                BlasLike.dset(1, 1.0, AA, M, i + M * ind);
                BlasLike.dset(1, 1.0, AA, M, i + M * (n + i - m));
            }
            for (var i = m + buysellI; i < m + buysellI + longshortI; ++i)
            {
                var ind = longshortIndex[i - m - buysellI];
                BlasLike.dset(1, 1.0, AA, M, i + M * ind);
                BlasLike.dset(1, 1.0, AA, M, i + M * (n + i - m));
            }
            var mult = kappa / (1.0 - kappa);
            BlasLike.dsccopyvec(n, -gamma / (1 - gamma), alpha, CC);
            BlasLike.daddvec(n, CC, fixedSecondOrder, CC);
            if (buysellvars && useCosts)
            {
                for (var i = 0; i < n; ++i)
                {
                    if (initial[i] <= L[i]) CC[i] += mult * buy[i];
                    else if (initial[i] >= U[i]) CC[i] -= mult * sell[i];
                    else if (initial[i] > L[i] && initial[i] < U[i]) CC[i] += mult * buy[i];
                }
                for (var i = 0; i < buysellI; ++i)
                {
                    var ind = buysellIndex[i];
                    CC[n + i] = mult * (buy[ind] + sell[ind]);
                }
            }
            else if (buysellvars)
            {
                BlasLike.dsetvec(buysellI, 0, CC, n);
            }
            var cnum = m + buysellI + longshortI;
            var cnumTurn = -1;
            var forcedTurn = 0.0;
            var forcedI = 0.0;
            var fixedTurn = 0.0;
            for (var i = 0; i < nfixed; ++i)
            {
                fixedTurn += Math.Abs(U[nfixed - i - 1 + n + m] - initial[i + n]) / 2.0;
            }
            for (var i = 0; i < n; ++i)
            {
                if (buysellIndex_inverse[i] == -1)
                {
                    if (initial[i] >= U[i])
                    {
                        forcedI += -initial[i];
                        forcedTurn += U[i] - initial[i];
                    }
                }
            }
            if (delta < 2.0)
            {
                cnumTurn = cnum;
                for (var i = 0; i < n; ++i)
                {
                    if (buysellIndex_inverse[i] == -1)
                    {
                        if (initial[i] <= L[i])
                            BlasLike.dset(1, 1.0, AA, M, cnum + M * i);
                        else if (initial[i] >= U[i])
                        {
                            BlasLike.dset(1, -1.0, AA, M, cnum + M * i);
                        }
                    }
                    else
                        BlasLike.dset(1, 1.0, AA, M, cnum + M * i);//sum w =sum buy+ initial + initial-sell
                }
                BlasLike.dset(buysellI, 2.0, AA, M, cnum + M * n);//2sum sell
                LL[N + cnum] = BlasLike.dsumvec(n, initial) + 2 * forcedI;
                UU[N + cnum] = 2.0 * (delta - fixedTurn) + BlasLike.dsumvec(n, initial) + 2 * forcedI;
                cnum++;
            }
            var extraLong = 0.0;
            var extraShort = 0.0;
            for (var i = nfixed - 1; i >= 0; --i)
            {
                extraLong += Math.Max(0, L[i + n + m]);
                extraShort += Math.Min(0, L[i + n + m]);
            }
            var cnumVal = -1;
            if (value > 0)
            {
                cnumVal = cnum;
                LL[N + cnum] = valuel - extraLong;
                UU[N + cnum] = value - extraLong;
                // BlasLike.dset(n, 1.0, AA, M, cnum);//L+S
                for (var i = 0; i < n; ++i)
                {
                    if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                    {
                        if (LL[i] >= 0) BlasLike.dset(1, 1, AA, M, cnum + i * M);
                        else if (UU[i] <= 0) BlasLike.dset(1, 0, AA, M, cnum + i * M);
                    }
                    else
                    {
                        BlasLike.dset(1, 1, AA, M, cnum + i * M);
                    }
                }
                BlasLike.dset(longshortI, 1.0, AA, M, cnum + M * (n + buysellI));//-S
                for (var i = 0; i < longshortbuysell; ++i)
                {
                    var ii = longshortbuysellIndex[i];
                    ii = buysellIndex_inverse[ii];
                    BlasLike.dset(1, 1, AA, M, cnum + M * (n + ii));
                }
                cnum++;
            }
            var cnumRminmax = -1;
            var cnumRmax = -1;
            if (rmax > 0 && rmax == rmin)
            {//rmax=-S/L i.e rmax*L+S=0
                cnumRminmax = cnum;
                LL[N + cnum] = -extraLong * rmax - extraShort;
                UU[N + cnum] = -extraLong * rmax - extraShort;
                //   BlasLike.dset(n, rmax, AA, M, cnum);//rmax*(L+S)
                for (var i = 0; i < n; ++i)
                {
                    if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                    {
                        if (LL[i] >= 0) BlasLike.dset(1, rmax, AA, M, cnum + i * M);
                        else if (UU[i] <= 0) BlasLike.dset(1, 1.0, AA, M, cnum + i * M);
                    }
                    else
                    {
                        BlasLike.dset(1, rmax, AA, M, cnum + i * M);
                    }
                }
                BlasLike.dset(longshortI, rmax - 1.0, AA, M, cnum + M * (n + buysellI));//-S*rmax+S
                for (var i = 0; i < longshortbuysell; ++i)
                {
                    var ii = longshortbuysellIndex[i];
                    ii = buysellIndex_inverse[ii];
                    BlasLike.dset(1, rmax - 1.0, AA, M, cnum + M * (n + ii));
                }
                cnum++;
            }
            else if (rmax > 0)
            {//-S/L<rmax rmax*L+S>0
                cnumRmax = cnum;
                LL[N + cnum] = -extraLong * rmax - extraShort;
                UU[N + cnum] = useIP ? BlasLike.lm_max : 10;
                //  BlasLike.dset(n, rmax, AA, M, cnum);
                for (var i = 0; i < n; ++i)
                {
                    if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                    {
                        if (LL[i] >= 0) BlasLike.dset(1, rmax, AA, M, cnum + i * M);
                        else if (UU[i] <= 0) BlasLike.dset(1, 1.0, AA, M, cnum + i * M);
                    }
                    else
                    {
                        BlasLike.dset(1, rmax, AA, M, cnum + i * M);
                    }
                }
                BlasLike.dset(longshortI, rmax - 1.0, AA, M, cnum + M * (n + buysellI));
                for (var i = 0; i < longshortbuysell; ++i)
                {
                    var ii = longshortbuysellIndex[i];
                    ii = buysellIndex_inverse[ii];
                    BlasLike.dset(1, rmax - 1.0, AA, M, cnum + M * (n + ii));
                }
                cnum++;
            }
            var cnumRmin = -1;
            if (rmin > 0 && rmin != rmax)
            {//-S/L>rmin rmin*L+S<0
                cnumRmin = cnum;
                LL[N + cnum] = extraLong * rmin + extraShort;
                UU[N + cnum] = useIP ? BlasLike.lm_max : 10;
                //  BlasLike.dset(n, -rmin, AA, M, cnum);
                for (var i = 0; i < n; ++i)
                {
                    if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                    {
                        if (LL[i] >= 0) BlasLike.dset(1, -rmin, AA, M, cnum + i * M);
                        else if (UU[i] <= 0) BlasLike.dset(1, -1.0, AA, M, cnum + i * M);
                    }
                    else
                    {
                        BlasLike.dset(1, -rmin, AA, M, cnum + i * M);
                    }
                }
                BlasLike.dset(longshortI, -(rmin - 1.0), AA, M, cnum + M * (n + buysellI));
                for (var i = 0; i < longshortbuysell; ++i)
                {
                    var ii = longshortbuysellIndex[i];
                    ii = buysellIndex_inverse[ii];
                    BlasLike.dset(1, -(rmin - 1.0), AA, M, cnum + M * (n + ii));
                }
                cnum++;
            }
            if (nabs > 0)
            {
                for (var con = 0; con < nabs; ++con)
                {
                    var Lfixed = 0.0;
                    for (var i = 0; i < nfixed; ++i)
                    {
                        Lfixed += Math.Abs(L[n + m + nfixed - i - 1] * A_abs[(n + i) * nabs + con]);
                    }
                    LL[N + cnum] = L_abs[con] - Lfixed;
                    UU[N + cnum] = U_abs[con] - Lfixed;
                    for (var i = 0; i < n; ++i)
                    {
                        if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                        {
                            if (LL[i] >= 0) BlasLike.dcopy(1, A_abs, nabs, AA, M, con + i * nabs, cnum + i * M);//L
                            else if (UU[i] <= 0) BlasLike.dsccopy(1, -1, A_abs, nabs, AA, M, con + i * nabs, cnum + i * M);//-S
                        }
                        else
                        {
                            BlasLike.dcopy(1, A_abs, nabs, AA, M, con + i * nabs, cnum + i * M);//L+S
                        }
                    }
                    for (var i = 0; i < longshortI; ++i)
                    {
                        var ind = longshortIndex[i];
                        BlasLike.dset(1, 2.0 * A_abs[ind * nabs + con], AA, M, cnum + M * (n + buysellI + i));//-2S
                    }
                    for (var i = 0; i < longshortbuysell; ++i)
                    {
                        var ii = longshortbuysellIndex[i];
                        ii = buysellIndex_inverse[ii];
                        BlasLike.dset(1, 2.0 * A_abs[ii * nabs + con], AA, M, cnum + M * (n + ii));//-2S
                    }
                    cnum++;
                }
            }
            var cnumGross = -1;
            if (mabs > 0)
            {
                for (var con = 0; con < mabs; ++con)
                {
                    var Lfixed = 0.0;
                    var basecon = I_a[con];
                    for (var i = 0; i < nfixed; ++i)
                    {
                        Lfixed += Math.Abs(L[n + m + nfixed - i - 1] * A[(n + i) * m + basecon]);
                    }
                    if (basecon == 0) cnumGross = cnum;
                    LL[N + cnum] = L_abs[con + nabs] - Lfixed;
                    UU[N + cnum] = U_abs[con + nabs] - Lfixed;
                    for (var i = 0; i < n; ++i)
                    {
                        if (longshortIndex_inverse[i] == -1 && longshortbuysellIndex_inverse[i] == -1)
                        {
                            if (LL[i] >= 0) BlasLike.dcopy(1, A, m, AA, M, basecon + i * m, cnum + i * M);//L
                            else if (UU[i] <= 0) BlasLike.dsccopy(1, -1, A, m, AA, M, basecon + i * m, cnum + i * M);//-S
                        }
                        else
                        {
                            BlasLike.dcopy(1, A, m, AA, M, basecon + i * m, cnum + i * M);//L+S
                        }
                    }
                    for (var i = 0; i < longshortI; ++i)
                    {
                        var ind = longshortIndex[i];
                        BlasLike.dset(1, 2.0 * A[ind * m + basecon], AA, M, cnum + M * (n + buysellI + i));//-2S
                    }
                    for (var i = 0; i < longshortbuysell; ++i)
                    {
                        var ii = longshortbuysellIndex[i];
                        ii = buysellIndex_inverse[ii];
                        BlasLike.dset(1, 2.0 * A[ii * m + basecon], AA, M, cnum + M * (n + ii));//-2S
                    }
                    cnum++;
                }
            }
            this.names = names;
            this.L = LL;
            this.U = UU;
            this.A = AA;
            this.n = N;
            this.m = M;
            this.initial = initial;
            this.gamma = gamma;
            this.c = CC;
            if (useIP)
            {
                var LLL = new double[N + M];
                back = InteriorOpt(1e-11, WW, LLL);
            }
            else
            {
                var LAMBDAS = new double[N + M];
                for (var i = 0; i < n; ++i)
                {
                    if (UU[i] > 0)
                        WW[i] = UU[i] / n;
                    else
                        WW[i] = 0.5 * (UU[i] + LL[i]) / n;
                }
                for (var i = 0; i < buysellI; ++i)
                {
                    var ind = buysellIndex[i];
                    WW[i + n] = Math.Abs(initial[ind]) < BlasLike.lm_eps ? 0 : Math.Max(0, (initial[ind] - WW[ind]));
                }
                for (var i = 0; i < longshortI; ++i)
                {
                    var ind = longshortIndex[i];
                    WW[i + n + buysellI] = 0;
                }
                this.w = WW;
                WriteInputs("./optinput2");
                back = ActiveOpt(0, WW, LAMBDAS);
                Console.WriteLine($"back = {back}");
            }
            ColourConsole.WriteLine("_______________________________________________________________________________________________________________________", ConsoleColor.Green);
            ColourConsole.WriteEmbeddedColourLine($"[yellow]{"Asset",12}[/yellow]\t[cyan]{"WEIGHT-INITIAL or WEIGHT",25}[/cyan]\t[red]{"SELL or SHORT",12}[/red]\t[darkcyan]{"BUY or LONG",12}[/darkcyan]\t[green]{"INITIAL or 0",12}[/green]\t\t[darkmagenta]{"LIMIT",12}[/darkmagenta]");
            for (var i = 0; i < n; ++i)
            {
                if (buysellvars && buysellIndex_inverse[i] != -1)
                {
                    var ind = buysellIndex_inverse[i];
                    var c1 = BlasLike.ddot(N, AA, M, WW, 1, ind + m);
                    if (Math.Abs(WW[i] - initial[i]) > 1e-6)
                    {
                        if (WW[i] <= initial[i]) ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i],12}[/yellow]\t[cyan]{(WW[i] - initial[i]),25:F8}[/cyan]\t[red]{WW[ind + n],12:F8}[/red]\t[darkcyan]{(c1 - initial[i]),12:F8}[/darkcyan]\t[green]{initial[i],12:F8}[/green]\t\t[darkmagenta]{(!useIP ? (UU[ind + N + m] - c1) : 10),12:f2}[/darkmagenta]");
                        else ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i],12}[/yellow]  [cyan]{(WW[i] - initial[i]),25:F8}[/cyan]\t[red]{WW[ind + n],12:F8}[/red]\t[darkcyan]{(c1 - initial[i]),12:F8}[/darkcyan]\t[green]{initial[i],12:F8}[/green]\t\t[darkmagenta]{(!useIP ? (UU[ind + N + m] - c1) : 10),12:f2}[/darkmagenta]");
                    }
                }
                if (longshortI > 0 && longshortIndex_inverse[i] != -1)
                {
                    var ind = longshortIndex_inverse[i];
                    var c1 = BlasLike.ddot(N, AA, M, WW, 1, ind + m + buysellI);
                    if (Math.Abs(WW[i]) > 1e-6)
                    {
                        if (WW[i] <= 0) ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i],12}[/yellow]\t[cyan]{(WW[i]),25:F8}[/cyan]\t[red]{WW[ind + n + buysellI],12:F8}[/red]\t[darkcyan]{(c1),12:F8}[/darkcyan]\t[green]{0,12:F8}[/green]\t\t[darkmagenta]{(!useIP ? (UU[ind + N + m + buysellI] - c1) : 10),12:f2}[/darkmagenta]");
                        else ColourConsole.WriteEmbeddedColourLine($"[yellow]{names[i],12}[/yellow] [cyan]{(WW[i]),25:F8}[/cyan]\t[red]{WW[ind + n + buysellI],12:F8}[/red]\t[darkcyan]{(c1),12:F8}[/darkcyan]\t[green]{0,12:F8}[/green]\t\t[darkmagenta]{(!useIP ? (UU[ind + N + m + buysellI] - c1) : 10),12:f2}[/darkmagenta]");
                    }
                }
            }
            ColourConsole.WriteLine("_______________________________________________________________________________________________________________________", ConsoleColor.Green);
            if (cnumTurn != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Turnover constraint:[/darkyellow]\t[red]{LL[N + cnumTurn],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumTurn),20:f16}[/cyan]\t[green]{UU[N + cnumTurn],20:f16}[/green]");
            if (cnumVal != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Long Value constraint:[/darkyellow]\t[red]{LL[N + cnumVal],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumVal),20:f16}[/cyan]\t[green]{UU[N + cnumVal],20:f16}[/green]");
            if (cnumGross != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Gross Value constraint:[/darkyellow]\t[red]{LL[N + cnumGross],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumGross),20:f16}[/cyan]\t[green]{UU[N + cnumGross],20:f16}[/green]");
            if (cnumRminmax != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Rminmax constraint:[/darkyellow]\t[red]{LL[N + cnumRminmax],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumRminmax),20:f16}[/cyan]\t[green]{UU[N + cnumRminmax],20:f16}[/green]");
            if (cnumRmin != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Rmin constraint:[/darkyellow]\t\t[red]{LL[N + cnumRmin],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumRmin),20:f16}[/cyan]\t[green]{UU[N + cnumRmin],20:f16}[/green]");
            if (cnumRmax != -1) ColourConsole.WriteEmbeddedColourLine($"[darkyellow]Test Rmax constraint:[/darkyellow]\t\t[red]{LL[N + cnumRmax],20:f16}[/red]\t[cyan]{BlasLike.ddot(N, AA, M, WW, 1, cnumRmax),20:f16}[/cyan]\t[green]{UU[N + cnumRmax],20:f16}[/green]");
            var turn2 = fixedTurn;
            if (buysellI > 0) turn2 += (-forcedTurn + BlasLike.dsumvec(buysellI, WW, n) + (BlasLike.dsumvec(n, WW) - BlasLike.dsumvec(n, initial)) * 0.5);
            var shortsideS = -extraShort;
            for (var i = 0; i < n; ++i)
            {
                var cc = 0;
                if ((longshortIndex_inverse[i] == -1 || longshortbuysellIndex_inverse[i] == -1) && UU[i] <= 0)
                    shortsideS -= WW[i];
                if ((cc = longshortbuysellIndex_inverse[i]) != -1)
                    shortsideS += WW[buysellIndex_inverse[i] + n];
            }
            if (longshortI > 0) shortsideS += BlasLike.dsumvec(longshortI, WW, n + buysellI);
            wback = new double[n + nfixed];
            BlasLike.dcopyvec(n, WW, wback);
            var alphaFixed = 0.0;
            if (nfixed > 0)
            {
                BlasLike.daddvec(m, L, boundLU, L, n, 0, n);
                BlasLike.daddvec(m, U, boundLU, U, n, 0, n);
                n += nfixed;
                if (debugLevel == 2)
                {
                    ActiveSet.Optimise.printV("L end before covert back", L, -1, n - nfixed);
                    ActiveSet.Optimise.printV("U end before covert back", U, -1, n - nfixed);
                }
                Order.bound_reorganise(0, n, n - nfixed, m, L);
                if (debugLevel == 2) ActiveSet.Optimise.printV("L end", L, -1, n - nfixed);
                Order.bound_reorganise(0, n, n - nfixed, m, U);
                if (debugLevel == 2) ActiveSet.Optimise.printV("U end", U, -1, n - nfixed);
                BlasLike.dcopyvec(nfixed, L, wback, n - nfixed, n - nfixed);
                alphaFixed = BlasLike.ddotvec(nfixed, alpha, wback, n - nfixed, n - nfixed);
                Order.Reorder(n, mainorderInverse, wback);
                Order.Reorder(n, mainorderInverse, L);
                Order.Reorder(n, mainorderInverse, U);
                Order.Reorder(n, mainorderInverse, alpha);
                Order.Reorder(n, mainorderInverse, initial);
                if (bench != null) Order.Reorder(n, mainorderInverse, bench);
                if (buy != null) Order.Reorder(n, mainorderInverse, buy);
                if (sell != null) Order.Reorder(n, mainorderInverse, sell);
                if (names != null) Order.Reorder(n, mainorderInverse, names);
                Order.Reorder_gen(n, mainorderInverse, A, m, 1, true);
                if (A_abs != null) Order.Reorder_gen(n, mainorderInverse, A_abs, nabs, 1, true);
                if (Q != null && nfac == -1) Order.ReorderSymm(n, mainorderInverse, Q);
                else if (Q != null && nfac >= 0)
                {
                    Order.Reorder(n, mainorderInverse, Q);
                    Order.Reorder_gen(n, mainorderInverse, Q, nfac, 1, true, n);
                }
            }
            var eret = BlasLike.ddotvec(n, alpha, wback);
            var printAlphas = false;
            for (var i = 0; printAlphas && i < n; ++i)
            {
                ColourConsole.WriteEmbeddedColourLine($"[magenta]{names[i]}[/magenta] [green]{alpha[i],12:f8}[/green] [red]{wback[i],12:f8}[/red]");
            }
            var nfixedold = nfixed;
            nfixed = 0;
            var variance = Variance(wback);
            var risk = 0.0;
            if (bench != null)
            {
                BlasLike.dsubvec(n, wback, bench, wback);
                risk = Math.Sqrt(Variance(wback));
                BlasLike.daddvec(n, wback, bench, wback);
            }
            else risk = Math.Sqrt(variance);
            nfixed = nfixedold;
            var costbase = 0.0;
            var initbase = 0.0;
            var longside = BlasLike.dsumvec(n, wback);
            var shortside = 0.0;
            for (var i = 0; i < n; ++i)
            {
                if (wback[i] < 0) shortside += wback[i];
            }
            longside -= shortside;
            if (buy != null && sell != null)
                for (var i = 0; i < n; ++i)
                {
                    if (U[i] == L[i]) continue;
                    if (initial[i] > U[i])
                    {
                        costbase -= sell[i] * wback[i];
                        initbase -= sell[i] * initial[i];
                    }
                    else
                    {
                        costbase += buy[i] * wback[i];
                        initbase += buy[i] * initial[i];
                    }
                }
            var eretA = alphaFixed * gamma / (1 - gamma) + BlasLike.ddotvec(n - nfixed, fixedSecondOrder, WW) - BlasLike.ddotvec(n - nfixed, CC, WW) + kappa / (1 - kappa) * costbase;
            if (gamma != 0) eretA /= gamma / (1 - gamma);
            var costA = 0.0;
            for (var i = 0; useCosts && i < buysellI; ++i)
            {
                var ind = buysellIndex[i];
                costA += WW[n - nfixed + i] * (buy[ind] + sell[ind]);
            }
            if (buysellI > 0) costA += costbase - initbase;
            ColourConsole.WriteEmbeddedColourLine($"[magenta]Gross Value:[/magenta]\t\t\t[darkyellow]{longside - shortside,20:f16}[/darkyellow]");
            ColourConsole.WriteEmbeddedColourLine($"[green]Longside={longside,20:f16}[/green]\t[red]Shortside={shortside,20:f16}[/red] [magenta]({-shortsideS,20:f16})[/magenta]");
            ColourConsole.WriteEmbeddedColourLine($"[magenta]-Short/Long:[/magenta]\t\t\t[red]{rmin,20:f16}[/red]\t[cyan]{-shortside / longside,20:f16}[/cyan]\t[green]{rmax,20:f16}[/green]");
            ColourConsole.WriteEmbeddedColourLine($"Variance:\t\t\t[green]{variance,20:f16}[/green]");
            ColourConsole.WriteEmbeddedColourLine($"Risk:\t\t\t\t[green]{risk,20:f16}[/green]");
            ColourConsole.WriteEmbeddedColourLine($"Return:\t\t\t\t[green]{eret,20:f16}:[/green]\t[cyan]{eretA,20:f16}[/cyan]");
            var benchmarkExtra = 0.0;
            nfixedold = nfixed;
            nfixed = 0;
            if (Q != null)
            {
                var implied = new double[ntrue];
                if (bench != null)
                {
                    hessmull(n, Q, bench, implied);
                    for (var i = 0; i < ntrue; ++i)
                    {
                        if (U[i] != L[i])
                            benchmarkExtra -= wback[i] * implied[i];
                    }
                }
            }
            nfixed = nfixedold;
            var turnover = 0.0;
            var cost = 0.0;
            var costFixed = 0.0;
            for (var i = 0; i < n; ++i)
            {
                turnover += Math.Abs(wback[i] - initial[i]);
                if ((buy != null) && (sell != null))
                {
                    var diff = (wback[i] - initial[i]);
                    cost += diff > 0 ? diff * buy[i] : -diff * sell[i];
                    if (U[i] == L[i])
                    {
                        costFixed += diff > 0 ? diff * buy[i] : -diff * sell[i];
                    }
                }
            }
            var utility = -eret * gamma / (1 - gamma) + kappa / (1 - kappa) * (cost + initbase) + 0.5 * variance + benchmarkExtra + alphaFixed * gamma / (1 - gamma) - costFixed * kappa / (1 - kappa) - fixedVariance;
            var utilityA = -BlasLike.ddotvec(n - nfixed, fixedSecondOrder, WW) + BlasLike.ddotvec(N, CC, WW) + 0.5 * variance + benchmarkExtra - fixedVariance;
            ColourConsole.WriteEmbeddedColourLine($"Utility:\t\t\t[green]{utility,20:f16}:[/green]\t[cyan] {utilityA,20:f16}[/cyan]");
            ColourConsole.WriteEmbeddedColourLine($"Turnover:\t\t\t[green]{turnover * 0.5,20:f16}:[/green]\t[cyan]{turn2,20:f16}[/cyan]");
            ColourConsole.WriteEmbeddedColourLine($"Cost:\t\t\t\t[green]{cost,20:f16}:[/green]\t[cyan]{costA + costFixed,20:f16}[/cyan]");
            for (var i = 0; i < m; ++i)
            {
                var ccval = BlasLike.ddot(n, A, m, wback, 1, i);
                ColourConsole.WriteEmbeddedColourLine($"[magenta]Portfolio constraint {(i + 1),3}:[/magenta]\t[cyan]{ccval,20:f16}[/cyan]\t([red]{L[i + n],20:f16},{U[i + n],20:f16}[/red])");
            }
            //            ActiveSet.Optimise.printV("optimal weights", WW, n);
            if (buysellI > 0 || longshortI > 0)
            {
                if (Math.Abs(shortside + shortsideS) > BlasLike.lm_rooteps * 2)
                    back = 6;
                if (Math.Abs(utility - utilityA) > BlasLike.lm_rooteps)
                    back = 6;
                if (kappa > 1e-14 && Math.Abs(cost - costA - costFixed) > BlasLike.lm_eps * 10)
                    back = 6;
            }
            BACK = back;
            return back;
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
                Q = new double[n * (n + 1) / 2];
                for (int i = 0, ij = 0; i < n; ++i, ij += i)
                {
                    Q[ij + i] = 10;
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
            this.names = names;
            this.L = LL;
            this.U = UU;
            this.A = AA;
            this.n = N;
            this.m = M;
            this.gamma = 0.5;
            this.c = cc;
            double alphamax = 0, alphamin = 0;
            BlasLike.dxminmax(alpha.Length, alpha, 1, ref alphamax, ref alphamin);
            Console.WriteLine($"c range ({alphamax},{alphamin}) ratio {(alphamin / alphamax):E8}");
            if (useIP)
            {
                this.ntrue = n;
                this.mtrue = m;
                var LLL = new double[N + M];
                var back = InteriorOpt(BlasLike.lm_eps2, ww, LLL);
            }
            else
            {
                BlasLike.dsetvec(n, 1.0 / n, ww);
                for (var i = 0; i < tlen; ++i)
                {
                    var LOSS = R - BlasLike.ddot(n, DATA, tlen, ww, 1, i);
                    ww[i + n] = Math.Max(0, LOSS);
                }
                this.w = ww;
                ntrue = n;
                mtrue = m;
                WriteInputs("./gainloss");
                var LAMBDAS = new double[N + M];
                var back = ActiveOpt(activeLp, ww, LAMBDAS);
            }
            var gain = 0.0;
            var loss = 0.0;
            for (var i = 0; i < tlen; ++i)
            {
                var GL = BlasLike.ddot(n, DATA, tlen, ww, 1, i) - R;
                gain += Math.Max(0, GL);
                loss += -Math.Min(0, GL);
            }
            var lossV = BlasLike.dsumvec(tlen, ww, n);
            var gainV = 0.0;
            for (var i = 0; i < tlen; ++i)
            {
                gainV += BlasLike.ddot(N, AA, M, ww, 1, i + m) - R;
            }
            ColourConsole.WriteEmbeddedColourLine($"[green]Total GAIN = \t\t\t\t{gain,12:F8}[/green]");
            ColourConsole.WriteEmbeddedColourLine($"[green]Total GAIN (check from opt variables) = {gainV,12:F8}[/green]");
            ColourConsole.WriteEmbeddedColourLine($"[red]Total LOSS = \t\t\t\t{loss,12:F8}[/red]");
            ColourConsole.WriteEmbeddedColourLine($"[red]Total LOSS (check from opt variables) = {lossV,12:F8}[/red]");
            var variance = Variance(ww);
            var expret = BlasLike.ddotvec(n, ww, alpha);
            Console.WriteLine($"Return {expret,12:F8}");
            Console.WriteLine($"Variance {variance,12:F8}");
            Console.WriteLine($"Utility {-expret + lossV * lambda / tlen + 0.5 * variance}");
            for (var i = 0; i < n; ++i)
            {
                Console.WriteLine($"{names[i]}\t{ww[i],12:F8}");
            }
            for (var i = 0; i < M; ++i)
            {
                var constraint = BlasLike.ddot(N, AA, M, ww, 1, i);
                if (i < m) ColourConsole.WriteEmbeddedColourLine($"[cyan]Constraint {i,3}[/cyan] = [magenta]{constraint,12:F8}[/magenta]");
                else if (useIP) ColourConsole.WriteEmbeddedColourLine($"[cyan]Constraint {i,3}[/cyan] = [magenta]{constraint,12:F8}[/magenta][red]  LOSS var {ww[n + i - 1],12:F8}[/red]");
                else ColourConsole.WriteEmbeddedColourLine($"[cyan]Constraint {i,3}[/cyan] = [magenta]{constraint,12:F8}[/magenta]  {UU[N + i]} [red]LOSS var {ww[n + i - 1],12:F8}[/red]  {UU[n + i - 1]}");
            }
            for (var i = 0; i < m; ++i)
            {
                var ccval = BlasLike.dsumvec(n, ww);
                ColourConsole.WriteEmbeddedColourLine($"Portfolio constraint {i,3}: [cyan]{ccval}[/cyan]");
            }

        }
        public virtual void OptimiseTest()
        {
            if (ntrue == 0) ntrue = n;
            if (mtrue == 0) mtrue = m;
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
                        AA[i * (m + 1) + 1] = -c[i];
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
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            WriteInputs("./basic");
            var llambda = new double[n + m];
            L[n + 1] = -10;
            var ok = ActiveOpt(0, w, llambda);
            ActiveSet.Optimise.printV("w from Active Set", w);
            Console.WriteLine($"Variance from Active Set:\t\t{Variance(w)}");
            Factorise.dmxmulv(m, n, A, w, Aw);
            ActiveSet.Optimise.printV("Constraints", Aw);
            L[n + 1] = -BlasLike.lm_max;
            var LLL = new double[n + m];
            var ww = new double[n];
            var ip = InteriorOpt(BlasLike.lm_eps2, ww, LLL);
            Console.WriteLine($"Variance from IP:\t\t{Variance(ww)}");
            Factorise.dmxmulv(m, n, A, ww, Aw);
            ActiveSet.Optimise.printV("Constraints", Aw);
        }
        public string inFile = "";
        public int n;
        ///<summary>Number of assets with lower==upper</summary>
        public int nfixed = 0;
        public int debugLevel = 1;
        public double[] fixedW = null;
        public double fixedVariance = 0;
        public int ntrue = 0;
        public int mtrue = 0;
        public int m;
        public double gamma;
        public double kappa;
        public double delta;
        public double[] w = null;
        public double[] wback = null;
        public double[] initial = null;
        public double[] bench = null;
        public double[] L = null;
        public double[] U = null;
        public double[] buy = null;
        public double[] sell = null;
        public double[] A = null;
        public double[] c = null;
        public double[] Q = null;
        public string[] names;
        public virtual int makeQ()
        {
            var nn = ntrue * (ntrue + 1) / 2;
            if (Q.Length == nn)
                return 0;
            else
                return -10;
        }
        public void hessmulltest(int nn, double[] QQ, double[] x, double[] hx)
        {
            BlasLike.dzerovec(nn, hx);
        }
        public void hessmulltest(int nn, int nrowh, int ncolh, int j, double[] QQ, double[] x, double[] hx)
        {
            BlasLike.dzerovec(nn, hx);
        }
        public virtual void hessmull(int nn, int nrowh, int ncolh, int j, double[] QQ, double[] x, double[] hx)
        {
            if (Q != null)
            {
                Factorise.CovMul(ntrue, Q, x, hx, 0, 0, 0, 'U', nfixed);
                BlasLike.dzerovec(nn - ntrue + nfixed, hx, ntrue - nfixed);
            }
            else BlasLike.dzerovec(nn, hx);
        }
        public virtual void hessmull(int nn, double[] QQ, double[] x, double[] hx)
        {
            if (Q != null)
            {
                Factorise.CovMul(ntrue, Q, x, hx, 0, 0, 0, 'U', nfixed);
                BlasLike.dzerovec(nn - ntrue + nfixed, hx, ntrue - nfixed);
            }
            else BlasLike.dzerovec(nn, hx);
        }
        ///<Summary>Set the upper and lower bounds to allow only the long/short and/or
        ///buy/sell value given by w</Summary>
        ///<param name="n">Number of assets</param>
        ///<param name="L">Lower bound array</param>
        ///<param name="U">Upper bound array</param>
        ///<param name="initial">Initial weight array</param>
        ///<param name="w">weight array that defines buy/sell and/or long/short</param>
        ///<param name="setFixed">If true, set the bounds to constrain the fixed wieghts</param>
        public void BoundsSetToSign(int n, double[] L, double[] U, double[] initial, double[] w, bool setFixed = false, double tol = 1e-14)
        {
            for (var i = 0; i < n; ++i)
            {
                if (setFixed && Math.Abs(w[i]) < tol)
                {
                    U[i] = 0; L[i] = 0;
                }
                else if (setFixed && Math.Abs(w[i] - initial[i]) < tol)
                {
                    U[i] = initial[i]; L[i] = initial[i];
                }
                else
                {
                    if (L[i] < 0 && U[i] > 0)
                    {
                        if (w[i] >= 0) L[i] = 0;
                        else if (w[i] < 0) U[i] = 0;
                    }
                    if (L[i] < initial[i] && U[i] > initial[i])
                    {
                        if (w[i] >= initial[i]) L[i] = initial[i];
                        else if (w[i] < initial[i]) U[i] = initial[i];
                    }
                }
            }
        }
        public double Variance(double[] w)
        {
            var Qx = new double[w.Length];
            hessmull(w.Length, Q, w, Qx);
            return BlasLike.ddotvec(w.Length, w, Qx);
        }
        public int ActiveOpt(int lp = 0, double[] www = null, double[] LAM = null)
        {
            if (ntrue == 0) ntrue = n;
            if (mtrue == 0) mtrue = m;
            var obj = 0.0;
            var iter = 10;
            var cextra = new double[n + nfixed];
            var opt = new ActiveSet.Optimise();
            if (lp == 0) opt.h = hessmull;
            if (bench != null && lp == 0)
            {
                var nfixedo = nfixed;
                nfixed = 0;
                opt.h(ntrue, 0, 0, 0, Q, bench, cextra);
                BlasLike.dnegvec(ntrue - nfixed, cextra);
                nfixed = nfixedo;
                BlasLike.dzerovec(nfixed, cextra, ntrue - nfixed);
            }
            BlasLike.daxpyvec(n, 1.0, c, cextra);
            if (www == null)
            {
                w = new double[n];
                BlasLike.dsetvec(n, 1.0 / n, w);
            }
            else
                w = www;
            var back = opt.QPopt(n, m, w, L, U, A, cextra, Q, ref obj, ref iter, lp, LAM);
            if (LAM != null) UtilityAnalysis(LAM, cextra);
            Console.WriteLine($"objective:\t\t{obj}; {iter} iterations");
            return back;
        }

        public int InteriorOpt(double conv = 1e-16, double[] wback = null, double[] LLL = null)
        {
            if (ntrue == 0) ntrue = n;
            if (mtrue == 0) mtrue = m;
            if (true)
            {
                for (var i = 0; i < n; ++i)
                {
                    if (U[i] == 1 && (L[i] == 0 || L[i] == -1)) U[i] = BlasLike.lm_max;
                    else if (L[i] == -1 && U[i] == 0) L[i] = -BlasLike.lm_max;
                }
            }

            var slacklarge = 0;
            var cextra = new double[n + nfixed];
            var CTEST = new double[n + nfixed];
            var slackb = 0;
            var slackL = 0;
            var slackU = 0;
            for (var i = 0; i < n; ++i)
            {
                if (U[i] != BlasLike.lm_max && U[i] != 0) slacklarge++;
                else if (U[i] != BlasLike.lm_max && L[i] != -BlasLike.lm_max && L[i] != 0) slacklarge++;
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
                else if (U[i] != BlasLike.lm_max && L[i] != -BlasLike.lm_max && L[i] != 0)
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
            if (bench != null && Q != null)
            {
                var nfixedo = nfixed;
                nfixed = 0;
                hessmull(ntrue, Q, bench, cextra);
                BlasLike.dnegvec(ntrue - nfixed, cextra);
                nfixed = nfixedo;
                BlasLike.dzerovec(nfixed, cextra, ntrue - nfixed);
            }
            BlasLike.daxpyvec(n, 1.0, c, cextra);
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
                else if (L[i] < 0 && U[i] == 0)
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
            /*     sign = null;
                 for (var i = 0; i < cextra.Length; ++i)
                 {
                      CTEST[i] = Math.Abs(cextra[i]);
                     //CTEST[i] = cextra[i];
                 }*/
            w = new double[n];
            BlasLike.dsetvec(n, 1.0 / n, w);
            var HH = new double[ntrue * (ntrue + 1) / 2];
            if (Q != null && Q.Length == (ntrue * (ntrue + 1) / 2)) HH = Q;
            else if (Q != null) Factorise.Fac2Cov(ntrue, (int)(Q.Length / ntrue) - 1, Q, HH);
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
            BlasLike.dcopyvec(n, ww, w);
            if (true && LLL != null)
            {
                BlasLike.dsccopyvec(n, 1, IOPT.z, LLL);
                BlasLike.dsccopyvec(m, 1, IOPT.y, LLL, 0, n);
                for (var i = 0; i < slackb; ++i)
                {
                    LLL[n + slackToConstraintBOTH[i]] += IOPT.y[m + slacklarge + i];
                }
                for (var i = 0; i < slacklarge; ++i)
                {
                    LLL[slacklargeConstraint[i]] -= IOPT.z[n + i];
                }
                var oldQ = Q;
                Q = null;
                UtilityAnalysis(LLL, CTEST);
                Q = oldQ;
            }
            if (back < -10) Console.WriteLine($"Failed -- too many iterations");
            if (back < 0) Console.WriteLine($"Normal Matrix became ill-conditioned");
            if (back == 6) Console.WriteLine("INFEASIBLE");
            else if (Q != null)
            {
                IOPT = new InteriorPoint.Optimise(slacklarge + n + totalConstraintslack, m + slacklarge + slackb, ww, null, bb, cextra, ntrue, HH);
                IOPT.h = hessmull;
                IOPT.nh = ntrue - nfixed;
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
                Console.WriteLine(BlasLike.ddotvec(n, w, testmul));
                var kk = new Portfolio("");
                kk.ntrue = ntrue;
                kk.Q = HH;
                kk.nfixed = nfixed;
                kk.hessmull(n, HH, w, testmul);
                Console.WriteLine(BlasLike.ddotvec(n, w, testmul));
                IOPT.alphamin = 1e-8;
                IOPT.conv = conv;
                IOPT.compConv = Math.Max(conv, IOPT.compConv);
                back = IOPT.Opt("QP", null, null, false, UL, sign);
                BlasLike.dcopyvec(n, ww, w);
                if (true && LLL != null)
                {
                    BlasLike.dsccopyvec(n, 1, IOPT.z, LLL);
                    BlasLike.dsccopyvec(m, 1, IOPT.y, LLL, 0, n);
                    for (var i = 0; i < slackb; ++i)
                    {
                        LLL[n + slackToConstraintBOTH[i]] += IOPT.y[m + slacklarge + i];
                    }
                    for (var i = 0; i < slacklarge; ++i)
                    {
                        LLL[slacklargeConstraint[i]] -= IOPT.z[n + i];
                    }
                    UtilityAnalysis(LLL, cextra);
                }
                if (back < -10) Console.WriteLine($"Failed -- too many iterations");
                else if (back < 0) Console.WriteLine($"Normal Matrix became ill-conditioned");
            }
            if (wback != null) BlasLike.dcopyvec(wback.Length, ww, wback);
            return back;
        }
        public static void printVector<T>(string name, T[] a, StreamWriter dave, int linelimit = 1)
        {
            if (a == null) return;
            if (linelimit == 1) linelimit = a.Length; //Hack so that line with one item is not treated as a scalar
            dave.WriteLine(name);
            for (int i = 0; i < a.Length; ++i)
            {
                var p = a[i].GetType();
                if (p.FullName == "System.Double")
                    dave.Write($"{a[i],20:F16} ");
                else
                    dave.Write($"{a[i]} ");
                if (i % linelimit == (linelimit - 1)) dave.Write("\n");
            }
        }

        public static void printVector<T>(string name, T[] a, StreamWriter dave, bool uplo)
        {
            if (a == null) return;
            var n = (int)((-1 + Math.Sqrt(1 + 8 * a.Length)) * 0.5);
            dave.WriteLine(name);
            if (uplo)
                for (int i = 0, ij = 1, ic = 0; ij <= n; ++i, ic++)
                {
                    var p = a[ic].GetType();
                    if (p.FullName == "System.Double")
                        dave.Write($"{a[ic],20:F16} ");
                    else
                        dave.Write($"{a[ic]} ");
                    if (i % (ij) == (ij - 1)) { dave.Write("\n"); ij++; i = -1; }
                }
            else
                for (int i = 0, ij = n, ic = 0; ic < a.Length; ++i, ic++)
                {
                    var p = a[ic].GetType();
                    if (p.FullName == "System.Double")
                        dave.Write($"{a[ic],12:F8} ");
                    else
                        dave.Write($"{a[ic]} ");
                    if (i % (ij) == (ij - 1)) { dave.Write("\n"); ij--; i = -1; }
                }
        }
    }
    public class FPortfolio : Portfolio
    {
        public override void WriteInputs(string filename)
        {
            double[] HH = null;
            if (Q != null)
            {
                HH = new double[ntrue * (ntrue + 1) / 2];
                Factorise.Fac2Cov(ntrue, nfac, Q, HH);
            }
            using (StreamWriter writer = new StreamWriter(filename))
            {
                writer.WriteLine("n");
                writer.WriteLine(n);
                writer.WriteLine("ntrue");
                writer.WriteLine(ntrue);
                writer.WriteLine("m");
                writer.WriteLine(m);
                writer.WriteLine("mtrue");
                writer.WriteLine(mtrue);
                printVector("initial", initial, writer);
                printVector("WW", w, writer);
                printVector("LL", L, writer);
                printVector("UU", U, writer);
                printVector("bench", bench, writer);
                printVector("QQ", HH, writer, true);
                printVector("AA", A, writer, m);
                printVector("CC", c, writer);
            }
        }
        public FPortfolio(string file) : base(file)
        {
            inFile = file;
            if (inFile != "")
                using (var OptData = new InputSomeData())
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
            if (Q != null)
            {
                if (nfixed > 0)
                {
                    Factorise.FacMul(ntrue, nfac, Q, x, hx, 0, 0, 0, nfixed);
                }
                else
                    Factorise.FacMul(ntrue, nfac, Q, x, hx);
                BlasLike.dzerovec(nn - ntrue + nfixed, hx, ntrue - nfixed);
            }
            else BlasLike.dzerovec(nn, hx);
        }
        public override void hessmull(int nn, double[] QQ, double[] x, double[] hx)
        {
            if (Q != null)
            {
                if (nfixed > 0)
                {
                    Factorise.FacMul(ntrue, nfac, Q, x, hx, 0, 0, 0, nfixed);
                }
                else
                    Factorise.FacMul(ntrue, nfac, Q, x, hx);
                BlasLike.dzerovec(nn - ntrue + nfixed, hx, ntrue - nfixed);
            }
            else BlasLike.dzerovec(nn, hx);
        }
        public int nfac;
        public double[] FL = null;
        public double[] SV = null;
        public double[] FC = null;
        public override int makeQ()
        {
            var nn = (nfac + 1) * ntrue;
            Q = new double[nn];
            return Factorise.FMP(ntrue, nfac, FC, SV, FL, Q);
        }
    }
}

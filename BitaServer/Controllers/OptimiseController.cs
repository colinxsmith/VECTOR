using Microsoft.AspNetCore.Mvc;
using Blas;
using Solver;
using DataFile;

namespace BitaServer.Controllers;

[ApiController]
[Route("[controller]")]
public class OptimiseController : ControllerBase
{

    private readonly ILogger<OptimiseController> _logger;

    public OptimiseController(ILogger<OptimiseController> logger)
    {
        _logger = logger;
    }
    [HttpPost("test")]
    public Optimise[] Post(Optimise op)
    {
        if (op.digit != null)
        {
            op.tdigit = op.digit != null ? Portfolio.Portfolio.check_digit((double)op.digit) : null;
        }
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return new[] { op };
    }
    [HttpGet()]
    public IEnumerable<Optimise> Blank()
    {
        return new[] { new Optimise() };
    }
    [HttpGet("LOSS")]
    public Optimise[] GetLOSS()
    {
        var op = new Optimise();
        var targetR = 0.0;
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "DATA L U A R";
            CVarData.intFields = "n tlen m";
            CVarData.stringFields = "names";
            try
            {
                CVarData.Read("./GLdist");
            }
            catch
            {
                CVarData.Read("../GLdist");
            }
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.A = CVarData.mapDouble["A"];
            op.DATA = CVarData.mapDouble["DATA"];
            op.names = CVarData.mapString["names"];
            op.n = CVarData.mapInt["n"][0];
            op.m = CVarData.mapInt["m"][0];
            op.tlen = CVarData.mapInt["tlen"][0];
            targetR = CVarData.mapDouble["R"][0];
        }
        op.TargetReturn = new double[op.tlen];
        BlasLike.dsetvec(op.tlen, targetR, op.TargetReturn);
        op.LOSSopt = false;
        op.Gstrength = 1;
        //     op.tail=0.05;
        op.LOSSmin = -1;
        op.LOSSmax = 10;
        op.gamma = 0;
        return new[] { op };
    }

    [HttpPost("LOSS")]
    public Optimise[] PostLOSS(Optimise op)
    {
        if (!op.gamma.HasValue) op.gamma = 0.5;
        if (!op.kappa.HasValue) op.kappa = -1;
        var TR = new double[op.tlen];
        if (op.TargetReturn == null)
        {
            op.TargetReturn = TR;
        }
        double[] Q;
        if (op.Q == null)
        {
            Q = new double[op.n.GetValueOrDefault() * (op.n.GetHashCode() + 1) / 2]; op.Q = Q;
            var ij = 0;
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                for (var j = 0; j <= i; ++j)
                {
                    Q[ij++] = Solver.Factorise.covariance(op.tlen, op.DATA, op.DATA, i * op.tlen, j * op.tlen);
                }
            }
        }
        else Q = op.Q;
        var opt = new Portfolio.Portfolio("");
        var ones = new double[op.tlen];
        op.alpha = new double[op.n.GetValueOrDefault()];
        opt.Q = Q;
        BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
        Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, 0, 0, 0, true);
        op.back = opt.BasicOptimisation(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), -1,
           op.A, op.L, op.U, op.gamma.GetValueOrDefault(), op.kappa.GetValueOrDefault(), -1, -1, -1, -1, -1, op.alpha, op.initial, null, null,
           op.names, false, 0, null, null, null, 0, null, op.tlen,
            op.Gstrength.GetValueOrDefault(), op.DATA, 0.0, op.TargetReturn, op.LOSSopt.GetValueOrDefault(), op.LOSSmin.GetValueOrDefault(), op.LOSSmax.GetValueOrDefault());
        op.w = opt.wback;
        op.CVARGLprob = opt.CVARGLprob;
        op.message = Portfolio.Portfolio.OptMessages(Math.Abs(op.back.GetValueOrDefault()));
        op.breakdown = new double[op.n.GetValueOrDefault()];
        op.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, op.breakdown); op.mctr = new double[op.n.GetValueOrDefault()];
        var lcheck = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.breakdown);
        opt.RiskBreakdown(op.w, null, op.mctr);
        op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
        op.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        return new[] { op };
    }

    [HttpGet("ETL")]
    public Optimise[] GetETL()
    {
        var op = new Optimise();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "DATA L U A";//The DATA is LOSS i.e. - returns
            CVarData.intFields = "n tlen m";
            CVarData.stringFields = "names";
            try
            {
                CVarData.Read("./GLdist");
            }
            catch
            {
                CVarData.Read("../GLdist");
            }
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.A = CVarData.mapDouble["A"];
            op.DATA = CVarData.mapDouble["DATA"];
            op.names = CVarData.mapString["names"];
            op.n = CVarData.mapInt["n"][0];
            op.m = CVarData.mapInt["m"][0];
            op.tlen = CVarData.mapInt["tlen"][0];
        }
        BlasLike.dnegvec(op.DATA.Length, op.DATA);
        op.ETLopt = false;
        op.Gstrength = 1;
        //     op.tail=0.05;
        op.ETLmin = -1;
        op.ETLmax = 1;
        op.gamma = 0;
        return new[] { op };
    }
    [HttpPost("ETL")]
    public Optimise[] PostETL(Optimise op)
    {
        if (!op.gamma.HasValue) op.gamma = 0.5;
        if (!op.kappa.HasValue) op.kappa = -1;
        double[] Q;
        if (op.Q == null)
        {
            Q = new double[op.n.GetValueOrDefault() * (op.n.GetHashCode() + 1) / 2]; op.Q = Q;
            var ij = 0;
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                for (var j = 0; j <= i; ++j)
                {
                    Q[ij++] = Solver.Factorise.covariance(op.tlen, op.DATA, op.DATA, i * op.tlen, j * op.tlen);
                }
            }
        }
        else Q = op.Q;
        var opt = new Portfolio.Portfolio("");
        var ones = new double[op.tlen];
        op.alpha = new double[op.n.GetValueOrDefault()];
        opt.Q = Q;
        BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
        Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, 0, 0, 0, true);
        BlasLike.dnegvec(op.n.GetValueOrDefault(), op.alpha);
        op.back = opt.BasicOptimisation(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), -1,
           op.A, op.L, op.U, op.gamma.GetValueOrDefault(), op.kappa.GetValueOrDefault(), -1, -1, -1, -1, -1, op.alpha, op.initial, null, null,
           op.names, false, 0, null, null, null, 0, null, op.tlen,
            op.Gstrength.GetValueOrDefault(), op.DATA, op.tail,
            null, op.ETLopt.GetValueOrDefault(), op.ETLmin.GetValueOrDefault(),
            op.ETLmax.GetValueOrDefault());
        op.w = opt.wback;
        op.CVARGLprob = opt.CVARGLprob;
        op.message = Portfolio.Portfolio.OptMessages(Math.Abs(op.back.GetValueOrDefault()));
        op.breakdown = new double[op.n.GetValueOrDefault()];
        double VAR = 0;
        int ind = -1;
        op.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref ind, op.breakdown);
        op.VAR = VAR;
        op.VARindex = ind;
        op.mctr = new double[op.n.GetValueOrDefault()];
        opt.RiskBreakdown(op.w, null, op.mctr);
        op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
        op.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        return new[] { op };
    }
    [HttpGet("test")]
    public Optimise[] Get(double? d1)
    {
        double testdigit;
        if (d1 != null) testdigit = d1.GetValueOrDefault();
        else testdigit = 123.9999999999;
        var op = new Optimise();
        op.digit = testdigit;
        op.tdigit = Portfolio.Portfolio.check_digit(testdigit);
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return new[] { op };
    }
    [HttpGet("test/n")]
    public Optimise[] Getn()
    {
        double testdigit = 11.000000000001;
        var op = new Optimise();
        op.n = 500;
        op.digit = testdigit;
        op.tdigit = Portfolio.Portfolio.check_digit(testdigit);
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return new[] { op };
    }
    int Optimise_internalCVPAFbl(int n, int nfac, string[] names, double[] w, int m,
                                    double[] A, double[] L, double[] U, double[] alpha,
                                    double[] benchmark, double[] Q, double gamma, double[] initial,
                                    double delta, double[] buy, double[] sell, double kappa, int basket,
                                    int trades, /*int revise, int costs,*/ double min_holding,
                                    double min_trade,
                                    /*int m_LS, int Fully_Invested, */double Rmin, double Rmax,
                                    int round, double[] min_lot, double[] size_lot, int[] shake,
                                    /*int ncomp, double[] Composite,*/ double LSValue,
                                    /*int npiece, double[] hpiece, double[] pgrad,*/
                                    int nabs, double[] Abs_A, int mabs, int[] I_A, double[] Abs_U,
                                    double[] FC, double[] FL, double[] SV, double minRisk, double maxRisk,
                                    ref double ogamma, double[] mask,/* int log, string logfile,*/
                                    /*int downrisk, double downfactor,*/
                                    int longbasket, int shortbasket,
                                    int tradebuy, int tradesell,/* double zetaS, double zetaF,*/
                                    /*double ShortCostScale,*/ double LSValuel, double[] Abs_L, double[] breakdown)
    {
        var back = -1;
        Portfolio.Portfolio op;
        if (nfac == -1)
        {
            var op1 = new Portfolio.Portfolio("");
            op = op1;
        }
        else
        {
            var op1 = new Portfolio.FPortfolio("");
            op1.nfac = nfac;
            op = op1;
        }
        op.n = n;
        op.m = m;
        op.A = A;
        op.L = L;
        op.U = U;
        op.bench = benchmark;
        op.buy = buy;
        op.sell = sell;
        op.delta = delta;
        op.gamma = gamma;
        op.delta = delta;
        op.initial = initial;
        op.kappa = kappa;
        op.names = names;
        op.Q = Q;
        if (initial == null) initial = new double[n];
        back = op.BasicOptimisation(n, m, nfac, A, L, U, gamma, kappa, delta, LSValue, LSValuel, Rmin, Rmax,
        alpha, initial, buy, sell, names, false, nabs, Abs_A, Abs_L, Abs_U, mabs, I_A);
        BlasLike.dcopyvec(n, op.wback, w);
        for (var i = 0; i < n; ++i)
        {
            w[i] = Portfolio.Portfolio.check_digit(1e2 * w[i]) * 1e-2;
        }
        if (breakdown != null) op.RiskBreakdown(w, op.bench, breakdown);
        var info = new Portfolio.Portfolio.INFO();
        info.A = A;
        info.alpha = alpha;
        info.bench = benchmark;
        info.buy = buy;
        info.delta = delta;
        info.initial = initial;
        info.kappa = kappa;
        info.L =(double[]) L.Clone();
        info.m = m;
        info.n = n;
        info.names = names;
        info.nfac = nfac;
        info.sell = sell;
        info.U =(double[]) U.Clone();
        if (basket > 0 || trades > 0)
        {
            op.DropRisk(basket, trades, -1, info);
            BlasLike.dcopyvec(n, op.wback, w);back=op.BACK;
        }
        if (back<=1&&  maxRisk > 0 && minRisk > 0)
        {
            var riskhere = BlasLike.ddotvec(n, w, breakdown);
            if (benchmark != null) riskhere -= BlasLike.ddotvec(n, breakdown, benchmark);
            if (!(riskhere <= maxRisk && riskhere >= minRisk))
            {
                if (maxRisk == minRisk) info.target = maxRisk;
                else if (maxRisk < riskhere) info.target = maxRisk;
                else if (minRisk > riskhere) info.target = minRisk;
                // op.CalcRisk(gamma, info);
                if(basket<0&&trades<0){
                if (info.target == minRisk)
                    ogamma = ActiveSet.Optimise.Solve1D(op.CalcRisk, gamma, 1 - BlasLike.lm_eps8, 0, info);
                else
                    ogamma = ActiveSet.Optimise.Solve1D(op.CalcRisk, 0, gamma, 0, info);
                }else {op.DropRisk(basket,trades,info.target,info);ogamma=gamma=op.gamma;}
                BlasLike.dcopyvec(n, op.wback, w);
                if (breakdown != null) op.RiskBreakdown(w, op.bench, breakdown);
            }
            else ogamma = gamma;
        }
        else ogamma = gamma;
        if (back <= 1 && ogamma > 1) back = 16;
        if (back > 1) return back;
        if (min_holding > 0 || min_trade > 0 || round == 1)
        {
            var Op = new Portfolio.Portfolio.OptParamRound();
            Op.basket = basket;
            Op.trades = trades;
            Op.lower = L;
            Op.m = m;
            Op.n = n;
            Op.upper = U;
            Op.minholdlot = null;
            Op.mintradelot = null;
            var mintrade = min_trade < 0 ? null : new double[n];
            if (mintrade != null) BlasLike.dsetvec(n, min_trade, mintrade);
            var minhold = min_holding < 0 ? null : new double[n];
            if (minhold != null) BlasLike.dsetvec(n, min_holding, minhold);
            var roundw = (double[])op.wback.Clone();
            Op.x = op.wback; Op.MoreInfo = info;
            for (var i = 0; i < n; ++i) shake[i] = -1;
            if (round == 1)
            {
                op.Rounding(basket, trades, initial, min_lot, size_lot, roundw, null, null, info);
                for (var i = 0; i < n; ++i)
                {
                    roundw[i] = Portfolio.Portfolio.check_digit(1e2 * roundw[i]) * 1e-2;
                }
                op.roundcheck(n, roundw, initial, min_lot, size_lot, shake);
            }
            else
            {
                op.Thresh(Op, mintrade == null ? null : initial, mintrade == null ? minhold : mintrade, roundw, mintrade == null ? null : minhold);
                for (var i = 0; i < n; ++i)
                {
                    roundw[i] = Portfolio.Portfolio.check_digit(1e2 * roundw[i]) * 1e-2;
                }
                op.thresh_check(n, roundw, mintrade == null ? null : initial, L, U, mintrade == null ? minhold : mintrade, mintrade == null ? null : minhold, BlasLike.lm_eps8, shake);
            }
            BlasLike.dcopyvec(n, roundw, op.wback);
            BlasLike.dcopyvec(n, op.wback, w);
            foreach (var i in shake)
            {
                if (i != -1) ColourConsole.WriteEmbeddedColourLine($"[green]{op.names[i]}[/green][red] was not rounded properly! {op.w[i],26:e16}[/red]");
            }
            if (breakdown != null) op.RiskBreakdown(w, op.bench, breakdown);
        }
        return back;
    }
    [HttpGet("general")]
    public Optimise[] GetGen(double? delta, double? gamma, double? maxRisk, double? minRisk, double? min_holding, double? min_trade, int? basket, int? trades)
    {
        var op = new Optimise();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "alpha bench gamma initial delta buy sell kappa min_holding min_trade minRisk maxRisk rmin rmax min_lot size_lot value valuel mask A L U Q A_abs Abs_U Abs_L SV FC FL";
            CVarData.intFields = "n nfac m basket longbasket shortbasket tradebuy tradesell tradenum nabs mabs I_A round";
            CVarData.stringFields = "names";
            try
            {
                CVarData.Read("./generalopt");
            }
            catch
            {
                CVarData.Read("../generalopt");
            }
            op.n = CVarData.mapInt["n"][0];
            op.nfac = CVarData.mapInt["nfac"][0];
            op.m = CVarData.mapInt["m"][0];
            op.names = CVarData.mapString["names"];
            op.A = CVarData.mapDouble["A"];
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.alpha = CVarData.mapDouble["alpha"];
            op.bench = CVarData.mapDouble["bench"];
            op.basket = CVarData.mapInt["basket"][0];
            op.trades = CVarData.mapInt["tradenum"][0];
            try { op.Q = CVarData.mapDouble["Q"]; } catch { op.Q = null; }
            try { op.FL = CVarData.mapDouble["FL"]; } catch { op.FL = null; }
            try { op.FC = CVarData.mapDouble["FC"]; } catch { op.FC = null; }
            try { op.SV = CVarData.mapDouble["SV"]; } catch { op.SV = null; }
            op.buy = CVarData.mapDouble["buy"];
            op.sell = CVarData.mapDouble["sell"];
            op.mask = CVarData.mapDouble["mask"];
            op.initial = CVarData.mapDouble["initial"];
            op.delta = CVarData.mapDouble["delta"][0];
            op.gamma = CVarData.mapDouble["gamma"][0];
            op.kappa = CVarData.mapDouble["kappa"][0];
            if (delta != null) op.delta = delta;
            if (gamma != null) op.gamma = gamma;
            if (maxRisk != null) op.maxRisk = maxRisk.GetValueOrDefault();
            if (minRisk != null) op.minRisk = minRisk.GetValueOrDefault();
            if (basket != null) op.basket = basket.GetValueOrDefault();
            if (trades != null) op.trades = trades.GetValueOrDefault();
            op.Abs_A = CVarData.mapDouble["A_abs"];
            op.Abs_L = CVarData.mapDouble["Abs_L"];
            op.Abs_U = CVarData.mapDouble["Abs_U"];
            op.nabs = CVarData.mapInt["nabs"][0];
            op.mabs = CVarData.mapInt["mabs"][0];
            op.I_A = CVarData.mapInt["I_A"];
            op.round = CVarData.mapInt["round"][0];
            op.min_lot = CVarData.mapDouble["min_lot"];
            op.size_lot = CVarData.mapDouble["size_lot"];
            op.min_holding = CVarData.mapDouble["min_holding"][0];
            op.min_trade = CVarData.mapDouble["min_trade"][0];
            if (min_holding != null) op.min_holding = min_holding;
            if (min_trade != null) op.min_trade = min_trade;
        }
        double ogamma = op.ogamma.GetValueOrDefault();
        op.shake = new int[op.n.GetValueOrDefault()];
        var breakdown = new double[op.n.GetValueOrDefault()];
        op.w = new double[op.n.GetValueOrDefault()];
        op.back = Optimise_internalCVPAFbl(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
        op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta.GetValueOrDefault(),
        op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket.GetValueOrDefault(), op.trades.GetValueOrDefault(), op.min_holding.GetValueOrDefault(),
        op.min_trade.GetValueOrDefault(), op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
        op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
        op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
        op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
        op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown);

        op.ogamma = ogamma;
        op.message = Portfolio.Portfolio.OptMessages(op.back.GetValueOrDefault());
        op.mctr = breakdown;
        op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
        if (op.bench != null) op.risk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.mctr);
        op.risk = Portfolio.Portfolio.check_digit(1e2 * op.risk.GetValueOrDefault()) * 1e-2;


        return new[] { op };
    }
}

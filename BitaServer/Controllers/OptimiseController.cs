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
    public Optimise Post(Optimise op)
    {
        if (op.digit != null)
        {
            op.tdigit = op.digit != null ? Portfolio.Portfolio.check_digit((double)op.digit) : null;
        }
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return op;
    }
    [HttpGet()]
    public Optimise Blank()
    {
        return new Optimise();
    }
    [HttpGet("LOSS")]
    public Optimise GetLOSS()
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
        return op;
    }

    [HttpPost("LOSS")]
    public Optimise PostLOSS(Optimise op)
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
            op.Gstrength, op.DATA, 0.0, op.TargetReturn, op.LOSSopt.GetValueOrDefault(), op.LOSSmin.GetValueOrDefault(), op.LOSSmax.GetValueOrDefault());
        op.w = opt.wback;
        op.CVARGLprob = opt.CVARGLprob;
        op.result = new Optimise.checkv();
        op.message = Portfolio.Portfolio.OptMessages(Math.Abs(op.back.GetValueOrDefault()));
        op.result.breakdown = new double[op.n.GetValueOrDefault()];
        op.result.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, op.result.breakdown); op.result.mctr = new double[op.n.GetValueOrDefault()];
        var lcheck = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.breakdown);
        opt.RiskBreakdown(op.w, null, op.result.mctr);
        op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
        op.result.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        return op;
    }

    [HttpGet("ETL")]
    public Optimise GetETL()
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
        return op;
    }
    [HttpPost("ETL")]
    public Optimise PostETL(Optimise op)
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
            op.Gstrength, op.DATA, op.tail,
            null, op.ETLopt.GetValueOrDefault(), op.ETLmin.GetValueOrDefault(),
            op.ETLmax.GetValueOrDefault());
        op.w = opt.wback;
        op.CVARGLprob = opt.CVARGLprob;
        op.result = new Optimise.checkv();
        op.message = Portfolio.Portfolio.OptMessages(Math.Abs(op.back.GetValueOrDefault()));
        op.result.breakdown = new double[op.n.GetValueOrDefault()];
        double VAR = 0;
        int ind = -1;
        op.result.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref ind, op.result.breakdown);
        op.result.VAR = VAR;
        op.result.VARindex = ind;
        op.result.mctr = new double[op.n.GetValueOrDefault()];
        opt.RiskBreakdown(op.w, null, op.result.mctr);
        op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
        op.result.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        return op;
    }
    [HttpGet("test")]
    public Optimise Get(double? d1)
    {
        double testdigit;
        if (d1 != null) testdigit = d1.GetValueOrDefault();
        else testdigit = 123.9999999999;
        var op = new Optimise();
        op.digit = testdigit;
        op.tdigit = Portfolio.Portfolio.check_digit(testdigit);
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return op;
    }
    [HttpGet("test/n")]
    public Optimise Getn()
    {
        double testdigit = 11.000000000001;
        var op = new Optimise();
        op.n = 500;
        op.digit = testdigit;
        op.tdigit = Portfolio.Portfolio.check_digit(testdigit);
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return op;
    }
    [HttpGet]
    [Route("general")]
    public Optimise GetGen(bool? doOpt, int? round, double? min_lot, double? size_lot,
     double? Gstrength, double? LOSSmax, double? LOSSmin, double? ETLmax, double? ETLmin,
      double? targetR, string? datafile, double? delta, double? gamma, double? maxRisk,
      double? minRisk, double? min_holding, double? min_trade, int? basket, int? trades,
      string? logfile)
    {
        var op = new Optimise();
        if (doOpt != null) op.doOpt = doOpt.GetValueOrDefault();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "alpha bench gamma initial delta buy sell kappa min_holding min_trade minRisk lambda maxRisk rmin rmax Rmin Rmax min_lot size_lot LSvalue LSvaluel value valuel mask A L U Q A_abs Abs_A Abs_U Abs_L SV FC FL DATA tail R";
            CVarData.intFields = "n nfac m basket longbasket shortbasket trades tradebuy tradesell tradenum nabs mabs I_A round tlen";
            CVarData.stringFields = "names logfile";
            op.datafile = "generalopt";
            if (datafile != null) op.datafile = datafile;
            try
            {
                CVarData.Read($"./{op.datafile}");
            }
            catch
            {
                try
                {
                    CVarData.Read($"../{op.datafile}");
                }
                catch { op.message = $"Input file error \"{op.datafile}\""; return op; }
            }
            op.n = CVarData.mapInt["n"][0];
            op.nfac = CVarData.mapInt["nfac"][0];
            op.m = CVarData.mapInt["m"][0];
            op.names = CVarData.mapString["names"];
            op.A = CVarData.mapDouble["A"];
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.alpha = CVarData.mapDouble["alpha"];
            try { op.logfile = CVarData.mapString["logfile"][0]; } catch {; }
            try { op.value = CVarData.mapDouble["value"][0]; } catch {; }
            try { op.value = CVarData.mapDouble["LSvalue"][0]; } catch {; }
            try { op.Gstrength = CVarData.mapDouble["lambda"][0]; } catch {; }
            try { op.rmin = CVarData.mapDouble["rmin"][0]; } catch {; }
            try { op.rmax = CVarData.mapDouble["rmax"][0]; } catch {; }
            try { op.rmin = CVarData.mapDouble["Rmin"][0]; } catch {; }
            try { op.rmax = CVarData.mapDouble["Rmax"][0]; } catch {; }
            try { op.valuel = CVarData.mapDouble["valuel"][0]; } catch {; }
            try { op.valuel = CVarData.mapDouble["LSvaluel"][0]; } catch {; }
            try { op.bench = CVarData.mapDouble["bench"]; } catch { op.bench = null; }
            op.basket = CVarData.mapInt["basket"][0];
            try { op.trades = CVarData.mapInt["tradenum"][0]; } catch {; }
            try { op.trades = CVarData.mapInt["trades"][0]; } catch {; }
            try { op.Q = CVarData.mapDouble["Q"]; } catch { op.Q = null; }
            try { op.FL = CVarData.mapDouble["FL"]; } catch { op.FL = null; }
            try { op.FC = CVarData.mapDouble["FC"]; } catch { op.FC = null; }
            try { op.SV = CVarData.mapDouble["SV"]; } catch { op.SV = null; }
            try { op.buy = CVarData.mapDouble["buy"]; } catch { op.buy = null; }
            try { op.sell = CVarData.mapDouble["sell"]; } catch { op.sell = null; }
            try { op.mask = CVarData.mapDouble["mask"]; } catch { op.mask = null; }
            try { op.initial = CVarData.mapDouble["initial"]; } catch { op.initial = null; }
            try { op.delta = CVarData.mapDouble["delta"][0]; } catch { op.delta = null; }
            op.gamma = CVarData.mapDouble["gamma"][0];
            op.kappa = CVarData.mapDouble["kappa"][0];
            if (delta != null) op.delta = delta;
            if (gamma != null) op.gamma = gamma;
            op.maxRisk = CVarData.mapDouble["maxRisk"][0];
            op.minRisk = CVarData.mapDouble["minRisk"][0];
            if (maxRisk != null) op.maxRisk = maxRisk.GetValueOrDefault();
            if (minRisk != null) op.minRisk = minRisk.GetValueOrDefault();
            if (basket != null) op.basket = basket.GetValueOrDefault();
            if (trades != null) op.trades = trades.GetValueOrDefault();
            try { op.Abs_A = CVarData.mapDouble["Abs_A"]; } catch {; }
            try { op.Abs_A = CVarData.mapDouble["A_abs"]; } catch {; }
            try { op.Abs_L = CVarData.mapDouble["Abs_L"]; } catch {; }
            try { op.Abs_U = CVarData.mapDouble["Abs_U"]; } catch {; }
            try { op.nabs = CVarData.mapInt["nabs"][0]; } catch { op.nabs = null; }
            try { op.mabs = CVarData.mapInt["mabs"][0]; } catch { op.mabs = null; }
            try { op.I_A = CVarData.mapInt["I_A"]; } catch { op.I_A = null; }
            try { op.round = CVarData.mapInt["round"][0]; } catch { op.round = null; }
            try { op.min_lot = CVarData.mapDouble["min_lot"]; } catch { op.min_lot = null; }
            if (min_lot != null && op.min_lot != null) op.min_lot[0] = min_lot.GetValueOrDefault();
            if (op.min_lot != null && op.min_lot.Length == 1)
            {
                var keep = op.min_lot[0];
                op.min_lot = new double[op.n.GetValueOrDefault()];
                BlasLike.dsetvec(op.n.GetValueOrDefault(), keep, op.min_lot);
            }
            try { op.size_lot = CVarData.mapDouble["size_lot"]; } catch { op.size_lot = null; }
            if (size_lot != null && op.size_lot != null) op.size_lot[0] = size_lot.GetValueOrDefault();
            if (op.size_lot != null && op.size_lot.Length == 1)
            {
                var keep = op.size_lot[0];
                op.size_lot = new double[op.n.GetValueOrDefault()];
                BlasLike.dsetvec(op.n.GetValueOrDefault(), keep, op.size_lot);
            }
            try { op.min_holding = CVarData.mapDouble["min_holding"][0]; } catch { op.min_holding = -1; }
            try { op.min_trade = CVarData.mapDouble["min_trade"][0]; } catch { op.min_trade = -1; }
            if (min_holding != null) op.min_holding = min_holding.GetValueOrDefault();
            if (min_trade != null) op.min_trade = min_trade.GetValueOrDefault();
            if (round != null) op.round = round.GetValueOrDefault();
            try { op.tlen = CVarData.mapInt["tlen"][0]; } catch { op.tlen = 0; }
            if (logfile != null) op.logfile = logfile;
            if (op.tlen > 0)
            {if(op.nfac==null)op.nfac=-1;
                op.DATA = CVarData.mapDouble["DATA"];
                try { op.tail = CVarData.mapDouble["tail]"][0]; } catch { op.tail = 0.05; }
                try { targetR = CVarData.mapDouble["R"][0]; } catch { targetR = null; }
                if (Gstrength != null) op.Gstrength = Gstrength.GetValueOrDefault();
                if (LOSSmin != null) op.LOSSmin = LOSSmin;
                if (LOSSmax != null) op.LOSSmax = LOSSmax;
                if (ETLmax != null) op.ETLmax = ETLmax;
                if (ETLmin != null) op.ETLmin = ETLmin;
                if (op.ETLmax != null && op.ETLmin != null) op.ETLopt = true;
                if (op.LOSSmax != null && op.LOSSmin != null) op.LOSSopt = true;
                if (op.alpha == null)
                {
                    var ones = (double[])new double[op.tlen];
                    op.alpha = new double[op.n.GetValueOrDefault()];
                    BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
                    Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, 0, 0, 0, true);
                }
                if (op.Q == null && op.nfac.GetValueOrDefault() < 0)
                {
                    var Q = new double[op.n.GetValueOrDefault() * (op.n.GetHashCode() + 1) / 2]; op.Q = Q;
                    var ij = 0;
                    for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
                    {
                        for (var j = 0; j <= i; ++j)
                        {
                            Q[ij++] = Solver.Factorise.covariance(op.tlen, op.DATA, op.DATA, i * op.tlen, j * op.tlen);
                        }
                    }
                    op.Q = Q;
                }

                if (targetR != null)
                {
                    if (op.LOSSmin == null) op.LOSSmin = -100;
                    if (op.LOSSmax == null) op.LOSSmax = 100;
                    op.TargetReturn = new double[op.tlen];
                    BlasLike.dsetvec(op.tlen, targetR.GetValueOrDefault(), op.TargetReturn);
                }
                else
                {
                    op.TargetReturn = null;  //probably redundant

                    if (op.ETLmin == null) op.ETLmin = -100;
                    if (op.ETLmax == null) op.ETLmax = 100;
                }
            }
        }

        double ogamma = op.ogamma.GetValueOrDefault();
        op.shake = new int[op.n.GetValueOrDefault()];
        var breakdown = new double[op.n.GetValueOrDefault()];
        op.result = new Optimise.checkv();
        if (op.doOpt)
        {
            op.w = new double[op.n.GetValueOrDefault()];
            bool CVARGLprob = false;
            op.back = Portfolio.Portfolio.OptimiseGeneral(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
            op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta.GetValueOrDefault(),
            op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket, op.trades, op.min_holding,
            op.min_trade, op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
            op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
            op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
            op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
            op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown, ref CVARGLprob, op.tlen, op.Gstrength,
            op.DATA, op.tail, op.TargetReturn, op.ETLopt.GetValueOrDefault() || op.LOSSopt.GetValueOrDefault(), op.TargetReturn == null ? op.ETLmin.GetValueOrDefault() : op.LOSSmin.GetValueOrDefault(),
            op.TargetReturn == null ? op.ETLmax.GetValueOrDefault() : op.LOSSmax.GetValueOrDefault(), op.logfile);

            op.ogamma = ogamma;
            op.CVARGLprob = CVARGLprob;
            op.message = Portfolio.Portfolio.OptMessages(op.back.GetValueOrDefault());
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                if (op.shake[i] == i) if (op.names == null) { op.message += $"\nAsset {i + 1} was not rounded properly"; }
                    else { op.message += $"\n{op.names[i]} was not rounded properly"; }
            }
            op.result.mctr = breakdown;
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            op.result.risk = Portfolio.Portfolio.check_digit(1e2 * op.result.risk.GetValueOrDefault()) * 1e-2;
        }
        else
        {
            if (op.w == null)
            {
                op.result = null;
                op.shake = null;
                op.ogamma = op.gamma;
                return op;
            }
        }
        if (op.nfac < 0)
        {
            var cov = new Portfolio.Portfolio("");
            cov.ntrue = op.n.GetValueOrDefault();
            cov.Q = op.Q;

            cov.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
        }
        else
        {
            var fac = new Portfolio.FPortfolio("");
            fac.ntrue = op.n.GetValueOrDefault();
            fac.nfac = op.nfac.GetValueOrDefault();
            if (op.SV != null && op.nfac > -1)
            {
                fac.SV = op.SV;
                fac.FC = op.FC;
                fac.FL = op.FL;
                fac.makeQ();
            }
            else fac.Q = op.Q;
            fac.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            if (op.SV != null && op.nfac > -1)
            {
                op.result.FX = new double[op.nfac.GetValueOrDefault()];
                op.result.Fmctr = new double[op.nfac.GetValueOrDefault()];
                op.result.SPmctr = new double[op.n.GetValueOrDefault()];
                fac.FactorRiskAttribution(op.w, op.bench, op.result.FX, op.result.Fmctr, op.result.SPmctr);
                op.result.facrisk = BlasLike.ddotvec(op.nfac.GetValueOrDefault(), op.result.FX, op.result.Fmctr);
                op.result.specrisk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.SPmctr);
                if (op.bench != null) op.result.specrisk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.SPmctr);
            }
        }
        if (true)
        {
            op.result.basket = 0;
            op.result.trades = 0;
            op.result.minhold = BlasLike.lm_max;
            op.result.gross = 0;
            op.result.longvalue = 0; op.result.shortvalue = 0;
            var iii=0;
            foreach (var ww in op.w)
            {
                if (ww > 0)
                {
                    op.result.longvalue += ww;
                    op.result.gross += ww;
                }
                else
                {
                    op.result.gross -= ww;
                    op.result.shortvalue -= ww;
                }
                if (op.result.longvalue > 0) op.result.shortoverlong = op.result.shortvalue / op.result.longvalue;
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.result.basket++;
                    if(op.L!=null&&op.U!=null&&op.L[iii]!=op.U[iii])                 op.result.minhold = Math.Min(Math.Abs(ww), op.result.minhold.GetValueOrDefault());
                }
                iii++;
            }
        }
        if (op.initial != null)
        {
            BlasLike.dsubvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
            op.result.trades = 0;
            op.result.mintrade = BlasLike.lm_max;
            op.result.turnover = 0;
            if (op.buy != null && op.sell != null)
            {
                op.result.cost = 0;
                for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
                {
                    op.result.cost += op.w[i] > 0 ? op.w[i] * op.buy[i] : -op.w[i] * op.sell[i];
                }
            }
            var iii=0;
            foreach (var ww in op.w)
            {
                op.result.turnover += Math.Abs(ww);
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.result.trades++;
                   if(op.L!=null&&op.U!=null&&op.L[iii]!=op.U[iii]) op.result.mintrade = Math.Min(Math.Abs(ww), op.result.mintrade.GetValueOrDefault());
                }
                iii++;
            }
            op.result.turnover *= 0.5;
            BlasLike.daddvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
        }

        op.result.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        if (op.tlen > 0)
        {
            var breakd = (double[])new double[op.n.GetValueOrDefault()];
            if (op.TargetReturn == null)
            {
                var VAR = 0.0;
                var VARindex = 0;
                op.result.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref VARindex, breakd);
                op.result.VAR = VAR;
                op.result.VARindex = VARindex;
                op.result.breakdown = breakd;
            }
            else
            {
                op.result.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, breakd);
                op.result.breakdown = breakd;
            }
        }

        return op;
    }
    [HttpPost("general")]
    public Optimise PostGen(Optimise op)
    {
        var breakdown = new double[op.n.GetValueOrDefault()];
        op.result = new Optimise.checkv();
        if(op.tlen>0){if(op.nfac==null)op.nfac=-1;
            if (op.alpha == null)
        {
            var ones = (double[])new double[op.tlen];
            op.alpha = new double[op.n.GetValueOrDefault()];
            BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
            Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, 0, 0, 0, true);
        }
        if (op.Q == null && op.nfac.GetValueOrDefault() < 0)
        {
            var Q = new double[op.n.GetValueOrDefault() * (op.n.GetHashCode() + 1) / 2]; op.Q = Q;
            var ij = 0;
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                for (var j = 0; j <= i; ++j)
                {
                    Q[ij++] = Solver.Factorise.covariance(op.tlen, op.DATA, op.DATA, i * op.tlen, j * op.tlen);
                }
            }
            op.Q = Q;
        }}

        if (op.doOpt)
        {
            op.w = new double[op.n.GetValueOrDefault()];
            op.shake = new int[op.n.GetValueOrDefault()];
            double ogamma = 0;
            bool CVARGLprob = false;
            op.back = Portfolio.Portfolio.OptimiseGeneral(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
            op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta.GetValueOrDefault(),
            op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket, op.trades, op.min_holding,
            op.min_trade, op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
            op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
            op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
            op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
            op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown, ref CVARGLprob, op.tlen, op.Gstrength,
            op.DATA, op.tail, op.TargetReturn, op.ETLopt.GetValueOrDefault() || op.LOSSopt.GetValueOrDefault(), op.TargetReturn == null ? op.ETLmin.GetValueOrDefault() : op.LOSSmin.GetValueOrDefault(),
            op.TargetReturn == null ? op.ETLmax.GetValueOrDefault() : op.LOSSmax.GetValueOrDefault(), op.logfile);


            op.CVARGLprob = CVARGLprob;
            op.ogamma = ogamma;
            op.message = Portfolio.Portfolio.OptMessages(op.back.GetValueOrDefault());
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                if (op.shake[i] == i) if (op.names == null) { op.message += $"\nAsset {i + 1} was not rounded properly"; }
                    else { op.message += $"\n{op.names[i]} was not rounded properly"; }
            }
        }
        op.result.mctr = new double[op.n.GetValueOrDefault()];
        if (op.nfac < 0)
        {
            var cov = new Portfolio.Portfolio("");
            cov.ntrue = op.n.GetValueOrDefault();
            cov.Q = op.Q;

            cov.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
        }
        else
        {
            var fac = new Portfolio.FPortfolio("");
            fac.ntrue = op.n.GetValueOrDefault();
            fac.nfac = op.nfac.GetValueOrDefault();
            if (op.SV != null && op.nfac > -1)
            {
                fac.SV = op.SV;
                fac.FC = op.FC;
                fac.FL = op.FL;
                fac.makeQ();
            }
            else fac.Q = op.Q;
            fac.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            if (op.SV != null && op.nfac > -1)
            {
                op.result.FX = new double[op.nfac.GetValueOrDefault()];
                op.result.Fmctr = new double[op.nfac.GetValueOrDefault()];
                op.result.SPmctr = new double[op.n.GetValueOrDefault()];
                fac.FactorRiskAttribution(op.w, op.bench, op.result.FX, op.result.Fmctr, op.result.SPmctr);
                op.result.facrisk = BlasLike.ddotvec(op.nfac.GetValueOrDefault(), op.result.FX, op.result.Fmctr);
                op.result.specrisk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.SPmctr);
                if (op.bench != null) op.result.specrisk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.SPmctr);
            }
        }
        if (true)
        {
            op.result.basket = 0;
            op.result.trades = 0;
            op.result.minhold = BlasLike.lm_max;
            op.result.gross = 0;
            op.result.longvalue = 0; op.result.shortvalue = 0;
            var iii=0;
            foreach (var ww in op.w)
            {
                if (ww > 0)
                {
                    op.result.longvalue += ww;
                    op.result.gross += ww;
                }
                else
                {
                    op.result.gross -= ww;
                    op.result.shortvalue -= ww;
                }
                if (op.result.longvalue > 0) op.result.shortoverlong = op.result.shortvalue / op.result.longvalue;
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.result.basket++;
                    if(op.L!=null&&op.U!=null&&op.L[iii]!=op.U[iii])op.result.minhold = Math.Min(Math.Abs(ww), op.result.minhold.GetValueOrDefault());
                }
                iii++;
            }
        }
        if (op.initial != null)
        {
            BlasLike.dsubvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
            op.result.trades = 0;
            op.result.mintrade = BlasLike.lm_max;
            op.result.turnover = 0;
            if (op.buy != null && op.sell != null)
            {
                op.result.cost = 0;
                for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
                {
                    op.result.cost += op.w[i] > 0 ? op.w[i] * op.buy[i] : -op.w[i] * op.sell[i];
                }
            }
            var iii=0;
            foreach (var ww in op.w)
            {
                op.result.turnover += Math.Abs(ww);
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.result.trades++;
                    if(op.L!=null&&op.U!=null&&op.L[iii]!=op.U[iii])op.result.mintrade = Math.Min(Math.Abs(ww), op.result.mintrade.GetValueOrDefault());
                }
                iii++;
            }
            op.result.turnover *= 0.5;
            BlasLike.daddvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
        }

        op.result.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        if (op.tlen > 0)
        {
            var breakd = (double[])new double[op.n.GetValueOrDefault()];
            if (op.TargetReturn == null)
            {
                var VAR = 0.0;
                var VARindex = 0;
                op.result.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref VARindex, breakd);
                op.result.VAR = VAR;
                op.result.VARindex = VARindex;
                op.result.breakdown = breakd;
            }
            else
            {
                op.result.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, breakd);
                op.result.breakdown = breakd;
            }
        }
        return op;
    }
}

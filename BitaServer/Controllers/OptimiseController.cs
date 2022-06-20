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
    [HttpGet("general")]
    public Optimise[] GetGen(bool?doOpt,  double? min_lot, double? size_lot, double? Gstrength, double? LOSSmax, double? LOSSmin, double? ETLmax, double? ETLmin, double? targetR, string? datafile, double? delta, double? gamma, double? maxRisk, double? minRisk, double? min_holding, double? min_trade, int? basket, int? trades)
    {
        var op = new Optimise();
        if(doOpt!=null)op.doOpt=doOpt.GetValueOrDefault();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "alpha bench gamma initial delta buy sell kappa min_holding min_trade minRisk maxRisk rmin rmax min_lot size_lot value valuel mask A L U Q A_abs Abs_U Abs_L SV FC FL DATA tail R";
            CVarData.intFields = "n nfac m basket longbasket shortbasket tradebuy tradesell tradenum nabs mabs I_A round tlen";
            CVarData.stringFields = "names";
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
                catch { op.message = $"No input file \"{op.datafile}\""; return new[] { op }; }
            }
            op.n = CVarData.mapInt["n"][0];
            op.nfac = CVarData.mapInt["nfac"][0];
            op.m = CVarData.mapInt["m"][0];
            op.names = CVarData.mapString["names"];
            op.A = CVarData.mapDouble["A"];
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.alpha = CVarData.mapDouble["alpha"];
            op.value=CVarData.mapDouble["value"][0];
            op.rmin=CVarData.mapDouble["rmin"][0];
            op.rmax=CVarData.mapDouble["rmax"][0];
            op.valuel=CVarData.mapDouble["valuel"][0];
            try { op.bench = CVarData.mapDouble["bench"]; } catch { op.bench = null; }
            op.basket = CVarData.mapInt["basket"][0];
            try { op.trades = CVarData.mapInt["tradenum"][0]; } catch { op.trades = -1; }
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
            op.maxRisk=CVarData.mapDouble["maxRisk"][0];
            op.minRisk=CVarData.mapDouble["minRisk"][0];
            if (maxRisk != null) op.maxRisk = maxRisk.GetValueOrDefault();
            if (minRisk != null) op.minRisk = minRisk.GetValueOrDefault();
            if (basket != null) op.basket = basket.GetValueOrDefault();
            if (trades != null) op.trades = trades.GetValueOrDefault();
            try { op.Abs_A = CVarData.mapDouble["A_abs"]; } catch { op.Abs_A = null; }
            try { op.Abs_L = CVarData.mapDouble["Abs_L"]; } catch { op.Abs_L = null; }
            try { op.Abs_U = CVarData.mapDouble["Abs_U"]; } catch { op.Abs_U = null; }
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
            try { op.tlen = CVarData.mapInt["tlen"][0]; } catch { op.tlen = 0; }
            if (targetR == null) targetR = 0.0;
            if (op.tlen > 0)
            {
                op.DATA = CVarData.mapDouble["DATA"];
                try { op.tail = CVarData.mapDouble["tail]"][0]; } catch { op.tail = 0.05; }
                targetR = CVarData.mapDouble["R"][0];
                if (Gstrength != null) op.Gstrength = Gstrength;
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
            }
            op.TargetReturn = new double[op.tlen];
            BlasLike.dsetvec(op.tlen, targetR.GetValueOrDefault(), op.TargetReturn);
        }

        double ogamma = op.ogamma.GetValueOrDefault();
        op.shake = new int[op.n.GetValueOrDefault()];
        var breakdown = new double[op.n.GetValueOrDefault()];
        if (op.doOpt)
        {
        op.w = new double[op.n.GetValueOrDefault()];
            op.back = Portfolio.Portfolio.OptimiseGeneral(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
            op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta.GetValueOrDefault(),
            op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket, op.trades, op.min_holding,
            op.min_trade, op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
            op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
            op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
            op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
            op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown, op.tlen, op.Gstrength.GetValueOrDefault(),
            op.DATA, op.tail, op.TargetReturn, op.ETLopt.GetValueOrDefault() || op.LOSSopt.GetValueOrDefault(), op.TargetReturn == null ? op.ETLmin.GetValueOrDefault() : op.LOSSmin.GetValueOrDefault(),
            op.TargetReturn == null ? op.ETLmax.GetValueOrDefault() : op.LOSSmax.GetValueOrDefault());

            op.ogamma = ogamma;
            op.message = Portfolio.Portfolio.OptMessages(op.back.GetValueOrDefault());
            op.mctr = breakdown;
            op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
            if (op.bench != null) op.risk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.mctr);
            op.risk = Portfolio.Portfolio.check_digit(1e2 * op.risk.GetValueOrDefault()) * 1e-2;
        }
        else{
            if(op.w==null)return new[] { op };
            if(op.mctr==null)op.mctr=new double[op.n.GetValueOrDefault()];
        }
        if (op.nfac < -1)
        {
            var cov = new Portfolio.Portfolio("");
            cov.ntrue = op.n.GetValueOrDefault();
            cov.Q = op.Q;

            cov.RiskBreakdown(op.w, op.bench, op.mctr);
            op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
            if (op.bench != null) op.risk -= op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.mctr);
        }
        else
        {
            var fac = new Portfolio.FPortfolio("");
            fac.ntrue = op.n.GetValueOrDefault();
            fac.nfac = op.nfac.GetValueOrDefault();
            if (op.SV != null)
            {
                fac.SV = op.SV;
                fac.FC = op.FC;
                fac.FL = op.FL;
                fac.makeQ();
            }
            else fac.Q = op.Q;
            fac.RiskBreakdown(op.w, op.bench, op.mctr);
            op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.mctr);
            if (op.bench != null) op.risk -= op.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.mctr);
            if (op.SV != null)
            {op.FX=new double[op.nfac.GetValueOrDefault()];
            op.Fmctr=new double[op.nfac.GetValueOrDefault()];
            op.SPmctr=new double[op.n.GetValueOrDefault()];
                fac.FactorRiskAttribution(op.w, op.bench, op.FX, op.Fmctr, op.SPmctr);
                op.facrisk = BlasLike.ddotvec(op.nfac.GetValueOrDefault(), op.FX, op.Fmctr);
                op.specrisk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.SPmctr);
                if(op.bench!=null)op.specrisk-=BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.SPmctr);
            }
        }
        if (true)
        {
            op.achievedbasket = 0;
            op.achievedminhold = BlasLike.lm_max;
            foreach (var ww in op.w)
            {
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.achievedbasket++;
                    op.achievedminhold = Math.Min(Math.Abs(ww), op.achievedminhold);
                }
            }
        }
        if (op.initial != null)
        {
            BlasLike.dsubvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
            op.achievedtrades = 0;
            op.achievedmintrade = BlasLike.lm_max;
            foreach (var ww in op.w)
            {
                if (Math.Abs(ww) > BlasLike.lm_eps8)
                {
                    op.achievedtrades++;
                    op.achievedmintrade = Math.Min(Math.Abs(ww), op.achievedmintrade);
                }
            }
            BlasLike.daddvec(op.n.GetValueOrDefault(), op.w, op.initial, op.w);
        }

        op.expreturn = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.alpha);
        if (op.tlen > 0)
        {
            var breakd = (double[])new double[op.n.GetValueOrDefault()];
            if (op.TargetReturn == null)
            {
                var VAR = 0.0;
                var VARindex = 0;
                op.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref VARindex, breakd);
                op.VAR = VAR;
                op.VARindex = VARindex;
                op.breakdown = breakd;
            }
            else
            {
                op.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, breakd);
                op.breakdown = breakd;
            }
        }

        return new[] { op };
    }
}

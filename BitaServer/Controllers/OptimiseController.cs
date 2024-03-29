using Microsoft.AspNetCore.Mvc;
using Blas;
using Solver;
using DataFile;
using Microsoft.Extensions.Hosting.WindowsServices;
namespace BitaServer.Controllers;

[ApiController]
[Route("[controller]")]
public class OptimiseController : ControllerBase
{
    private readonly ILogger<OptimiseController> _logger;
    private static double rounder(double k) => k;//Not sure whether to round the risks given back. Set code here if wanted.
    public OptimiseController(ILogger<OptimiseController> logger)
    {
        _logger = logger;
    }
    [HttpPost("test")]
    public Test Post(Test op)
    {
        if (op.digit != null)
        {
            op.tdigit = op.digit != null ? Portfolio.Portfolio.check_digit((double)op.digit) : null;
        }
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green]");
        _logger.LogInformation($"POST test at {DateTimeOffset.Now}");
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
        op.delta = -1;
        op.TargetReturn = new double[op.tlen];
        BlasLike.dsetvec(op.tlen, targetR, op.TargetReturn);
        op.LOSSopt = false;
        op.Gstrength = 1;
        //     op.tail=0.05;
        op.LOSSmin = -1;
        op.LOSSmax = 10;
        op.gamma = 0;
        _logger.LogInformation($"GET LOSS at {DateTimeOffset.Now}");
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
            Q = new double[op.n.GetValueOrDefault() * (op.n.GetValueOrDefault() + 1) / 2]; op.Q = Q;
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
        Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, atran: true);
        op.back = opt.BasicOptimisation(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), -1,
           op.A, op.L, op.U, op.gamma.GetValueOrDefault(), op.kappa.GetValueOrDefault(), -1, -1, -1, -1, -1, op.alpha, op.initial, null, null,
           op.names, false, 0, null, null, null, 0, null, op.tlen,
            op.Gstrength.GetValueOrDefault(), op.DATA, 0.0, op.TargetReturn, op.LOSSopt.GetValueOrDefault(), op.LOSSmin.GetValueOrDefault(), op.LOSSmax.GetValueOrDefault());
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
        _logger.LogInformation($"POST LOSS at {DateTimeOffset.Now}");
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
        op.delta = -1;
        BlasLike.dnegvec(op.DATA.Length, op.DATA);
        op.ETLopt = false;
        op.Gstrength = 1;
        //     op.tail=0.05;
        op.ETLmin = -1;
        op.ETLmax = 1;
        op.gamma = 0;
        _logger.LogInformation($"GET ETL at {DateTimeOffset.Now}");
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
            Q = new double[op.n.GetValueOrDefault() * (op.n.GetValueOrDefault() + 1) / 2]; op.Q = Q;
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
        Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, atran: true);
        BlasLike.dnegvec(op.n.GetValueOrDefault(), op.alpha);
        op.back = opt.BasicOptimisation(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), -1,
           op.A, op.L, op.U, op.gamma.GetValueOrDefault(), op.kappa.GetValueOrDefault(), -1, -1, -1, -1, -1, op.alpha, op.initial, null, null,
           op.names, false, 0, null, null, null, 0, null, op.tlen,
            op.Gstrength.GetValueOrDefault(), op.DATA, op.tail,
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
        _logger.LogInformation($"POST ETL at {DateTimeOffset.Now}");
        return op;
    }
    [HttpGet("test")]
    public Test Get(double? d1)
    {
        double testdigit;
        if (d1 != null) testdigit = d1.GetValueOrDefault();
        else testdigit = 123.9999999999;
        var op = new Test();
        op.digit = testdigit;
        op.tdigit = Portfolio.Portfolio.check_digit(testdigit);
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green]");
        _logger.LogInformation($"GET test at {DateTimeOffset.Now}");
        return op;
    }
    [HttpGet("test/n")]
    public Optimise Getn()
    {
        var op = new Optimise();
        op.n = 500;
        ColourConsole.WriteEmbeddedColourLine($"[yellow]{op.n}[/yellow]");
        return op;
    }
    [HttpGet]
    [Route("general")]
    public Object GetGen(bool? doOpt, int? round, double? min_lot, double? size_lot,double?value,double?valuel,
     double? Gstrength, double? LOSSmax, double? LOSSmin, double? ETLmax, double? ETLmin,double?tail,
      double? targetR, string? datafile, double ?delta, double? gamma, double? maxRisk,
      double? minRisk, double? min_holding, double? min_trade, int? basket, int? trades,
      string? logfile,int? ncomp, bool negdata = false)
    {
        var op = new Optimise();
        var lic = new Licensing.Licence();
        var ok = lic.CheckLicence();
        op.VersionString = lic.VersionString;
        op.isLicensed = ok;
        if (!ok) return Problem(title: "Bad licence", detail: lic.VersionString,type:lic.connectedNames);
        if (!ok) return op;
        if (doOpt != null) op.doOpt = doOpt.GetValueOrDefault();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "alpha bench gamma initial delta buy sell kappa min_holding min_trade minRisk lambda maxRisk rmin rmax Rmin Rmax min_lot size_lot LSValue LSValuel value valuel mask A L U Q A_abs Abs_A Abs_U Abs_L SV FC FL DATA tail R CVarMin CVarMax CVar_averse ETLorLOSSmin ETLorLOSSmax Composites";
            CVarData.intFields = "n nfac m basket longbasket shortbasket trades tradebuy tradesell tradenum nabs mabs I_A round tlen number_included CVar_constraint ETLorLOSSconstraint ncomp";
            CVarData.stringFields = "names logfile";
            op.datafile = "generalopt";
            if (datafile != null) op.datafile = datafile;
            var ContentRootPath = AppContext.BaseDirectory;
#if DEBUG
            ContentRootPath = "./";
#endif
            try
            {
                CVarData.Read($"{ContentRootPath}{op.datafile}");
            }
            catch
            {
                try
                {
                    CVarData.Read($"../{op.datafile}");
                }
                catch (Exception ppp) { op.message = $"{ppp.Message} Input file error \"{op.datafile}\""; return Problem(title: "No DATA file", detail: lic.VersionString,type:op.message); }
            }
            op.n = CVarData.mapInt["n"][0];
            try { op.nfac = CVarData.mapInt["nfac"][0]; } catch { op.nfac = -1; }
            op.m = CVarData.mapInt["m"][0];
            op.names = CVarData.mapString["names"];
            op.A = CVarData.mapDouble["A"];
            if(op.Aas2D==null)op.Aas2D = Portfolio.Portfolio.oneD2twoD(op.m.GetValueOrDefault(), op.n.GetValueOrDefault(), op.A, transpose: false);
            op.L = CVarData.mapDouble["L"];
            op.U = CVarData.mapDouble["U"];
            op.alpha = CVarData.mapDouble["alpha"];
            try { op.ncomp = CVarData.mapInt["ncomp"][0]; } catch {op.ncomp=null; }
            try { op.composites = CVarData.mapDouble["Composites"]; } catch {op.composites=null; }
            try { op.logfile = CVarData.mapString["logfile"][0]; } catch {; }
            try { op.value = CVarData.mapDouble["value"][0]; } catch {; }
            try { op.value = CVarData.mapDouble["LSValue"][0]; } catch {; }
            try { op.Gstrength = CVarData.mapDouble["lambda"][0]; } catch {; }
            try { op.rmin = CVarData.mapDouble["rmin"][0]; } catch {; }
            try { op.rmax = CVarData.mapDouble["rmax"][0]; } catch {; }
            try { op.rmin = CVarData.mapDouble["Rmin"][0]; } catch {; }
            try { op.rmax = CVarData.mapDouble["Rmax"][0]; } catch {; }
            try { op.valuel = CVarData.mapDouble["valuel"][0]; } catch {; }
            try { op.valuel = CVarData.mapDouble["LSValuel"][0]; } catch {; }
            try { op.bench = CVarData.mapDouble["bench"]; } catch { op.bench = null; }
            op.basket = CVarData.mapInt["basket"][0];
            try { op.trades = CVarData.mapInt["tradenum"][0]; } catch {; }
            try { op.trades = CVarData.mapInt["trades"][0]; } catch {; }
            try { op.Q = CVarData.mapDouble["Q"]; } catch { op.Q = null; }
            try { op.FL = CVarData.mapDouble["FL"]; } catch { op.FL = null; }
            if(op.FL!=null)op.FLas2D=Portfolio.Portfolio.oneD2twoD(op.n.GetValueOrDefault(),op.nfac.GetValueOrDefault(),op.FL);
            try { op.FC = CVarData.mapDouble["FC"]; } catch { op.FC = null; }
            try { op.SV = CVarData.mapDouble["SV"]; } catch { op.SV = null; }
            try { op.buy = CVarData.mapDouble["buy"]; } catch { op.buy = null; }
            try { op.sell = CVarData.mapDouble["sell"]; } catch { op.sell = null; }
            try { op.mask = CVarData.mapDouble["mask"]; } catch { op.mask = null; }
            try { op.initial = CVarData.mapDouble["initial"]; } catch { op.initial = null; }
            try { op.delta = CVarData.mapDouble["delta"][0]; } catch { op.delta = -1; }
            op.gamma = CVarData.mapDouble["gamma"][0];
            try { op.kappa = CVarData.mapDouble["kappa"][0]; } catch { op.kappa = -1; }
            if (delta != null) op.delta = delta.GetValueOrDefault();
            if (gamma != null) op.gamma = gamma;
            try { op.maxRisk = CVarData.mapDouble["maxRisk"][0]; } catch {; }
            try { op.minRisk = CVarData.mapDouble["minRisk"][0]; } catch {; }
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
            if(round.GetValueOrDefault()==1){op.round=round;}

            if (min_lot != null) { op.min_lot = new double[1]; op.min_lot[0] = min_lot.GetValueOrDefault(); }
            if (op.min_lot != null && op.min_lot.Length == 1)            op.min_lot=Portfolio.Portfolio.one2many(op.n.GetValueOrDefault(),op.min_lot[0]);
            
            try { op.size_lot = CVarData.mapDouble["size_lot"]; } catch { op.size_lot = null; }
            if (size_lot != null) { op.size_lot = new double[1]; op.size_lot[0] = size_lot.GetValueOrDefault(); }
            if (op.size_lot != null && op.size_lot.Length == 1)op.size_lot=Portfolio.Portfolio.one2many(op.n.GetValueOrDefault(),op.size_lot[0]);

            try { op.min_holding = CVarData.mapDouble["min_holding"][0]; } catch { op.min_holding = -1; }
            try { op.min_trade = CVarData.mapDouble["min_trade"][0]; } catch { op.min_trade = -1; }
            if (min_holding != null) op.min_holding = min_holding.GetValueOrDefault();
            if (min_trade != null) op.min_trade = min_trade.GetValueOrDefault();
            if (round != null) op.round = round.GetValueOrDefault();
            try { op.tlen = CVarData.mapInt["tlen"][0]; } catch { op.tlen = 0; }
            if (logfile != null) op.logfile = logfile;
            if(ncomp>0)op.ncomp=ncomp;
            if(value!=null)op.value=value.GetValueOrDefault();
            if(valuel!=null)op.valuel=valuel.GetValueOrDefault();
            if (op.tlen > 0)
            {
                if (op.nfac == null) op.nfac = -1;
                op.DATA = CVarData.mapDouble["DATA"];
                if (negdata) BlasLike.dnegvec(op.DATA.Length, op.DATA);
                try { Gstrength = CVarData.mapDouble["CVar_averse"][0]; } catch {; }
                try { op.tail = CVarData.mapDouble["tail"][0]; }
                catch
                {
                    try
                    {
                        var number_included = CVarData.mapInt["number_included"][0];
                        op.tail = 1 - (double)((number_included)) / op.tlen;
                    }
                    catch { op.tail = 0.02; }
                }
            if(tail!=null)op.tail=tail.GetValueOrDefault();
                if (targetR == null) try { targetR = CVarData.mapDouble["R"][0]; } catch { targetR = null; }
                if (Gstrength != null) op.Gstrength = Gstrength.GetValueOrDefault();
                if (LOSSmin != null) op.LOSSmin = LOSSmin;
                if (LOSSmax != null) op.LOSSmax = LOSSmax;
                if (ETLmax != null) op.ETLmax = ETLmax;
                if (ETLmin != null) op.ETLmin = ETLmin;
                if (op.ETLmax != null || op.ETLmin != null) op.ETLopt = true;
                if(op.ETLopt.GetValueOrDefault() && (op.ETLmax==null||op.ETLmin==null))return Problem(detail:$"ETL max is {op.ETLmax}. ETL min is {op.ETLmin}",title:"Set both bounds for ETL");
                var ETLopt = false;
                try { op.ETLmin = CVarData.mapDouble["CVarMin"][0]; } catch {; }
                try { op.ETLmax = CVarData.mapDouble["CVarMax"][0]; } catch {; }
                try { ETLopt = CVarData.mapInt["CVar_constraint"][0] != 0 ? true : false; } catch {; }
                if (ETLopt) op.ETLopt = ETLopt;
                if (op.LOSSmax != null || op.LOSSmin != null) op.LOSSopt = true;
                if(op.LOSSopt.GetValueOrDefault() && (op.LOSSmax==null||op.LOSSmin==null))return Problem(detail:$"Loss max is {op.LOSSmax}. Loss min is {op.LOSSmin}",title:"Set both bounds for loss");
                if (op.alpha == null)
                {
                    var ones = (double[])new double[op.tlen];
                    op.alpha = new double[op.n.GetValueOrDefault()];
                    BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
                    Factorise.dmxmulv(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, atran: true);
                for(var i=0;i<op.ncomp.GetValueOrDefault();++i){
                    op.alpha[i+op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault()]=BlasLike.ddotvec(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault(),op.composites,op.alpha,astart:i*(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault()));
                }
                }
                if (op.Q == null && op.nfac.GetValueOrDefault() < 0)
                {var nhere=op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                    var Q = new double[nhere * (nhere + 1) / 2]; op.Q = Q;
                    var ij = 0;
                    for (var i = 0; i < nhere; ++i)
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
                    int ETLorLOSSconstraint;
                    try
                    {
                        ETLorLOSSconstraint = CVarData.mapInt["ETLorLOSSconstraint"][0];
                        op.LOSSopt = ETLorLOSSconstraint == 1;
                    }
                    catch {; }
                    double ETLorLOSSmin;
                    try
                    {
                        ETLorLOSSmin = CVarData.mapDouble["ETLorLOSSmin"][0];
                        op.LOSSmin = ETLorLOSSmin;
                    }
                    catch {; }
                    double ETLorLOSSmax;
                    try
                    {
                        ETLorLOSSmax = CVarData.mapDouble["ETLorLOSSmax"][0];
                        op.LOSSmax = ETLorLOSSmax;
                    }
                    catch {; }
                    if (op.LOSSmin == null) op.LOSSmin = -100;
                    if (op.LOSSmax == null) op.LOSSmax = 100;
                    op.TargetReturn=Portfolio.Portfolio.one2many(op.tlen,targetR.GetValueOrDefault());
                }
                else
                {
                    op.TargetReturn = null;  //probably redundant

                    int ETLorLOSSconstraint;
                    try
                    {
                        ETLorLOSSconstraint = CVarData.mapInt["ETLorLOSSconstraint"][0];
                        op.ETLopt = ETLorLOSSconstraint == 1;
                    }
                    catch {; }
                    double ETLorLOSSmin;
                    try
                    {
                        ETLorLOSSmin = CVarData.mapDouble["ETLorLOSSmin"][0];
                        op.ETLmin = ETLorLOSSmin;
                    }
                    catch {; }
                    double ETLorLOSSmax;
                    try
                    {
                        ETLorLOSSmax = CVarData.mapDouble["ETLorLOSSmax"][0];
                        op.ETLmax = ETLorLOSSmax;
                    }
                    catch {; }
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
            if (op.names != null && op.names.Length < op.n.GetValueOrDefault()) op.names = null;
            op.w = new double[op.n.GetValueOrDefault()];
            bool CVARGLprob = false;
#if !NOTUSETRY //This is defined in the .csproj file
            try
            {
#endif
                op.back = Portfolio.Portfolio.OptimiseGeneral(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
                op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta,
                op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket, op.trades, op.min_holding,
                op.min_trade, op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
                op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
                op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
                op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
                op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown, ref CVARGLprob, op.tlen, op.Gstrength.GetValueOrDefault(),
                op.DATA, op.tail, op.TargetReturn, op.ETLopt.GetValueOrDefault() || op.LOSSopt.GetValueOrDefault(), op.TargetReturn == null ? op.ETLmin.GetValueOrDefault() : op.LOSSmin.GetValueOrDefault(),
                op.TargetReturn == null ? op.ETLmax.GetValueOrDefault() : op.LOSSmax.GetValueOrDefault(), op.logfile, ncomp: op.ncomp.GetValueOrDefault(), compw: op.composites);
#if !NOTUSETRY
            }
            catch (Exception popt) { op.message = $"{popt.Message}"; return Problem(title: "Optimiser input scalar variable error", detail: lic.VersionString, type: op.message); }
#endif

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
            op.result.risk = rounder(op.result.risk.GetValueOrDefault());
        }
        else
        {
            if (op.w == null)
            {
                op.result = null;
                op.shake = null;
                op.ogamma = null;
                return op;
            }
        }
        if (op.m.GetValueOrDefault() > 0 && op.A != null)
        {
            op.result.cval = new double[op.m.GetValueOrDefault()];
            for (var i = 0; i < op.m; ++i)
            {
                op.result.cval[i] = BlasLike.ddot(op.n.GetValueOrDefault(), op.A, op.m.GetValueOrDefault(), op.w, 1, i);
            }
        }
        if (op.nfac < 0)
        {
            var cov = new Portfolio.Portfolio("");
            cov.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
            cov.Q = op.Q;
            cov.ncomp=op.ncomp.GetValueOrDefault();
            cov.compw=op.composites;
            cov.makeQ();

            cov.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            op.result.risk = rounder(op.result.risk.GetValueOrDefault());
        }
        else
        {
            var fac = new Portfolio.FPortfolio("");
            fac.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
            fac.nfac = op.nfac.GetValueOrDefault();
            fac.ncomp=op.ncomp.GetValueOrDefault();
            fac.compw=op.composites;
            if (op.SV != null && op.nfac > -1)
            {
                fac.SV = op.SV;
                fac.FC = op.FC;
                fac.FL = op.FL;
            }
            else fac.Q = op.Q;
            fac.makeQ();
            fac.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            op.result.risk = rounder(op.result.risk.GetValueOrDefault());
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
            var iii = 0;
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
                    if (op.L != null && op.U != null && op.L[iii] != op.U[iii]) op.result.minhold = Math.Min(Math.Abs(ww), op.result.minhold.GetValueOrDefault());
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
            var iii = 0;
            foreach (var ww in op.w)
            {
                op.result.turnover += Math.Abs(ww);
                if (Math.Abs(ww) > 1e-8)
                {
                    op.result.trades++;
                    if (op.L != null && op.U != null && op.L[iii] != op.U[iii]) op.result.mintrade = Math.Min(Math.Abs(ww), op.result.mintrade.GetValueOrDefault());
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
                op.result.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref VARindex, breakd,ncomp:op.ncomp.GetValueOrDefault(),compw:op.composites);
                op.result.VAR = VAR;
                op.result.VARindex = VARindex;
                op.result.breakdown = breakd;
            }
            else
            {
                op.result.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, breakd,ncomp:op.ncomp.GetValueOrDefault(),compw:op.composites);
                op.result.breakdown = breakd;
            }
        }

        _logger.LogInformation($"GET general at {DateTimeOffset.Now}");
        return op;
    }
    [HttpPost("general")]
    [RequestSizeLimit(bytes: 2_000_000_000)]
    public Object PostGen(Optimise op)
    {if(op.getmethod!=null&&op.getmethod.ToLower().StartsWith("optimiseinternal"))op.doOpt=true;
        var breakdown = new double[op.n.GetValueOrDefault()];
        var lic = new Licensing.Licence();
        var ok = lic.CheckLicence();
        op.VersionString = lic.VersionString;
        op.isLicensed = ok;
        if (!ok) return Problem(title: "Bad licence", detail: lic.VersionString,type:lic.connectedNames);
        if (!(op.A == null && op.Aas2D == null))
        {
            if (op.Aas2D != null)
            {
                op.A = Portfolio.Portfolio.twoD2oneD(op.m.GetValueOrDefault(), op.n.GetValueOrDefault(), op.Aas2D, transpose: false);
            }
            else if (op.transposeLinearConstraintArray)
            {
                if (op.A != null)
                {
                    op.Aas2D = Portfolio.Portfolio.oneD2twoD(op.m.GetValueOrDefault(), op.n.GetValueOrDefault(), op.A, transpose: true);
                    Factorise.dmx_transpose(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), op.A, op.A);
                }

            }
        }
        if (!ok) return op;
        op.result = new Optimise.checkv();
        if (op.tlen > 0)
        {
            if (op.nfac == null) op.nfac = -1;
            if (op.alpha == null)
            {
                var ones = (double[])new double[op.tlen];
                op.alpha = new double[op.n.GetValueOrDefault()];
                BlasLike.dsetvec(op.tlen, 1.0 / op.tlen, ones);
                Factorise.dmxmulv(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault(), op.tlen, op.DATA, ones, op.alpha, atran: true);
                for(var i=0;i<op.ncomp.GetValueOrDefault();++i){
                    op.alpha[i+op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault()]=BlasLike.ddotvec(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault(),op.composites,op.alpha,astart:i*(op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault()));
                }
            }
            if (op.Q == null && op.nfac.GetValueOrDefault() < 0)
            {var nhere=op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                var Q = new double[nhere * (nhere + 1) / 2]; op.Q = Q;
                var ij = 0;
                for (var i = 0; i < nhere; ++i)
                {
                    for (var j = 0; j <= i; ++j)
                    {
                        Q[ij++] = Solver.Factorise.covariance(op.tlen, op.DATA, op.DATA, i * op.tlen, j * op.tlen);
                    }
                }
                op.Q = Q;
            }
        }
        if (op.names != null && op.names.Length < op.n.GetValueOrDefault()) op.names = null;
        if(op.round>0){
            if( op.min_lot!=null&&  op.min_lot.Length==1){
                op.min_lot=Portfolio.Portfolio.one2many(op.n.GetValueOrDefault(),op.min_lot[0]);
            }
            if( op.size_lot!=null&&  op.size_lot.Length==1){
                op.size_lot=Portfolio.Portfolio.one2many(op.n.GetValueOrDefault(),op.size_lot[0]);
            }
        }
        if (op.doOpt)
        {
            if (op.SV != null && op.nfac > -1)
            {
                if (op.FLas2D != null)
                    op.FL = Portfolio.Portfolio.twoD2oneD(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.FLas2D, transpose: false);
            }
            op.w = new double[op.n.GetValueOrDefault()];
            op.shake = new int[op.n.GetValueOrDefault()];
            double ogamma = 0;
            bool CVARGLprob = false;
            op.back = Portfolio.Portfolio.OptimiseGeneral(op.n.GetValueOrDefault(), op.nfac.GetValueOrDefault(), op.names,
            op.w, op.m.GetValueOrDefault(), op.A, op.L, op.U, op.alpha, op.bench, op.Q, op.gamma.GetValueOrDefault(), op.initial, op.delta,
            op.buy, op.sell, op.kappa.GetValueOrDefault(), op.basket, op.trades, op.min_holding,
            op.min_trade, op.rmin, op.rmax, op.round.GetValueOrDefault(), op.min_lot, op.size_lot, op.shake, op.value,
            op.nabs.GetValueOrDefault(), op.Abs_A, op.mabs.GetValueOrDefault(), op.I_A, op.Abs_U,
            op.FC, op.FL, op.SV, op.minRisk, op.maxRisk, ref ogamma,
            op.mask, op.longbasket.GetValueOrDefault(), op.shortbasket.GetValueOrDefault(), op.tradebuy.GetValueOrDefault(),
            op.tradesell.GetValueOrDefault(), op.valuel, op.Abs_L, breakdown, ref CVARGLprob, op.tlen, op.Gstrength.GetValueOrDefault(),
            op.DATA, op.tail, op.TargetReturn, op.ETLopt.GetValueOrDefault() || op.LOSSopt.GetValueOrDefault(), op.TargetReturn == null ? op.ETLmin.GetValueOrDefault() : op.LOSSmin.GetValueOrDefault(),
            op.TargetReturn == null ? op.ETLmax.GetValueOrDefault() : op.LOSSmax.GetValueOrDefault(), op.logfile, 0,ncomp:op.ncomp.GetValueOrDefault(),compw:op.composites);


            op.CVARGLprob = CVARGLprob;
            op.ogamma = ogamma;
            op.message = Portfolio.Portfolio.OptMessages(op.back.GetValueOrDefault());
            for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
            {
                if (op.shake[i] == i) if (op.names == null) { op.message += $"\nAsset {i + 1} was not rounded properly"; }
                    else { op.message += $"\n{op.names[i]} was not rounded properly"; }
            }
            if(op.getmethod!=null&&op.getmethod.ToLower().StartsWith("optimiseinternal")){
            var backop=new OptimiseInternal();
            backop.VersionString=op.VersionString;
            backop.isLicensed=op.isLicensed;
            backop.back=op.back;
            backop.w=op.w;
            backop.shake=op.shake;
            backop.message=op.message;
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
            return backop;}
        }
        else if (op.w == null)
        {
            op.result = null;
            op.message = null;
            op.ogamma = null;
            op.shake = null;
            op.CVARGLprob = null;
            if (op.nfac > -1 && op.SV != null && op.getmethod != null)
            {//Could use Contains() instead of StartsWidth()
                if (op.getmethod.ToLower().StartsWith("factormodelprocess"))
                {
                    var fac = new Portfolio.FPortfolio("");
                    fac.ntrue = op.n.GetValueOrDefault();
                    fac.nfac = op.nfac.GetValueOrDefault();
                    fac.SV = op.SV;
                    fac.FC = op.FC;
                    if (op.FLas2D != null)
                        fac.FL = Portfolio.Portfolio.twoD2oneD(fac.ntrue, fac.nfac, op.FLas2D, transpose: false);
                    else
                    {
                        fac.FL = op.FL;
                        if (op.FL == null) return Problem(title: "Factor Loadings not set", detail: "Both FL and FLas2D may not be null!");
                    }
                    fac.makeQ();
                    var Qinv=(double[]?)null;
                    var withinverse=false;
                    if (op.getmethod.ToLower().Contains("inverse")){withinverse=true;
                    Qinv=(double[])fac.Q.Clone();
                    var back=Factorise.FMPinverse(fac.ntrue,fac.nfac,Qinv);
                    if (back!=0)withinverse=false;}
                    if(!withinverse){
                    var opnew = new FactorModelProcess();
                    //opnew.FLbacktest=Portfolio.Portfolio.oneD2twoD(fac.ntrue,fac.nfac,op.FL);
                    opnew.VersionString = op.VersionString;
                    opnew.isLicensed = op.isLicensed;
                    opnew.QMATRIX = fac.Q;
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                    return opnew;}else{
                    var opnew = new FactorModelProcessInverse();
                    opnew.VersionString = op.VersionString;
                    opnew.isLicensed = op.isLicensed;
                    opnew.QMATRIX = fac.Q;
                    opnew.QMATRIXinverse=Qinv;
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                    return opnew;}
                }

                else if (op.getmethod.ToLower().StartsWith("factor2cov"))
                {
                    var fac = new Portfolio.FPortfolio("");
                    fac.ntrue = op.n.GetValueOrDefault();
                    fac.nfac = op.nfac.GetValueOrDefault();
                    fac.SV = op.SV;
                    fac.FC = op.FC;
                    if (op.FLas2D != null)
                        fac.FL = Portfolio.Portfolio.twoD2oneD(fac.ntrue, fac.nfac, op.FLas2D, transpose: false);
                    else
                        fac.FL = op.FL;
                    fac.makeQ();
                    var opnew = new Factor2COV();
                    opnew.VersionString = op.VersionString;
                    opnew.isLicensed = op.isLicensed;
                    opnew.COV = fac.Factor2Cov();
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                    return opnew;
                }

                else if (op.getmethod.ToLower().StartsWith("factor2var"))
                {
                    var fac = new Portfolio.FPortfolio("");
                    fac.ntrue = op.n.GetValueOrDefault();
                    fac.nfac = op.nfac.GetValueOrDefault();
                    fac.SV = op.SV;
                    fac.FC = op.FC;
                    if (op.FLas2D != null)
                        fac.FL = Portfolio.Portfolio.twoD2oneD(fac.ntrue, fac.nfac, op.FLas2D, transpose: false);
                    else
                        fac.FL = op.FL;
                    fac.makeQ();
                    var opnew = new Factor2VAR();
                    opnew.VersionString = op.VersionString;
                    opnew.isLicensed = op.isLicensed;
                    opnew.VAR = fac.Factor2Var();
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                    return opnew;
                }
            }
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
            return op;
        }
        else
        {
            if (op.alpha == null) op.alpha = new double[op.n.GetValueOrDefault()];
        }
        if (op.getmethod != null)
        {
            if (op.getmethod.ToLower().StartsWith("riskproperties"))
            {
                var props = new RiskProperties();
                props.VersionString = op.VersionString;
                props.isLicensed = op.isLicensed;
                Portfolio.Portfolio riskprop;
                if (op.nfac > -1)
                {
                    var fac = new Portfolio.FPortfolio("");
                    riskprop = fac;
                    fac.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                    fac.nfac = op.nfac.GetValueOrDefault();
                    if (op.SV != null && op.nfac > -1)
                    {
                        fac.SV = op.SV;
                        fac.FC = op.FC;
                        if (op.FLas2D != null)
                            fac.FL = Portfolio.Portfolio.twoD2oneD(fac.ntrue, fac.nfac, op.FLas2D, transpose: false);
                        else fac.FL = op.FL;
                    }
                    else fac.Q = op.Q;
                }
                else
                {
                    var cov = new Portfolio.Portfolio("");
                    riskprop = cov;
                    cov.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                    cov.Q = op.Q;
                }
                riskprop.n=op.n.GetValueOrDefault();
                riskprop.ncomp=op.ncomp.GetValueOrDefault();
                riskprop.compw=op.composites;
                riskprop.makeQ();
                props.marginalContributionToTotalRisk = new double[riskprop.n];
                riskprop.RiskBreakdown(op.w,breakdown: props.marginalContributionToTotalRisk);
                props.totalRisk = BlasLike.ddotvec(riskprop.n, op.w, props.marginalContributionToTotalRisk);
                props.totalRisk = rounder(props.totalRisk.GetValueOrDefault());
                if (op.bench != null)
                {
                    props.assetBetas = new double[riskprop.n];
                    props.marginalContributionToActiveRisk = new double[riskprop.n];
                    props.marginalContributionToBenchmarkRisk = new double[riskprop.n];
                    props.marginalContributionToResidualRisk = new double[riskprop.n];
                    riskprop.RiskBreakdown(op.w, op.bench, props.marginalContributionToActiveRisk);
                    riskprop.RiskBreakdown(op.w, op.bench, props.marginalContributionToResidualRisk, props.assetBetas);
                    props.portfolioBeta = BlasLike.ddotvec(riskprop.ntrue, op.w, props.assetBetas);
                    riskprop.RiskBreakdown(op.bench, null, props.marginalContributionToBenchmarkRisk);
                    props.activeRisk = BlasLike.ddotvec(riskprop.ntrue, op.w, props.marginalContributionToActiveRisk) - BlasLike.ddotvec(riskprop.ntrue, op.bench, props.marginalContributionToActiveRisk);
                    props.activeRisk = rounder(props.activeRisk.GetValueOrDefault());
                    props.benchmarkRisk = BlasLike.ddotvec(riskprop.ntrue, op.bench, props.marginalContributionToBenchmarkRisk);
                    props.benchmarkRisk = rounder(props.benchmarkRisk.GetValueOrDefault());
                    BlasLike.daxpyvec(riskprop.ntrue, -props.portfolioBeta.GetValueOrDefault(), op.bench, op.w);
                    props.residualRisk = BlasLike.ddotvec(riskprop.ntrue, op.w, props.marginalContributionToResidualRisk);
                    props.residualRisk = rounder(props.residualRisk.GetValueOrDefault());
                    BlasLike.daxpyvec(riskprop.ntrue, props.portfolioBeta.GetValueOrDefault(), op.bench, op.w);
                }
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                return props;
            }

            else if (op.getmethod.ToLower().StartsWith("risks"))
            {
                var props = new Risks();
                props.VersionString = op.VersionString;
                props.isLicensed = op.isLicensed;
                Portfolio.Portfolio risks;
                if (op.nfac > -1)
                {
                    var fac = new Portfolio.FPortfolio("");
                    risks = fac;
                    fac.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                    fac.nfac = op.nfac.GetValueOrDefault();
                    if (op.SV != null && op.nfac > -1)
                    {
                        fac.SV = op.SV;
                        fac.FC = op.FC;
                        if (op.FLas2D != null)
                            fac.FL = Portfolio.Portfolio.twoD2oneD(fac.ntrue, fac.nfac, op.FLas2D, transpose: false);
                        else fac.FL = op.FL;
                    }
                    else fac.Q = op.Q;
                }
                else
                {
                    var cov = new Portfolio.Portfolio("");
                    risks = cov;
                    cov.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
                    cov.Q = op.Q;
                }
                risks.n=op.n.GetValueOrDefault();
                risks.ncomp=op.ncomp.GetValueOrDefault();
                risks.compw=op.composites;
                risks.makeQ();
                var mctr = new double[risks.n];
                risks.RiskBreakdown(op.w, breakdown: mctr);
                props.totalRisk = BlasLike.ddotvec(risks.ntrue, op.w, mctr);
                props.totalRisk = rounder(props.totalRisk.GetValueOrDefault());
                if (op.bench != null)
                {
                    var betas = new double[risks.n];
                    risks.RiskBreakdown(op.w, op.bench, mctr);
                    props.activeRisk = BlasLike.ddotvec(risks.ntrue, op.w, mctr) - BlasLike.ddotvec(risks.ntrue, op.bench, mctr);
                    props.activeRisk = rounder(props.activeRisk.GetValueOrDefault());
                    risks.RiskBreakdown(op.w, op.bench, mctr, betas);
                    props.portfolioBeta = BlasLike.ddotvec(risks.ntrue, op.w, betas);
                    BlasLike.daxpyvec(risks.ntrue, -props.portfolioBeta.GetValueOrDefault(), op.bench, op.w);
                    props.residualRisk = BlasLike.ddotvec(risks.ntrue, op.w, mctr);
                    props.residualRisk = rounder(props.residualRisk.GetValueOrDefault());
                    BlasLike.daxpyvec(risks.ntrue, props.portfolioBeta.GetValueOrDefault(), op.bench, op.w);
                    risks.RiskBreakdown(op.bench, breakdown: mctr);
                    props.benchmarkRisk = BlasLike.ddotvec(risks.ntrue, op.bench, mctr);
                    props.benchmarkRisk = rounder(props.benchmarkRisk.GetValueOrDefault());
                }
             _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
                return props;
            }
        }
        op.result.mctr = new double[op.n.GetValueOrDefault()];


        if (op.nfac < 0)
        {
            var cov = new Portfolio.Portfolio("");
            cov.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
            cov.Q = op.Q;
            if (op.bench != null)
            {
                op.result.BETA = new double[cov.ntrue];
                cov.RiskBreakdown(op.w, op.bench, op.result.mctr, op.result.BETA);//breakdown is for residual position. Probably don't want this so call again below
                op.result.portBETA = BlasLike.ddotvec(cov.n, op.w, op.result.BETA);
            }
            cov.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            op.result.risk = rounder(op.result.risk.GetValueOrDefault());
        }
        else
        {
            var fac = new Portfolio.FPortfolio("");
            fac.ntrue = op.n.GetValueOrDefault()-op.ncomp.GetValueOrDefault();
            fac.nfac = op.nfac.GetValueOrDefault();
            if (op.SV != null && op.nfac > -1)
            {
                fac.SV = op.SV;
                fac.FC = op.FC;
                fac.FL = op.FL;
                fac.makeQ();
            }
            else fac.Q = op.Q;
            if (op.bench != null)
            {
                op.result.BETA = new double[op.n.GetValueOrDefault()];
                fac.RiskBreakdown(op.w, op.bench, op.result.mctr, op.result.BETA);//breakdown is for residual position. Probably don't want this so call again below
                op.result.portBETA = BlasLike.ddotvec(fac.ntrue, op.w, op.result.BETA);
            }
            fac.RiskBreakdown(op.w, op.bench, op.result.mctr);
            op.result.risk = BlasLike.ddotvec(op.n.GetValueOrDefault(), op.w, op.result.mctr);
            if (op.bench != null) op.result.risk -= BlasLike.ddotvec(op.n.GetValueOrDefault(), op.bench, op.result.mctr);
            op.result.risk = rounder(op.result.risk.GetValueOrDefault());
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
            if (op.m.GetValueOrDefault() > 0 && op.A != null)
            {
                op.result.cval = new double[op.m.GetValueOrDefault()];
                for (var i = 0; i < op.m; ++i)
                {
                    op.result.cval[i] = BlasLike.ddot(op.n.GetValueOrDefault(), op.A, op.m.GetValueOrDefault(), op.w, 1, i);
                }
            }
            op.result.basket = 0;
            op.result.trades = 0;
            op.result.minhold = BlasLike.lm_max;
            op.result.gross = 0;
            op.result.longvalue = 0; op.result.shortvalue = 0;
            var iii = 0;
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
                    if (op.L != null && op.U != null && op.L[iii] != op.U[iii]) op.result.minhold = Math.Min(Math.Abs(ww), op.result.minhold.GetValueOrDefault());
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
            var iii = 0;
            foreach (var ww in op.w)
            {
                op.result.turnover += Math.Abs(ww);
                if (Math.Abs(ww) > 1e-8)
                {
                    op.result.trades++;
                    if (op.L != null && op.U != null && op.L[iii] != op.U[iii]) op.result.mintrade = Math.Min(Math.Abs(ww), op.result.mintrade.GetValueOrDefault());
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
                op.result.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail, ref VAR, ref VARindex, breakd,ncomp:op.ncomp.GetValueOrDefault(),compw:op.composites);
                op.result.VAR = VAR;
                op.result.VARindex = VARindex;
                op.result.breakdown = breakd;
            }
            else
            {
                op.result.LOSS = Portfolio.Portfolio.LOSS(op.n.GetValueOrDefault(), op.w, op.DATA, op.TargetReturn, breakd,ncomp:op.ncomp.GetValueOrDefault(),compw:op.composites);
                op.result.breakdown = breakd;
            }
        }
        _logger.LogInformation($"POST general at {DateTimeOffset.Now}");
        return op;
    }
}

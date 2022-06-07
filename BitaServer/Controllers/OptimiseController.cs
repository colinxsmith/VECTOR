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
    [HttpGet("ETL")]
    public Optimise[] GetETL()
    {
        var op = new Optimise();
        using (var CVarData = new InputSomeData())
        {
            CVarData.doubleFields = "DATA";//The DATA is LOSS i.e. - returns
            CVarData.intFields = "n tlen";
            CVarData.stringFields = "names";
            try
            {
                CVarData.Read("./GLdist");
            }
            catch
            {
                CVarData.Read("../GLdist");
            }
            op.DATA = CVarData.mapDouble["DATA"];
            op.names = CVarData.mapString["names"];
            op.n = CVarData.mapInt["n"][0];
            op.tlen = CVarData.mapInt["tlen"][0];
        }
        op.m = 1;
        op.L = new double[op.n.Value + op.m.Value];
        op.U = new double[op.n.Value + op.m.Value];
        op.A = new double[op.n.Value * op.m.Value];
        BlasLike.dsetvec(op.n.Value, 1, op.U);
        BlasLike.dsetvec(op.n.Value, 1, op.A);
        op.L[op.n.Value] = op.U[op.n.Value] = 1.0;
        op.ETLopt = false;
        op.Gstrength = 1;
        return new[] { op };
    }
    [HttpPost("ETL")]
    public Optimise[] PostETL(Optimise op)
    {
        if (!op.gamma.HasValue) op.gamma = 0.5;
        if (!op.kappa.HasValue) op.kappa = op.gamma;
        double[] Q;
        if(op.Q==null){Q = new double[op.n.GetValueOrDefault() * (op.n.GetHashCode() + 1) / 2];op.Q=Q;
                var ij = 0;
        for (var i = 0; i < op.n.GetValueOrDefault(); ++i)
        {
            for (var j = 0; j <= i; ++j)
            {
                Q[ij++] = Solver.Factorise.covariance(op.tlen.GetValueOrDefault(), op.DATA, op.DATA, i * op.tlen.GetValueOrDefault(), j * op.tlen.GetValueOrDefault());
            }
        }}
        else Q=op.Q;
        var opt = new Portfolio.Portfolio("");
        var ones = new double[op.tlen.GetValueOrDefault()];
        var alpha = new double[op.n.GetValueOrDefault()];
        opt.Q = Q;
        BlasLike.dsetvec(op.tlen.GetValueOrDefault(), 1.0 / op.tlen.GetValueOrDefault(), ones);
        Factorise.dmxmulv(op.n.GetValueOrDefault(), op.tlen.GetValueOrDefault(), op.DATA, ones, alpha, 0, 0, 0, true);
        op.back = opt.BasicOptimisation(op.n.GetValueOrDefault(), op.m.GetValueOrDefault(), -1,
           op.A, op.L, op.U, op.gamma.GetValueOrDefault(), op.kappa.GetValueOrDefault(), -1, -1, -1, -1, -1, alpha, op.initial, null, null,
           op.names, false, 0, null, null, null, 0, null, op.tlen.GetValueOrDefault(),
            op.Gstrength.GetValueOrDefault(), op.DATA, op.tail.GetValueOrDefault(),
            null, op.ETLopt.GetValueOrDefault(), op.ETLmin.GetValueOrDefault(),
            op.ETLmax.GetValueOrDefault());
        op.w = opt.wback;
        op.breakdown = new double[op.n.GetValueOrDefault()];
        double VAR = 0;
        int ind = -1;
        op.ETL = Portfolio.Portfolio.ETL(op.n.GetValueOrDefault(), op.w, op.DATA, op.tail.GetValueOrDefault(), ref VAR, ref ind, op.breakdown);
        op.VAR = VAR;
        op.VARindex = ind;
        return new[] { op };
    }
    [HttpGet("test")]
    public Optimise[] Get()
    {
        double testdigit = 123.9999999999;
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
}

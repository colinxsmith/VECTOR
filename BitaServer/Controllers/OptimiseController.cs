using Microsoft.AspNetCore.Mvc;

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
    public IEnumerable<Optimise> Post(Optimise op)
    {
        if (op.digit != null)
        {
            op.tdigit = op.digit != null ? Portfolio.Portfolio.check_digit((double)op.digit) : null;
        }
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        var back = (Optimise[])new Optimise[1];
        back[0] = op;
        return back;
    }
    [HttpGet("test")]
    public IEnumerable<Optimise> Get()
    {
        double testdigit = 123.9999999999;
        var back = (Optimise[])new Optimise[1];
        back[0] = new Optimise();
        back[0].digit = testdigit;
        back[0].tdigit = Portfolio.Portfolio.check_digit(testdigit);
        var op=back[0];
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return back;
    }
    [HttpGet("test/n")]
    public IEnumerable<Optimise> Getn()
    {
        double testdigit = 11.000000000001;
        var back = (Optimise[])new Optimise[1];
        back[0] = new Optimise();
        back[0].n = 500;
        back[0].digit = testdigit;
        back[0].tdigit = Portfolio.Portfolio.check_digit(testdigit);
        var op = back[0];
        ColourConsole.WriteEmbeddedColourLine($"[red]{op.digit}[/red] [green]{op.tdigit}[/green] [yellow]{op.n}[/yellow]");
        return back;
    }
}

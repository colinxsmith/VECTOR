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
    [HttpPost]
    public IEnumerable<Optimise> Post(Optimise op)
    {
        var dd = Portfolio.Portfolio.check_digit(op.digit);
        op.tdigit = dd;
        var back = Enumerable.Range(0, 1).Select(aa => op).ToArray();
        return back;
    }
    [HttpGet(Name = "GetOptimise")]
    public IEnumerable<Optimise> Get()
    {
        double testdigit = 2.9999999999;
        var back = Enumerable.Range(0, 1).Select(aa => new Optimise
        {
            digit = testdigit,
            tdigit = Portfolio.Portfolio.check_digit(testdigit)
        }).ToArray();
        return back;
    }
}

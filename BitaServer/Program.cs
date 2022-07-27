using Microsoft.Extensions.Hosting.WindowsServices;
using Microsoft.Extensions.Logging.EventLog;
using Microsoft.Extensions.Hosting.Systemd;
var options = new WebApplicationOptions
{
    Args = args,
    ContentRootPath = (WindowsServiceHelpers.IsWindowsService() || SystemdHelpers.IsSystemdService()) ? AppContext.BaseDirectory : default
};
var builder = WebApplication.CreateBuilder(options);
// Add services to the container.

builder.Services.AddControllers();
// Learn more about configuring Swagger/OpenAPI at https://aka.ms/aspnetcore/swashbuckle
builder.Services.AddEndpointsApiExplorer();
if (OperatingSystem.IsWindows())
{
    builder.Host
    .UseWindowsService()
    .ConfigureLogging((_, logging) => logging.AddEventLog())
    .ConfigureLogging(logging =>
    {
        logging.ClearProviders();
        logging.AddConsole();
    });//Seems not to work for getting log info
    builder.Services.Configure<EventLogSettings>(conf =>
    {
        conf.LogName = "BitaServer";
        conf.SourceName = "Optimiser Server";
        conf.MachineName = null;
    });
}
else if (OperatingSystem.IsLinux())
{
    builder.Host.UseSystemd();
}
//builder.WebHost.UseUrls("http://*:7779"); //Don't do this
//builder.Services.AddSwaggerGen();
var app = builder.Build();

// Configure the HTTP request pipeline.
/*if (app.Environment.IsDevelopment())
{
    app.UseSwagger();
    app.UseSwaggerUI();
}*/

//app.UseHttpsRedirection();
app.UseStaticFiles();
//app.UseAuthorization();

app.MapControllers();

//app.Run();

await app.RunAsync();
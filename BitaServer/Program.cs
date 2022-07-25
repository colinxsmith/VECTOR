using Microsoft.Extensions.Hosting.WindowsServices;
using Microsoft.Extensions.Logging.EventLog;
var options = new WebApplicationOptions
{
    Args = args,
    ContentRootPath = WindowsServiceHelpers.IsWindowsService() ? AppContext.BaseDirectory : default
};
var builder = WebApplication.CreateBuilder(options);

// Add services to the container.

builder.Services.AddControllers();
// Learn more about configuring Swagger/OpenAPI at https://aka.ms/aspnetcore/swashbuckle
builder.Services.AddEndpointsApiExplorer();
builder.Host.UseWindowsService();
//builder.Services.AddSwaggerGen();
if (OperatingSystem.IsWindows())
{//Seems not to work for getting log info
    builder.Services.Configure<EventLogSettings>(conf =>
    {
        conf.LogName = "BitaServer";
        conf.SourceName = "BitaServer";
    });
}
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
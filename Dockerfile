FROM mcr.microsoft.com/dotnet/sdk:8.0 AS build
WORKDIR /source

COPY BitaServer/*.csproj ./BitaServer/
COPY BlasLike/*.csproj ./BlasLike/

COPY BitaServer/. ./BitaServer/
COPY BlasLike/. ./BlasLike/
WORKDIR /source/BitaServer
RUN dotnet publish -c release -o /app 
COPY UseBlas/bin/Debug/net8.0/licence /app
COPY generalopt /app

FROM mcr.microsoft.com/dotnet/aspnet:8.0 AS runtime
WORKDIR /app
COPY --from=build /app ./
EXPOSE 7778
ENV PATH=/usr/bin:/bin:/app
ENTRYPOINT ["BitaServer"]
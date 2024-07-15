FROM mcr.microsoft.com/dotnet/sdk:8.0 AS build
WORKDIR /source

COPY BitaServer/*.csproj ./BitaServer/
COPY BlasLike/*.csproj ./BlasLike/

COPY BitaServer/. ./BitaServer/
COPY BlasLike/. ./BlasLike/
WORKDIR /source/BitaServer
RUN dotnet publish -c release -o /app 
#COPY UseBlas/bin/Debug/net8.0/licence /app
#I ran  future -b 13101D54 10/07/2024 23/12/2024 to get licence_base_n
COPY licence_base_n /app/licence
COPY generalopt /app

FROM alpine AS runtime
RUN apk add aspnetcore-runtime-8.0
WORKDIR /app
COPY --from=build /app ./
EXPOSE 7778
#ENV PATH=/usr/bin:/bin:/app
ENTRYPOINT ["dotnet","BitaServer.dll"]
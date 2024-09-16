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

# Showing how to run as colin and set up password copied from COR stuff. Edit in properly sometime
FROM mcr.microsoft.com/dotnet/aspnet:8.0 AS runtime
#curl wget and ping are not necessary but helpful when testing
RUN apt update && apt upgrade -y && apt install -y curl wget inetutils-ping sudo
#Show how to run as user. Set up the the usual sudo that would be available on linux debian-like
#run passwd from the linux prompt to set user password
#RUN useradd -m -N -s/bin/bash -u 1000 -p "" colin && usermod -aG sudo colin
#Or if you have one, set an encrypted password
#I generated this on linux using perl command crypt $password, $salt
RUN useradd -m -N -s/bin/bash -u 1000 -p "aaTUg7pzsFzvk" colin && usermod -aG sudo colin
#BitaServer is configured to serve on port 7777 so we must expose it here,
#then when running use -p 1234:7777 to map 7777 onto an external port eg. 1234
EXPOSE 7777

WORKDIR /home/colin
COPY --from=build /app ./
#give user access to run programs just copied from /app
RUN chown -R colin /home/colin
USER colin
WORKDIR /home/colin
ENV PATH=/usr/bin:/bin:/home/colin
ENTRYPOINT ["BitaServer"]

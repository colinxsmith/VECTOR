FROM mcr.microsoft.com/dotnet/aspnet:8.0 AS runtime
WORKDIR /top
COPY pusher/ /top
EXPOSE 7778
ENV PATH=$PATH:/top
ENTRYPOINT ["BitaServer"]
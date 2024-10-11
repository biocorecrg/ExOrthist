FROM rocker/r-ver:3.6.3

RUN apt-get -y update && apt-get install -y libnlopt-dev libxml2-dev libcurl4-openssl-dev \
  && rm -rf /var/lib/apt/lists/*

COPY deps.R /tmp
RUN Rscript /tmp/deps.R > /tmp/deps.log

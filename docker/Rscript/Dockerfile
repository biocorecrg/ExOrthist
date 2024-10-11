FROM rocker/r-ver:3.6.3

RUN apt-get -y update && apt-get install -y libnlopt-dev \
  && rm -rf /var/lib/apt/lists/*

COPY deps.R /tmp
RUN Rscript /tmp/deps.R > /tmp/deps.log

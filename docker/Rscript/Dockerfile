FROM rocker/r-ver:4.4.1

RUN apt-get -y update && apt-get install -y libnlopt-dev zlib1g-dev \
  libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libtiff5-dev libjpeg-dev build-essential libcurl4-gnutls-dev libxml2-dev \ 
  libssl-dev libfontconfig1-dev cmake \
  && rm -rf /var/lib/apt/lists/*

ARG REPO=http://cran.us.r-project.org

RUN install2.r -s -e -r ${REPO} remotes
RUN Rscript -e "remotes::install_github('nathan-russell/hashmap')"
RUN install2.r -s -e -r ${REPO} igraph ggplot2
RUN install2.r -s -e -r ${REPO} pbkrtest lme4 pbkrtest lme4
RUN install2.r -s -e -r ${REPO} ggpubr cowplot reshape2 nloptr

RUN rm -rf /tmp/downloaded_packages \
  && strip /usr/local/lib/R/site-library/*/libs/*.so

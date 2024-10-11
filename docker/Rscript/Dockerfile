FROM rocker/r-ver:4

RUN apt-get -y update && apt-get install -y libnlopt-dev zlib1g-dev \
  libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev \ 
  libssl-dev libfontconfig1-dev \
  && rm -rf /var/lib/apt/lists/*

RUN install2.r --error devtools 
COPY deps.R /tmp
RUN Rscript /tmp/deps.R > /tmp/deps.log
RUN install2.r --error igraph ggplot2 nplotr 
# RUN Rscript -e "install.packages(c('pbkrtest', 'lme4', 'pbkrtest', 'lme4'), repos='http://cran.us.r-project.org');"
# RUN Rscript -e "install.packages(c('ggpubr', 'cowplot', 'rshape2'), repos='http://cran.us.r-project.org');"
# devtools::install_github('nathan-russell/hashmap')


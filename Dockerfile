FROM biocorecrg/centos-perlbrew-pyenv3-java:centos7

# File Author / Maintainer
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

ARG R_VERSION=3.6.0
ARG MAFT_VERSION=x.x.x

#upgrade pip
RUN pip install --upgrade pip


#Installing R
RUN yum install -y epel-release libxml2-devel libcurl-devel 
RUN yum install R-${R_VERSION} -y
RUN mkdir -p /usr/share/doc/R-${R_VERSION}/html
RUN Rscript -e "install.packages(c('data.table','flexdashboard','dplyr','plyr','ggExtra','ggplot2','hexbin','knitr','optparse','RColorBrewer','reshape2'), repos='http://cran.us.r-project.org')"


# Clean cache
RUN yum clean all 

#cleaning
RUN rm -fr *.tar.gz; rm -fr *.bz2
RUN rm -rf /var/cache/yum
ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8

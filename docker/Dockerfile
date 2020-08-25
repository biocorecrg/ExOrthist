FROM biocorecrg/centos-perlbrew-pyenv3-java:centos7

# File Author / Maintainer
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

ARG R_VERSION=3.6.0
ARG MAFFT_VERSION=7.429

#upgrade pip
RUN pip install --upgrade pip

#Installing MAFFT
RUN bash -c 'curl -k -L https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-with-extensions-src.tgz > mafft.tgz'
RUN tar -zvxf mafft.tgz
RUN cd mafft-${MAFFT_VERSION}-with-extensions/core;make clean; make; make install; cd ../../; rm mafft.tgz
 
#Installing R
RUN yum install -y epel-release libxml2-devel libcurl-devel 
RUN yum install R-${R_VERSION} -y
RUN mkdir -p /usr/share/doc/R-${R_VERSION}/html
RUN Rscript -e "install.packages(c('igraph'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('ggplot2'), repos='http://cran.us.r-project.org')"
RUN yum install -y NLopt-devel.x86_64
RUN Rscript -e "install.packages(c('nloptr'), repos='http://cran.us.r-project.org')" 
RUN Rscript -e "install.packages(c('pbkrtest', 'lme4', 'pbkrtest', 'lme4'), repos='http://cran.us.r-project.org')"
#RUN echo ciao
RUN Rscript -e "install.packages(c('ggpubr', 'cowplot'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('hashmap', 'reshape'), repos='http://cran.us.r-project.org')"
#hashmap
#reshape
#Installing python modules
RUN yum install -y openssl
RUN pip install --upgrade virtualenv
RUN pip install pandas

# Clean cache
RUN yum clean all 

#cleaning
RUN rm -fr *.tar.gz; rm -fr *.bz2
RUN rm -rf /var/cache/yum
ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8

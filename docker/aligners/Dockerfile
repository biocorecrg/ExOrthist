FROM perl:5.24-threaded-buster

ARG MAFFT_VERSION=7.429

RUN apt-get update -y && apt-get install -y build-essential clustalw && rm -rf /var/lib/apt/lists/*

#Installing MAFFT
WORKDIR /tmp
RUN curl -k -L https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-with-extensions-src.tgz > mafft.tgz
RUN tar -zvxf mafft.tgz
RUN cd mafft-${MAFFT_VERSION}-with-extensions/core;make clean; make; make install; cd ../../; rm mafft.tgz

WORKDIR /

# ENV LC_ALL=en_US.UTF8
# ENV LANG=en_US.UTF8

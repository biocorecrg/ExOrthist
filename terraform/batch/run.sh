#!/usr/bin/env bash

git clone https://github.com/biocorecrg/ExOrthist
cd ExOrthist
git checkout aws

nextflow run main.nf -with-singularity -bg -config params.config.test -profile awsbatch -with-report -with-trace > /mnt/data-nf/log.txt

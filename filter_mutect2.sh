#!/bin/bash

for f in */*snpEff.ann.vcf.gz; \
sudo docker run -v ./:/data/ biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools view -f "PASS" $f | gzip >  "${$f%%.*}.ann.passOnly.vcf.gz";done

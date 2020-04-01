#!/bin/bash

while read x
do
bwa mem -P genomic_hg38/GCF_000001405.38_GRCh38.p12_genomic.fna ${x}_1.paired.fastq ${x}_2.paired.fastq > align_hg38_${x}.sam
done<$1

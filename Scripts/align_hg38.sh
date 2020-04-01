#!/bin/bash

mkdir genomic_hg38
cd genomic_hg38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz
gunzip GCF_000001405.38_GRCh38.p12_genomic.fna.gz
bwa index GCF_000001405.38_GRCh38.p12_genomic.fna
cd ../

while read x
do
bwa mem -P genomic_hg38/GCF_000001405.38_GRCh38.p12_genomic.fna ${x}_1.paired.fastq ${x}_2.paired.fastq > align_hg38_${x}.sam
done<$1

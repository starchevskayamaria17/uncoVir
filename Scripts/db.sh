#/bin/bash

mkdir db
cd db

#Скачивание и индексирование генома hg38

mkdir genomic_hg38
cd genomic_hg38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz
gunzip GCF_000001405.38_GRCh38.p12_genomic.fna.gz
bwa index GCF_000001405.38_GRCh38.p12_genomic.fna
cd ../


#Создание базы данных UNIVEC

mkdir UNIVEC
cd UNIVEC
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
makeblastdb -in UniVec_Core -dbtype nucl -out UniVec_Core.db

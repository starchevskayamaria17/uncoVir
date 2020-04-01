#!/bin/bash

while read x
do
{ПУТЬ ДО ДИРЕКТОРИИ С BBMAP}/bbmap/repair.sh in=${x}_1.fastq in2=${x}_2.fastq out=${x}_1.paired.fastq out2=${x}_2.paired.fastq outs=${x}.unpaired.fastq
mv ${x}_1.pair.fastq ${x}_1.fastq
mv ${x}_2.pair.fastq ${x}_2.fastq
done<$1

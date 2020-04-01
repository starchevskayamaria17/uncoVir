#!/bin/bash

while read x
do
spades.py --meta -1 ${x}_1.fastq -2 ${x}_2.fastq -o AssemblySpades_${x}
done<$1

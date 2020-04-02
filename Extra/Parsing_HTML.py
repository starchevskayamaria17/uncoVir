#Скрипт создает файл Overseq_del.sh, при запуске которого все перепредставленные последовательности будут удалены из FATSTQ-файла.
#Для выполнения необходим файл fastqc

import pandas as pd
table = pd.read_html('fastqc.html')
k=len(table)
if k == 2:
    string = table [1]
    my_list = string['Sequence'].tolist()
    print(my_list)
    f=open('seq.txt', 'w')
    f.writelines(["%s\n" % item  for item in my_list])
    f.close()
    f=open('seq.txt', 'r').read().splitlines()
    f1=open('Overseq_del.sh', 'w')
    print('#!/bin/bash', file=f1)
    for i in range(len(f)):
        print("""awk 'BEGIN {OFS = "\\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if ( seq !~ /""", f[i], "/ )  {print header, seq, qheader, qseq}}' < $1.fastq > $1.",i,".fastq;",sep='', file=f1)
        print('rm $1.fastq;', file=f1)
        break
    for i in range(1, len(f)):
        print("""awk 'BEGIN {OFS = "\\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if ( seq !~ /""", f[i], "/ )  {print header, seq, qheader, qseq}}' < $1.", i-1, ".fastq > $1.",i,".fastq;",sep='', file=f1)
        print('rm $1.', i-1, '.fastq;', sep='', file=f1)
else:
    print('No overrepresented sequences')

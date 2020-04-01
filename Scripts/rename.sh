#!bin/bash/

#Переименовать ID fastq-файлов, добавить в начало каждого ID название файла
#На вход подаем файл с SRR (без указания forward/reverse)

while read x
do

awk -v x="$x" '{print (NR%4 == 1) ? "@" x "_" ++i "/1" : $0}' ${x}_1.fastq > ${x}_1_rename.fastq
mv ${x}_1_rename.fastq ${x}_1.fastq

awk -v x="$x" '{print (NR%4 == 1) ? "@" x "_" ++i "/2" : $0}' ${x}_2.fastq > ${x}_2_rename.fastq
mv ${x}_2_rename.fastq ${x}_2.fastq

done<$1

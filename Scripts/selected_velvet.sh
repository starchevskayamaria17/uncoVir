#!bin/bash 

#Отбор контигов после сборки Velvet с покрытием более 10 и длиной более 300

#Получим список всех id контигов

grep '>' contigs.fa > all_name_contigs.txt

#В id заключена информация о длине контига, а также покрытии (cov_ЧИСЛО), получим список id контигов с покрытием менее 10

grep  'cov_0\.' all_name_contigs.txt > name_bad_coverage_contigs.txt
grep  'cov_1\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_2\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_3\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_4\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_5\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_6\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_7\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_8\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_9\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt

#Преобразуем ID контига в релевантный для seqfilter формат (уберем '>')

cat name_bad_coverage_contigs.txt | cut -c 2- > name_bad_coverage_contigs_1.txt
mv name_bad_coverage_contigs_1.txt name_bad_coverage_contigs.txt

#На основании списка ID, отфильтруем исходный fasta-файл с контигами (опция "n" для инфертированного вывода), а также отфильтруем контиги по длине (опция "m"- минимальная длина контига)

/seqfilter \
	-i contigs.fasta \
	-l name_bad_coverage_contigs.txt \
	-o contigs_coverage_more_10.fasta \
	-m 300 \
	-n \

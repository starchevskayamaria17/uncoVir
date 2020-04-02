#!bin/bash 

#Отбор контигов после сборки Spades с покрытием более 30 и длиной контига более 500
#Получим список всех id контигов

grep '>' contigs.fasta > all_name_contigs.txt

#В id заключена информация о длине контига, а также покрытии (cov_ЧИСЛО), получим список id контигов с покрытием менее 30

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
grep  'cov_10\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_11\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_12\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_13\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_14\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_15\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_16\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_17\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_18\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_19\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_20\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_21\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_22\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_23\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_24\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_25\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_26\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_27\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_28\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_29\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt
grep  'cov_30\.' all_name_contigs.txt >> name_bad_coverage_contigs.txt

#Преобразуем ID контига в релевантный для seqfilter формат (уберем '>')

cat name_bad_coverage_contigs.txt | cut -c 2- > name_bad_coverage_contigs_1.txt
mv name_bad_coverage_contigs_1.txt name_bad_coverage_contigs.txt

#На основании списка ID, отфильтруем исходный fasta-файл с контигами (опция "n" для инфертированного вывода), а также отфильтруем контиги по длине (опция "m"- минимальная длина контига)

/seqfilter \
	-i contigs.fasta \
	-l name_bad_coverage_contigs.txt \
	-o contigs_coverage_more_30.fasta \
	-m 500 \
	-n \

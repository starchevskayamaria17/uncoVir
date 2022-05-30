# VirMercenary

Pipeline for viruses detection in NGS data. 

```
conda env create snakemake.yml
```

Creating databases for classification using BLAST

Download databases of nucleotide and amino acid sequences. Attention, the size of the database NCBI nt in compressed format is 179G, and the database NCBI nr is 128G. After unzipping, you will need about 500G of free disk space. You can use other smaller databases for your needs.

```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz





Предполагается, что данные секвенирования представлены в виде forward/reverse файлов в формате FASTQ. 



Для single прочтений и других форматов смотрите в директории Extra. 
Также предполагается, что проведен контроль качетсва образцов (FastQC) и данные предобработаны (удалены адаптеры, низкокачественные и короткие последовательности). Примеры обработки можно найти в директории Extra.

До начала работы запустите: source export_path.sh.

Подготовьте файл db.sh, добавив все необходимые базы данных и референсные геномы. Запустите.

Последовательность запуска скриптов:
1. rename.sh 

Необязательно. Переименовывает id fastq-файла в соответствие с названием образца. Удобно в случае последующего пулирования образцов.

2. Filtration_UNIVEC.sh 

Контроль контаминации: элиминирование синтетических последовательностей (векторов, линкеров, адаптеров и др.)

3. repair.sh

Позволяет сохранить парность прочтений.

4. align_hg38.sh

Выравнивание прочтений на геном человека. Вы можете также указать любой другой геном. 

5. filtration_hg38.sh

Фильтрация выравненных на геном человека прочтений. Указав другой геном в п.4 прочтения будут отфильтрованы соответственно на представленный геном. 

6. Assembly_spades.sh

Сборка фильтрованных прочтений с использованим SPAdes.

7. Assembly_velvet.sh 

Сборка фильтрованных прочтения с использованием VELVET.

8. selected_spades.sh

Отбор контигов, полученных в результате сборки SPAdes.

9. selected_velvet.sh

Отбор контигов, полученных в результате сборки VELVET.

# VirMercenary

Pipeline for viruses detection in NGS data. 

```
conda env create snakemake.yml
```

1) Reference genome indexing

Download human genome and host genome, and unzip archives. The genome of the Colorado potato beetle was used as the host organism. You can specify the genome of any organism you are interested in.

```
wget 
```

2) Creating databases for classification using BLAST

Download databases of nucleotide and amino acid sequences and unzip archives. Warning, the size of the database NCBI nt in compressed format is 179G, and the database NCBI nr is 128G. After unzipping, you will need about 500G of free disk space. You can use other smaller databases for your needs.

```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
unzip nt.gz nr.gz
```

Create BLAST databases. Warning, the BLAST program is required. Work in snakemake.yml environment.

```
makeblastdb -in nt -dbtype nucl -out [:path:]/nt
makeblastdb -in nr -dbtype prot -out [:path:]/nr
```

3) Creating databases for classification using Kraken2

As an example, the viral nucleotide sequences from NCBI RefSeq database are used. You can create a custom database for any nucleotide sequences in the format fasta.

```
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
unzip viral.?.1.genomic.fna.gz
```

Get the NCBI taxonomy files
```
kraken2-build --download-taxonomy --db viral_custom
```

Add custom reference data
```
kraken2-build --add-to-library viral.1.1.genomic.fna --db viral_custom
kraken2-build --add-to-library viral.2.1.genomic.fna --db viral_custom
kraken2-build --add-to-library viral.3.1.genomic.fna --db viral_custom
```
Finalize the database
```
kraken2-build --build --db $DBNAME
```

4) 
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz



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

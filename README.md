# VirMercenary
Pipeline for viruses detection 

Сейчас pipeline не может работать автоматически. Здесь приведена последовательность запуска скриптов из директории Scripts. Пути для запуска находятся в файле export_path.sh для horse. 

Предполагается, что данные секвенирования представлены в виде forward/reverse файлов в формате FASTQ. Для single прочтений и других форматов смотрите в директории Extra. 
Также предполагается, что проведен контроль качетсва образцов (FastQC) и данные предобработаны (удалены адаптеры, низкокачественные и короткие последовательности). 
Пример команды для предобработки данных с использованием Trimmomatic:
java -jar Trimmomatic-0.36/trimmomatic-0.36.jar PE ${x}_1.fastq ${x}_2.fastq ${x}_1.paired.fastq ${x}_1.unpaired.fastq ${x}_2.paired.fastq ${x}_2.unpaired.fastq 
\ILLUMINACLIP:/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true SLIDINGWINDOW:6:10 LEADING:13 TRAILING:13 MINLEN:36
Далее образцы следует передать repair.sh для соблюдния парности прочтений. 

До начала работы запустите: source export_path.sh

Последовательность запуска скриптов:
1. rename.sh 
2. Filtration_UNIVEC.sh
3. repair.sh
4. align_hg38.sh
5. filtration_hg38.sh
6. Assembly_spades.sh
7. Assembly_velvet.sh 

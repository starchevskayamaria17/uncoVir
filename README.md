# Pipeline for viruses detection in NGS data. 

## Installation and launch


```
conda env create snakemake.yml
```

## 2. Preparation of reference genomes and databases. Creating a configuration file. 

### 1) Reference genome indexing

Download human genome and host genome, and unzip archives. The genome of the Colorado potato beetle was used as the host organism. You can specify the genome of any organism you are interested in.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Leptinotarsa_decemlineata/latest_assembly_versions/GCA_000500325.2_Ldec_2.0/GCA_000500325.2_Ldec_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
unzip *.gz
```

Index genomes for subsequent alignment using the program bwa-mem
```
bwa index GCA_000500325.2_Ldec_2.0_genomic.fna
bwa index GCA_000001405.28_GRCh38.p13_genomic.fna
```

### 2) Creating databases for remove contamination and classification contigs using BLAST

Download databases of synthetic, nucleotide and amino acid sequences and unzip archives. Warning, the size of the database NCBI nt in compressed format is 179G, and the database NCBI nr is 128G. After unzipping, you will need about 500G of free disk space. You can use other smaller databases for your needs.

```
wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
unzip nt.gz nr.gz
```

Create BLAST databases. Warning, the BLAST program is required. Work in snakemake.yml environment.

```
makeblastdb -in UniVec -dbtype nucl -out UniVec
makeblastdb -in nt -dbtype nucl -out nt
makeblastdb -in nr -dbtype prot -out nr
```

### 3) Creating databases for classification using Kraken2

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

### 4) 
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

## 

As input files, the pipeline receives paired-end sequencing data in fastq-format. The work of the pipeline can be divided into the following stages (Fig. 1):
#### 1. Data preprocessing. 
The trimmomatic tool performs preliminary cleaning of data from adapters (the path to the file with adapters must be specified in the config file) and low-quality reads with standard parameters for most data: SLIDINGWINDOW:6:10 LEADING:13 TRAILING:13 MINLEN:50.
#### 2. Removal of contamination. 
Filtering against the human genome using the bwa-mem program with parameters: -B 10 -O 17,17 -E 17,17'. Filtering using the BLASTn program with specific parameters (-reward 1 -penalty -5 -gapopen 3 -gapextend 3  -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000) against a database of synthetic sequences NCBI UniVec.
#### 3. Selection of non-aligned reads for the host genome.
Filtering against the host genome using the bwa-mem program. Further, **pipeline works only with reads that are non-aligned to the host genome**.
#### 4. Assembly of the metagenome into contigs.
Non-aligned reads to the host genome are collected using SPADes with the option --meta. Contigs are filtered by length over 500 and coverage over 5.
#### 5. Classification of contigs. 
The contigs filtered by length and coverage are classified by the program BLASTn against the database NCBI nt. Unclassified contigs are selected and classified using BLASTx program against the database NCBI nr.
#### 6. Analysis of the spectrum of k-mers.
Non-aligned reads to the host genome are classified by the program Kraken2 against Kraken2 databases (see p.2.3). 
#### 7. Alignment of reads on virus genomes. 
Non-aligned reads to the host genome are aligned to viral genomes using a program Bowtie2. This method of classification makes it possible to estimate how evenly the viral genome is covered by reads and to take into account reads not collected into contigs.

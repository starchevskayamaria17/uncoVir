# Pipeline for viruses detection in NGS data 

The pipeline was designed to classify the viral sequences from raw NGS paired-end reads obtained on the Illumina platform. The main pipeline strategy is the assembly of contamination-free reads into contigs with their subsequent classification by BLAST tools. Additionally the pipeline makes the k-mer analysis with Kraken2 and aligns the reads against the database of viral genomes. The pipeline was created with Snakemake.

## 1. Dependence

Anaconda (version > 4.12)

Python (version > 3.9.9) 

Installation of additional programs is not required.

## 2. Setting the environment

After downloading the directory, you need to set up a working environment in which all the necessary software packages will be installed for the pipeline to work. To do this, run:

```
conda env create snakemake.yml 
conda activate snakemake
```
## 3. Pipeline launch

Before starting the pipeline, refer to paragraph 4 to create a configuration file and to build the databases. When everything is ready do: 

```
snakemake -s snakefile -j 2                                     #j - number of threads
```

## 4. Preparation of reference genomes and databases. Creating a configuration file. 

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

### 4) Creating databases for classification using alignments by Bowtie2

Merge all viral reference genomes into one file (viral.fasta) and create bowtie2 index database (database name: viral_bowtie2):
```
bowtie2-build viral.fasta viral_bowtie2
```
We recommend using a variety of virus databases with a uniform distribution of virus strains, removing duplicate sequences first.

### 5) Download taxonomy archives

To get taxonomy-marked BLAST and alignment reporting tables, you need archives that store NCBI Taxonomy. Please note that accession ID databases
sequences (especially amino acids) depend on which database you are using. Make sure accession ID match.

Download and unzip taxonomy archive:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```
Download and unzip archive accession ID for database of NCBI nt:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
```

Download and unzip archives accession ID and merge its for database of NCBI nr:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot_db.accession2taxid.gz
gunzip dead_prot.accession2taxid.gz pdb.accession2taxid.gz prot_db.accession2taxid.gz
cat dead_prot.accession2taxid pdb.accession2taxid prot.accession2taxid > prot_db.accession2taxid 
```

## 5. Description of the pipeline

As input files, the pipeline receives paired-end sequencing data in fastq-format. The work of the pipeline can be divided into the following stages (Fig. 1):
#### 1. Data preprocessing. 
The trimmomatic tool performs preliminary cleaning of data from adapters (the path to the file with adapters must be specified in the config file) and low-quality reads with standard parameters for most data: SLIDINGWINDOW:6:10 LEADING:13 TRAILING:13 MINLEN:50.
#### 2. Removal of contamination. 
Filtering against the human genome using the bwa-mem program with parameters: -B 10 -O 17,17 -E 17,17'. Filtering using the BLASTn program against a database of synthetic sequences NCBI UniVec with specific parameters: -reward 1 -penalty -5 -gapopen 3 -gapextend 3  -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000.
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

## 6. Results
#### 1. Statistics. 

| Name | Number of reads or contigs | Total length of reads or contigs, bp | Minimum length of reads or contigs | Average length of reads or contigs | Maximum length of reads or contigs | Q1 | Q2 | Q3 | N50 | Q20(%) | Q30(%) |
|---------------------------|----|----|----|----|----|----|----|----|----|----|----|
| raw forward/reverse read | 4629691 | 420444472 | 80 | 90.8 | 91 | 91.0 | 91.0 | 91.0 | 91 | 99.97 | 99.11 |
| clipped forward/reverse read | 4629691 | 419997824 | 76 | 90.7 | 91 | 91.0 | 91.0 | 91.0 | 91 | 99.97 | 99.11 |
| clipped unpair forward read | 0 | 0 | 0 | 0.0 | 0 | 0.0 | 0.0 | 0.0 | 0 | 0.0 | 0.0 |
| clipped unpair reverse read | 0 | 0 | 0 | 0.0 | 0 | 0.0 | 0.0 | 0.0 | 0 | 0.0 | 0.0 |
| free contamination forward/reverse read | 2001870 | 181674797 | 78 | 90.8 | 91 | 91.0 | 91.0 | 91.0 | 91 | 99.97 | 99.06 |
| non-aligned Ldec forward/reverse read | 40508 | 3675362 | 78 | 90.7 | 91 | 91.0 | 91.0 | 91.0 | 91 | 99.96 | 99.01 |
| raw contigs ||||||||||| 1211 | 422400 | 59 | 348.8 | 9215 | 235.0 | 283.0 | 379.0 | 344 | 0.0 | 0.0 |
| length and coverage contigs over 500 and 20 respectively ||||||||||| 47 | 48991 | 500 | 1042.4 | 9215 | 607.5 | 716.0 | 838.0 | 888 | 0.0 | 0.0 |
| contigs unclassified by the BLASTn | 21 | 17751 | 515 | 845.3 | 1654 | 712.0 | 745.0 | 921.0 | 776 | 0.0 | 0.0 |
| contigs classified as viral by the BLASTn | 2 | 15041 | 5826 | 7520.5 | 9215 | 5826.0 | 7520.5 | 9215.0 | 9215 | 0.0 | 0.0 |
| contigs unclassified by the BLASTx | 6 | 4255 | 615 | 709.2 | 783 | 659.0 | 723.0 | 752.0 | 745 | 0.0 | 0.0 |

Completeness of assembly of unfiltered contigs: 84.4623 


#### 2. Classification table based on the BLASTn (NCBI nt).
  
  
| qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore | accession.version | taxid | kindom | phylum | class | order | family | genus | species |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| NODE_1_length_9215_cov_18.424017 | EF026075.1 | 99.826 | 9210 | 16 | 0 | 1 | 9210 | 484 | 9693 | 0.0 | 16920.0 | EF026075.1 | 12216 | Viruses | Pisuviricota | Stelpaviricetes | Patatavirales | Potyviridae | Potyvirus | Potato virus Y |
| NODE_1_length_9215_cov_18.424017 | MN216357.1 | 99.815 | 9210 | 17 | 0 | 1 | 9210 | 439 | 9648 | 0.0 | 16914.0 | MN216357.1 | 12216 | Viruses | Pisuviricota | Stelpaviricetes | Patatavirales | Potyviridae | Potyvirus | Potato virus Y |
| NODE_2_length_5826_cov_22.717207 | AB364946.1 | 98.681 | 1289 | 17 | 0 | 4479 | 5767 | 1 | 1289 | 0.0 | 2287.0 | AB364946.1 | 12169 | Viruses | Kitrinoviricota | Alsuviricetes | Tymovirales | Betaflexiviridae | Carlavirus | Potato virus S |
| NODE_2_length_5826_cov_22.717207 | AY512653.1 | 96.082 | 1276 | 50 | 0 | 4551 | 5826 | 1 | 1276 | 0.0 | 2080.0 | AY512653.1 | 12169 | Viruses | Kitrinoviricota | Alsuviricetes | Tymovirales | Betaflexiviridae | Carlavirus | Potato virus S |
| NODE_2_length_5826_cov_22.717207 | S45593.1 | 94.973 | 1313 | 57 | 2 | 4523 | 5826 | 1 | 1313 | 0.0 | 2050.0 | S45593.1 | 12169 | Viruses | Kitrinoviricota | Alsuviricetes | Tymovirales | Betaflexiviridae | Carlavirus | Potato virus S |
| NODE_9_length_1079_cov_11.799805 | DQ631667.1 | 98.072 | 830 | 16 | 0 | 61 | 890 | 1 | 830 | 0.0 | 1445.0 | DQ631667.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_9_length_1079_cov_11.799805 | DQ631668.1 | 98.777 | 736 | 9 | 0 | 61 | 796 | 1 | 736 | 0.0 | 1310.0 | DQ631668.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_9_length_1079_cov_11.799805 | DQ631669.1 | 100.0 | 664 | 0 | 0 | 61 | 724 | 1 | 664 | 0.0 | 1227.0 | DQ631669.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_17_length_978_cov_12.575298 | HM175847.1 | 99.284 | 978 | 7 | 0 | 1 | 978 | 64 | 1041 | 0.0 | 1768.0 | HM175847.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_26_length_888_cov_13.318127 | XM_023170369.1 | 99.81 | 527 | 1 | 0 | 362 | 888 | 813 | 287 | 0.0 | 968.0 | XM_023170369.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_30_length_843_cov_23.440355 | XM_023172897.1 | 99.524 | 840 | 4 | 0 | 4 | 843 | 3601 | 2762 | 0.0 | 1530.0 | XM_023172897.1 | 7539 | Eukaryota | Arthropoda | Insecta | Coleoptera | Chrysomelidae | Leptinotarsa | Leptinotarsa decemlineata |
| NODE_32_length_833_cov_24.933162 | XM_011070431.1 | 71.935 | 367 | 87 | 12 | 209 | 567 | 2246 | 1888 | 2.49e-14 | 93.5 | XM_011070431.1 | 103372 | Eukaryota | Arthropoda | Insecta | Hymenoptera | Formicidae | Acromyrmex | Acromyrmex echinatior |
| NODE_32_length_833_cov_24.933162 | XM_011070430.1 | 71.935 | 367 | 87 | 12 | 209 | 567 | 2303 | 1945 | 2.49e-14 | 93.5 | XM_011070430.1 | 103372 | Eukaryota | Arthropoda | Insecta | Hymenoptera | Formicidae | Acromyrmex | Acromyrmex echinatior |
  
  
#### 3. Classification table based on the BLASTx (NCBI nr). 
#### 4. Classification table based on the BLASTx (NCBI nr). 

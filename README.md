# Pipeline for viruses detection in NGS data 

The pipeline was designed to classify the viral sequences from raw NGS paired-end reads obtained on the Illumina platform. The main pipeline strategy is the assembly of contamination-free reads into contigs with their subsequent classification by BLAST tools. Additionally the pipeline makes the k-mer analysis with Kraken2 and aligns the reads against the database of viral genomes. The pipeline was created with Snakemake.

## 1. Dependences

Anaconda (version > 4.12)

Python (version > 3.9.9) 

## 2. Setting the environment

After downloading the directory, you need to set up a working environment in which all the necessary software packages will be installed. To do this, run:

```
conda env create uncovir.yml 
conda activate uncovir
```
## 3. Pipeline launch

Before starting the pipeline, refer to paragraph 4 to create a configuration file and to build the databases. When everything is ready do: 

```
snakemake -s snakefile -j 2                                     #j - number of threads
```

## 4. Preparing the reference genomes and databases, creating the configuration file 

### 1) Reference genome indexing

Download human genome and host genome, and unzip archives. The genome of the Colorado potato beetle was used as the host organism. You can specify the genome of any organism you are interested in.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Leptinotarsa_decemlineata/latest_assembly_versions/GCA_000500325.2_Ldec_2.0/GCA_000500325.2_Ldec_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
unzip *.gz
```

Index the genomes for subsequent alignment using the program Bowtie2
```
bowtie2-build GCA_000500325.2_Ldec_2.0_genomic.fna GCF_000500325.1_Ldec_2.0_genomic
bowtie2-build GCA_000001405.28_GRCh38.p13_genomic.fna GCF_000001405.38_GRCh38.p13_genomic
```

### 2) Creating the databases for contamination removal and contigs classification with BLAST tools

Download the databases of synthetic, nucleotide and amino acid sequences and unzip the archives. **NB!** The size of the NCBI nt database in compressed format is 179G, and the NCBI nr database &mdash; 128G. After unzipping, you will need about 500G of free disk space. You can use other smaller databases for your needs.

```
wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
unzip nt.gz nr.gz
```

Create BLAST databases. Warning, the BLAST program is required. **NB!** Work in *uncovir* environment.

```
makeblastdb -in UniVec -dbtype nucl -out UniVec
makeblastdb -in nt -dbtype nucl -out nt
makeblastdb -in nr -dbtype prot -out nr
```

### 3) Creating the databases for classification with Kraken2

As an example, the viral nucleotide sequences from NCBI RefSeq database are used. You can create a custom database for any nucleotide sequences in fasta format.

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

### 4) Creating the databases for classification with Bowtie2 mapping

Merge all viral reference genomes into one file (viral.fasta) and create bowtie2 index database (database name: viral_bowtie2):
```
bowtie2-build viral.fasta viral_bowtie2
```
We recommend using a variety of virus databases with a uniform distribution of virus strains, removing duplicate sequences first.

### 5) Download the taxonomy archives

To get the taxonomy-marked BLAST and alignment reporting tables, you'll need the NCBI Taxonomy archives. Please note that sequences accession IDs (and especially for amino acid sequences) strongly depend on which database you are using. Make sure the accession IDs match between the taxonomy tables and your sequence databases.

Download and unzip taxonomy archive:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```
Download and unzip accession IDs for NCBI nt database:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
```

Download and unzip the archives with amino acid sequences accession IDs and merge them into prot_db.accession2taxid:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot_db.accession2taxid.gz
gunzip dead_prot.accession2taxid.gz pdb.accession2taxid.gz prot_db.accession2taxid.gz
cat dead_prot.accession2taxid pdb.accession2taxid prot.accession2taxid > prot_db.accession2taxid 
```

## 5. Description of the pipeline

As input files, the pipeline receives the paired-end sequencing data in fastq format. The pipeline can be divided into the following stages:
#### 1. Data preprocessing. 
The *trimmomatic* tool performs preliminary cleaning of the data from adapters (the path to the file with adapters must be specified in the config file) and low-quality reads with standard parameters: SLIDINGWINDOW:6:10 LEADING:13 TRAILING:13 MINLEN:50.
#### 2. Contamination removal. 
Filtering against the human genome using the *Bowtie2* program. Filtering using the *BLASTn* program against a database of synthetic sequences NCBI UniVec with specific parameters: -reward 1 -penalty -5 -gapopen 3 -gapextend 3  -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000.
#### 3. Selection of non-aligned reads for the host genome.
Filtering against the host genome using the *Bowtie2* program. Further, **pipeline works only with reads that are non-aligned to the host genome**.
#### 4. Contigs assembly.
The reads non-aligned to the host genome are assembled using *SPADes* with the option --meta. Be default the contigs are filtered by length over 500 and coverage over 5.
#### 5. Contigs classification. 
The remaining contigs (selected by length and coverage) are classified by *BLASTn* against the NCBI nt database. The unclassified contigs are subsequently classified with *BLASTx* against the NCBI nr database.
#### 6. K-mer spectrum analysis.
The reads non-aligned to the host genome are classified with *Kraken2* against Kraken2 databases (see p.2.3). 
#### 7. Alignment of reads on virus genomes. 
The reads non-aligned to the host genome are aligned to viral genomes with *Bowtie2*. This allows to estimate how evenly the viral genome is covered by the matched reads and to take into account the reads that were not assembled into contigs.

![alt "Pipeline scheme"](https://github.com/starchevskayamaria17/uncoVir/blob/master/misc/pipeline.png?raw=true)


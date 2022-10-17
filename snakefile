configfile: "config.yml"
SAMPLES, = glob_wildcards("samples/{sample}_1.fastq")
READS=["1", "2"]

rule all: 
	input: 
		expand(config['tmp'] + '{sample}_1.fastq', sample=SAMPLES), 
		expand(config['tmp'] + '{sample}_2.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_1.paired.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_1.unpaired.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.paired.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.unpaired.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.unpaired.fastq', sample=SAMPLES),
		expand(config['tmp'] + 'align_hg38_{sample}.sam', sample=SAMPLES),
		expand(config['tmp'] + 'align_hg38_{sample}_unmapped.bam', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_1.Nohg38.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.Nohg38.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_1.Nohg38.fasta', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.Nohg38.fasta', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_f_blast_univec_result.csv', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_r_blast_univec_result.csv', sample=SAMPLES),		
		expand(config['tmp'] + '{sample}_1.NoUnivec.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_2.NoUnivec.fastq', sample=SAMPLES),
		expand(config['results'] + '{sample}_1.filtered.fastq', sample=SAMPLES),
		expand(config['results'] + '{sample}_2.filtered.fastq', sample=SAMPLES),
		expand(config['tmp'] + '{sample}_1.filtered.unpaired.fastq', sample=SAMPLES),         
		expand(config['tmp'] + 'align_Ldec_{sample}.sam', sample=SAMPLES),
		expand(config['tmp'] + 'align_Ldec_{sample}_unmapped.bam', sample=SAMPLES),
		expand(config['results'] + '{sample}_1.NoLdec.fastq', sample=SAMPLES),
		expand(config['results'] + '{sample}_2.NoLdec.fastq', sample=SAMPLES),
		expand(config['results'] + '{sample}_assembly/contigs.fasta', sample=SAMPLES),
		expand(config['results'] + '{sample}_assembly', sample=SAMPLES),
		expand(config['results'] + '{sample}_contigs_5_more_coverage_500_more_len.fasta', sample=SAMPLES),
		expand(config['results'] + '{sample}_blastn_result.csv', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.taxonomy.blastn.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_contigs_no_blastn.fasta', sample=SAMPLES),
		expand(config['results'] + '{sample}_contigs_viral_blastn.fasta', sample=SAMPLES),
		expand(config['results'] + '{sample}_blastx_result.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_contigs_no_blastx.fasta', sample=SAMPLES),
		expand(config['results'] + '{sample}_kraken_result.txt', sample=SAMPLES),
		expand(config['results'] + '{sample}.stats.csv', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.bowtie.sam', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.bowtie.mapped.bam', sample=SAMPLES), 
		expand(config['tmp'] + '{sample}.bowtie.mapped.sam', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.bowtie.mapped.NoHeader.sam', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.accession.txt', sample=SAMPLES),
		expand(config['tmp'] + '{sample}.taxonomy.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_taxonomy_blastn.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_taxonomy_blastx.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_bowtie2_stats.txt', sample=SAMPLES), 
		expand(config['tmp'] + '{sample}.bowtie.coverage.trim.csv', sample=SAMPLES),
		expand(config['results'] + '{sample}_coverage_taxonomy.csv', sample=SAMPLES)

rule rename:
	input:
		forward_reads = config['samples'] + '{sample}_1.fastq', 
		reverse_reads = config['samples'] + '{sample}_2.fastq'
	output:
		rename_reads_f = config['tmp'] + '{sample}_1.fastq', 
		rename_reads_r = config['tmp'] + '{sample}_2.fastq'
	run:
		shell(""" awk  '{{print (NR%4 == 1) ? "@" ++i "/1" : $0}}' {input.forward_reads} > {output.rename_reads_f} """)
		shell("""awk  '{{print (NR%4 == 1) ? "@" ++i "/2" : $0}}' {input.reverse_reads} > {output.rename_reads_r}""")
		
rule trimmomatic:
	input: 
		forward_reads = config['tmp'] + '{sample}_1.fastq', 
		reverse_reads = config['tmp'] + '{sample}_2.fastq',
		adapters = config['adapters']
		
	output: 
		f_trim_reads = config['tmp'] + '{sample}_1.paired.fastq', 
		f_trim_un_reads = config['tmp'] + '{sample}_1.unpaired.fastq', 
		r_trim_reads = config['tmp'] + '{sample}_2.paired.fastq', 
		r_trim_un_reads = config['tmp'] + '{sample}_2.unpaired.fastq'
	shell: 
		"trimmomatic PE {input.forward_reads} {input.reverse_reads} {output.f_trim_reads} {output.f_trim_un_reads} {output.r_trim_reads} {output.r_trim_un_reads} ILLUMINACLIP:{input.adapters}:2:30:10:1:true SLIDINGWINDOW:6:10 LEADING:13 TRAILING:13 MINLEN:50"	

### Align to the human genome to remove contamination in the future
rule bwa_hg38:
	input: 
		reference = config['ref_hg38'], 
		forward_reads = config['tmp'] + '{sample}_1.paired.fastq', 
		reverse_reads = config['tmp'] + '{sample}_2.paired.fastq'
	output: 
		config['tmp'] + 'align_hg38_{sample}.sam'
	params:
		threads = config['threads_bwa'],
		name_db = config['ref_hg38_name'] 
	shell:
		"cd {input.reference} && bowtie2 -x {params.name_db} -1 {input.forward_reads} -2 {input.reverse_reads}  -p {params.threads} -S {output}"

rule bwa_unmapped_hg38_samtools:
	input: 
		config['tmp'] + 'align_hg38_{sample}.sam'
	output: 
		config['tmp'] + 'align_hg38_{sample}_unmapped.bam'
	params:
		config['threads_samtools']
	shell: 
		"samtools view -@ {params} -bS {input} | samtools view -@ {params} -b -F 2 | samtools sort -n > {output}"

rule filtration_unmapped_hg38_reads:		
	input:
		config['tmp'] + 'align_hg38_{sample}_unmapped.bam'
	output:
		forward_reads = config['tmp'] + '{sample}_1.Nohg38.fastq',
		reverse_reads = config['tmp'] + '{sample}_2.Nohg38.fastq'
	shell:
		"bedtools bamtofastq -i {input} -fq {output.forward_reads} -fq2 {output.reverse_reads}"

### Search syntetic sequences in UniVec database 
rule fq2fa:
		input:
			forward_reads = config['tmp'] + '{sample}_1.Nohg38.fastq',
			reverse_reads = config['tmp'] + '{sample}_2.Nohg38.fastq'
		output:
			forward_reads_fa = config['tmp'] + '{sample}_1.Nohg38.fasta',
			reverse_reads_fa = config['tmp'] + '{sample}_2.Nohg38.fasta'
		run:
			shell("""seqtk seq -a {input.forward_reads} > {output.forward_reads_fa}""")
			shell("""seqtk seq -a {input.reverse_reads} > {output.reverse_reads_fa}""")
	
rule blastn_univec:
	input:
		forward_reads_fa = config['tmp'] + '{sample}_1.Nohg38.fasta',
		reverse_reads_fa = config['tmp'] + '{sample}_2.Nohg38.fasta',
		univec_db = config['db_univec_path']
	output: 
		report_f = config['tmp'] + '{sample}_f_blast_univec_result.csv',
		report_r = config['tmp'] + '{sample}_r_blast_univec_result.csv'
	params:
		threads = config['threads_blast'],
		name = config['db_univec_name']	
	run: 
		shell("""cd {input.univec_db} && blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3  -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -num_threads {params.threads} \
		-db {params.name} -query {input.forward_reads_fa}  -outfmt 6 -out {output.report_f}""")
		shell("""cd {input.univec_db} && blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3  -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -num_threads {params.threads} \
		-db {params.name} -query {input.reverse_reads_fa}  -outfmt 6 -out {output.report_r}""")

rule fastq_no_univec:
	input:
		report_f = config['tmp'] + '{sample}_f_blast_univec_result.csv',
		forward_reads_fa = config['tmp'] + '{sample}_1.Nohg38.fasta',
		forward_reads = config['tmp'] + '{sample}_1.Nohg38.fastq',
		report_r = config['tmp'] + '{sample}_r_blast_univec_result.csv',
		reverse_reads_fa = config['tmp'] + '{sample}_2.Nohg38.fasta',
		reverse_reads = config['tmp'] + '{sample}_2.Nohg38.fastq'
	output:
		id_fastq_f = config['tmp'] + 'id_{sample}_1.Nohg38.txt',
		id_univec_f = config['tmp'] + 'id_univec_{sample}_1.Nohg38.txt',
		id_no_univec_f = config['tmp'] + 'id_no_univec_{sample}_1.Nohg38.txt',
		forward_no_univec = config['tmp'] + '{sample}_1.NoUnivec.fastq',
		id_fastq_r = config['tmp'] + 'id_{sample}_2.Nohg38.txt',
		id_univec_r = config['tmp'] + 'id_univec_{sample}_2.Nohg38.txt',
		id_no_univec_r = config['tmp'] + 'id_no_univec_{sample}_2.Nohg38.txt',
		reverse_no_univec = config['tmp'] + '{sample}_2.NoUnivec.fastq'
	run:
		shell (""" cut -f 1 {input.report_f} | sort | uniq | sort > {output.id_univec_f}; \
		grep '^>' {input.forward_reads_fa} | cut -d'>' -f 2- | sort > {output.id_fastq_f}; \
		comm -13 {output.id_univec_f} {output.id_fastq_f} > {output.id_no_univec_f}; \ 
		seqtk subseq {input.forward_reads} {output.id_no_univec_f} > {output.forward_no_univec}; \
		cut -f 1 {input.report_r} | sort | uniq | sort > {output.id_univec_r}; \
		grep '^>' {input.reverse_reads_fa} | cut -d'>' -f 2- | sort > {output.id_fastq_r}; \
		comm -13 {output.id_univec_r} {output.id_fastq_r} > {output.id_no_univec_r}; \ 
		seqtk subseq {input.reverse_reads} {output.id_no_univec_r} > {output.reverse_no_univec} """)

rule repair: 
	input:
		forward_no_univec = config['tmp'] + '{sample}_1.NoUnivec.fastq',
		reverse_no_univec = config['tmp'] + '{sample}_2.NoUnivec.fastq'
	output:
		forward_reads = config['results'] + '{sample}_1.filtered.fastq',
		reverse_reads = config['results'] + '{sample}_2.filtered.fastq',
		unpaired_reads = config['tmp'] + '{sample}_1.filtered.unpaired.fastq'
	shell:
		"repair.sh in={input.forward_no_univec} in2={input.reverse_no_univec} out={output.forward_reads} out2={output.reverse_reads} outs={output.unpaired_reads}"
	
rule bwa_Ldec:
		input: 
			reference = config['ref_ldec'],
			forward_reads = config['results'] + '{sample}_1.filtered.fastq',
			reverse_reads = config['results'] + '{sample}_2.filtered.fastq'
		output: 
			config['tmp'] + 'align_Ldec_{sample}.sam'
		params:
			t = config['threads_bowtie2'],
			name = config['ref_ldec_name'] 
		shell:
			"cd {input.reference} && bowtie2 -x {params.name} -1 {input.forward_reads} -2 {input.reverse_reads} -p {params.t} -S {output}"

rule bwa_unmapped_Ldec_samtools:
	input: 
		config['tmp'] + 'align_Ldec_{sample}.sam'
	output: 
		config['tmp'] + 'align_Ldec_{sample}_unmapped.bam'
	params:
		config['threads_samtools']
	shell: 
		"samtools view -@ {params} -bS {input} | samtools view -@ {params} -b -f 4  | samtools sort -n  > {output}"


rule filtration_unmapped_Ldec_reads:		
	input:
		config['tmp'] + 'align_Ldec_{sample}_unmapped.bam'
	output:
		forward_reads = config['results'] + '{sample}_1.NoLdec.fastq',
		reverse_reads = config['results'] + '{sample}_2.NoLdec.fastq'
	shell:
		"bedtools bamtofastq -i {input} -fq {output.forward_reads} -fq2 {output.reverse_reads}"

rule spades:
	input: 
		forward_reads = config['results'] + '{sample}_1.NoLdec.fastq',
		reverse_reads = config['results'] + '{sample}_2.NoLdec.fastq'
	output:
		dir_assembly = directory(config['results'] + '{sample}_assembly'),
		contigs = config['results'] + '{sample}_assembly/contigs.fasta'
	params:
		threads = config['threads_spades'],
		memory = config['memory_spades']
	shell:
		"spades.py --meta -1 {input.forward_reads} -2 {input.reverse_reads} -t {params.threads} -m {params.memory} -o {output.dir_assembly}"


rule selection_of_contigs:
	input:
		contigs = config['results'] + '{sample}_assembly/contigs.fasta'
	output:
		id_bad_coverage = config['tmp'] + '{sample}_id_bad_coverage.txt',
		id_contigs = config['tmp'] + '{sample}_id_contigs.txt',
		id_5_more_coverage = config['tmp'] + '{sample}_id_5_more_coverage.txt',
		contigs_5_more_coverage = config['tmp'] + '{sample}_contigs_5_more_coverage.fasta',
		contigs_5_more_coverage_500_more_len = config['results'] + '{sample}_contigs_5_more_coverage_500_more_len.fasta'
	shell: 
		"""
		grep -E 'cov_0\.|cov_1\.|cov_2\.|cov_3\.|cov_4\.' {input.contigs} | cut -c 2- | sort > {output.id_bad_coverage};
		grep '>' {input.contigs} | cut -c 2- | sort > {output.id_contigs};
		comm -13 {output.id_bad_coverage} {output.id_contigs} > {output.id_5_more_coverage};
		seqtk subseq {input.contigs} {output.id_5_more_coverage} > {output.contigs_5_more_coverage};
		seqkit seq -m 500 {output.contigs_5_more_coverage} > {output.contigs_5_more_coverage_500_more_len}
		"""
		
rule blastn:
	input:  
		contigs = config['results'] + '{sample}_contigs_5_more_coverage_500_more_len.fasta',
		nt_db = config['db_nt']		
	output: 
		report = config['results'] + '{sample}_blastn_result.csv'
	params:
		threads = config['threads_blast'],
		name = config['db_nt_name']
	shell: 
		"cd {input.nt_db} &&  blastn -query {input.contigs} -db {params.name} -num_threads {params.threads} -outfmt 6 -out {output.report}"


rule taxonkit_blastn: 
	input:
		report = config['results'] + '{sample}_blastn_result.csv',
		db_taxid = config['db_taxid'],
		db_accession = config['db_nucl_accession']
	output:
		accession_uniq = config['tmp'] + '{sample}.accession.blastn.uniq.txt',
		acc2taxid = config['tmp'] + '{sample}.acc2taxid.blastn.txt',
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.blastn.csv'
	shell:
		"""
		cut -f 2 {input.report}  | cut -d. -f 1 | sort | uniq > {output.accession_uniq};
		csvtk grep -t -f accession {input.db_accession} -P {output.accession_uniq}  | csvtk cut -t -f accession.version,taxid | sed 1d > {output.acc2taxid};    
		taxonkit lineage --data-dir {input.db_taxid} -i 2 {output.acc2taxid} | taxonkit reformat --data-dir {input.db_taxid} -f '{{k}}|{{p}}|{{c}}|{{o}}|{{f}}|{{g}}|{{s}}' -t --fill-miss-rank -i 3 | csvtk cut -t -f 1,2,4 | csvtk -H -t sep -f 3 -s '|' -R | csvtk add-header -t -n accession.version,taxid,kindom,phylum,class,order,family,genus,species > {output.taxonomy_csv} 
		"""

rule taxonomy_blastn:
	input:
		report = config['results'] + '{sample}_blastn_result.csv', 
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.blastn.csv',
		blast_taxonomy = config['blast_taxonomy']
	output:
		blastn_taxonomy = config['results'] + '{sample}_taxonomy_blastn.csv'
	shell:
		"python {input.blast_taxonomy} {input.report} {input.taxonomy_csv} {output.blastn_taxonomy}"


rule no_blastn:
	input:
		contigs = config['results'] + '{sample}_contigs_5_more_coverage_500_more_len.fasta',
		report = config['results'] + '{sample}_blastn_result.csv',
		blastn_taxonomy = config['results'] + '{sample}_taxonomy_blastn.csv'
	output:
		id_good_contigs = config['tmp'] + '{sample}_id_good_contigs.txt',
		id_blastn = config['tmp'] + '{sample}_id_blastn.txt',
		id_no_blastn = config['tmp'] + '{sample}_id_no_blastn.txt',
		contigs_no_blastn = config['results'] + '{sample}_contigs_no_blastn.fasta',
		id_viral_blastn = config['tmp'] + '{sample}_id_viral_blastn.txt',
		contigs_viral_blastn = config['results'] + '{sample}_contigs_viral_blastn.fasta',
		contigs_no_blastn_and_viral = config['results'] + '{sample}_contigs_no_blastn_and_viral.fasta'
	shell:
		"""
		grep '>' {input.contigs} | cut -d'>' -f 2 | sort > {output.id_good_contigs}	
		cut -f 1 {input.report} | sort | uniq | sort > {output.id_blastn};
		comm -13 {output.id_blastn} {output.id_good_contigs} > {output.id_no_blastn};
		seqtk subseq {input.contigs} {output.id_no_blastn} > {output.contigs_no_blastn};
		grep -i 'Virus' {input.blastn_taxonomy} |  cut -f 2 | sort | uniq > {output.id_viral_blastn};
		seqtk subseq {input.contigs} {output.id_viral_blastn} > {output.contigs_viral_blastn};
		cat {output.contigs_viral_blastn} {output.contigs_no_blastn} > {output.contigs_no_blastn_and_viral}
		"""

rule blastx:
	input:  
		contigs = config['results'] + '{sample}_contigs_no_blastn_and_viral.fasta',
		nr_db = config['db_nr']
	output: 
		report = config['results'] + '{sample}_blastx_result.csv'
	params:
		threads = config['threads_blast'],
		name = config['db_nr_name']
	shell: 
		"cd {input.nr_db} &&  blastx -query {input.contigs} -db {params.name} -num_threads {params.threads} -outfmt 6 -out {output.report}"			

rule no_blastx:
	input:
		contigs = config['results'] + '{sample}_contigs_no_blastn_and_viral.fasta',
		report = config['results'] + '{sample}_blastx_result.csv'
	output:
		id_good_contigs = config['tmp'] + '{sample}_id_good_contigs.txt',
		id_blastx = config['tmp'] + '{sample}_id_blastx.txt',
		id_no_blastx = config['tmp'] + '{sample}_id_no_blastx.txt',
		contigs_no_blastx = config['results'] + '{sample}_contigs_no_blastx.fasta'
	shell:
		"""
		grep '>' {input.contigs} | cut -d'>' -f 2 | sort > {output.id_good_contigs}	
		cut -f 1 {input.report} | sort | uniq | sort > {output.id_blastx};
		comm -13 {output.id_blastx} {output.id_good_contigs} > {output.id_no_blastx};
		seqtk subseq {input.contigs} {output.id_no_blastx} > {output.contigs_no_blastx}
		"""

rule taxonkit_blastx: 
	input:
		report = config['results'] + '{sample}_blastx_result.csv',
		db_taxid = config['db_taxid'],
		db_accession = config['db_prot_accession']
	output:
		accession_uniq = config['tmp'] + '{sample}.accession.blastx.uniq.txt',
		acc2taxid = config['tmp'] + '{sample}.acc2taxid.blastx.txt',
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.blastx.csv'
	shell:
		"""
		cut -f 2 {input.report}  |  sort | uniq > {output.accession_uniq};
		csvtk grep -t -f accession {input.db_accession} -P {output.accession_uniq}  | csvtk cut -t -f accession.version,taxid | sed 1d > {output.acc2taxid};    
		taxonkit lineage --data-dir {input.db_taxid} -i 2 {output.acc2taxid} | taxonkit reformat --data-dir {input.db_taxid} -f '{{k}}|{{p}}|{{c}}|{{o}}|{{f}}|{{g}}|{{s}}' -t --fill-miss-rank -i 3 | csvtk cut -t -f 1,2,4 | csvtk -H -t sep -f 3 -s '|' -R | csvtk add-header -t -n accession.version,taxid,kindom,phylum,class,order,family,genus,species > {output.taxonomy_csv} 
		"""

rule taxonomy_blastx:
	input:
		report = config['results'] + '{sample}_blastx_result.csv', 
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.blastx.csv',
		blast_taxonomy = config['blast_taxonomy']
	output:
		blastx_taxonomy = config['results'] + '{sample}_taxonomy_blastx.csv'
	shell:
		"python {input.blast_taxonomy} {input.report} {input.taxonomy_csv} {output.blastx_taxonomy}"


rule reads_stat:
	input: 
		raw_reads = config['samples'] + '{sample}_1.fastq',
		trimmomatic_pair = config['tmp'] + '{sample}_1.paired.fastq',
		trimmomatic_unpair_forward = config['tmp'] + '{sample}_1.unpaired.fastq',
		trimmomatic_unpair_reverse = config['tmp'] + '{sample}_2.unpaired.fastq',
		without_contamination_reads = config['results'] + '{sample}_1.filtered.fastq',
		nonalign_host_reads_1 = config['results'] + '{sample}_1.NoLdec.fastq',
		nonalign_host_reads_2 = config['results'] + '{sample}_2.NoLdec.fastq',
		contigs_raw = config['results'] + '{sample}_assembly/contigs.fasta',
		contigs_filtr = config['results'] + '{sample}_contigs_5_more_coverage_500_more_len.fasta',
		contig_no_blastn = config['results'] + '{sample}_contigs_no_blastn.fasta',
		contigs_viral_blastn = config['results'] + '{sample}_contigs_viral_blastn.fasta',
		contig_no_blastx = config['results'] + '{sample}_contigs_no_blastx.fasta',
		rename_script = config['rename']
	output:
		list_file = config['tmp'] + '{sample}.list.txt',
		alignment = config['tmp'] + '{sample}_reads_on_contigs.bam',
		stats = config['tmp'] + '{sample}.stats.csv',
		stats_rename = config['results'] + '{sample}.stats.csv'
	params:
		threads_seqkit = config['seqkit'],
		threads_bwa = config['threads_bwa'],
		threads_samtools = config['threads_samtools'] 
	shell:
		"""
		echo "{input.raw_reads}" > {output.list_file};
		echo "{input.trimmomatic_pair}" >> {output.list_file};
		echo "{input.trimmomatic_unpair_forward}" >> {output.list_file};
		echo "{input.trimmomatic_unpair_reverse}" >> {output.list_file};
		echo "{input.without_contamination_reads}" >> {output.list_file};
		echo "{input.nonalign_host_reads_1}" >> {output.list_file};
		echo "{input.contigs_raw}" >> {output.list_file};
		echo "{input.contigs_filtr}" >> {output.list_file};
		echo "{input.contig_no_blastn}" >> {output.list_file};
		echo "{input.contigs_viral_blastn}" >> {output.list_file};
		echo "{input.contig_no_blastx}" >> {output.list_file};
		seqkit stat -a -e -j {params.threads_seqkit} -T --infile-list {output.list_file} > {output.stats};
		python {input.rename_script} {output.stats} {output.stats_rename};
		bwa index {input.contigs_raw};
		bwa mem -t {params.threads_bwa} {input.contigs_raw} {input.nonalign_host_reads_1} {input.nonalign_host_reads_2} | samtools sort -@ {params.threads_samtools} > {output.alignment};
		echo "Completeness of assembly of unfiltered contigs:" >> {output.stats_rename};
		samtools stats {output.alignment} | grep 'raw total sequences:\|reads mapped:' | cut -f 3 | tr '\n' ' ' | awk '{{print $2/$1*100}}' >> {output.stats_rename}  
		"""

rule kraken2:
	input:
		forward_reads = config['results'] + '{sample}_1.filtered.fastq',
		reverse_reads = config['results'] + '{sample}_2.filtered.fastq',
		db_kraken = config['db_kraken']
	output:
		report = config['results'] + '{sample}_kraken_result.txt',
		out = config['tmp'] + '{sample}_output.kraken.txt'
	params:
		threads = config['threads_kraken']
	shell:
		"kraken2 --paired {input.forward_reads} {input.reverse_reads} --db {input.db_kraken} --report {output.report} --output {output.out} --threads {params.threads}"
	


rule bowtie2:
	input: 
		reference = config['db_bowtie2'], 
		forward_reads = config['results'] + '{sample}_1.NoLdec.fastq', 
		reverse_reads = config['results'] + '{sample}_2.NoLdec.fastq'
	output: 
		config['tmp'] + '{sample}.bowtie.sam'
	params:
		config['threads_bowtie2']
	log:
		bowtie_log = config['results'] + '{sample}_bowtie2_stats.txt'
	shell:
		"""
		bowtie2 -x {input.reference} -1 {input.forward_reads} -2 {input.reverse_reads} -a --no-unal -p {params} -S {output};
		2>{log.bowtie_log}
		"""


rule mapped_reads:
	input: 
		config['tmp'] + '{sample}.bowtie.sam'
	output: 
		mapped_bam = config['tmp'] + '{sample}.bowtie.mapped.bam',
		mapped_sam = config['tmp'] + '{sample}.bowtie.mapped.sam',
		mapped_sam_noheader = config['tmp'] + '{sample}.bowtie.mapped.NoHeader.sam',
		mapped_accession = config['tmp'] + '{sample}.accession.txt'
	params:
		config['threads_samtools']
	shell: 
		"""
		samtools sort -@ {params} {input} | samtools view -@ {params} -b -F 4 > {output.mapped_bam};
		samtools view -@ {params} -h -o {output.mapped_sam} {output.mapped_bam};
		grep -v '^@' {output.mapped_sam} > {output.mapped_sam_noheader};
		cut -f 3 {output.mapped_sam_noheader} > {output.mapped_accession} 
		"""

rule taxonkit: 
	input:
		mapped_accession = config['tmp'] + '{sample}.accession.txt',
		db_taxid = config['db_taxid'],
		db_accession = config['db_nucl_accession']
	output:
		accession_normal = config['tmp'] + '{sample}.accession.normal.txt',
		accession_uniq = config['tmp'] + '{sample}.accession.uniq.txt',
		acc2taxid = config['tmp'] + '{sample}.acc2taxid.txt',
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.csv'
	shell:
		"""
		cut -d' ' -f 1 {input.mapped_accession}  | cut -d: -f 1 | cut -d. -f 1 > {output.accession_normal}; 
		cut -d' ' -f 1 {input.mapped_accession}  | cut -d: -f 1 | cut -d. -f 1 | sort | uniq > {output.accession_uniq};
		csvtk grep -t -f accession {input.db_accession} -P {output.accession_uniq}  | csvtk cut -t -f accession,taxid | sed 1d > {output.acc2taxid};    
		taxonkit lineage --data-dir {input.db_taxid} -i 2 {output.acc2taxid} | taxonkit reformat --data-dir {input.db_taxid} -f '{{k}}|{{p}}|{{c}}|{{o}}|{{f}}|{{g}}|{{s}}' -t --fill-miss-rank -i 3 | csvtk cut -t -f 1,2,4 | csvtk -H -t sep -f 3 -s '|' -R | csvtk add-header -t -n accession,taxid,kindom,phylum,class,order,family,genus,species > {output.taxonomy_csv} 
		"""

rule classification_alignment:
	input: 
		mapped_bam = config['tmp'] + '{sample}.bowtie.mapped.bam',
		alignment_taxonomy = config['alignment_taxonomy'],
		taxonomy_csv = config['tmp'] + '{sample}.taxonomy.csv'
	output:
		coverage = config['tmp'] + '{sample}.bowtie.coverage.csv',
		coverage_trim = config['tmp'] + '{sample}.bowtie.coverage.trim.csv',
		coverage_taxonomy = config['results'] + '{sample}_coverage_taxonomy.csv' 
	shell:
		"""
		samtools coverage  {input.mapped_bam} > {output.coverage};
		awk 'BEGIN{{FS=OFS="\t"}} {{sub(/\..*/,"",$1)}} 1' {output.coverage} | awk 'BEGIN{{FS=OFS="\t"}} {{sub(/:.*/,"",$1)}} 1' > {output.coverage_trim};
		python {input.alignment_taxonomy} {output.coverage_trim} {input.taxonomy_csv} {output.coverage_taxonomy}
		"""

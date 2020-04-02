#!/bin/bash

#!/bin/bash

PATH=/samtools-1.7:$PATH
PATH=/seqfilter:$PATH
PATH=/sratoolkit.2.9.6-1-centos_linux64/bin:$PATH
PATH=/SPAdes-3.13.1-Linux/bin:$PATH
PATH=/ncbi-blast-2.9.0+/bin:$PATH
PATH=/Trimmomatic-0.36:$PATH
PATH=/FastQC:$PATH
PATH=/velvet_1.2.10:$PATH
PATH=/SPAdes-3.13.1-Linux:$PATH
PATH=/quast-4.6.3:$PATH
PATH=/ncbi-blast-2.9.0+:$PATH
PATH=/bwa-0.7.15:$PATH
PATH=/bbmap:$PATH
exportPATH
echo$PATH

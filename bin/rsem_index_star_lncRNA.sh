#! /bin/sh
genome_fa=/data/public_data/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa

rsem-prepare-reference -p 25 --star --gtf lncRNA/lncRNA.final.v2.gtf  $genome_fa lncRNA/lncRNA_RSEM_index/lncRNA.final.v2.RSEM.star

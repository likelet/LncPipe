#! /bin/sh
genome_fa=/data/public_data/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa

rsem-prepare-reference -p 25 --star --gtf lncRNA/protein_coding.final.gtf  $genome_fa lncRNA/pc_RSEM_index/protein_coding.final.RSEM.star

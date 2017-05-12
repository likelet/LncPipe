#! /bin/sh
annotation_dir=/HOME/nsfc2015_555/NSFC/WORKSPACE/database/GRCH38/Homo_sapiens/UCSC/hg38/Annotation/Genes
genome_fa=/HOME/nsfc2015_555/NSFC/WORKSPACE/database/GRCH38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
name=/HOME/nsfc2015_555/NSFC/WORKSPACE/database/GRCH38/Homo_sapiens/UCSC/hg38/Sequence/RSEMIndex/star_knownGene
rsem-prepare-reference -p 25 --star --gtf $annotation_dir/knowngene.gtf --transcript-to-gene-map $annotation_dir/knownIsoforms_symbol.txt $genome_fa $name

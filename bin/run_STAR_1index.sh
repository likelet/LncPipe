#!/bin/sh
cat /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/*_1pass/SJ.out.tab > /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/all.SJ.out.tab
mkdir -pv /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index
cd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index
STAR --runThreadN 25 --runMode genomeGenerate --genomeDir /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index  --genomeFastaFiles /HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf --sjdbFileChrStartEnd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/all.SJ.out.tab --sjdbOverhang 149 --limitSjdbInsertNsj 1600000

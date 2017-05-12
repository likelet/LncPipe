#!/usr/bin/perl -w
use strict;

open FH,">run_STAR_1index.sh" or die;

print FH "#!/bin/sh\n";
print FH "cat /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/*_1pass/SJ.out.tab > /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/all.SJ.out.tab\n";
print FH "mkdir -pv /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index\n";
print FH "cd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index\n";
print FH "STAR --runThreadN 25 --runMode genomeGenerate --genomeDir /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index  --genomeFastaFiles /HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf --sjdbFileChrStartEnd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/all.SJ.out.tab --sjdbOverhang 149 --limitSjdbInsertNsj 1600000\n";

`chmod +x run_STAR_1index.sh`;
`yhbatch -p nsfc3 run_STAR_1index.sh`;

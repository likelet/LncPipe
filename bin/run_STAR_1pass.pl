#!/usr/bin/perl -w
use strict;

`mkdir -pv STAR`;

open FH,"CRC_datainfo.txt" or die;

while(<FH>){
	chomp;
	my ($id,$seq1,$seq2)=split "\t";
	open OUT,">STAR/$id\_STAR_1pass.sh" or die;
	print OUT "#!/bin/sh\n";
	print OUT "mkdir -pv /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/$id\_1pass\n";
	print OUT "cd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/$id\_1pass\n";
	print OUT "gzip -dc $seq1 > $id.1.fq\n";
	print OUT "gzip -dc $seq2 > $id.2.fq\n";
	print OUT "STAR --runThreadN 25 --genomeDir /HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/STARIndex --readFilesIn $id.1.fq $id.2.fq --chimSegmentMin 20 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\n";
	`chmod +x STAR/$id\_STAR_1pass.sh`;
	`yhbatch -p nsfc3 STAR/$id\_STAR_1pass.sh`;
}

#!/usr/bin/perl -w
use strict;


open FH,"CRC_datainfo.txt" or die;

while(<FH>){
	chomp;
	my ($id,$seq1,$seq2)=split "\t";
	open OUT,">STAR/$id\_STAR_2pass.sh" or die;
	print OUT "#!/bin/sh\n";
	print OUT "mkdir /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/$id\_2pass\n";
	print OUT "cd /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/$id\_2pass\n";
	print OUT "STAR --runThreadN 25 --genomeDir /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/1pass_index --readFilesIn ../$id\_1pass/$id.1.fq ../$id\_1pass/$id.2.fq  --outSAMtype BAM SortedByCoordinate  --chimSegmentMin 20 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\n";
	`chmod +x STAR/$id\_STAR_2pass.sh`;
	`yhbatch -p nsfc3 STAR/$id\_STAR_2pass.sh`;
}

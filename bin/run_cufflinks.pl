#!/usr/bin/perl -w
use strict;

open FH,"CRC_datainfo.txt" or die;
my $gtf="/HOME/nsfc2015_555/NSFC/WORKSPACE/database/genecode_hg38/gencode.v24.annotation.gtf";
my $ref="/HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa";
my $rRNAmask="/HOME/nsfc2015_555/NSFC/WORKSPACE/database/UCSC_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/rRNA.gtf";
while(<FH>){
	chomp;
	my ($id,$seq1,$seq2)=split "\t";
	open OUT,">cufflinks/$id.cufflinks.sh" or die;
	print OUT "#!/bin/sh\n";
	print OUT "cufflinks -g $gtf -b $ref -M $rRNAmask --library-type fr-firststrand --max-multiread-fraction 0.25 --3-overhang-tolerance 2000 -o /NSFCGZ/nsfc2015_555/WORKSPACE/data/zuo/project/CRC_RNASeq/cufflinks/$id\_cufflinks -p 25 /HOME/nsfc2015_555/NSFC/WORKSPACE/data/zuo/project/CRC_RNASeq/STAR/$id\_2pass/Aligned.sortedByCoord.out.bam\n";
	`chmod +x cufflinks/$id.cufflinks.sh`;
	`yhbatch -p nsfc3 cufflinks/$id.cufflinks.sh`;
}


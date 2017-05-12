#!/usr/bin/perl -w
use strict;

open FH,"lncRNA/gencode.v24.protein_coding.gtf" or die;

while(<FH>){
	$_=~/gene_name "(.+?)";/;
	my $genename=$1;
	$_=~s/gene_id ".+?"/gene_id "$genename"/;
	print $_;
}

#!/usr/bin/perl -w
use strict;

open FH,shift or die;

while(<FH>){
	$_=~/gene_name "(.+?)";/;
	my $genename=$1;
	$_=~s/gene_id ".+?"/gene_id "$genename"/;
	print $_;
}

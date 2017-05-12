#!/usr/bin/perl -w
use strict;

my $gtf=shift;

open FH,$gtf or die;
my %count;
while(<FH>){
	chomp;
	my @field=split "\t";
	$field[8]=~/transcript_id \"(.+?)\"/;
			
	my $tid=$1;
	if(!exists $count{$tid}){
		$count{$tid}=1;
	}else{	
		$count{$tid}++;
	}
}

foreach my $k (keys %count){
	print $k."\t".$count{$k}."\n";
}

#!/usr/bin/perl -w
use strict;

open FH,"lncRNA/novel.longRNA.CPAT.out" or die;

<FH>;
my %info;
while(<FH>){
	chomp;
	my @field=split "\t";
	if($field[5]>=0.364){
		$info{$field[0]}{CPAT}="TUCP";
	}else{
		$info{$field[0]}{CPAT}="lncRNA";
	}
}

open FH,"lncRNA/novel.longRNA.PLEK.out" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	$field[2]=~/>(.+?) /;
	my $id=$1;
	if($field[0] eq "Coding"){
		$info{$id}{PLEK}="TUCP";
	}else{
		$info{$id}{PLEK}="lncRNA";
	}
}

open FH,"lncRNA/novel.longRNA.exoncount.txt" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	$info{$field[0]}{EXONCOUNT}=$field[1];
}

foreach my $id (sort keys %info){
	print $id."\t".$info{$id}{CPAT}."\t".$info{$id}{PLEK}."\t".$info{$id}{EXONCOUNT}."\t";
	if($info{$id}{CPAT} eq "lncRNA" && $info{$id}{PLEK} eq "lncRNA"){
		print "lncRNA\n";
	}else{
		print "TUCP\n";
	}
	
}

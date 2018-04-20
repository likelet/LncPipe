#!/usr/bin/perl -w
use strict;

my $gtf=shift;

open FH,$gtf or die;

my %map;
while(<FH>){
	chomp;
	if (/^#/){
	next;
}
	my @field=split "\t";
	#print $field[8]."\n";
	$field[8]=~/transcript_id \"(.+?)\"/;
	
	my $tid=$1;
	if(!exists $map{$tid}){
		$map{$tid}=$_;
	}else{
		$map{$tid}.="\n".$_;
	}
}

my $idfile=shift;

open FH,$idfile or die;

while(<FH>){
	chomp;
	print $map{$_}."\n";
}

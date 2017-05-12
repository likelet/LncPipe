#!/usr/bin/perl -w
use strict;

open FH,"lncRNA/lncRNA.final.gtf" or die;

my %class;
my %g2t;
my %trans_len;
my %exon_num;
while(<FH>){
	chomp;
	my @field=split "\t";
	$_=~/gene_id "(.+?)"/;
	my $gid=$1;
	$_=~/transcript_id "(.+?)"/;
	my $tid=$1;
	$class{$tid}=$field[1];
	$g2t{$tid}=$gid;
	my $len=$field[4]-$field[3];
	$trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
	$exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
}
open FH,"lncRNA/protein_coding.final.gtf" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	$_=~/gene_id "(.+?)"/;
	my $gid=$1;
	$_=~/transcript_id "(.+?)"/;
	my $tid=$1;
	$class{$tid}="protein_coding";
	$g2t{$tid}=$gid;
	my $len=$field[4]-$field[3];
	$trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
	$exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
}

open FH,"lncRNA/lncRNA.final.CPAT.out" or die;

<FH>;

while(<FH>){
	chomp;
	my @field=split "\t";
	my $tid=$field[0];
	print $g2t{$tid}."\t".$tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\n";
}

open FH,"lncRNA/protein_coding.final.CPAT.out" or die;

<FH>;

while(<FH>){
	chomp;
	my @field=split "\t";
	my $tid=$field[0];
	print $g2t{$tid}."\t".$tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\n";
}

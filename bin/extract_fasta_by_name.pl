#!/usr/bin/perl

use strict;
use warnings;

$ARGV[1] or die "use extractSeq.pl LIST FASTA > OUT\n";

my $list = shift @ARGV;
my $fasta = shift @ARGV;

my %select;
open L, "$list" or die;
while (<L>) {
    chomp;
    
    $select{$_} = 1;
}
close L;



open F, "$fasta" or die;
my $id;
while (<F>) {
	
	if($_=~/>(.+?):/){
		$id = $1;
		
	}

	
    print  "$_" if (exists $select{$id});
}
close F;
close L;

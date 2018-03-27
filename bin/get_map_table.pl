#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Data::Dumper;

my $gtf_file = '';
GetOptions('gtf_file=s'=>\$gtf_file);

unless($gtf_file){
  print "\n\nRequired parameters missing\n\n";
  print "Usage:  get_map_table.pl --gtf_file=final_all.gtf \n\n";
  exit();
}
die "\n\nCould not find GTF file: $gtf_file\n\n" unless(-e $gtf_file);

#Build a map of transcript to gene IDs
my %trans;
my $gtf_fh = IO::File->new($gtf_file, 'r');
while (my $gtf_line = $gtf_fh->getline) {
  chomp($gtf_line);
  my @gtf_entry = split("\t", $gtf_line);
  next unless $gtf_entry[2] eq 'exon';
  my $g_id = '';
  my $t_id = '';
  my $type = '';
  if ($gtf_entry[8] =~ /gene_id\s\"(\S+)\";\s\S+/){
    $g_id = $1;
	if($gtf_entry[8] =~ /protein_coding/){
		$type = "protein_coding";
	}else{	$type = $gtf_entry[1]; }
  }
  
  if ($gtf_entry[8] =~ /transcript_id\s\"(\S+)\";\s\S+/){
    $t_id = $1;
	if($gtf_entry[8] =~ /protein_coding/){
		$type = "protein_coding";
	}else{	$type = $gtf_entry[1]; }
  }
  my @array = '' ;
  $array[0] = $g_id;
  $array[1] = $t_id;
  $array[2] = $type;
  die "\n\nCould not identify gene and transcript id in GTF transcript line:\n$gtf_line\n\n" unless ($g_id && $t_id);
  $trans{$t_id} = \@array;
  
}
$gtf_fh->close;

foreach my $k (keys %trans) {
#foreach my $k (sort {$a<=>$b} keys %trans) {  
   print join "\t",@{$trans{$k}};  ## no sorted!!!
   print "\n";
}
exit;




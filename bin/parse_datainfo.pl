#!/usr/bin/perl -w
use strict;
use File::Basename;


open FH,"datainfo1.txt" or die;
my %datainfo;
while(<FH>){
	chomp;
	my ($fn,$dir,$suffix)=fileparse $_;
	my $id=(split "\_",$fn)[0];
	if($fn=~/1.clean.fq.gz/){
		$datainfo{$id}{1}=$_;
	}
	if($fn=~/2.clean.fq.gz/){
            $datainfo{$id}{2}=$_;
        }
	
}
close FH;
open FH,"datainfo2.txt" or die;

while(<FH>){
        chomp;
        my ($fn,$dir,$suffix)=fileparse $_;
         my $id=(split "\_",$fn)[0];
         if($fn=~/1.clean.fq.gz/){
                 $datainfo{$id}{1}=$_;
         }
         if($fn=~/2.clean.fq.gz/){
             $datainfo{$id}{2}=$_;
         }
         
  }

foreach my $id (sort keys %datainfo){
	print $id."\t".$datainfo{$id}{1}."\t".$datainfo{$id}{2}."\n";
}


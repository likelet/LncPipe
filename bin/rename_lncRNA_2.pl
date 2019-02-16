#!/usr/bin/perl -w
use strict;
#die ("usage: <Genecode gtf file> <lncipedia gtf file>") unless @ARGV > 2;
#print "#Query file ".$ARGV[0]." with file_number ".$ARGV[1]."\n";

my %know_lnc;
open FH,"known.lncRNA.bed" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	if ($field[7] eq "exon"){
		$know_lnc{$field[0].'\t'.$field[1].'\t'.$field[5]} = $field[3];
		$know_lnc{$field[0].'\t'.$field[2].'\t'.$field[5]} = $field[3];
	}
}


my %genecode;my %lncpedia;
if (@ARGV == 2){
open FH,"$ARGV[0]" or die;
while(<FH>){
	chomp;
	if ($_ =~ /^#/){
		next;
	}
	my @field=split "\t";
	if ($field[2] ne "exon"){
		next;
	}
	$_=~/gene_id "(.+?)"/;
	my $gene_name=$1;
	my $loc1 = $field[0].'\t'.($field[3]-1).'\t'.$field[6];
	my $loc2 = $field[0].'\t'.$field[4].'\t'.$field[6];
	if (defined($know_lnc{$loc1})){
			$genecode{$know_lnc{$loc1}} = $gene_name;
	}
	if (defined($know_lnc{$loc2})){
			$genecode{$know_lnc{$loc2}} = $gene_name;
	}
}
open FH,"$ARGV[1]" or die;
while(<FH>){
	chomp;
	if ($_ =~ /^#/){
		next;
	}
	my @field=split "\t";
	if ($field[2] ne "exon"){
		next;
	}
	$_=~/gene_id "(.+?)"/;
	my $gene_name=$1;
	my $loc1 = $field[0].'\t'.($field[3]-1).'\t'.$field[6];
	my $loc2 = $field[0].'\t'.$field[4].'\t'.$field[6];
	if (defined($know_lnc{$loc1})){
			$lncpedia{$know_lnc{$loc1}} = $gene_name;
	}
	if (defined($know_lnc{$loc2})){
			$lncpedia{$know_lnc{$loc2}} = $gene_name;
	}
}
}elsif (@ARGV == 1){
open FH,"$ARGV[0]" or die;
while(<FH>){
	chomp;
	if ($_ =~ /^#/){
		next;
	}
	my @field=split "\t";
	if ($field[2] ne "exon"){
		next;
	}
	$_=~/gene_id "(.+?)"/;
	my $gene_name=$1;
	my $loc1 = $field[0].'\t'.($field[3]-1).'\t'.$field[6];
	my $loc2 = $field[0].'\t'.$field[4].'\t'.$field[6];
	if (defined($know_lnc{$loc1})){
			$genecode{$know_lnc{$loc1}} = $gene_name;
	}
	if (defined($know_lnc{$loc2})){
			$genecode{$know_lnc{$loc2}} = $gene_name;
	}
}
}else{
	die ("usage: at least one gtf file is needed!!!")
}
my %exon;
my %gene;

open FH,"known.lncRNA.bed" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	my $chr=$field[0];
	my $genename="known:".$field[3];
	my $start=$field[1];
	my $end=$field[2];
	my $strand=$field[5];
	$_=~/transcript_id "(.+?)"/;
	my $transid=$1;
	my $loc=$chr."\t".($start+1)."\t".$end."\t".$strand;
	$exon{$genename}{$transid}{$loc}=1;
	$gene{$genename}{CHR}=$chr;
	$gene{$genename}{STRAND}=$strand;
	if(exists $gene{$genename}{START}){
		$gene{$genename}{START}=$gene{$genename}{START}<$start?$gene{$genename}{START}:$start;
	}else{
		$gene{$genename}{START}=$start;
	}
	if(exists $gene{$genename}{END}){
		$gene{$genename}{END}=$gene{$genename}{END}>$end?$gene{$genename}{END}:$end;
	}else{
		$gene{$genename}{END}=$end;
	}

}

open FH,"novel.lncRNA.stringent.filter.bed" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	my $chr=$field[0];
	my $genename="novel:".$field[3];
	my $start=$field[1];
	my $end=$field[2];
	my $strand=$field[5];
	$_=~/transcript_id "(.+?)"/;
	my $transid=$1;
	my $loc=$chr."\t".($start+1)."\t".$end."\t".$strand;
	$exon{$genename}{$transid}{$loc}=1;
	$gene{$genename}{CHR}=$chr;
	$gene{$genename}{STRAND}=$strand;
	if(exists $gene{$genename}{START}){
		$gene{$genename}{START}=$gene{$genename}{START}<$start?$gene{$genename}{START}:$start;
	}else{
		$gene{$genename}{START}=$start;
	}
	if(exists $gene{$genename}{END}){
		$gene{$genename}{END}=$gene{$genename}{END}>$end?$gene{$genename}{END}:$end;
	}else{
		$gene{$genename}{END}=$end;
	}

}
open OUT,">lncRNA.for_anno.bed" or die;
foreach my $k (keys %gene){
	print OUT $gene{$k}{CHR}."\t".$gene{$k}{START}."\t".$gene{$k}{END}."\t".$k."\t.\t".$gene{$k}{STRAND}."\n";
}

`sort-bed lncRNA.for_anno.bed > lncRNA.for_anno.srt.bed`;

`closest-features --dist lncRNA.for_anno.srt.bed  gencode.protein_coding.gene.bed > lncRNA.for_anno.srt.neighbour.txt`;

open FH,"lncRNA.for_anno.srt.neighbour.txt" or die;
my %map;
my %genename2gene;
my $naidx=0;
while(<FH>){
	chomp;
	my @field=split /\|/;
	my ($chr,$start,$end,$geneid,$tmp2,$strand)=split "\t",$field[0];
	if($strand ne "+" and $strand ne "-"){
		#print $field[0]."\n";
		next;
	}
	my $up_gene="";
	my $up_dist=999999;
	my $up_strand;
	my $down_gene="";
	my $down_dist=999999;
	my $down_strand;
	my $closet_genename="";
	if($field[1] ne "NA"){
		$field[1]=~/gene_name "(.+?)"/;
		$up_gene=$1;
		$up_dist=abs($field[2]);
		my @tmp=split "\t",$field[1];
		$up_strand=$tmp[5];
	}
	if($field[3] ne "NA"){
		$field[3]=~/gene_name "(.+?)"/;
		$down_gene=$1;
		$down_dist=abs($field[4]);
		my @tmp=split "\t",$field[3];
		$down_strand=$tmp[5];
	}
	if($field[1] eq "NA" and $field[3] eq "NA"){
		$naidx++;
		my $genename="NA-$naidx";
		$map{$genename}{$geneid}="NA";
		#print $genename."\n";
	}else{
	if($up_dist < $down_dist){
		my $genename;
		if($up_dist==0){
			if($strand ne $up_strand){
				$genename=$up_gene."-AS";
				$genename2gene{$genename} = $up_gene;
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $up_dist?$up_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$up_dist;
				}

			}else{
				#print "LALALALALA"."\n";
			}
		}else{
			$genename="LINC-".$up_gene;
			$genename2gene{$genename} = $up_gene;
			if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $up_dist?$up_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$up_dist;
				}
		}
	}else{
		my $genename;
		if($down_dist==0){
			if($strand ne $down_strand){
				$genename=$down_gene."-AS";
				$genename2gene{$genename} = $down_gene;
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
			}else{
				#print "HAHAHAHAHAHAA"."\n";
			}
		}else{
			$genename="LINC-".$down_gene;
			$genename2gene{$genename} = $down_gene;
			if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
		}
	}
	}
}
my %MSTRG2genename;
open OUT1,">lncRNA.final.v2.gtf" or die;
open OUT2,">lncRNA.final.v2.map" or die;
foreach my $genename (keys %map){
	my $tmp1=$map{$genename};
	my %tmp1=%$tmp1;
	my $gindex=0;
	my $out_gene = $genename2gene{$genename};
	foreach my $id (sort {$tmp1{$a}<=>$tmp1{$b}} keys %tmp1){
		my ($class,$cuffid)=split ":",$id;
		my $dist=$tmp1{$id};
		my $tmp2=$exon{$id};
		my %tmp2=%$tmp2;
		$gindex++;
		my $geneid=$genename."-".$gindex;
		$MSTRG2genename{$cuffid} = $geneid;
		my $tindex=0;
		print OUT2 $geneid."\t".$cuffid."\n";
		foreach my $tid (keys %tmp2){
			my $tmp3=$tmp2{$tid};
			my %tmp3=%$tmp3;
			$tindex++;
			my $transid=$geneid.":".$tindex;
			foreach my $exon (keys %tmp3){
				my @loc=split "\t",$exon;
				print OUT1 $loc[0]."\t$class\texon\t".$loc[1]."\t".$loc[2]."\t.\t".$loc[3]."\t.\tgene_id \"$geneid\"; transcript_id \"$transid\"; dist \"$dist\"; "."closet_gene \"".$out_gene."\";"."\n";
			}

		}
	}
}
my %all_data;
foreach my $mstr(sort(keys %genecode)){
	$all_data{$mstr} = 1;
}
foreach my $mstr(sort(keys %lncpedia)){
	$all_data{$mstr} = 1;
}
open OUT3,">lncRNA.mapping.file" or die;
foreach my $mstr(sort(keys %all_data)){
	if(defined($MSTRG2genename{$mstr}) && defined($lncpedia{$mstr}) && defined($genecode{$mstr})){
		print OUT3 $mstr."\t".$MSTRG2genename{$mstr}."\t".$genecode{$mstr}."\t".$lncpedia{$mstr}."\n";
	}elsif (defined($MSTRG2genename{$mstr}) && defined($genecode{$mstr})){
		print OUT3 $mstr."\t".$MSTRG2genename{$mstr}."\t".$genecode{$mstr}."\t\n"
	}elsif (defined($MSTRG2genename{$mstr}) && defined($lncpedia{$mstr})){
		print OUT3 $mstr."\t".$MSTRG2genename{$mstr}."\t\t".$lncpedia{$mstr}."\n"
	}else{
		next;
	}
}

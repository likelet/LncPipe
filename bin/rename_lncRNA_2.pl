#!/usr/bin/perl -w
use strict;



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
		my @tmp=split "\t",$field[1];
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
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $up_dist?$up_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$up_dist;
				}
				
			}
		}else{
			$genename="LINC-".$up_gene;
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
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
			}
		}else{
			$genename="LINC-".$down_gene;
			if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
		}
	}
	}

}

open OUT1,">lncRNA.final.v2.gtf" or die;
open OUT2,">lncRNA.final.v2.map" or die;
foreach my $genename (keys %map){
	my $tmp1=$map{$genename};
	my %tmp1=%$tmp1;
	my $gindex=0;
	foreach my $id (sort {$tmp1{$a}<=>$tmp1{$b}} keys %tmp1){
		my ($class,$cuffid)=split ":",$id;
		my $dist=$tmp1{$id};
		my $tmp2=$exon{$id};
		my %tmp2=%$tmp2;
		$gindex++;
		my $geneid=$genename."-".$gindex;
		my $tindex=0;
		print OUT2 $geneid."\t".$cuffid."\n";
		foreach my $tid (keys %tmp2){
			my $tmp3=$tmp2{$tid};
			my %tmp3=%$tmp3;
			$tindex++;
			my $transid=$geneid.":".$tindex;
			foreach my $exon (keys %tmp3){
				my @loc=split "\t",$exon;
				print OUT1 $loc[0]."\t$class\texon\t".$loc[1]."\t".$loc[2]."\t.\t".$loc[3]."\t.\tgene_id \"$geneid\"; transcript_id \"$transid\"; dist \"$dist\";\n";
			}

		}
	}	
	

	
}

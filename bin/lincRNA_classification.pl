#!perl
use strict;
use warnings;

#input format
die "USAGE: lncRNA_file genecode_file OUT_file" unless @ARGV == 3;
open IN1,"$ARGV[0]" || die "$!";
open IN2,"$ARGV[1]" || die "$!";
open OUT,"> $ARGV[2]" || die "$!";
my %gene_data;
my %gene_dir;
#read genecode file and store the data
while(<IN2>){
	chomp;
	my @data = split /\t/,$_;
	if ($data[2] ne 'exon'){
		next;
	}
	$data[8] =~ /gene_name\s\"(.*)\";\stranscript_type/;
	my $gene_name = $1;
	my @arr = [$data[3],$data[4]];
	push(@{$gene_data{$gene_name}},@arr);
	$gene_dir{$gene_name} = $data[6];
}

#read the lncRNA file per line and classify
my @results;
while(<IN1>){
	chomp;
	my $result;
	my $smname;
	my @data = split /\t/,$_;
	#if ($data[8] !~ /transcript_id\s\"LINC/){        #remove not LINC
	#	next;
	#}
	my @rna = split /;/,$data[8];
	$rna[1] =~ /\"(.*)\"/;
	my $linc = $1;
	$rna[2] =~ /\"(\d+)\"/;
	my $distance = $1;
	$_ =~ /closet_gene "(.*)"/;
	$smname = $1;
	if ($distance > 1000){                             #c1
		$result = $linc."\t"."Intergenic";
		push(@results,$result);
	}elsif($distance > 0){                             #c1 and c2      
		if($gene_dir{$smname} eq $data[6]){
			$result = $linc."\t"."Intergenic";
			push(@results,$result);
		}else{
			$result = $linc."\t"."Bidirectional";
			push(@results,$result);
		}
	}elsif($distance == 0){                         #c3 and c4 and c5
		if ($gene_dir{$smname} ne $data[6]){          #c3
			$result = $linc."\t"."Antisense";
			push(@results,$result);
		}else{                                        #c4 and c5
			my $m = 0;
			foreach my $gename(keys %gene_data){
				if($gename eq $smname){
					my @num = @{$gene_data{$gename}};
					for (my $i =0;$i <= $#num;$i++){
						if(($num[$i][0] < $data[4])&&($num[$i][1] > $data[3])){
							$m++;
						}
					}
				}
			}
			if ($m > 0){
				$result = $linc."\t"."Exonic Sense";
				push(@results,$result);
			}else{
				$result = $linc."\t"."Intronic Sense";
				push(@results,$result);
			}
		}
	}
}
my %hash;@results = grep { ! $hash{$_} ++ } @results;
foreach my $result(@results){
	print OUT $result."\n";
}
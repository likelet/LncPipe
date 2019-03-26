#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Design.pl input.txt > Design.txt
#
#  DESCRIPTION: Generate an experiment design file
#
#  INPUT FILES:
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: //2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use strict;

open IN, shift or die $!;
print "SampleID\tcondition\n";

while(<IN>){
 chomp;
 my @F = split /:/;
 my @G = split /,/, $F[1];
 for (@G){
  print "$_\t$F[0]\n";
 }
}
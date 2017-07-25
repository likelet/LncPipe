#!/bin/sh
# required software 
#cuffmergepath=cuffmerge
#cuffcomparepath=cuffcompare
#gtf2bedpath=gtf2bed
#sort-bed



#################
#run on tianhe#
################ 
#Step 1 reads mapping by STAR
perl run_STAR_1pass.pl
perl run_STAR_1index.pl
perl run_STAR_2pass.pl
#Step 2 assembling with cufflinks 
perl run_cufflinks.pl
perl run_cuffmerge.pl
sh run_cuffcompare.sh
#Step 3 filtering transcripts
#   counting exon number based on gtf
perl get_exoncount.pl cuffmerge/merged.gtf > cuffmerge/transcript_exoncount.txt
#   filtering novel lncRNA based on cuffmerged trascripts 
awk '$3=="x"||$3=="u"||$3=="i"{print $0}' cuffmerge/merged.merged.gtf.tmap > lncRNA/novel.gtf.tmap
#   excluding length smaller than 200 nt
awk '$11>200{print}' lncRNA/novel.gtf.tmap > lncRNA/novel.longRNA.gtf.tmap
#   extract gtf
awk '{print $5}' lncRNA/novel.longRNA.gtf.tmap |perl extract_gtf_by_name.pl cuffmerge/merged.gtf - >lncRNA/novel.longRNA.gtf
perl get_exoncount.pl lncRNA/novel.longRNA.gtf > lncRNA/novel.longRNA.exoncount.txt
# get fasta sequencing from gtf
sh get_fasta_from_gff3.sh
#Step 4 predicting  coding potential
# coding poteintial prediction 
sh run_PLEK.sh
sh run_CPAT.sh
sh run_CNCI.sh


################
# run locally###
################
#merged transcripts 
perl integrate_novel_transcripts.pl > lncRNA/novel.longRNA.txt

awk '$4>1{print $1}' lncRNA/novel.longRNA.txt|perl extract_gtf_by_name.pl cuffmerge/merged.gtf - >lncRNA/novel.longRNA.stringent.gtf
# retain lncRNA only by coding ability
awk '$4>1&&$5=="lncRNA"{print $1}' lncRNA/novel.longRNA.txt|perl extract_gtf_by_name.pl cuffmerge/merged.gtf - >lncRNA/novel.lncRNA.stringent.gtf
awk '$4>1&&$5=="TUCP"{print $1}' lncRNA/novel.longRNA.txt|perl extract_gtf_by_name.pl cuffmerge/merged.gtf - >lncRNA/novel.TUCP.stringent.gtf
#########################################################
#compile public lncRNA annotation (gencode and lncipedia)#
########################################################
#lncRNA.gtflist 
#gencode.v24.long_noncoding_RNAs.gtf
#lncipedia_4_0_hg38.gtf
cd /data/public_data/hg38

cuffmerge -o /data/public_data/hg38/merged_lncRNA /data/public_data/hg38/lncRNA.gtflist 

cat /data/public_data/hg38/gencode.v24.annotation.gtf |grep "protein_coding" > /data/public_data/hg38/gencode.v24.protein_coding.gtf

#run cuffcompare  merged gtf with gencode annotation
cuffcompare -o /data/public_data/hg38/merged_lncRNA -r /data/public_data/hg38/gencode.v24.protein_coding.gtf -p 12 /data/public_data/hg38/merged_lncRNA/merged.gtf

awk '$3=="u"||$3=="x"{print $5}' /data/public_data/hg38/merged_lncRNA/merged_lncRNA.merged.gtf.tmap |sort|uniq|perl extract_gtf_by_name.pl /data/public_data/hg38/merged_lncRNA/merged.gtf - > /data/public_data/hg38/merged_lncRNA/merged.filter.gtf

cp /data/public_data/hg38/gencode.v24.protein_coding.gtf lncRNA/gencode.v24.protein_coding.gtf
cp /data/public_data/hg38/merged_lncRNA/merged.filter.gtf lncRNA/known.lncRNA.gtf

#further filter novel lncRNA
 cuffcompare -o lncRNA/filter -r  lncRNA/known.lncRNA.gtf -p 12 lncRNA/novel.lncRNA.stringent.gtf
awk '$3=="u"||$3=="x"{print $5}' lncRNA/filter.novel.lncRNA.stringent.gtf.tmap |sort|uniq|perl extract_gtf_by_name.pl lncRNA/novel.lncRNA.stringent.gtf - > lncRNA/novel.lncRNA.stringent.filter.gtf

#rename lncRNAs according to neighbouring protein coding genes

 awk '$3=="gene"{print }' lncRNA/gencode.v24.protein_coding.gtf > lncRNA/gencode.v24.protein_coding.gene.gtf

 #############################

gtf2bed < lncRNA/gencode.v24.protein_coding.gene.gtf |sort-bed - > lncRNA/gencode.v24.protein_coding.gene.bed
gtf2bed < lncRNA/novel.lncRNA.stringent.filter.gtf |sort-bed - > lncRNA/novel.lncRNA.stringent.filter.bed
gtf2bed < lncRNA/known.lncRNA.gtf |sort-bed - > lncRNA/known.lncRNA.bed


perl rename_lncRNA_2.pl

perl rename_proteincoding.pl > lncRNA/protein_coding.final.gtf



#use RSEM to do quantification
rsem_index_star_lncRNA.sh
rsem_index_star_pc.sh
perl run_RSEM_STAR.pl
Rscript get_rsem_matrix.R

#coding potential (CPAT) for final lncRNA set and protein_coding set
sh run_CPAT_for_lncRNA_final_gtf.sh
sh run_CPAT_for_pc_final_gtf.sh


# some statistics

#compare coding poteintial between novel lncRNA, annotated lncRNA and proteins

perl compare_basic_charac.pl > lncRNA/basic_charac.txt



#DE analysis
cd lncRNA
rsem-run-ebseq lncRNA.rsem.count.txt 20,20 lncRNA.rsem.EBSeq.out.v2.txt
rsem-run-ebseq pc.rsem.count.txt 20,20 pc.rsem.EBSeq.out.txt

cut -f1,12-21,32-41 lncRNA.rsem.count.txt > lncRNA.rsem.count.noR.txt
rsem-run-ebseq lncRNA.rsem.count.noR.txt 20,10 lncRNA.rsem.EBSeq.TvN.out.txt

cut -f1,12-21,32-41 pc.rsem.count.txt > pc.rsem.count.noR.txt
rsem-run-ebseq pc.rsem.count.noR.txt 20,10 pc.rsem.EBSeq.TvN.out.txt

cut -f1,22-41 lncRNA.rsem.count.txt > lncRNA.rsem.count.tumor.txt
rsem-run-ebseq lncRNA.rsem.count.tumor.txt 10,10 lncRNA.rsem.EBSeq.RvT.out.txt

cut -f1,22-41 pc.rsem.count.txt > pc.rsem.count.tumor.txt
rsem-run-ebseq pc.rsem.count.tumor.txt 10,10 pc.rsem.EBSeq.RvT.out.txt




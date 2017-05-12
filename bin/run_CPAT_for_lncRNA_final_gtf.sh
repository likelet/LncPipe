#!/bin/sh
#gffread lncRNA/lncRNA.final.gtf -o-|gffread -  -g /data/public_data/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -w lncRNA/lncRNA.final.fa -W
cpat.py -g lncRNA/lncRNA.final.fa -x /data/program/CPAT/prebuilt_models/Human_Hexamer.tab -d /data/program/CPAT/prebuilt_models/Human_logitModel.RData -o lncRNA/lncRNA.final.CPAT.out

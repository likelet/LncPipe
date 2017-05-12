#!/bin/sh
gffread lncRNA/protein_coding.final.gtf -o-|gffread -  -g /data/public_data/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -w lncRNA/protein_coding.final.fa -W
cpat.py -g lncRNA/protein_coding.final.fa -x /data/program/CPAT/prebuilt_models/Human_Hexamer.tab -d /data/program/CPAT/prebuilt_models/Human_logitModel.RData -o lncRNA/protein_coding.final.CPAT.out

#!/usr/bin/env bash
for file in *.hisat2_summary.txt
do
grep -H '' $file | perl -F':|\s+|\(|\)|(>?\d\stimes?)' -lanE '/%/ and /aligned/ ? say join qq{\t}, @F : next' | tr -s '\t' '\t' > ${file%%.hisat2_summary.txt}.hisat2_format.tsv
done

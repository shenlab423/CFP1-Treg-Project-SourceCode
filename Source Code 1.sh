#!/bin/bash

for IP_sample in *.bam  
do
bamToBed -i ${IP_sample} |~/tools/extend_single ~/tools/mm9.chrom.sizes 300|sort -k1,1 -k2,2n > ${IP_sample/.bam/}.bed
IP_mappable_read_count=$(samtools view -F 0x0004 ${IP_sample} | wc -l)
for peakfile in mm9_promoter.bed  
do
cut -f 1,2,3 ${peakfile}|coverageBed -b ${IP_sample/.bam/}.bed -a - | awk -v OFS='\t' -v SIZE=$IP_mappable_read_count '{print $1,$2,$3,($4*1000000/SIZE)*1000/($3-$2)}' > ${IP_sample/.bam/}-${peakfile/.bed/.rpkm}
done
done

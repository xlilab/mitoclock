#!/usr/bin/env bash
#Description: 
#Author: Wang Zhenguo
#Created Time: 2022/01/07 10:29

echo "Usage: sh $0 [ TISSUE ]"

#set -e -u -o pipefail

TISSUE=$1

echo -e "Indiv\tCPN" > $TISSUE.copy_number.txt
for i in `ls ~/data/datasets/gtex/RNA_bam/v8/RNAseq_bam/$TISSUE/*.bam`
do
	ID=`echo $i |sed 's/.*\(GTEX-[0-9A-Z]\{4,5\}\)-.*/\1/'`
	samtools idxstats $i |\
		head -n 25 |\
		awk 'BEGIN{OFS="\t"} NR <= 24 {nuclear = nuclear + $3} NR == 25 {chrM = $3} END{print "'$ID'",chrM / (nuclear + chrM)}' >> $TISSUE.copy_number.txt
done

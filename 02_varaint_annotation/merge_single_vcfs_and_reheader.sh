#!/usr/bin/env bash

echo "Usage: $0 [ tissue_file out_dir ]"

set -e -u -o pipefail

tissue_file=$1
out_dir=`readlink -f $2`

for tissue in `cat $tissue_file`
do
	echo "processing $tissue..."
	cd ../01_variant_calling/$tissue/output_vcf/
	
	#################################################

	echo "merging original vcfs..."
	for i in `ls *.original.vcf`
	do
		bcftools view $i -Oz -o `basename $i .original.vcf`.original.vcf.gz
		tabix -p vcf `basename $i .original.vcf`.original.vcf.gz
	done
	
	bcftools merge -0 -m none --threads 10 *.original.vcf.gz |\
		bcftools annotate --set-id "%CHROM\_%POS\_%REF\_%ALT" -Oz -o $out_dir/$tissue.final.split.varID.vcf.gz
	bcftools query -l $out_dir/$tissue.final.split.varID.vcf.gz >$out_dir/$tissue.indivID.id
	sed -i -r 's/(GTEX-[0-9A-Z]{4,5}).+/\1/' $out_dir/$tissue.indivID.id
	bcftools reheader -s $out_dir/$tissue.indivID.id $out_dir/$tissue.final.split.varID.vcf.gz -o $out_dir/$tissue.original.vcf.gz
	tabix -p vcf $out_dir/$tissue.original.vcf.gz
	
	rm $out_dir/$tissue.final.split.varID.vcf.gz
	echo "merging original vcf done!"
	
	#################################################
	
	#echo "merging consensus vcfs..."
	#for i in `ls *.consensus.vcf`
	#do
	#	bcftools view -i 'FILTER="PASS" & TYPE="snp"' $i -Oz -o `basename $i .consensus.vcf`.consensus.PASS.vcf.gz
	#	tabix -p vcf `basename $i .consensus.vcf`.consensus.PASS.vcf.gz
	#done
	#
	#bcftools merge -0 -m none --threads 10 *.original.PASS.vcf.gz |\
	#	bcftools annotate --set-id "%CHROM\_%POS\_%REF\_%ALT" -Oz -o $out_dir/$tissue.final.split.varID.vcf.gz
	#bcftools query -l $out_dir/$tissue.final.split.varID.vcf.gz >$out_dir/$tissue.indivID.id
	#sed -i -r 's/(GTEX-[0-9A-Z]{4,5}).+/\1/' $out_dir/$tissue.indivID.id
	#bcftools reheader -s $out_dir/$tissue.indivID.id $out_dir/$tissue.final.split.varID.vcf.gz -o $out_dir/$tissue.original.PASS.vcf.gz
	#tabix -p vcf $out_dir/$tissue.original.PASS.vcf.gz
	#
	#rm $out_dir/$tissue.final.split.varID.vcf.gz
	#echo "merging original vcf done!"
	
	#################################################
	
	cd ../../
	echo "$tissue finished!"
done

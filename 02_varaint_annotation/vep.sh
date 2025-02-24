#!/usr/bin/env bash
#Description: 
#Author: Wang Zhenguo
#Created Time: 2021/01/05 20:02

set -e -u -o pipefail

echo "Usage: sh $0 [VCF]"

VCF=$1
REF=/path_to/gatk/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta

vep --species homo_sapiens \
	--assembly GRCh38 \
	--input_file $VCF \
	--format vcf \
	--output_file `basename $VCF .vcf.gz`.vep_distance0.everything.vcf \
	--cache \
	--dir_cache /path_to/VEP/cache_101 \
	--fasta $REF \
	--cache_version 101 \
	--distance 0 \
	--allele_number \
	--vcf \
	--offline \
	--force_overwrite \
	--everything \
	--distance 0 \
	--plugin LoF,loftee_path:/path_to/VEP/LOFTEE/loftee,gerp_bigwig:/path_to/VEP/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/path_to/VEP/LOFTEE/human_ancestor.fa.gz,conservation_file:/path_to/VEP/LOFTEE/loftee.sql,phylocsf_data:/path_to/VEP/LOFTEE/phylocsf_gerp.sql \
	--custom /path_to/clinvar/GRCh38/clinvar_20230121.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
	--custom /path_to/gnomad/gnomAD/r3.x/gnomad.genomes.v3.1.sites.chrM.vcf.gz,gnomAD,vcf,exact,0,AN,AC_hom,AF_hom,AC_het,AF_het,max_hl 

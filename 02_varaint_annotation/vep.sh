#!/usr/bin/env bash
#Description: 
#Author: Wang Zhenguo
#Created Time: 2021/01/05 20:02

set -e -u -o pipefail

echo "Usage: sh $0 [VCF]"

VCF=$1
REF=/home/wangzhenguo/data/datasets/ref/gatk/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta

alias vep=/picb/lilab/tools/ensembl-vep-release-101/vep
vep --species homo_sapiens \
	--assembly GRCh38 \
	--input_file $VCF \
	--format vcf \
	--output_file `basename $VCF .vcf.gz`.vep_distance0.everything.vcf \
	--cache \
	--dir_cache /picb/lilab/annotation/VEP/cache_101 \
	--fasta $REF \
	--cache_version 101 \
	--distance 0 \
	--allele_number \
	--vcf \
	--offline \
	--force_overwrite \
	--everything \
	--distance 0 \
	--plugin LoF,loftee_path:/picb/lilab/annotation/VEP/LOFTEE/loftee,gerp_bigwig:/picb/lilab/annotation/VEP/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/picb/lilab/annotation/VEP/LOFTEE/human_ancestor.fa.gz,conservation_file:/picb/lilab/annotation/VEP/LOFTEE/loftee.sql,phylocsf_data:/picb/lilab/annotation/VEP/LOFTEE/phylocsf_gerp.sql \
	--custom /picb/lilab/annotation/clinvar/GRCh38/clinvar_20230121.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
	--custom /picb/lilab/annotation/gnomad/gnomAD/r3.x/gnomad.genomes.v3.1.sites.chrM.vcf.gz,gnomAD,vcf,exact,0,AN,AC_hom,AF_hom,AC_het,AF_het,max_hl 

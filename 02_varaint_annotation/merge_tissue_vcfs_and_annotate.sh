#!/usr/bin/env bash

out_dir=`readlink -f $1`

bcftools merge -m none --force-samples $out_dir/*original.vcf.gz |\
   bcftools view --drop-genotypes -Oz -o $out_dir/All_Tissues.only_variant.vcf.gz
sh vep.sh $out_dir/All_Tissues.only_variant.vcf.gz
python2 read_vep_vcf.py --vcf $out_dir/All_Tissues.only_variant.vep_distance0.everything.vcf >$out_dir/All_Tissues.only_variant.vep_distance0.everything.txt


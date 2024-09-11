#!/usr/bin/env bash

bcftools merge -m none --force-samples *original.vcf.gz |\
   bcftools view --drop-genotypes -Oz -o All_Tissues.only_variant.vcf.gz
sh vep.sh All_Tissues.only_variant.vcf.gz
python2 read_vep_vcf.py --vcf All_Tissues.only_variant.vep_distance0.everything.vcf >All_Tissues.only_variant.vep_distance0.everything.txt


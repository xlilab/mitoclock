# Mitochondrial Clock of Aging
Variant calling, filtering, plotting codes and supplementary tables for the manuscript "Tissue-specific mitochondrial clonal mosaicism as a molecular clock of aging".

<img width="571" alt="image" src="https://github.com/xlilab/mitoclock/assets/7442902/c7c6fa6e-4645-42de-91bf-be169917be8d">

# To run the code
## Install the dependencies
External software
- cromwell-52.jar
- GATK 4.2.0.0
- picard 2.23.3
- bwa
- haplocheckCLI
- bcftools
- samtools
- VEP

Unix packages
- GNU parallel
- tabix
- bgzip

## Download required files
Download from https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0 :
- files included in 01_variant_calling/mitochondria_m2_RNA_wdl/MitochondriaPipeline.inputs.json

Download from dbGaP:
- GTEx RNA bam files, these files should be organized by tissue

## Setup
Change the software path by running:
```
sed -i s#/path_to#/path_to_your_software_directory#g 01_variant_calling/mitochondria_m2_RNA_wdl/*wdl
sed -i s#/path_to#/path_to_your_software_directory#g 01_variant_calling/*sh
```
Change the reference file path by running:
```
sed -i s#/path_to#/path_to_your_reference_file_directory#g 01_variant_calling/mitochondria_m2_RNA_wdl/MitochondriaPipeline.inputs.json
sed -i s#/path_to#/path_to_your_reference_file_directory#g 02_varaint_annotation/vep.sh
```

# Pipeline
## Variant calling
```
cd 01_variant_calling
bash run_gatk_mitochondria_rna_pipeline_with_multi_tissues.sh tissue.list bam_dir
```
## Variant annotation
```
cd 02_varaint_annotations
bash merge_single_vcfs_and_reheader.sh tissue.list out_dir
bash merge_tissue_vcfs_and_annotate.sh out_dir
```
## Copy number calculation
```
cd 03_copy_number_calculation
bash calc_gtex_mt_chr_coverage.sh
bash calc_gtex_tissue_mt_copy_number.sh
bash estimate_RNA_copy_number.sh
```
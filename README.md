# Mitochondrial Clock of Aging
Variant calling, filtering, plotting codes and supplementary tables for the manuscript **"Tissue-specific mitochondrial clonal mosaicism as a molecular clock of aging"**.

<img width="571" alt="image" src="https://github.com/xlilab/mitoclock/assets/7442902/c7c6fa6e-4645-42de-91bf-be169917be8d">

# To run the code
## Install the dependencies
External software
- cromwell-52.jar
- GATK 4.2.0.0
- picard 2.23.3
- bwa 0.7.17
- haplocheckCLI
- bcftools 1.12
- samtools 1.13
- ensembl-vep 101.0

Unix packages
- GNU parallel
- tabix
- bgzip

R packages
- tidyverse
- data.table
- VGAM
- vcfR
- roll
- ggpubr
- boot
- ComplexHeatmap
- ggbeeswarm
- aplot
- ggtree
- ggthemes
- ggpie
- circlize

## Download required files
Download from https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0 :
- files included in `01_variant_calling/mitochondria_m2_RNA_wdl/MitochondriaPipeline.inputs.json`

Download from dbGaP :
- GTEx RNA bam files, these files should be organized by tissues, for example `Adipose-Subcutaneous/GTEX-1117F-0226-SM-5GZZ7.Aligned.sortedByCoord.out.patched.md.bam`

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
## 1. Variant calling
```
cd 01_variant_calling
bash run_gatk_mitochondria_rna_pipeline_with_multi_tissues.sh tissue.list bam_dir
```
## 2. Variant annotation
```
cd 02_varaint_annotations
bash merge_single_vcfs_and_reheader.sh tissue.list out_dir
bash merge_tissue_vcfs_and_annotate.sh out_dir
```
## 3. Copy number calculation
```
cd 03_copy_number_calculation
bash estimate_RNA_copy_number.sh
```
## 4. Tables formation
```
cd 04_tables_formation
Rscript estimate_sequencing_error_beta-binomial.R
Rscript generate_VxI_tables.R
Rscript processing_Pmtrna_Igtex_T_coverage.R
Rscript generate_Igtex_T_CPN.R
Rscript read_fixed_files.R
Rscript filter_variants.R
```
## 5. Figure plotting
Enter Rstudio, load data by runnig `load("04_tables_formation/analysis_base.Rdata")`, then execute codes in `05_figure_plotting`.

# Run example
To get unfiltered VCF files:
```
bash run_gatk_mitochondria_rna_pipeline_with_multi_tissues.sh 06_sample_input/tissue.list 06_sample_input
```
To plot figures (execute in Rstudio):  
Enter Rstudio, load data by runnig `load("06_sample_input/analysis_base.Rdata")`, then execute codes in `05_figure_plotting`.

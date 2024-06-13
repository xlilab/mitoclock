library(circlize)
library(ggpie)
library(ggthemes)
library(ggtree)
library(aplot)
library(ggbeeswarm)
library(ComplexHeatmap)
library(VGAM)
library(boot)
library(data.table)
library(ggpubr)
library(tidyverse)

# Functions ---------------------------------------------------------------

bootstrap_median <- function(data){
  f_median <- function(x,indices){
    d <- x[indices]
    return(median(d))
  }
  
  res.boot <- boot(data,statistic = f_median,R = 1000)
  res.ci <- boot.ci(res.boot,type = "perc")
  return(c(res.boot$t0,res.ci$percent[4:5]))
}

bootstrap_mean <- function(data){
  f_mean <- function(x,indices){
    d <- x[indices]
    return(mean(d))
  }
  
  res.boot <- boot(data,statistic = f_mean,R = 1000)
  res.ci <- boot.ci(res.boot,type = "perc")
  return(c(res.boot$t0,res.ci$percent[4:5]))
}

bootstrap_corr <- function(data){
  f_corr <- function(x,indices){
    d <- x[indices,]
    return(cor(d)[1,2])
  }
  
  res.boot <- boot(data,statistic = f_corr,R = 1000)
  res.ci <- boot.ci(res.boot,type = "perc")
  return(c(res.boot$t0,res.ci$percent[4:5]))
}

bootstrap_slope <- function(data){
  f_slope <- function(x,indices){
    d <- x[indices,]
    fit <- lm(HFN ~ AGE,d)
    return(summary(fit)$coefficients[2,1])
  }
  
  res.boot <- boot(data,statistic = f_slope,R = 1000)
  res.ci <- boot.ci(res.boot,type = "perc")
  return(c(res.boot$t0,res.ci$percent[4:5]))
}

bootstrap_multi_tissue <- function(data,tissues,var_name,boot_function){
  res <- sapply(tissues,function(x){
    data <- data %>% 
      filter(tissue == x) %>% 
      pull(var_name)
    if(boot_function == "mean"){
      res <- bootstrap_mean(data)
    }else if(boot_function == "median"){
      res <- bootstrap_median(data)
    }
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "tissue")
  return(res)
}

# Main --------------------------------------------------------------------

# chrM gtf
chrM_gtf <- rtracklayer::import("~/data/scientific_project/mitochondrial_genome/annotation/gencode.v30.primary_assembly.chrM.annotation.gtf") %>% 
  as.data.frame() %>%
  filter(type == "gene") %>% 
  mutate(color = ifelse(str_detect(gene_type,"Mt_rRNA"),"#c9a2cf",
                        ifelse(str_detect(gene_type,"Mt_tRNA"),"#fcac55","#9dbfdc")),
         gene_type = ifelse(str_detect(gene_type,"tRNA"),"tRNA",
                            ifelse(str_detect(gene_type,"rRNA"),"rRNA","Protein coding"))) %>% 
  select(seqnames:width,gene_type,color) %>% 
  rbind(data.frame(seqnames = rep("chrM",2),
                   start = c(1,16024),
                   end = c(576,16569),
                   width = c(576,546),
                   gene_type = "D-loop",
                   color = rep("#ffdc57",2)))

# Color information
gtex_color <- fread("~/data/datasets/gtex/v8/gtex_tissue_colors_correct.txt")
gtex_49_tissues <- read_tsv("~/data/scientific_project/mitochondrial_genome/copy_number/RNA_estimate/gtex_49_tissues.list",col_names = F) %>% 
  filter(!X1 %in% c("Cells-Culturedfibroblasts","Cells-EBV-transformedlymphocytes")) %>% 
  pull(X1)
gtex_color <- gtex_color %>% 
  filter(tissue %in% gtex_49_tissues)

# Germ layer definition
endoderm <- c("Colon-Sigmoid","Colon-Transverse","Esophagus-GastroesophagealJunction",
              "Esophagus-Mucosa","Esophagus-Muscularis","Liver","Lung","Pancreas",
              "Stomach","SmallIntestine-TerminalIleum","Thyroid")
mesoderm <- c("Adipose-Subcutaneous","Adipose-Visceral_Omentum","Artery-Aorta","Artery-Coronary","Artery-Tibial",
              "Heart-AtrialAppendage","Heart-LeftVentricle","Kidney-Cortex","Muscle-Skeletal","Ovary","Prostate",
              "Spleen","Testis","Uterus","Vagina","WholeBlood")
ectoderm <- c("AdrenalGland","Brain","Breast-MammaryTissue","MinorSalivaryGland","Nerve-Tibial","Pituitary",
              "Skin-NotSunExposed_Suprapubic","Skin-SunExposed_Lowerleg")
gtex_color <- gtex_color %>% 
  mutate(germ_layers = ifelse(tissue %in% endoderm,"endoderm",
                              ifelse(tissue %in% mesoderm,"mesoderm","ectoderm")))

# Indiv information
Igtex <- fread("~/data/scientific_project/mitochondrial_genome/tables/20220311_Igtex.txt") %>% 
  rename(Indiv = SUBJID) %>% 
  mutate(SEX = ifelse(SEX == 1,"M","F"))
gtex_haplogroup <- fread("~/data/scientific_project/mitochondrial_genome/variant_calling/gtex_838_individuals.haplogroup",select = 1:2) %>% 
  mutate(Macrogroup = ifelse(str_sub(Haplogroup,1,1) == "L",str_sub(Haplogroup,1,2),str_sub(Haplogroup,1,1))) %>% 
  mutate(Macrogroup = ifelse(Macrogroup == "L0","L1",Macrogroup))
Igtex <- Igtex %>% 
  left_join(gtex_haplogroup,by = c("Indiv" = "SampleID"))

# Tissue information
Igtex_T <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230312_Igtex_T.txt") %>% 
  mutate(Indiv = str_extract(SAMPID,"GTEX-\\w+")) %>% 
  rename(tissue_site_detail = SMTSD) %>% 
  filter(SMAFRZE == "RNASEQ") %>% 
  left_join(gtex_color) %>% 
  filter(!is.na(tissue))

# Coverage
filt_coverage_list <- readRDS("20230530_Pmtrna_Igtex_T_coverage_stat.rds")
filt_sample_coverage <- do.call(rbind,lapply(filt_coverage_list,"[[",1))
filt_pos_10bp <- do.call(c,lapply(filt_coverage_list,"[[",2)) %>% unique()
Igtex_T_med_coverage <- do.call(rbind,lapply(filt_coverage_list,"[[",3))
Pmtrna_T_coverage <- do.call(rbind,lapply(filt_coverage_list,"[[",4))

# Contamination
Igtex_T_haplocheck <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230529_Igtex_T_haplocheck.txt")
filt_contamination <- Igtex_T_haplocheck %>% 
  filter(Contamination == "YES") %>% 
  select(Indiv,tissue)

# Sites
NUMT_FP <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230331_NUMT_FP.txt") %>% 
  mutate(ID = str_c(CHROM,POS,REF,ALT,sep = "_"))
TRN <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230401_TRN.bed")
pos_tRNA <- TRN %>% 
  select(V2,V3) %>% 
  mutate(V2 = V2 + 1) %>% 
  apply(1,function(x) seq(x[1],x[2]))
pos_tRNA <- do.call(c,pos_tRNA)
blacklist <- c(66:71,300:316,513:526,3106:3107,12418:12425,16182:16194,566:573,295,2617,13710)
# haplogroup_marker <- fread("~/data/scientific_project/mitochondrial_genome/annotation/MITOMAP/20221204_Markers_Found_at_more_than_80%_in_Macrogroup.txt")

region_gtf <- apply(chrM_gtf,1,function(x){
  seq(x[2],x[3])
})
region_gtf <- do.call(c,region_gtf)
pos_nc <- setdiff(c(1:16569),region_gtf)

# Variant annotation
Vmtrna_Gmt <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230529_Vmtrna_Gmt.txt")
MLC <- read_tsv("~/data/scientific_project/mitochondrial_genome/annotation/MLC_mitochondrial_local_constraint_score.txt")

Vmtrna_Gmt$var_type <- sapply(Vmtrna_Gmt$Consequence,function(x){
  if(str_detect(x,"synonymous|stop_retained")){
    "Synonymous,#2D898B"
  }else if(str_detect(x,"missense")){
    "Missense,#FF7D00"
  }else if(str_detect(x,"intergenic")){
    "D-loop,#ADB5BD"
  }else if(str_detect(x,"non_coding")){
    "rRNA,#5A189A"
  }else if(str_detect(x,"stop_gained|stop_lost|start_lost|frameshift|ablation")){
    "LoF,#E5383B"
  }else{
    "Others,#e8e8e8"
  }
})
Vmtrna_Gmt <- Vmtrna_Gmt %>% 
  mutate(var_color = str_split(var_type,",",2,simplify = T) %>% .[,2]) %>% 
  mutate(var_type = str_split(var_type,",",2,simplify = T) %>% .[,1]) %>% 
  mutate(ID = str_c(CHROM,POS,REF,ALT,sep = "_")) %>% 
  mutate(gnomAD_AF_hom = ifelse(is.na(gnomAD_AF_hom),0,gnomAD_AF_hom),
         gnomAD_AF_het = ifelse(is.na(gnomAD_AF_het),0,gnomAD_AF_het)) %>% 
  mutate(gnomAD_AF_max = ifelse(gnomAD_AF_hom > gnomAD_AF_het,gnomAD_AF_hom,gnomAD_AF_het)) %>% 
  left_join(select(MLC,-Consequence),by = c("POS" = "Position","REF" = "Reference","ALT" = "Alternate")) %>% 
  mutate(var_type = ifelse(var_type == "Missense",
                           ifelse(MLC_score < .7,"Missense\nnon-conserved","Missense\nconserved"),var_type)) %>% 
  mutate(var_type = ifelse(var_type == "rRNA",
                           ifelse(MLC_score < .7,"rRNA\nnon-conserved","rRNA\nconserved"),var_type)) %>% 
  mutate(var_color = ifelse(var_type == "Missense\nnon-conserved","#ffb703",
                            ifelse(var_type == "rRNA\nnon-conserved","#E4C1F9",var_color))) %>% 
  mutate(var_type = factor(var_type,levels = c("D-loop","Synonymous","Missense\nnon-conserved","Missense\nconserved",
                                               "LoF","rRNA\nnon-conserved","rRNA\nconserved","Others")))

mitomap <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230202_MITOMAP_Disease_Mutations.txt")
Vmtrna_Gmt <- Vmtrna_Gmt %>% 
  left_join(mitomap) %>% 
  mutate(mitomap_cfrm = ifelse(str_detect(Status,"Cfrm"),TRUE,FALSE))

# gnomAD
gnomAD <- fread("~/data/scientific_project/mitochondrial_genome/annotation/gnomad.genomes.v3.1.sites.chrM.vcf.vep_distance0.everything.txt")
gnomAD$var_type <- sapply(gnomAD$Consequence,function(x){
  if(str_detect(x,"synonymous|stop_retained")){
    "Synonymous"
  }else if(str_detect(x,"missense")){
    "Missense"
  }else if(str_detect(x,"intergenic")){
    "D-loop"
  }else if(str_detect(x,"stop_gained|stop_lost|start_lost|frameshift|ablation")){
    "LoF"
  }else{
    "Others"
  }
})
gnomAD$var_type <- apply(gnomAD[,c("REF","ALT","SYMBOL","var_type")],1,function(x){
  if(is.na(x[3])){
    x[4]
  }else if(str_detect(x[3],"^MT-T") & str_length(x[1]) == str_length(x[2])){
    "tRNA"
  }else if(str_detect(x[3],"MT-R") & str_length(x[1]) == str_length(x[2])){
    "rRNA"
  }else{
    x[4]
  }
})
gnomAD <- gnomAD %>% 
  mutate(gnomAD_AF_hom = ifelse(is.na(gnomAD_AF_hom),0,gnomAD_AF_hom),
         gnomAD_AF_het = ifelse(is.na(gnomAD_AF_het),0,gnomAD_AF_het)) %>% 
  mutate(gnomAD_AF_max = ifelse(gnomAD_AF_hom > gnomAD_AF_het,gnomAD_AF_hom,gnomAD_AF_het)) %>% 
  left_join(select(MLC,-Consequence),by = c("POS" = "Position","REF" = "Reference","ALT" = "Alternate")) %>% 
  mutate(var_type = ifelse(var_type == "Missense",
                           ifelse(MLC_score < .7,"Missense\nnon-conserved","Missense\nconserved"),var_type)) %>% 
  mutate(var_type = ifelse(var_type == "rRNA",
                           ifelse(MLC_score < .7,"rRNA\nnon-conserved","rRNA\nconserved"),var_type)) %>% 
  mutate(var_type = ifelse(var_type == "tRNA",
                           ifelse(MLC_score < .7,"tRNA\nnon-conserved","tRNA\nconserved"),var_type))

# xCell score
Igtex_T_xCell_score <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230712_Igtex_T_xCell_score_10type.txt") %>% 
  filter(!tissue %in% c("Cells-Culturedfibroblasts","Cells-EBV-transformedlymphocytes"))

# Igtex_T_CPN
Igtex_T_CPN <- fread("~/data/scientific_project/mitochondrial_genome/tables/20221209_Igtex_T_CPN.txt")

# Vmtrna_Ipancreas8_C
Vmtrna_Ipancreas8_C <- fread("~/data/scientific_project/mitochondrial_genome/tables/20231123_Vmtrna_Ipancreas8_C_filtered.txt")

# Vmtrna_Igtex_T table
Vmtrna_Igtex_T_raw <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230530_Vmtrna_Igtex_T_original_fisher_binomp_betabinomp.txt") %>%
  filter(!tissue %in% c("Cells-Culturedfibroblasts","Cells-EBV-transformedlymphocytes"))

# Vmtrna_Igtex_T full table
Vmtrna_Igtex_T_raw_include_indel <- fread("~/data/scientific_project/mitochondrial_genome/tables/20230529_Vmtrna_Igtex_T_original.txt") %>%
  filter(!tissue %in% c("Cells-Culturedfibroblasts","Cells-EBV-transformedlymphocytes"))

# Rename tissue name
gtex_color <- gtex_color %>% 
  mutate(tissue_site_detail = ifelse(str_detect(tissue_site_detail,"BA9|BA24|cervical c-1"),
                                     str_remove(tissue_site_detail," \\(BA9\\)| \\(BA24\\)| \\(cervical c-1\\)"),tissue_site_detail))

library(VGAM)
library(vcfR)
library(data.table)
library(tidyverse)

# GTEx mitochondrial tables
# Vmtrna_Igtex ------------------------------------------------------------

process_single_vcf <- function(file){
  data <- read.vcfR(file)
  data <- vcfR2tidy(data,single_frame = T,info_fields = c("AS_SB_TABLE","MPOS","TLOD"),format_fields = c("AD","AF","DP")) %>%
    .$dat %>%
    mutate(ID = str_c(CHROM,POS,REF,ALT,sep = "_"),
           ref_AD_F = str_split(AS_SB_TABLE,"\\||,",simplify = T) %>%
             .[,1],
           ref_AD_R = str_split(AS_SB_TABLE,"\\||,",simplify = T) %>%
             .[,2],
           alt_AD_F = str_split(AS_SB_TABLE,"\\||,",simplify = T) %>%
             .[,3],
           alt_AD_R = str_split(AS_SB_TABLE,"\\||,",simplify = T) %>%
             .[,4],
           Indiv = str_extract(Indiv,"GTEX-\\w{4,5}"),
           gt_AD_REF = str_split(gt_AD,pattern = ",",simplify = T)[,1],
           gt_AD_ALT = str_split(gt_AD,pattern = ",",simplify = T)[,2],
           tissue = str_split(file,"/",simplify = T) %>%
             .[1,8]) %>%
    # group_by(ID) %>%
    # filter(str_length(REF) == 1,
    #        str_length(ALT) == 1,
    #        min(str_split(alt_AD,",",simplify = T) %>%
    #              as.numeric()) >= 2) %>%
    # ungroup() %>%
    select(POS:ALT,FILTER,Indiv,gt_AF,gt_DP,gt_AD_REF,gt_AD_ALT,ref_AD_F,ref_AD_R,alt_AD_F,alt_AD_R,MPOS,TLOD,tissue)
  return(data)
}

tissues <- read.table("../01_variant_calling/gtex_tissues.list",header = F) %>%
  .$V1
# tissues <- c("Brain-Amygdala","Cells-Culturedfibroblasts","Cells-EBV-transformedlymphocytes","Liver","Testis")
files <- lapply(tissues,function(x){
  dir(str_c("../01_variant_calling/",x,"/output_vcf/"),pattern = "original.vcf$",full.names = T)
})
files <- do.call(c,files)

data_list <- lapply(files,process_single_vcf)
data <- do.call(rbind,data_list)
write_tsv(data,"Vmtrna_Igtex_T_original.txt")

# Fisher  -----------------------------------------------------------------

Vmtrna_Igtex_T_raw <- fread("Vmtrna_Igtex_T_original.txt")
Vmtrna_Igtex_T_raw <- Vmtrna_Igtex_T_raw %>%
  group_by(ID,Indiv,tissue) %>%
  mutate(log_OR_L = fisher.test(matrix(c(ref_AD_F,ref_AD_R,alt_AD_F,alt_AD_R),nrow = 2))[[2]][1] %>%
           log2(),
         log_OR_H = fisher.test(matrix(c(ref_AD_F,ref_AD_R,alt_AD_F,alt_AD_R),nrow = 2))[[2]][2] %>%
           log2()) %>%
  ungroup()
fwrite(Vmtrna_Igtex_T_raw,"Vmtrna_Igtex_T_original_fisher.txt",quote = F,sep = "\t")

# Binomial test -----------------------------------------------------------

calc_sequencing_error_p <- function(data){
  data_sub <- data %>% 
    select(gt_AD_ALT,gt_DP)
  
  data_sub$gt_AD_ALT <- apply(data_sub,1,function(x){
    set.seed(1234)
    res <- ifelse(x[2] >= 1000,rhyper(1,x[1],x[2] - x[1],1000),x[1])
    return(res)
  })
  data_sub <- data_sub %>% 
    mutate(gt_DP = ifelse(gt_DP >= 1000,1000,gt_DP))
  
  data$binom.p <- apply(data_sub,1,function(x){
    res <- binom.test(x[1],x[2],.001,alternative = "greater")
    return(res$p.value)
  })
  
  return(data)
}

Vmtrna_Igtex_T_raw <- fread("Vmtrna_Igtex_T_original_fisher.txt")
Vmtrna_Igtex_T_raw <- calc_sequencing_error_p(Vmtrna_Igtex_T_raw)
fwrite(Vmtrna_Igtex_T_raw,"Vmtrna_Igtex_T_original_fisher_binomp.txt",quote = F,sep = "\t")

# Beta-binomial test ------------------------------------------------------

add_extra_VxIxT <- function(data){
  
  # base substitution type
  data <- data %>% 
    mutate(basesub_type = paste(REF,ALT,sep = ">"))
  data$basesub_type <- sapply(data$basesub_type,function(x){
    case_when(x %in% c("A><*>","T><*>") ~ "T>all",
              x %in% c("C><*>","G><*>") ~ "C>all",
              x %in% c("A>G","T>C") ~ "T>C",
              x %in% c("A>C","T>G") ~ "T>G",
              x %in% c("A>T","T>A") ~ "T>A",
              x %in% c("G>A","C>T") ~ "C>T",
              x %in% c("G>C","C>G") ~ "C>G",
              x %in% c("G>T","C>A") ~ "C>A")
  })
  
  # variant side & minor allele count
  data <- data %>% 
    mutate(var_side = ifelse(gt_AD_REF >= gt_AD_ALT,"ref","alt"),
           minor_count = ifelse(gt_AD_REF >= gt_AD_ALT,gt_AD_ALT,gt_AD_REF),
           major_count = ifelse(gt_AD_REF >= gt_AD_ALT,gt_AD_REF,gt_AD_ALT))
  
  return(data)
}

calc_sequencing_error_p <- function(data,beta_a_b){
  data <- data %>% 
    filter(str_length(REF) == 1 & str_length(ALT) == 1) %>% 
    add_extra_VxIxT()
  
  data_sub <- data %>% 
    inner_join(beta_a_b) %>% 
    select(minor_count:b)
  
  data$beta_binom.p <- apply(data_sub,1,function(x) pzoibetabinom.ab(x[1],x[1] + x[2],x[3],x[4],lower.tail = F))
  
  return(data)
}

beta_a_b <- fread("beta_binomial_a_b_rna.txt")
Vmtrna_Igtex_T_raw <- fread("Vmtrna_Igtex_T_original_fisher_binomp.txt")
Vmtrna_Igtex_T_raw <- calc_sequencing_error_p(Vmtrna_Igtex_T_raw,beta_a_b)
fwrite(Vmtrna_Igtex_T_raw,"Vmtrna_Igtex_T_original_fisher_binomp_betabinomp.txt",quote = F,sep = "\t")

# Beta-binomial test (downsample) -----------------------------------------

calc_sequencing_error_p <- function(data,beta_a_b){
  set.seed(1234)
  data <- data %>% 
    group_by(Indiv,tissue,ID) %>% 
    mutate(gt_DP_down = ifelse(minor_count + major_count >= 1000,1000,gt_DP)) %>% 
    mutate(minor_count_down = ifelse(gt_DP_down == 1000,rhyper(1,minor_count,major_count,1000),minor_count)) %>% 
    mutate(major_count_down = gt_DP_down - minor_count_down) %>% 
    ungroup()
  
  data_sub <- data %>% 
    inner_join(beta_a_b) %>% 
    select(minor_count_down:b)
  
  data$beta_binom.p_down <- apply(data_sub,1,function(x) pzoibetabinom.ab(x[1],x[1] + x[2],x[3],x[4],lower.tail = F))
  return(data)
}

beta_a_b <- fread("beta_binomial_a_b_rna.txt")
Vmtrna_Igtex_T_raw <- fread("Vmtrna_Igtex_T_original_fisher_binomp_betabinomp.txt")
Vmtrna_Igtex_T_raw <- calc_sequencing_error_p(Vmtrna_Igtex_T_raw,beta_a_b)
fwrite(Vmtrna_Igtex_T_raw,"Vmtrna_Igtex_T_original_fisher_binomp_betabinomp_downsample.txt",quote = F,sep = "\t")

# Pmtrna_Igtex_T ----------------------------------------------------------

process_single_coverage <- function(file){
  data <- fread(file)
  data <- data %>%
    select(pos,coverage) %>%
    mutate(Indiv = str_split(file,"/",simplify = T) %>%
             .[1,11] %>%
             str_extract("GTEX-\\w+"),
           tissue = str_split(file,"/",simplify = T) %>%
             .[1,8]) %>%
    rename(POS = pos)
  return(data)
}

tissues <- read.table("../01_variant_calling/gtex_tissues.list",header = F) %>%
  .$V1
for(i in list(1:10,11:20,21:30,31:40,41:49)){
  files <- lapply(tissues[i],function(x){
    dir(str_c("../01_variant_calling/gtex_tissues.list/",x,"/output_coverage/"),pattern = "coverage.tsv$",full.names = T)
  })
  files <- do.call(c,files)

  data_list <- lapply(files,process_single_coverage)
  data <- do.call(rbind,data_list)
  print(i)
  fwrite(data,str_c("Pmtrna_Igtex_T_coverage_",str_c(range(i),collapse = "to"),".txt"),quote = F,sep = "\t")
}

# haplocheck --------------------------------------------------------------

process_single_haplocheck <- function(file){
  data <- read_tsv(file) %>%
    mutate(tissue = str_split(file,"/",simplify = T) %>%
             .[1,8]) %>%
    mutate(Indiv = str_extract(SampleID,"GTEX-[0-9A-Z]{4,5}")) %>% 
    select(-SampleID) %>% 
    select(Indiv,everything())
  return(data)
}

tissues <- read.table("../01_variant_calling/gtex_tissues.list",header = F) %>%
  .$V1
files <- lapply(tissues,function(x){
  dir(str_c("../01_variant_calling/",x,"/"),pattern = "haplocheck.txt",full.names = T)
})
files <- do.call(c,files)

data_list <- lapply(files,process_single_haplocheck)
data <- do.call(rbind,data_list)
write_tsv(data,"Igtex_T_haplocheck.txt")

source("R_scripts_clonal_expansion_paper/formal_analysis/01_read_fixed_files.R")

# Functions ---------------------------------------------------------------

filt_variant_hard_cutoff <- function(data){
  filt_indels <- data %>% 
    filter(!(str_length(REF) == 1 & str_length(ALT) == 1)) %>% 
    select(Indiv,ID,tissue)
  
  filt_multi_allelic <- data %>% 
    group_by(Indiv,tissue,POS) %>% 
    count() %>% 
    filter(n > 1) %>% 
    ungroup() %>% 
    select(Indiv,POS)
  
  return(list(filt_indels = filt_indels,
              filt_multi_allelic = filt_multi_allelic))
}

filt_sample_abnormal_hetero <- function(data,Igtex,Igtex_T){
  data <- data %>% 
    # filter(!(var_side == "ref" & p.adj >= .05),gt_AF < .95)
    filter(p.adj < .05)
  het_count <- data %>% 
    group_by(Indiv,tissue) %>% 
    count()
  covar <- Igtex %>% 
    left_join(Igtex_T) %>% 
    select(Indiv,tissue,SEX,AGE,Macrogroup,SMRIN,SMTSISCH,COHORT) %>% 
    mutate(SMTSISCH = ifelse(is.na(SMTSISCH),0,SMTSISCH))
  data_reg <- het_count %>% 
    left_join(covar)
  res.reg <- MASS::glm.nb(n ~ tissue + AGE + SEX + Macrogroup + SMTSISCH + COHORT + tissue:AGE,data_reg)
  data_reg$residuals <- res.reg$residuals
  
  opar <- par(no.readonly = T)
  par(mfrow = c(2,2))
  plot(res.reg)
  
  abs_resid_max <- data_reg$residuals %>% 
    min() %>% 
    abs()
  data_reg <- data_reg %>% 
    filter(abs(residuals) <= abs_resid_max)
  res.reg <- MASS::glm.nb(n ~ tissue + AGE + SEX + Macrogroup + SMTSISCH + COHORT + tissue:AGE,data_reg)
  
  plot(res.reg)
  par(opar)
  
  remain_sample <- data_reg %>% 
    select(Indiv,tissue)
  return(remain_sample)
}

get_RNA_het_variants <- function(data){
  # Samples
  data <- data %>% 
    anti_join(filt_sample_coverage) %>%
    anti_join(filt_contamination)
  
  # Sites
  data <- data %>% 
    filter(!POS %in% blacklist,
           !POS %in% pos_tRNA,
           !POS %in% filt_pos_10bp,
           !POS %in% NUMT_FP$POS,
           !POS %in% pos_nc)
  
  # Variants
  data <- data %>% 
    mutate(strand_bias = ifelse(gt_AF < .5 & !(log_OR_L < 0 & log_OR_H > 0) & (log_OR_H < -.5 | log_OR_L > .5),T,F)) %>%
    # mutate(strand_bias = ifelse(gt_AF < .01 & !(log_OR_L < 0 & log_OR_H > 0) & (log_OR_H < -1 | log_OR_L > 1),T,F)) %>%
    filter(MPOS >=  4) %>% 
    filter(strand_bias == F)
  
  filt_variant <- filt_variant_hard_cutoff(data)
  data <- data %>% 
    # anti_join(filt_variant$filt_indels) %>% 
    anti_join(filt_variant$filt_multi_allelic)
  
  data <- data %>% 
    mutate(p.adj = p.adjust(beta_binom.p))
  
  # Contamination (negative binomial)
  remain_sample <- filt_sample_abnormal_hetero(data,Igtex,Igtex_T)
  data <- data %>%
    inner_join(remain_sample)
  
  return(data)
}

# Main --------------------------------------------------------------------

Vmtrna_Igtex_T_flag <- get_RNA_het_variants(Vmtrna_Igtex_T_raw)

Vmtrna_Igtex_T <- Vmtrna_Igtex_T_flag %>% 
  filter(!(var_side == "ref" & p.adj >= .05)) %>% 
  mutate(basesub_type = factor(basesub_type,levels = c("T>C","C>T","C>A","C>G","T>G","T>A")))
  # mutate(basesub_type = factor(basesub_type,levels = c("T>C","C>T","C>A","C>G","T>G","T>A"),
  #                              labels = c("A:T>G:C","C:G>T:A","G:C>T:A","G:C>C:G","T:A>G:C","T:A>A:T")))

tissues_per_indiv <- Vmtrna_Igtex_T %>% 
  as_tibble() %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  group_by(Indiv) %>% 
  summarise(n_tissue = n())
indiv_id <- Vmtrna_Igtex_T %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  left_join(tissues_per_indiv) %>% 
  # filter(n < n_tissue)
  filter(n < n_tissue * .8)

Vmtrna_Igtex_T_het_var <- Vmtrna_Igtex_T %>% 
  # filter(gt_AF < .5)
  inner_join(indiv_id) %>% 
  filter(!(Indiv %in% c("GTEX-Y8E5","GTEX-13JUV") & POS == 12354),
         !(Indiv == "GTEX-113IC" & POS == 150))
Vmtrna_Igtex_T <- Vmtrna_Igtex_T %>% 
  filter(!(Indiv %in% c("GTEX-Y8E5","GTEX-13JUV") & POS == 12354),
         !(Indiv == "GTEX-113IC" & POS == 150))

Vmtrna_Igtex_T_HF <- Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,Indiv) %>%
  summarise(HF = sum(gt_AF)) %>% 
  ungroup()

Vmtrna_Igtex_T_HFN <-  Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,Indiv) %>% 
  summarise(HFN = n()) %>% 
  ungroup()

data_reg <- Vmtrna_Igtex_T_HFN %>% 
  left_join(Igtex) %>% 
  left_join(Igtex_T) %>% 
  select(Indiv,tissue,HFN,SEX,AGE,Macrogroup,SMRIN,SMTSISCH,COHORT) %>% 
  mutate(SMTSISCH = ifelse(is.na(SMTSISCH),0,SMTSISCH))

data_reg_HF <- Vmtrna_Igtex_T_HF %>% 
  left_join(Igtex) %>% 
  left_join(Igtex_T) %>% 
  select(Indiv,tissue,HF,SEX,AGE,Macrogroup,SMRIN,SMTSISCH,COHORT) %>% 
  mutate(SMTSISCH = ifelse(is.na(SMTSISCH),0,SMTSISCH))

tissue_n <- Vmtrna_Igtex_T %>% 
  as_tibble() %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  group_by(tissue) %>% 
  summarise(n_tissue = n())

# save.image("20240423.Rdata")

library(VGAM)
library(data.table)
library(tidyverse)

mpileup_Vmtdna_Igtex_T_raw <- fread("~/data/scientific_project/mitochondrial_genome/tables/20220515_mpileup_Vmtdna_Igtex_T.txt")

data <- mpileup_Vmtdna_Igtex_T_raw %>% 
  filter(Indiv == "GTEX-111YS") %>%
  add_extra_VxIxT()

data_sub <- data %>% 
  filter(basesub_type == "C>A",
         var_side == "ref")
reg_res <- vglm(cbind(gt_AD_ALT,gt_AD_REF) ~ 1,betabinomialff,data_sub,trace = T)
a <- Coef(reg_res)[1]
b <- Coef(reg_res)[2]
pzoibetabinom.ab(5,2000,a,b,lower.tail = F)


a_b <- read_rds("20230206_beta_binomial_a_b_dna.rds")

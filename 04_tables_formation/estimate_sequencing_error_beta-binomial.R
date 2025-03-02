library(VGAM)
library(data.table)
library(tidyverse)

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

estimate_beta_binomial_a_b <- function(data){

  # choose max alt counts
  index <- data %>%
    group_by(POS,Indiv,tissue,gt_AD_REF) %>%
    summarise(gt_AD_ALT = max(gt_AD_ALT))
  data <- data %>%
    inner_join(index)

  data <- add_extra_VxIxT(data)

  # sample variants grouped by var_side and basesub_type
  set.seed(1234)
  data <- data %>%
    filter(!is.na(basesub_type)) %>%
    group_by(var_side,basesub_type) %>%
    slice_sample(n = 10000)

  side_vartype_table <- data %>%
    select(var_side,basesub_type) %>%
    unique() %>%
    arrange(var_side,basesub_type) %>%
    filter(!str_detect(basesub_type,"all")) %>%
    na.omit()
  print(side_vartype_table)

  # choose var_side and basesub_type and estimate a,b
  res <- apply(side_vartype_table,1,function(x){
    if(str_starts(x[2],"T|A")){
      data_sub <- data %>%
        filter(var_side == x[1],
               basesub_type %in% c(x[2],"T>all"))
    }else{
      data_sub <- data %>%
        filter(var_side == x[1],
               basesub_type %in% c(x[2],"C>all"))
    }

    reg_res <- vglm(cbind(minor_count,major_count) ~ 1,betabinomialff,data_sub,trace = T)
    a <- Coef(reg_res)[1]
    b <- Coef(reg_res)[2]
    final_res <- tibble(var_side = x[1],
                        basesub_type = x[2],
                        a = a,
                        b = b)
    return(final_res)
  })
  res <- do.call(rbind,res)

  return(res)
}

mpileup_Vmtrna_Igtex_T_36indiv <- fread("mpileup_Vmtrna_Igtex_T_36_subset_sample.txt")

### choose one tissue for each individual
set.seed(1234)
indiv_tissue_36 <- mpileup_Vmtrna_Igtex_T_36indiv %>%
  select(Indiv,tissue) %>%
  unique() %>%
  group_by(Indiv) %>%
  slice_sample()
mpileup_Vmtrna_Igtex_T_36sample <- mpileup_Vmtrna_Igtex_T_36indiv %>%
  inner_join(indiv_tissue_36)

beta_binomial_a_b_rna <- estimate_beta_binomial_a_b(mpileup_Vmtrna_Igtex_T_36sample)
write_tsv(beta_binomial_a_b_rna,"beta_binomial_a_b_rna.txt")

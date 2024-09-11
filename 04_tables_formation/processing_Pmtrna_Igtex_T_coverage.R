library(roll)
library(data.table)
library(tidyverse)

# Low coverage samples, sites and median coverage
filt_sample_hard_cutoff <- function(data){
  filt_sample <- data %>%
    mutate(flag = ifelse(coverage < 200,T,F)) %>%
    group_by(Indiv,tissue) %>%
    summarise(sites_prop = sum(flag) / length(flag)) %>%
    filter(sites_prop >= .1) %>%
    select(Indiv,tissue)
  return(filt_sample)
}

filt_10bp_pos_dp <- function(tissue_name,data){
  # single tissue all individuals median depth
  pos_median_dp <- data %>%
    filter(tissue == tissue_name) %>%
    group_by(POS) %>%
    summarise(depth = median(coverage)) %>%
    mutate(index = 1:nrow(.))

  # filter 10bp window median depth < 200 pos
  pos_median_dp$roll_DP <- roll_median(pos_median_dp[,2,drop = T],width = 10)
  filt_pos <- pos_median_dp %>%
    filter(roll_DP < 200) %>%
    mutate(index = index - 9) %>%
    left_join(pos_median_dp,by = "index") %>%
    select(POS.y,POS.x) %>%
    apply(1,function(x) x[1]:x[2]) %>%
    unlist() %>%
    unique()

  return(filt_pos)
}

filt_10bp_multi_tissue <- function(data,gtex_tissues){
  filt_pos <- sapply(gtex_tissues,filt_10bp_pos_dp,data)
  filt_pos <- do.call(c,filt_pos) %>%
    unique() %>%
    sort()
  return(filt_pos)
}

files <- dir("~/data/scientific_project/mitochondrial_genome/tables/","20230529_Pmtrna_Igtex_T_coverage",full.names = T)
filt_coverage_list <- lapply(files,function(file){
  Pmtrna_Igtex_T_coverage <- fread(file)

  filt_sample <- filt_sample_hard_cutoff(Pmtrna_Igtex_T_coverage)

  gtex_tissues <- unique(Pmtrna_Igtex_T_coverage$tissue)
  pos_10bp_filt <- filt_10bp_multi_tissue(Pmtrna_Igtex_T_coverage %>%
                           anti_join(filt_sample),gtex_tissues)

  med_coverage <- Pmtrna_Igtex_T_coverage %>%
    group_by(Indiv,tissue) %>%
    summarise(med_coverage = median(coverage))

  mean_sd_coverage <- Pmtrna_Igtex_T_coverage %>%
    anti_join(filt_sample) %>%
    group_by(POS,tissue) %>%
    summarise(mean_coverage = mean(coverage),
              sd_coverage = sd(coverage))

  return(list(filt_sample,
              pos_10bp_filt,
              med_coverage,
              mean_sd_coverage))
})
saveRDS(filt_coverage_list,"20230530_Pmtrna_Igtex_T_coverage_stat.rds")

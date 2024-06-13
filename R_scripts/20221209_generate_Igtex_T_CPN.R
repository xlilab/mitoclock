library(tidyverse)

cpn_files <- dir("~/data/scientific_project/mitochondrial_genome/copy_number/RNA_estimate/",pattern = "copy_number.txt",full.names = T)
cpn_full <- data.frame()
for(file in cpn_files){
  data <- read_tsv(file)
  tissue <- str_split(file,"[/|.]",simplify = T) %>% 
    .[1,10]
  data <- data %>% 
    mutate(tissue = str_replace_all(tissue,"\\(.+\\)",""))
  cpn_full <- rbind(cpn_full,data)
}
cpn_full %>% 
  write_tsv("~/data/scientific_project/mitochondrial_genome/tables/20221209_Igtex_T_CPN.txt")

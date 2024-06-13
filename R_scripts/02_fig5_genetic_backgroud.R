
# Functions ---------------------------------------------------------------

plot_no.het_tissue_sex <- function(){
  data_plot <- Vmtrna_Igtex_T_HFN %>% 
    left_join(Igtex) %>% 
    left_join(gtex_color) %>% 
    filter(str_detect(tissue,"Breast|Liver|Muscle")|tissue == "Brain-Cerebellum") %>% 
    mutate(SEX = ifelse(SEX == "F","Female","Male"))
  
  p <- data_plot %>% 
    ggviolin("tissue_site_detail","HFN",fill = "SEX",legend = "right",palette = c("#b5d1ff","#657691"),
             add = "boxplot",ylab = "# Heteroplasmies",legend.title = "")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))+
    stat_compare_means(aes(group = SEX), label = "p.signif",hide.ns = T,color = "red",label.y.npc = .38)
  return(p)
}

plot_no.het_age_breast <- function(sex){
  # 20~50岁分10个bin
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == "Breast-MammaryTissue") %>% 
    left_join(Igtex) %>% 
    filter(SEX == sex) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),
                             labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),ordered_result = T)) %>%
    mutate(HF_interval = fct_relevel(HF_interval,rev),
           AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T)) %>% 
    group_by(Indiv,AGE,HF_interval) %>%
    count() %>% 
    pivot_wider(id_cols = c(Indiv,AGE),names_from = HF_interval,values_from = n) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  
  x <- ifelse(sex == "F",4,3)
  p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                  palette = rev(c("#43AA8B","#90BE6D","#F9C74F","#F8961E")[1:x]),
                  xlab = "Age",ylab = "# Heteroplasmies",legend = "none")+
    labs(fill = "Heteroplasmy",color = "Heteroplasmy")+
    ylim(0,5)+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  
  # 增加总克隆丰度图
  data_plot <- Vmtrna_Igtex_T_HF %>% 
    filter(tissue == "Breast-MammaryTissue") %>% 
    left_join(Igtex) %>% 
    filter(SEX == sex) %>% 
    group_by(Indiv,AGE) %>%
    summarise(HF = sum(HF)) %>% 
    mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T))
  
  data_simple <- sapply(as.character(sort(unique(data_plot$AGE))),function(x){
    data_plot_sub <- data_plot %>% 
      filter(AGE == x) %>% 
      pull(HF)
    res <- bootstrap_median(data_plot_sub)
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "AGE") %>% 
    mutate(AGE = factor(AGE))
  
  p2 <- ggplot(data_simple,aes(AGE,V1))+
    geom_bar(aes(fill = AGE),stat = "identity")+
    geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
    theme_pubr()+
    scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                           "469d89","358f80","248277","14746f","036666")))+
    xlab("Age")+
    ylab("Sum of Heteroplasmy")+
    ylim(0,.15)+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(fill = "none")
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

plot_no.het_age_testis <- function(){
  # 20~50岁分10个bin
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == "Testis") %>% 
    left_join(Igtex) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),
                             labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),ordered_result = T)) %>%
    mutate(HF_interval = fct_relevel(HF_interval,rev),
           AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T)) %>% 
    group_by(Indiv,AGE,HF_interval) %>%
    count() %>% 
    pivot_wider(id_cols = c(Indiv,AGE),names_from = HF_interval,values_from = n) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  
  p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                  palette = rev(c("#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                  xlab = "Age",ylab = "# Heteroplasmies",legend = "right")+
    labs(fill = "Heteroplasmy",color = "Heteroplasmy")+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  
  # 增加总克隆丰度图
  data_plot <- Vmtrna_Igtex_T_HF %>% 
    filter(tissue == "Testis") %>% 
    left_join(Igtex) %>% 
    group_by(Indiv,AGE) %>%
    summarise(HF = sum(HF)) %>% 
    mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T))
  
  data_simple <- sapply(as.character(sort(unique(data_plot$AGE))),function(x){
    data_plot_sub <- data_plot %>% 
      filter(AGE == x) %>% 
      pull(HF)
    res <- bootstrap_median(data_plot_sub)
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "AGE") %>% 
    mutate(AGE = factor(AGE))
  
  p2 <- ggplot(data_simple,aes(AGE,V1))+
    geom_bar(aes(fill = AGE),stat = "identity")+
    geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
    theme_pubr()+
    scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                           "469d89","358f80","248277","14746f","036666")))+
    xlab("Age")+
    ylab("Sum of Heteroplasmy")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(fill = "none")
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

plot_no.het_age_haplogroup <- function(){
  data <- Vmtrna_Igtex_T_HFN %>% 
    filter(tissue == "Skin-SunExposed_Lowerleg") %>%
    left_join(Igtex) %>% 
    mutate(Macrogroup = ifelse(str_detect(Macrogroup,"L"),"L","non L"))
  
  res_lm <- lapply(c("L","non L"),function(x){
    data_sub <- data %>% 
      filter(Macrogroup == x)
    res <- lm(HFN ~ AGE,data_sub) %>% 
      summary()
    return(tibble(tissue = x,
                  intercept = res$coefficients[1,1],
                  beta = res$coefficients[2,1],
                  p = res$coefficients[2,4]))
  })
  res_lm <- do.call(rbind,res_lm)
  formula <- res_lm %>% 
    mutate(across(where(is.numeric),~ round(.,2))) %>% 
    mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
    pull(formula)
  label <- str_c(c("L","Non L"),formula)
  p <- ggplot(data,aes(AGE,HFN,color = Macrogroup,fill = Macrogroup))+
    geom_smooth(method = "lm")+
    scale_fill_manual(values = c("#1f77b4","#ff7f0e")[c(2,1)],labels = label)+
    scale_color_manual(values = c("#1f77b4","#ff7f0e")[c(2,1)],labels = label)+
    xlab("Age")+
    ylab("# Heteroplasmies")+
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL))+
    theme_pubr()+
    theme(legend.position =c(.3,.9),)
  return(p)
}

plot_no.het_haplogroup <- function(){
  Igtex_macrogrp <- Igtex %>% 
    filter(between(AGE,50,70)) %>% 
    mutate(Macrogroup = ifelse(Macrogroup %in% c("B","R","J","T"),"R*",
                               ifelse(Macrogroup %in% c("X","A","N","W","I"),"N*",
                                      ifelse(Macrogroup %in% c("D","M","C","Z","L3"),"L3*",Macrogroup)))) %>%
    mutate(Macrogroup = ifelse(Macrogroup == "V","H",
                               ifelse(Macrogroup == "K","U",
                                      ifelse(Macrogroup %in% c("L1","L2"),"L0/1/2",Macrogroup))))
  Vmtrna_Igtex_T_HFN %>% 
    select(Indiv) %>% 
    unique() %>% 
    inner_join(Igtex_macrogrp) %>% 
    group_by(Macrogroup) %>% 
    count()
  p1 <- Vmtrna_Igtex_T_HFN %>% 
    filter(str_detect(tissue,"Skin")) %>%
    inner_join(Igtex_macrogrp) %>% 
    left_join(gtex_color) %>% 
    ggviolin("Macrogroup","HFN",fill = "tissue_site_detail",legend = "right",palette = c("#0000FF","#7777FF"),
             add = "boxplot",ylab = "# Heteroplasmies",legend.title = "")+
    ylim(0,18)+
    theme(axis.title.x = element_blank(),)+
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank())+
    stat_compare_means(aes(group = tissue), label = "p.signif",hide.ns = T,color = "red",label.y = 17)
  
  nwk <- "(L1:4,(L3:3,(N:2,(R:1,U:1,H:1):1[&&NHX:H=R]):1[&&NHX:H=N]):1[&&NHX:H=L3]);"
  # tree <- read.tree(text = nwk)
  tree <- treeio::read.nhx(textConnection(nwk))
  tree@phylo$tip.label <- c("L0/1/2","L3*","N*","R*","U","H")
  p4 <- ggtree(tree)+
    geom_label(aes(label=H),label.padding = unit(0.15,"lines"))+
    # geom_tiplab()+
    coord_flip()
  
  p <- p1 %>% insert_bottom(p4,height = .4)
  return(p)
}

plot_gwas_example <- function(){
  data <- Vmtrna_Igtex_T_raw_include_indel %>% 
    filter(ID == "chrM_16183_A_AC") %>% 
    group_by(Indiv) %>% 
    summarise(gt_AF = median(gt_AF)) %>% 
    mutate(HF = gt_AF)
  
  gtex_indivs <- read_tsv("~/data/scientific_project/mitochondrial_genome/gwas/20230526_GTEx_838_Indivs.txt",col_names = "Indiv")
  genotype_plog2 <- read_tsv("~/data/scientific_project/mitochondrial_genome/gwas/20230526_GTEx_genotype_chr17:64480334.txt",col_names = "GT")
  genotype <- cbind(gtex_indivs,genotype_plog2) %>%
    group_by(Indiv) %>% 
    mutate(DS = str_extract_all(GT,pattern = "\\d",simplify = T) %>% 
             as.numeric() %>% 
             sum()) %>%
    mutate(GT = cut(DS,c(0,0.1,1,2),include.lowest = T,labels = c("CC","CG","GG"))) %>% 
    ungroup()
  
  data_reg <- data %>% 
    left_join(genotype) %>% 
    left_join(Igtex)
  lm(HF ~ GT,data_reg) %>% 
    summary()
  
  p <- ggplot(data_reg)+
    geom_beeswarm(aes(GT,HF),corral = "wrap",size = .8)+
    geom_smooth(aes(DS+1,HF),method = "lm")+
    ylab("Heteroplasmy of chrM_16183_A_AC")+
    xlab("chr17:64480334\nPOLG2")+
    annotate("text",1,1,label = str_c("β = 0.18, p = 8.6e-04"))+
    theme_pubr()
  return(p)
}

plot_loeuf_decile <- function(){
  # read files
  mito_all <- fread("/picb/lilab5/wangzhenguo/data/scientific_project/mitochondrial_genome/annotation/MitoCarta3.0/Human.MitoCarta3.0_MitoPathways.txt")
  mito <- fread("/picb/lilab5/dongdanyue/dongdanyue/mito/output/MitoCarta_KEGG.txt")
  KEGG <- fread("/picb/lilab5/annotation/KEGG/Table/KEGG_LOEUF_Halflife_chinese.txt")
  LOEUF <- fread("/picb/lilab/annotation/gnomad/gnomAD/r2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz")
  aFC <- fread("~/data/scientific_project/mitochondrial_genome/tables/DongDanyue/GTEx_aFC.txt")
  gtf <- fread("/picb/lilab5/dongdanyue/lilab5/annotation/preprocessing/output/gencode_v26_geneTSS.txt") %>% 
    select(ID,short_name) %>% 
    mutate(ID = str_extract(ID,"[0-9,a-z,A-Z]*"))
  
  # prepare table
  oePLI <- LOEUF %>%
    select(gene,oe_lof,oe_lof_lower,oe_lof_upper,exp_lof) %>%
    filter(!duplicated(gene))
  oePLI_noNA <- oePLI %>%
    filter(!is.na(oe_lof_upper))
  
  aFC <- aFC %>%
    pivot_longer(cols = where(is.numeric),names_to = "tissue",values_to = "aFC")
  aFC_avg <- aFC %>%
    mutate(aFC_abs = abs(aFC)) %>% 
    group_by(gene) %>%
    summarise(aFC_abs = mean(aFC_abs,na.rm = T))
  
  oePLI_noNA <- oePLI_noNA %>% 
    inner_join(aFC_avg)
  
  type1 <- KEGG %>% 
    filter(Merge_TYPE == "DNA_replication_repair" | Merge_TYPE == "RNA_Protein_Process")
  type1.df <- type1 %>% 
    select(Gene,gene_name) %>% 
    unique() 
  
  mito_all_select <- mito_all %>% 
    select(Genes, Ensembl_gene_identifier) %>% 
    rename(Gene = Ensembl_gene_identifier,gene_name = Genes) %>% 
    unique()
  
  type2_1 <- mito_all %>% 
    filter(MitoPathway %in% c("mtDNA maintenance","mtRNA metabolism","Translation"))
  type2_1.df <- type2_1 %>% 
    select(Genes,Ensembl_gene_identifier) %>% 
    rename(Gene = Ensembl_gene_identifier,gene_name = Genes) %>%
    unique()
  
  type2_2 <- mito %>% 
    filter(KEGG_level2 %in% c("Mitochondrial transcription and translation factors",
                              "Mitochondrial DNA replication factors"))
  type2_2.df <- type2_2 %>% 
    select(Symbol) %>% 
    rename(gene_name = Symbol) %>%
    unique() %>% 
    left_join(gtf,by = c("gene_name" = "short_name")) %>% 
    left_join(mito_all_select) %>% 
    unique() %>% 
    mutate(Gene2 = case_when(!is.na(ID) ~ ID,
                             is.na(ID) ~ Gene)) %>% 
    select(gene_name,Gene2) %>%
    rename(Gene = Gene2) 
  
  type3_1 <- mito_all %>% 
    filter(`MitoPathways Hierarchy` %in% c("Mitochondrial dynamics and surveillance > Fission" ,
                                           "Mitochondrial dynamics and surveillance > Fusion",
                                           "Mitochondrial dynamics and surveillance > Mitophagy",
                                           "Mitochondrial dynamics and surveillance > Autophagy"))
  type3_1.df <- type3_1 %>% 
    select(Genes,Ensembl_gene_identifier) %>% 
    rename(Gene = Ensembl_gene_identifier,gene_name = Genes) %>%
    unique()
  
  type3_2 <- mito %>% 
    filter(KEGG_level2 %in% c("Mitochondrial dynamics","Mitophagy factors"))
  type3_2.df <- type3_2 %>% 
    select(Symbol) %>% 
    rename(gene_name = Symbol) %>%
    unique()
  
  t1 <- unique(setdiff(setdiff(type1.df$gene_name,type2_1.df$gene_name),type2_2.df$gene_name))
  t2 <- unique(type2_1.df$gene_name,type2_2.df$gene_name)
  
  selecttype <- data.frame(gene_name = c(t1,t2),
                           type_class = c(rep("Nuclear cellular machinery",length(t1)),
                                          rep("Mitochondrial central dogma",length(t2))))
  
  KEGG_select <- KEGG %>% 
    select(gene_name,gtex_v26_gene_id,Unique_KEGG_type,Inheritance) %>% 
    inner_join(oePLI_noNA,by = c("gene_name" = "gene"))
  
  KEGG_select$Unique_KEGG_type[KEGG_select$Unique_KEGG_type == "Enzyme_energy"] <- "Enzyme"
  KEGG_select$Unique_KEGG_type[KEGG_select$Unique_KEGG_type == "Transporter_energy"] <- "Transporter"
  
  alltype.df <- selecttype %>% 
    inner_join(KEGG_select)
  alltype.df$LOEUF_group <- cut(alltype.df$oe_lof_upper,breaks = round(quantile(oePLI_noNA$oe_lof_upper,seq(0,1,0.1),na.rm = T),
                                                                       digits = 2),
                                include.lowest = T)
  alltype.df$max_aFC_group <- cut(alltype.df$aFC_abs,breaks = round(quantile(oePLI_noNA$aFC_abs,seq(0,1,0.1),na.rm = T),digits = 2),
                                  include.lowest = T)
  LOEUF_type <- table(alltype.df$LOEUF_group)
  LOEUF_type.df <- data.frame(LOEUF_group = names(LOEUF_type),
                              DNA_decile = c(1:10))
  max_aFC_type <- table(alltype.df$max_aFC_group)
  max_aFC_type.df <- data.frame(max_aFC_group = names(max_aFC_type),
                                RNA_decile = c(1:10))
  
  alltype.df <- alltype.df %>% 
    left_join(LOEUF_type.df) %>% 
    left_join(max_aFC_type.df)
  alltype.df$Inheritance[alltype.df$Inheritance %in% c("undefined","other")] <- "Not reported"
  
  ## DNA
  data_up.df <- alltype.df %>% 
    count(type_class,DNA_decile,LOEUF_group,Inheritance) %>%
    group_by(type_class) %>% 
    mutate(freq = n / sum(n)) %>% 
    filter(!is.na(LOEUF_group)) %>% 
    arrange(DNA_decile)
  data_up.df$Inheritance <- factor(data_up.df$Inheritance ,levels = rev(c("Dominant","Ambiguous","Recessive","Not reported")))
  
  p1 <- ggbarplot(data_up.df,"LOEUF_group","n",fill = "Inheritance",color = "Inheritance",
                  facet.by = "type_class",scales = "free_y",legend = "top",
                  palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                  xlab = "LOEUF Decile",ylab = "# Genes")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(),
          axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
  
  ## RNA
  data_up.df <- alltype.df %>% 
    count(type_class,RNA_decile,max_aFC_group,Inheritance) %>%
    group_by(type_class) %>% 
    mutate(freq = n / sum(n)) %>% 
    filter(!is.na(max_aFC_group)) %>% 
    arrange(RNA_decile)
  data_up.df$Inheritance <- factor(data_up.df$Inheritance ,levels = rev(c("Dominant","Ambiguous","Recessive","Not reported")))
  
  p2 <- ggbarplot(data_up.df,"max_aFC_group","n",fill = "Inheritance",color = "Inheritance",
                  facet.by = "type_class",scales = "free_y",legend = "",
                  palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                  xlab = "Mean log2aFC Decile",ylab = "# Genes")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(),
          axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
  
  p <- plot_list(gglist = list(p1,p2),ncol = 1)
  return(p)
}

# Main --------------------------------------------------------------------

plot_no.het_tissue_sex()

ps <- lapply(list("F","M"),plot_no.het_age_breast)
aplot::plot_list(gglist = ps,nrow = 1,widths = c(1,1))

plot_no.het_age_testis()
plot_no.het_age_haplogroup()
plot_no.het_haplogroup()
plot_gwas_example()
plot_loeuf_decile()

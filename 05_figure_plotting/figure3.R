
# Figure 3A ---------------------------------------------------------------

indiv_tissue <- Vmtrna_Igtex_T_het_var %>% 
  as_tibble() %>% 
  select(Indiv,tissue) %>% unique() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE))
data_plot <- Vmtrna_Igtex_T_het_var %>% 
  group_by(Indiv,tissue,basesub_type) %>% 
  summarise(HFN = n()) %>% 
  ungroup() %>% 
  left_join(Igtex) %>% 
  select(Indiv,tissue,basesub_type,HFN,AGE)

basesub_types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
res_lm_list <- lapply(basesub_types,function(x){
  data_sub <- data_plot %>% 
    filter(basesub_type == x) %>% 
    right_join(indiv_tissue) %>% 
    mutate(HFN = ifelse(is.na(HFN),0,HFN),
           basesub_type = x)
  
  ## 各组织截距
  res <- lm(HFN ~ AGE + tissue,data_sub) %>% 
    summary()
  data_intercept <- res$coefficients %>% 
    as_tibble(rownames = "tissue") %>% 
    .[-2,] %>% 
    {.[1,1] <- "Adipose-Subcutaneous";.} %>% 
    mutate(tissue = str_remove_all(tissue,"tissue")) %>% 
    select(tissue,Estimate)
  xx <- data_intercept[1,2,drop = T]
  data_intercept <- data_intercept %>% 
    mutate(Estimate = Estimate + xx) %>% 
    {.[1,2] <- xx;.}
  
  ## 平均斜率和平均截距
  data_sub <- data_sub %>% 
    left_join(data_intercept) %>% 
    mutate(HFN = HFN - Estimate) %>% 
    mutate(HFN = HFN + mean(data_intercept$Estimate))
  res <- lm(HFN ~ AGE,data_sub) %>% 
    summary()
  res_lm <- tibble(basesub_type = x,
                   intercept = res$coefficients[1,1],
                   beta = res$coefficients[2,1],
                   p = res$coefficients[2,4])
  return(list(res_lm = res_lm,
              data_sub = data_sub))
})
res_lm <- tibble()
for(i in 1:length(basesub_types)){
  res_lm <- rbind(res_lm,res_lm_list[[i]][[1]])
}
print(res_lm)

formula <- res_lm %>% 
  mutate(across(where(is.numeric),~ round(.,3))) %>% 
  mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
  pull(formula)
label <- str_c(basesub_types,formula)

data_plot <- tibble()
for(i in 1:length(basesub_types)){
  data_plot <- rbind(data_plot,res_lm_list[[i]][[2]])
}
data_plot %>%
  mutate(basesub_type = factor(basesub_type,levels = c("T>C","C>T","C>A","C>G","T>G","T>A"))) %>% 
  ggplot(aes(AGE,HFN,color = basesub_type,fill = basesub_type))+
  geom_smooth(method = "lm")+
  scale_fill_manual(values = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"),labels = label[c(5,3,1,2,6,4)])+
  scale_color_manual(values = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"),labels = label[c(5,3,1,2,6,4)])+
  xlab("Age")+
  ylab("Number of mutations")+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL))+
  theme_pubr()+
  theme(legend.position = c(.3,.84),)

# Figure 3B ---------------------------------------------------------------

# left
ids <- Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,ID) %>% 
  count() %>% 
  group_by(tissue) %>% 
  arrange(desc(n)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  select(ID) %>% 
  unique() %>% 
  mutate(hotspot = T)

## slope of each variant
conditions <- crossing(ID = ids$ID,tissue = gtex_49_tissues)
p_res <- apply(conditions,1,function(x){
  indiv_tissue <- Vmtrna_Igtex_T_het_var %>% 
    as_tibble() %>% 
    select(Indiv,tissue) %>% 
    unique() %>% 
    filter(tissue == x[2])
  data <- Vmtrna_Igtex_T_het_var %>% 
    filter(ID == x[1]) %>% 
    right_join(indiv_tissue) %>% 
    left_join(Igtex) %>% 
    mutate(HFN = ifelse(is.na(gt_AF),0,gt_AF))
  res <- lm(HFN ~ AGE,data) %>% 
    summary() %>% 
    .$coefficients
  return(tibble(tissue = x[2],
                ID = x[1],
                beta = res[2,1],
                p = res[2,4]))
})
p_res <- do.call(rbind,p_res)
sig_ids <- p_res %>% 
  mutate(q = p.adjust(p,method = "fdr")) %>% 
  filter(q < .05) %>% 
  select(tissue,ID,beta)

data_plot <- sig_ids[c(4:5,7:20,22:30,32:33,35,44:52),] %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," ")) %>% 
  mutate(ID = str_c(ID,var_type,sep = " - ")) %>% 
  mutate(POS = str_extract(ID,"\\d+") %>% as.numeric()) %>% 
  mutate(POS = ifelse(POS >= 16024,POS - 20000,POS)) %>% 
  arrange(POS) %>% 
  left_join(gtex_color %>% 
              select(tissue,tissue_site_detail)) %>% 
  pivot_wider(id_cols = c(tissue_site_detail),names_from = ID,values_from = beta) %>% 
  mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
  as.data.frame() %>% 
  {rownames(.) <- .[,1];.} %>% 
  .[,-1]
Heatmap(data_plot,col = colorRamp2(c(0,1e-4,.001), c("white","#ff9a7d","red")),name = "beta",
        cluster_columns = F,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),
        show_row_dend = F,show_heatmap_legend = T)

# right
tissues <- c("Heart-LeftVentricle","Muscle-Skeletal","Testis")
ids <- c("chrM_16035_G_A","chrM_16327_C_T","chrM_13369_T_C")
conditions <- crossing(tissues,ids) %>% 
  mutate(ids = factor(ids,levels = c("chrM_16035_G_A","chrM_16327_C_T","chrM_13369_T_C"))) %>% 
  arrange(ids)
ps <- apply(conditions,1,function(x){
  indiv_tissue_sub <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == x[1]) %>% 
    select(Indiv,tissue) %>% 
    unique() %>% 
    left_join(Igtex %>% 
                select(Indiv,AGE))
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(ID == x[2],
           tissue == x[1]) %>% 
    right_join(indiv_tissue_sub) %>%
    mutate(HF = ifelse(is.na(gt_AF),0,gt_AF),
           ID = x[2],
           AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T))
  data_simple <- sapply(as.character(sort(unique(data_plot$AGE))),function(x){
    data_plot_sub <- data_plot %>% 
      filter(AGE == x) %>% 
      pull(HF)
    res <- bootstrap_mean(data_plot_sub)
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "AGE")
  
  if(x[2] == "chrM_16035_G_A"){
    p <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      ylim(0,.03)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            legend.position = "none")+
      guides(fill = "none")
  }else if(x[2] == "chrM_16327_C_T"){
    p <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      ylim(0,.013)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            legend.position = "none")+
      guides(fill = "none")
  }else{
    p <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      ylim(0,.013)+
      theme(legend.position = "none",
            axis.title = element_blank())+
      guides(fill = "none")
  }
  return(p)
})
aplot::plot_list(gglist = ps)

# Figure 3C ---------------------------------------------------------------

# row 0
data_slope <- Vmtrna_Igtex_T_het_var %>%
  mutate(nonROS = ifelse(basesub_type %in% c("T>C","C>T"),T,F)) %>%
  group_by(tissue,Indiv) %>%
  summarise(HFN = sum(nonROS)) %>%
  ungroup() %>% 
  left_join(Igtex) %>% 
  select(Indiv,tissue,AGE,HFN)

res_slope <- lapply(gtex_49_tissues,function(x){
  res_slope <- data_slope %>% 
    filter(tissue == x) %>% 
    select(AGE,HFN) %>% 
    bootstrap_slope()
  return(tibble(tissue = x,
                slope = res_slope[1] * 10,
                lower = res_slope[2] * 10,
                upper = res_slope[3] * 10))
})
res_slope <- do.call(rbind,res_slope)
data_plot <- res_slope %>% 
  left_join(gtex_color)
xlab_sort <- data_plot %>% 
  group_by(tissue_site_detail) %>% 
  summarise(slope = median(slope)) %>% 
  arrange(desc(slope)) %>% 
  pull(tissue_site_detail)

p0 <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  geom_hline(yintercept = 0,linetype = 2)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("# mutations per 10y")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = "none")

# row 1
ids <- sig_ids[c(4:5,7:20,22:30,32:33,35,41:42,44:52),] %>% 
  mutate(hotspot = T) %>% 
  select(ID,hotspot) %>% 
  unique() %>% 
  filter(ID != "chrM_64_C_A")

data_slope <- Vmtrna_Igtex_T_het_var %>% 
  left_join(ids) %>% 
  mutate(hotspot = ifelse(is.na(hotspot),F,hotspot)) %>% 
  group_by(Indiv,tissue) %>% 
  summarise(HFN = sum(hotspot)) %>% 
  ungroup() %>% 
  left_join(Igtex) %>% 
  select(Indiv,tissue,AGE,HFN)

res_slope <- lapply(gtex_49_tissues,function(x){
  res_slope <- data_slope %>% 
    filter(tissue == x) %>% 
    select(AGE,HFN) %>% 
    bootstrap_slope()
  return(tibble(tissue = x,
                slope = res_slope[1] * 10,
                lower = res_slope[2] * 10,
                upper = res_slope[3] * 10))
})
res_slope <- do.call(rbind,res_slope)
data_plot <- res_slope %>% 
  left_join(gtex_color)
xlab_sort <- data_plot %>%
  group_by(tissue_site_detail) %>%
  summarise(slope = median(slope)) %>%
  arrange(desc(slope)) %>%
  pull(tissue_site_detail)

p1 <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  geom_hline(yintercept = 0,linetype = 2)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("# mutations per 10y\nhotspots")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = "none")

# row 2
data_slope <- Vmtrna_Igtex_T_het_var %>% 
  left_join(ids) %>% 
  mutate(hotspot = ifelse(is.na(hotspot) & basesub_type %in% c("T>C","C>T"),T,F)) %>% 
  group_by(Indiv,tissue) %>% 
  summarise(HFN = sum(hotspot)) %>% 
  ungroup() %>% 
  left_join(Igtex) %>% 
  select(Indiv,tissue,AGE,HFN)

res_slope <- lapply(gtex_49_tissues,function(x){
  res_slope <- data_slope %>% 
    filter(tissue == x) %>% 
    select(AGE,HFN) %>% 
    bootstrap_slope()
  return(tibble(tissue = x,
                slope = res_slope[1] * 10,
                lower = res_slope[2] * 10,
                upper = res_slope[3] * 10))
})
res_slope <- do.call(rbind,res_slope)
data_plot <- res_slope %>% 
  left_join(gtex_color)
xlab_sort <- data_plot %>% 
  group_by(tissue_site_detail) %>% 
  summarise(slope = median(slope)) %>% 
  arrange(desc(slope)) %>% 
  pull(tissue_site_detail)

p2 <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  geom_hline(yintercept = 0,linetype = 2)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("# mutations per 10y\nsporadic")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = "none")

# row 3
Vmtrna_Igtex_T_HF_nonROS_OldAge <- Vmtrna_Igtex_T_het_var %>% 
  mutate(nonROS = ifelse(basesub_type %in% c("T>C","C>T"),gt_AF,0)) %>% 
  group_by(tissue,Indiv) %>% 
  summarise(HF = sum(nonROS)) %>% 
  ungroup() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70))

data_simple <- bootstrap_multi_tissue(Vmtrna_Igtex_T_HF_nonROS_OldAge,gtex_49_tissues,"HF","median")
data_simple <- data_simple %>% 
  left_join(gtex_color) %>% 
  arrange(tissue_site_detail)

xlab_sort <- data_simple %>% 
  arrange(desc(V1)) %>% 
  pull(tissue_site_detail)
xlab_color <- gtex_color %>% 
  filter(tissue_site_detail %in% xlab_sort) %>% 
  arrange(tissue_site_detail) %>% 
  pull(tissue_color_hex)

p3 <- ggplot(data_simple,aes(tissue_site_detail,V1))+
  geom_bar(aes(fill = tissue_site_detail),stat = "identity")+
  geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
  theme_pubr()+
  scale_fill_manual(values = xlab_color)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("Total clonal abundance")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11))+
  guides(fill = "none")

# row 4
data <- Vmtrna_Igtex_T_het_var %>% 
  filter(basesub_type %in% c("T>C","C>T")) %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70)) %>% 
  mutate(gt_AF = ifelse(is.na(gt_AF),0,gt_AF)) %>% 
  mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),
                           ordered_result = T,include.lowest = T)) %>%
  mutate(HFN = ifelse(gt_AF != 0,1,0)) %>% 
  group_by(Indiv,tissue,HF_interval) %>%
  summarise(n = sum(HFN)) %>% 
  pivot_wider(id_cols = c(Indiv,tissue),names_from = HF_interval,values_from = n) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
  pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
  group_by(tissue,HF_interval) %>% 
  summarise(n = mean(n)) %>% 
  mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1")))) %>% 
  left_join(gtex_color)

p4 <- ggplot(data,aes(tissue_site_detail,n,fill = HF_interval))+
  geom_bar(stat = "identity")+
  ylab("Number of mutations")+
  scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
  labs(fill = "Clonal abundance")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.position = "right")

# row 5
data <- Igtex_T_xCell_score %>% 
  left_join(gtex_color[,c("tissue","tissue_site_detail")]) %>% 
  pivot_longer(names_to = "Celltype",values_to = "score",cols = where(is.numeric)) %>% 
  group_by(tissue_site_detail,Celltype) %>% 
  summarise(score = mean(score)) %>% 
  mutate(Celltype = ifelse(Celltype %in% c("Epithelial cells","Hepatocytes","Keratinocytes","Sebocytes"),"Epithelial cells",Celltype)) %>% 
  mutate(Celltype = ifelse(Celltype == "Mesenchymal stem cells","MSC",Celltype),
         Celltype = factor(Celltype,levels = rev(c("Epithelial cells","MSC","Fibroblasts",
                                                   "Endothelial cells","Adipocytes","Neutrophils","Myocytes")))) %>% 
  group_by(tissue_site_detail,Celltype) %>% 
  summarise(score = max(score))

p5 <- ggplot(data,aes(tissue_site_detail,Celltype,fill = score))+
  geom_tile()+
  scale_fill_gradient(low = "white",high = "red")+
  scale_x_discrete(limits = xlab_sort)+
  # ylab("Cell type")+
  labs(fill = "Cell type\nenrichment")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank(),
        legend.position = "right")

# row 6
load("estgtex.rda")
p6 <- est.gtex$CIB$MB %>% 
  as_tibble(rownames = "SAMPID") %>% 
  left_join(Igtex_T) %>% 
  select(Indiv,tissue_site_detail,Neurons:Microglia) %>% 
  pivot_longer(names_to = "Celltype",values_to = "score",cols = where(is.numeric)) %>% 
  group_by(tissue_site_detail,Celltype) %>% 
  summarise(score = mean(score)) %>% 
  mutate(Celltype = ifelse(Celltype %in% c("Oligodendrocytes","Astrocytes","Microglia"),"Gliocyte",Celltype)) %>% 
  group_by(tissue_site_detail,Celltype) %>% 
  summarise(score = sum(score)) %>% 
  mutate(Celltype = factor(Celltype,levels = rev(c("Gliocyte","Neurons")))) %>% 
  ggplot(aes(tissue_site_detail,Celltype,fill = score))+
  geom_tile()+
  scale_fill_gradient(low = "white",high = "blue")+
  scale_x_discrete(limits = xlab_sort[str_starts(xlab_sort,"Brain")])+
  ylab("")+
  labs(fill = "Fraction")+
  theme_pubr()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60,size = 11),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank(),
        legend.position = "right")

p0 %>% insert_bottom(p1) %>% insert_bottom(p2) %>% insert_bottom(p3) %>% insert_bottom(p4) %>% insert_bottom(p5,height = .8) %>% insert_bottom(p6,height = .25)

# Figure 3D ---------------------------------------------------------------

# 4种组织
plot_het_age_tissue <- function(tissue_var,x){
  ## 变异数量图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == tissue_var) %>%
    left_join(Vmtrna_Gmt) %>% 
    mutate(HFN = ifelse(basesub_type %in% x,T,F)) %>%
    left_join(Igtex) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),ordered_result = T)) %>%
    mutate(HF_interval = fct_relevel(HF_interval,rev),
           AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T)) %>% 
    group_by(Indiv,AGE,HF_interval) %>%
    summarise(n = sum(HFN)) %>% 
    pivot_wider(names_from = HF_interval,values_from = n) %>% 
    mutate(across(where(is.numeric),~ replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  if(all(x == c("T>C","C>T"))){
    p1 <- ggplot(data_plot,aes(AGE,n,fill = HF_interval))+
      geom_bar(stat = "identity")+
      ylab("")+
      xlab("")+
      ylim(0,10)+
      scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
      theme_pubr()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            legend.position = "none")
  }else{
    p1 <- ggplot(data_plot,aes(AGE,n,fill = HF_interval))+
      geom_bar(stat = "identity")+
      ylab("")+
      xlab("")+
      ylim(0,2)+
      scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
      theme_pubr()+
      theme(legend.position = "none",
            axis.title = element_blank())
  }
  
  ## 总克隆丰度图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == tissue_var) %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(gt_AF = ifelse(basesub_type %in% x,gt_AF,0)) %>%
    left_join(Igtex) %>% 
    group_by(Indiv,AGE) %>%
    summarise(HF = sum(gt_AF)) %>% 
    mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T))
  data_simple <- sapply(as.character(sort(unique(data_plot$AGE))),function(x){
    data_plot_sub <- data_plot %>% 
      filter(AGE == x) %>% 
      pull(HF)
    res <- bootstrap_mean(data_plot_sub)
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "AGE")
  
  if(all(x == c("T>C","C>T"))){
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("")+
      ylab("")+
      # ggtitle("non ROS")+
      ylim(0,.18)+
      # ylim(0,.15)+
      # ylim(0,.05)+
      theme(axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = .5))+
      guides(fill = "none")
  }else{
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("")+
      ylab("")+
      # ggtitle("ROS")+
      # ylim(0,.18)+
      # ylim(0,.15)+
      ylim(0,.026)+
      theme(axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = .5))+
      guides(fill = "none")
  }
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

basesub_types <- list(c("T>C","C>T"),c("C>A","C>G"))
conditions <- crossing(gtex_49_tissues[c(27,24,1)],basesub_types)
ps <- map2(conditions$`gtex_49_tissues[c(27, 24, 1)]`,conditions$basesub_types,plot_het_age_tissue)

## 加入脑
plot_het_age_brain <- function(x){
  ## 变异数量图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(str_detect(tissue,"Brain")) %>%
    left_join(Vmtrna_Gmt) %>% 
    mutate(HFN = ifelse(basesub_type %in% x,T,F)) %>%
    left_join(Igtex) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),ordered_result = T)) %>%
    mutate(HF_interval = fct_relevel(HF_interval,rev),
           AGE = cut(AGE,breaks = seq(20,70,5),
                     # labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T)) %>% 
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T)) %>% 
    group_by(Indiv,tissue,AGE,HF_interval) %>%
    summarise(n = sum(HFN)) %>% 
    pivot_wider(names_from = HF_interval,values_from = n) %>% 
    mutate(across(where(is.numeric),~ replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  if(all(x == c("T>C","C>T"))){
    p1 <- ggplot(data_plot,aes(AGE,n,fill = HF_interval))+
      geom_bar(stat = "identity")+
      ylab("")+
      xlab("")+
      ylim(0,10)+
      scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
      theme_pubr()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            legend.position = "none")
  }else{
    p1 <- ggplot(data_plot,aes(AGE,n,fill = HF_interval))+
      geom_bar(stat = "identity")+
      ylab("")+
      xlab("")+
      ylim(0,2)+
      scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
      theme_pubr()+
      theme(legend.position = "none",
            axis.title = element_blank())
  }
  
  ## 总克隆丰度图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(str_detect(tissue,"Brain")) %>%
    left_join(Vmtrna_Gmt) %>% 
    mutate(gt_AF = ifelse(basesub_type %in% x,gt_AF,0)) %>%
    left_join(Igtex) %>% 
    group_by(Indiv,tissue,AGE) %>%
    summarise(HF = sum(gt_AF)) %>% 
    group_by(Indiv,AGE) %>% 
    summarise(HF = median(HF)) %>% 
    mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                     # labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T))
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T))
  data_simple <- sapply(as.character(sort(unique(data_plot$AGE))),function(x){
    data_plot_sub <- data_plot %>% 
      filter(AGE == x) %>% 
      pull(HF)
    res <- bootstrap_mean(data_plot_sub)
    return(res)
  }) %>% 
    t() %>% 
    as_tibble(rownames = "AGE")
  if(all(x == c("T>C","C>T"))){
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("")+
      ylab("")+
      # ggtitle("non ROS")+
      ylim(0,.18)+
      # ylim(0,.15)+
      # ylim(0,.05)+
      theme(axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = .5))+
      guides(fill = "none")
  }else{
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("")+
      ylab("")+
      # ggtitle("ROS")+
      # ylim(0,.18)+
      # ylim(0,.15)+
      ylim(0,.026)+
      theme(axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = .5))+
      guides(fill = "none")
  }
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

ps2 <- lapply(basesub_types,plot_het_age_brain)

ps <- c(ps,ps2)

## 合并non ROS与ROS
pss <- list()
i <- 1
j <- 1
while(j < 8){
  pss[[i]] <- aplot::plot_list(gglist = c(ps[j:(j+1)]),ncol = 1,heights = c(2,2))
  i <- i + 1
  j <- j + 2
}

## 添加组织名称
label_title <- lapply(c(gtex_color$tissue_site_detail[c(1,24,27)],"Brain"),function(x){
  res <- grid.text(x, x = 0.5, y = 0.5, just = "center", gp = gpar(fontsize = 14))
  return(res)
})

## 合并图像与组织
psss <- list()
for(i in 1:4){
  psss[[i]] <- aplot::plot_list(gglist = list(label_title[[i]],pss[[i]]),ncol = 1,heights = c(.045,1))
}

aplot::plot_list(gglist = psss[c(3,2,4,1)],ncol = 4)


# Figure 3E ---------------------------------------------------------------

# row 1
data_plot <- Vmtrna_Igtex_T_het_var %>% 
  mutate(ROS = ifelse(basesub_type %in% c("C>A","C>G"),gt_AF,0)) %>% 
  select(Indiv,tissue,ID,ROS) %>% 
  group_by(Indiv,tissue) %>% 
  summarise(HF = sum(ROS))

data_simple <- bootstrap_multi_tissue(data_plot,gtex_49_tissues,"HF","mean")
data_simple <- data_simple %>% 
  left_join(gtex_color) %>% 
  arrange(tissue_site_detail)

xlab_sort <- data_simple %>% 
  arrange(desc(V1)) %>% 
  pull(tissue_site_detail)
xlab_color <- gtex_color %>% 
  filter(tissue_site_detail %in% xlab_sort) %>% 
  arrange(tissue_site_detail) %>% 
  pull(tissue_color_hex)

p1 <- ggplot(data_simple,aes(tissue_site_detail,V1))+
  geom_bar(aes(fill = tissue_site_detail),stat = "identity")+
  geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
  theme_pubr()+
  scale_fill_manual(values = xlab_color)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("Total clonal abundance")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill = "none")

# row 2
p2 <- Igtex_T_CPN %>% 
  inner_join(gtex_color) %>% 
  ggplot(aes(tissue_site_detail,CPN,color = tissue_site_detail))+
  geom_boxplot()+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  ylab("Relative Copy Number")+
  theme_pubr()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
        axis.title.x = element_blank())+
  guides(color = "none")

# row 3
ids <- Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,ID,POS,basesub_type) %>% 
  count() %>% 
  group_by(tissue) %>% 
  arrange(desc(n)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," ")) %>% 
  mutate(ID = str_c(ID,var_type,sep = "\n")) %>% 
  filter(basesub_type %in% c("C>A","C>G")) %>% 
  select(ID,POS) %>% 
  unique() %>% 
  mutate(hotspot = T) %>% 
  arrange(POS)
p3 <- Vmtrna_Igtex_T_het_var %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," ")) %>% 
  mutate(ID = str_c(ID,var_type,sep = "\n")) %>% 
  mutate(ROS = ifelse(basesub_type %in% c("C>A","C>G"),T,F)) %>% 
  mutate(hotspot = ifelse(ID %in% c(ids$ID),ID,"Others")) %>% 
  group_by(Indiv,tissue,hotspot) %>% 
  summarise(HFN = sum(ROS)) %>% 
  pivot_wider(id_cols = c(Indiv,tissue),names_from = hotspot,values_from = HFN) %>% 
  mutate(across(where(is.numeric),~ replace_na(.,0))) %>% 
  pivot_longer(names_to = "hotspot",values_to = "HFN",cols = where(is.numeric)) %>% 
  group_by(tissue,hotspot) %>% 
  summarise(HFN = mean(HFN)) %>% 
  left_join(gtex_color) %>% 
  mutate(hotspot = factor(hotspot,levels = c(ids$ID,"Others"))) %>% 
  ggplot(aes(tissue_site_detail,HFN,fill = hotspot))+
  geom_bar(stat = "identity")+
  # scale_fill_manual(values = str_c("#",c("89c2d9","61a5c2","468faf","2c7da0","2a6f97","014f86","013a63","b3b3b3")))+
  # scale_fill_manual(values = str_c("#",c("656d4a","a4ac86","c2c5aa","b6ad90","a68a64","936639","7f4f24","b3b3b3")))+
  scale_fill_manual(values = str_c("#",c("89c2d9","468faf","a4ac86","656d4a","a68a64","936639","7f4f24","b3b3b3")))+
  ylab("Number of mutations")+
  labs(fill = "ROS Hotspots")+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

p1 %>% insert_bottom(p2) %>% insert_bottom(p3)



# Figure 4A ---------------------------------------------------------------

vars <- Vmtrna_Igtex_T_het_var %>% 
  as_tibble() %>% 
  select(Indiv,ID) %>% 
  unique()
indivs <- Vmtrna_Igtex_T %>% 
  inner_join(vars) %>% 
  left_join(gtex_color) %>% 
  select(Indiv,tissue,germ_layers) %>% 
  unique() %>% 
  group_by(Indiv,germ_layers) %>% 
  count() %>% 
  mutate(n = ifelse(n <= 3,n,3)) %>% 
  group_by(Indiv) %>% 
  summarise(n = sum(n)) %>% 
  filter(n == 9) %>% 
  pull(Indiv)
n_1tissue <- Vmtrna_Igtex_T %>% 
  inner_join(vars) %>% 
  filter(Indiv %in% indivs) %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  filter(n == 1)

df1 <- Vmtrna_Igtex_T %>% 
  as_tibble() %>% 
  inner_join(vars) %>% 
  filter(Indiv %in% indivs) %>% 
  anti_join(n_1tissue) %>% 
  left_join(gtex_color) %>%
  select(Indiv,ID,germ_layers) %>%
  unique() %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  rbind(n_1tissue %>% 
          mutate(n = 0)) %>% 
  mutate(n_germ = case_when(n == 0 ~ "Only 1 tissue",
                            n == 1 ~ "1 germ layer",
                            n == 2 ~ "2 germ layers",
                            n == 3 ~ "3 germ layers"))
df2 <- Vmtrna_Igtex_T %>% 
  anti_join(Vmtrna_Igtex_T_het_var %>% 
              select(Indiv,ID)) %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  mutate(n_germ = "Germline")

data_plot <- df1 %>% 
  rbind(df2) %>% 
  mutate(n_germ = factor(n_germ,levels = rev(c("Only 1 tissue","1 germ layer","2 germ layers","3 germ layers","Germline")))) %>% 
  left_join(Vmtrna_Gmt) %>% 
  group_by(n_germ,var_type) %>%
  count() %>% 
  group_by(n_germ) %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ungroup()

ggbarplot(data_plot,"n_germ","prop",fill = "var_type",color = "var_type",palette = unique(Vmtrna_Gmt$var_color)[c(1,7,8,9,5,3,4)],
          xlab = "",ylab = "Mutations %",legend = "right",legend.title = "Variant type",
          width = c(rep(1,7),rep(.06,7),rep(.11,7),rep(.05,7),rep(3.12,7)))+
  coord_flip()+
  theme(axis.ticks.y = element_blank())

# Figure 4B ---------------------------------------------------------------

plot_het_age_vartype <- function(tissue_sel){
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == tissue_sel) %>% 
    left_join(Igtex) %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T)) %>% 
    group_by(Indiv,AGE,var_type) %>%
    count() %>% 
    pivot_wider(id_cols = c(Indiv,AGE),names_from = var_type,values_from = n) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "var_type",values_to = "n") %>% 
    group_by(AGE,var_type) %>% 
    summarise(n = mean(n)) %>% 
    mutate(var_type = factor(var_type,levels = c("D-loop","Synonymous","Missense\nnon-conserved","Missense\nconserved",
                                                 "LoF","rRNA\nnon-conserved","rRNA\nconserved","Others")))
  fill_color <- Vmtrna_Gmt %>% 
    select(var_type,var_color) %>% 
    na.omit() %>% 
    unique() %>% 
    arrange(var_type)
  
  p <- ggplot(data_plot,aes(AGE,n,fill = var_type))+
    geom_bar(stat = "identity")+
    ylab("Number of mutations")+
    xlab("Age")+
    # ylim(0,10)+
    scale_fill_manual(values = fill_color$var_color)+
    ggtitle(tissue_sel)+
    theme_pubr()+
    theme(axis.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = .5))
  return(p)
}

ps <- lapply(c("Liver","AdrenalGland","Testis","Heart-LeftVentricle"),plot_het_age_vartype)
aplot::plot_list(gglist = ps,nrow = 2)

# Figure 4C ---------------------------------------------------------------

# row 1
data_plot <- Vmtrna_Igtex_T_het_var %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70)) %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(HF = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),gt_AF,0)) %>%
  group_by(tissue,Indiv) %>%
  summarise(HF = sum(HF))

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
  ylab("Total clonal abundance\ndeleterious mutations at 60y")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill = "none")

# row 2
data_slope <- Vmtrna_Igtex_T_het_var %>%
  left_join(Vmtrna_Gmt) %>%
  mutate(conserve = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),T,F)) %>%
  group_by(tissue,Indiv) %>%
  summarise(HFN = sum(conserve)) %>%
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
  ylab("# deleterious mutations\nper decade")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = "none")

# row 3
data_plot <- Vmtrna_Igtex_T_het_var %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70)) %>% 
  left_join(Vmtrna_Gmt) %>% 
  # mutate(HFN = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),T,F)) %>%
  group_by(Indiv,tissue,var_type) %>%
  # summarise(n = sum(HFN)) %>% 
  count() %>% 
  pivot_wider(id_cols = c(Indiv,tissue),names_from = var_type,values_from = n) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
  pivot_longer(cols = where(is.numeric),names_to = "var_type",values_to = "n") %>% 
  group_by(tissue,var_type) %>% 
  summarise(n = mean(n)) %>% 
  mutate(var_type = factor(var_type,levels = c("D-loop","Synonymous","Missense\nnon-conserved","Missense\nconserved",
                                               "LoF","rRNA\nnon-conserved","rRNA\nconserved","Others"))) %>% 
  left_join(gtex_color)
fill_color <- Vmtrna_Gmt %>% 
  select(var_type,var_color) %>% 
  na.omit() %>% 
  unique() %>% 
  arrange(var_type)
p3 <- ggplot(data_plot,aes(tissue_site_detail,n))+
  geom_bar(aes(fill = var_type),stat = "identity")+
  theme_pubr()+
  scale_fill_manual(values = fill_color$var_color)+
  ylab("Number of mutations\nat 60y")+
  labs(fill = "Variant type")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(),
        legend.position = "right")

# row 4
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

p4 <- ggplot(data,aes(tissue_site_detail,Celltype,fill = score))+
  geom_tile()+
  scale_fill_gradient(low = "white",high = "red")+
  ylab("Cell type")+
  labs(fill = "Cell type\nenrichment")+
  theme_pubr()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.position = "right")

p1 %>% insert_bottom(p2) %>% insert_bottom(p3) %>% insert_bottom(p4)

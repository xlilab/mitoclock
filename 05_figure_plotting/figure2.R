
# Figure 2A ---------------------------------------------------------------

vars <- Vmtrna_Igtex_T_het_var %>% 
  as_tibble() %>% 
  select(Indiv,ID) %>% 
  unique()

## choose individuals >=3 tissues per germ layer
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
indiv_tissue <- Vmtrna_Igtex_T %>% 
  inner_join(vars) %>% 
  filter(Indiv %in% indivs) %>% 
  as_tibble() %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  filter(!str_detect(tissue,"Brain")) %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                   labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T))

df1 <- Vmtrna_Igtex_T %>% 
  inner_join(n_1tissue) %>% 
  filter(!str_detect(tissue,"Brain")) %>% 
  left_join(gtex_color) %>% 
  select(Indiv,ID,germ_layers,tissue) %>% 
  group_by(Indiv,tissue) %>% 
  count() %>% 
  right_join(indiv_tissue) %>% 
  mutate(n = ifelse(is.na(n),0,n)) %>% 
  group_by(AGE,tissue) %>% 
  summarise(n = mean(n)) %>% 
  mutate(tissue = case_when(tissue %in% c("Esophagus-GastroesophagealJunction","Colon-Sigmoid","Colon-Transverse",
                                          "SmallIntestine-TerminalIleum","Pancreas","Vagina") ~ "Other Endoderm",
                            tissue %in% c("Esophagus-Muscularis","Heart-AtrialAppendage","Kidney-Cortex","Spleen",
                                          "Uterus","Ovary") ~ "Other Mesoderm",
                            tissue %in% c("Pituitary","Skin-NotSunExposed_Suprapubic","Breast-MammaryTissue",
                                          "Nerve-Tibial") ~ "Other Ectoderm",
                            tissue %in% c("Adipose-Subcutaneous","Adipose-Visceral_Omentum") ~ "Adipose",
                            tissue %in% c("Artery-Aorta","Artery-Coronary","Artery-Tibial") ~ "Artery",
                            TRUE ~ tissue)) %>% 
  group_by(AGE,tissue) %>% 
  summarise(n_var = mean(n))

df2 <- Vmtrna_Igtex_T %>% 
  as_tibble() %>% 
  inner_join(vars) %>% 
  filter(Indiv %in% indivs) %>% 
  anti_join(n_1tissue) %>% 
  filter(!str_detect(tissue,"Brain")) %>% 
  left_join(gtex_color) %>%
  select(Indiv,ID,germ_layers) %>%
  unique() %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  mutate(germ_layers = case_when(n == 0 ~ "Only 1 Tissue",
                                 n == 1 ~ "1 Germ Layer",
                                 n == 2 ~ "2 Germ Layer",
                                 n == 3 ~ "3 Germ Layer")) %>% 
  mutate(tissue = germ_layers) %>% 
  select(-n) %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  mutate(AGE = cut(AGE,breaks = seq(20,70,5),
                   labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T)) %>% 
  group_by(Indiv,AGE,tissue) %>% 
  count() %>% 
  group_by(AGE,tissue) %>% 
  summarise(n_var = mean(n))

data_plot <- df1 %>% 
  rbind(df2)
data_plot <- data_plot %>%
  left_join(gtex_color %>% 
              select(tissue,tissue_site_detail)) %>% 
  mutate(tissue_site_detail = ifelse(is.na(tissue_site_detail),tissue,tissue_site_detail)) %>% 
  mutate(tissue_site_detail = factor(tissue_site_detail,levels = c("Esophagus - Mucosa","Liver","Lung","Thyroid","Stomach",
                                                                   "Prostate","Other Endoderm",
                                                                   "Heart - Left Ventricle","Testis","Adipose","Muscle - Skeletal",
                                                                   "Artery","Whole Blood","Other Mesoderm",
                                                                   "Adrenal Gland","Skin - Sun Exposed (Lower leg)","Minor Salivary Gland",
                                                                   "Other Ectoderm",
                                                                   "1 germ layer","2 germ layer","3 germ layer")))
stack_color <- gtex_color %>% 
  mutate(tissue_site_detail = factor(tissue_site_detail,levels = levels(data_plot$tissue_site_detail))) %>% 
  arrange(tissue_site_detail) %>% 
  na.omit() %>% 
  pull(tissue_color_hex)
stack_color <- c(stack_color[1:6],"#ffa6c1",stack_color[7:8],"#FF6600",stack_color[9],"#FF5555",stack_color[10],
                 "#a1cca5",stack_color[11:13],"#468faf",str_c("#",c("f7e1d7","dedbd2","b0c4b1")))

ggbarplot(data_plot,"AGE","n_var",fill = "tissue_site_detail",color = "tissue_site_detail",palette = stack_color,
          xlab = "Age",ylab = "Number of mutations",legend = "right",legend.title = "")+
  theme(axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10))

# Figure 2B ---------------------------------------------------------------

data_plot <- Vmtrna_Igtex_T_het_var %>% 
  filter(tissue == "Liver") %>% 
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
                palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                xlab = "Age",ylab = "Number of mutations",legend = "right")+
  labs(fill = "Clonal abundance",color = "Clonal abundance")+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)

data_plot <- Vmtrna_Igtex_T_HF %>% 
  filter(tissue == "Liver") %>% 
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
  ylab("Total clonal abundance")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill = "none")

p2 %>% insert_bottom(p1)

# Figure 2C ---------------------------------------------------------------

data_plot <- Vmtrna_Igtex_T_HFN %>% 
  left_join(Igtex) %>% 
  left_join(gtex_color)

res_lm <- lapply(gtex_49_tissues,function(x){
  data_sub <- data_plot %>% 
    filter(tissue == x)
  res <- lm(HFN ~ AGE,data_sub) %>% 
    summary()
  return(tibble(tissue = x,
                intercept = res$coefficients[1,1],
                beta = res$coefficients[2,1],
                p = res$coefficients[2,4]))
})
res_lm <- do.call(rbind,res_lm) %>% 
  arrange(desc(beta))
print(res_lm)

tissues <- c("ADRNLG","ESPMCS","ESPMSL","HRTLV","SKINNS","SKINS")
data_plot <- data_plot %>% 
  filter(tissue_site_detail_abbr %in% tissues)

formula <- res_lm %>% 
  left_join(gtex_color) %>% 
  filter(tissue_site_detail_abbr %in% tissues) %>% 
  arrange(tissue_site_detail_abbr) %>%
  group_by(tissue) %>% 
  mutate(across(where(is.numeric),~ round(.,2))) %>% 
  mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
  pull(formula)
label <- str_c(tissues,formula)

ggplot(data_plot,aes(AGE,HFN,color = tissue_site_detail_abbr,fill = tissue_site_detail_abbr))+
  geom_smooth(method = "lm")+
  scale_fill_manual(values = c("#33DD33","#552200","#BB9988","#660099","#0000FF","#7777FF"),labels = label)+
  scale_color_manual(values = c("#33DD33","#552200","#BB9988","#660099","#0000FF","#7777FF"),labels = label)+
  xlab("Age")+
  ylab("Number of mutations")+
  coord_cartesian(ylim = c(0,11.5))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL))+
  theme_pubr()+
  # theme(legend.position = "right")
  theme(legend.position =c(.3,.85),)

# Figure 2D ---------------------------------------------------------------

data_slope <- Vmtrna_Igtex_T_HFN %>% 
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
  ylab("# mutations per 10y")+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = "none")

I_T_het_var <- Vmtrna_Igtex_T_het_var %>% 
  as_tibble() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70)) %>% 
  select(Indiv,tissue) %>% 
  unique()

data_simple <- bootstrap_multi_tissue(Vmtrna_Igtex_T_HF %>% 
                                        inner_join(I_T_het_var),gtex_49_tissues,"HF","median")
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

p2 <- ggplot(data_simple,aes(tissue_site_detail,V1))+
  geom_bar(aes(fill = tissue_site_detail),stat = "identity")+
  geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
  theme_pubr()+
  scale_fill_manual(values = xlab_color)+
  scale_x_discrete(limits = xlab_sort)+
  ylab("Total clonal abundance\nat 60y")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11))+
  guides(fill = "none")

data <- Vmtrna_Igtex_T_het_var %>%
  inner_join(I_T_het_var) %>% 
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

p3 <- ggplot(data,aes(tissue_site_detail,n,fill = HF_interval))+
  geom_bar(stat = "identity")+
  ylab("Number of mutations\nat 60y")+
  scale_x_discrete(limits = xlab_sort)+
  scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
  ylim(0,11.5)+
  labs(fill = "Clonal abundance")+
  theme_pubr()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.position = "right")

p1 %>% insert_bottom(p2) %>% insert_bottom(p3)

# Figure 2E ---------------------------------------------------------------

hfn_at_60 <- Vmtrna_Igtex_T_HFN %>% 
  left_join(Igtex) %>% 
  filter(between(AGE,50,70)) %>% 
  group_by(tissue) %>% 
  summarise(HFN = mean(HFN),n_tissue = n())
hf_at_60 <- Vmtrna_Igtex_T_HF %>% 
  left_join(Igtex) %>% 
  filter(between(AGE,50,70)) %>% 
  group_by(tissue) %>% 
  summarise(HF = mean(HF))

res_lm <- lapply(gtex_49_tissues,function(x){
  data_sub <- data_reg %>% 
    filter(tissue == x)
  res <- lm(HFN ~ AGE,data_sub) %>% 
    summary()
  return(tibble(tissue = x,
                intercept = res$coefficients[1,1],
                beta = res$coefficients[2,1],
                p = res$coefficients[2,4]))
})
res_lm <- do.call(rbind,res_lm) %>% 
  mutate(beta = beta * 10)

data_plot <- hfn_at_60 %>% 
  left_join(hf_at_60) %>% 
  left_join(res_lm) %>% 
  left_join(gtex_color)
ggplot(data_plot,aes(beta,HFN))+
  geom_point(aes(size = HF,color = tissue))+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  scale_size_continuous(breaks = c(.02,.05,.10))+
  ylim(0,11.5)+
  guides(color = "none")+
  ylab("Mutations at 60y")+
  xlab("Mutations gained per 10y")+
  labs(size = "Total clonal abundance")+
  theme_pubr()+
  stat_cor()

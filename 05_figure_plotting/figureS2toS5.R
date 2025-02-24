
# Figure S2 ---------------------------------------------------------------

# col 0
data_plot <- Vmtrna_Igtex_T %>% 
  select(tissue) %>% 
  unique() %>% 
  left_join(gtex_color) %>% 
  arrange(desc(tissue))
p0 <- data_plot %>% 
  ggplot(aes(tissue,1,fill = tissue))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = gtex_color$tissue_color_hex)+
  scale_y_continuous(breaks = c(0,1))+
  coord_flip()+
  theme_few()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

# col 1
p1 <- data_plot %>% 
  ggplot(aes(y = tissue,x = 0))+
  geom_text(label = unique(data_plot$tissue_site_detail),hjust = 0)+
  scale_y_discrete(limit = data_plot$tissue)+
  xlim(0,.2)+
  xlab("Tissue")+
  theme_few()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

# col 2
data_plot <- Vmtrna_Igtex_T %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  group_by(tissue) %>% 
  count() %>% 
  mutate(n = as.character(n))
p2 <- data_plot %>% 
  ggplot(aes(y = tissue,x = 0))+
  geom_text(label = data_plot$n)+
  # xlim(0,.2)+
  xlab("# Individuals")+
  theme_few()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

# col 3
data_plot <- Vmtrna_Igtex_T %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  left_join(Igtex_T) %>% 
  mutate(SMTSISCH = SMTSISCH / 60)

p3 <- data_plot %>% 
  ggviolin("tissue","SMTSISCH",fill = "tissue",color = "tissue",
           add = "median",add.params = list(color = "black",size = .25),
           legend = "none",palette = gtex_color$tissue_color_hex,
           ylab = "Ischemic Time (h)",)+
  scale_y_continuous(breaks = c(-6,0,6,12,18,24),limits = c(-6,24))+
  # ylim(-6,24)+
  coord_flip()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        panel.grid.major.x = element_line(colour = "gray",linewidth = .25))

# col 4
data_plot <- Vmtrna_Igtex_T %>% 
  select(Indiv,tissue) %>% 
  unique() %>% 
  left_join(Igtex)
p4 <- data_plot %>% 
  group_by(tissue,COHORT) %>% 
  count() %>% 
  group_by(tissue) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(COHORT = factor(COHORT,levels = rev(c("Postmortem","Organ Donor (OPO)","Surgical")),
                         labels = rev(c("Postmortem","Organ Donor","Surgical")))) %>% 
  ggplot(aes(tissue,prop,fill = COHORT))+
  geom_bar(stat = "identity")+
  ylab("Cohort Ratio")+
  labs(fill = "Cohort")+
  scale_fill_manual(values = c("#ef476f","#ffd166","#06d6a0"))+
  scale_y_continuous(breaks = c(0,1))+
  coord_flip()+
  theme_pubr()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "right")

# col 5
p5 <- data_plot %>% 
  ggviolin("tissue","AGE",fill = "tissue",color = "tissue",
           add = "median",add.params = list(color = "black",size = .25),
           legend = "none",palette = gtex_color$tissue_color_hex,
           ylab = "Age (years)",)+
  scale_y_continuous(breaks = seq(20,70,10),limits = c(20,70))+
  # ylim(20,70)+
  coord_flip()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        panel.grid.major.x = element_line(colour = "gray",linewidth = .25))

# col 6
p6 <- data_plot %>% 
  group_by(tissue,SEX) %>% 
  count() %>% 
  group_by(tissue) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(SEX = factor(SEX,labels = rev(c("Male","Female")),levels = rev(c("M","F")))) %>% 
  ggplot(aes(tissue,prop,fill = SEX))+
  geom_bar(stat = "identity")+
  ylab("Sex Ratio")+
  labs(fill = "Sex")+
  scale_fill_manual(values = rev(c("#267abf","#bf2626")))+
  scale_y_continuous(breaks = c(0,1))+
  coord_flip()+
  theme_pubr()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 10))

p1 %>% insert_left(p0,width = .05) %>% insert_right(p2,width = .4) %>% insert_right(p5,width = .4) %>% insert_right(p6,width = .4) %>% insert_right(p3,width = .4) %>% insert_right(p4,width = .4) # 12*10

# Figure S3 ---------------------------------------------------------------

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

# plot preparation
data <- Vmtrna_Igtex_T %>% 
  as_tibble() %>% 
  inner_join(vars) %>% 
  filter(Indiv %in% indivs)
Vmtrna_Igtex_N <- data %>% 
  group_by(Indiv,ID) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(quartile = cut(n,breaks = c(0,1,2,3,4,5,6,7,34),
                        labels = c("1","2","3","4","5","6","7",">=8"),include.lowest = T))
data <- data %>% 
  left_join(gtex_color) %>%
  select(Indiv,ID,germ_layers) %>%
  unique()

## reality
df1 <- data %>% 
  inner_join(Vmtrna_Igtex_N) %>% 
  select(Indiv,ID,germ_layers,quartile) %>% 
  unique() %>% 
  group_by(Indiv,ID,quartile) %>% 
  summarise(n_layer = n()) %>% 
  group_by(quartile,n_layer) %>% 
  count() %>% 
  group_by(quartile) %>% 
  mutate(percent = n / sum(n) * 100,
         n_layer = factor(n_layer,levels = c(1,2,3),labels = c("1 germ layer","2 germ layers","3 germ layers"))) %>% 
  mutate(n_layer = fct_rev(n_layer)) %>%
  mutate(source = "GTEx")

## permutation
# permute_list <- lapply(1:1000,function(x){
#   set.seed(x)
#   data <- Vmtrna_Igtex_T %>% 
#     as_tibble() %>% 
#     inner_join(vars) %>% 
#     filter(Indiv %in% indivs)
#   permute_label <- data %>%
#     select(Indiv,tissue) %>%
#     unique() %>%
#     group_by(Indiv) %>%
#     mutate(permute = sample(tissue))
#   data <- data %>%
#     left_join(permute_label) %>%
#     inner_join(Vmtrna_Igtex_N) %>%
#     inner_join(gtex_color,by = c("permute" = "tissue")) %>%
#     group_by(Indiv,ID,germ_layers) %>%
#     summarise(n_layer = n())
#   data_permute <- data %>%
#     inner_join(Vmtrna_Igtex_N) %>%
#     mutate(source = ifelse(n_layer == n,T,F)) %>%
#     select(Indiv,ID,quartile,source) %>%
#     unique() %>%
#     group_by(quartile,source) %>%
#     count() %>%
#     group_by(quartile) %>%
#     mutate(n_total = sum(n)) %>%
#     filter(source == F) %>%
#     mutate(n_same = n_total - n) %>%
#     select(-c(source,n)) %>%
#     mutate(percent = n_same / n_total,)
#   # source = "GTEx permutation")
#   return(data_permute)
# })
# permute_res <- do.call(rbind,permute_list) %>%
#   mutate(source = "permutation")
# saveRDS(permute_res,"20240407_germ_layer_permute1000_results_merge_tissue_1to7and8+.rds")

# permute_list <- lapply(1:1000,function(x){
#   set.seed(x)
#   data <- Vmtrna_Igtex_T %>%
#     as_tibble() %>%
#     inner_join(vars) %>%
#     filter(Indiv %in% indivs)
#   permute_label <- data %>%
#     select(Indiv,tissue) %>%
#     unique() %>%
#     group_by(Indiv) %>%
#     mutate(permute = sample(tissue))
#   data <- data %>%
#     left_join(permute_label) %>%
#     inner_join(Vmtrna_Igtex_N) %>%
#     inner_join(gtex_color,by = c("permute" = "tissue")) %>%
#     group_by(Indiv,ID,germ_layers) %>%
#     summarise(n_layer = n())
#   data_permute <- data %>%
#     inner_join(Vmtrna_Igtex_N) %>%
#     select(Indiv,ID,germ_layers,quartile) %>% 
#     unique() %>% 
#     group_by(Indiv,ID,quartile) %>% 
#     summarise(n_layer = n()) %>% 
#     group_by(quartile,n_layer) %>% 
#     count() %>% 
#     filter(!quartile %in% c(1,2)) %>% 
#     mutate(n_layer_m = ifelse(n_layer %in% c(1,2),"1~2","3")) %>% 
#     group_by(quartile,n_layer_m) %>% 
#     summarise(n = sum(n)) %>% 
#     group_by(quartile) %>% 
#     mutate(percent = n / sum(n) * 100,
#            n_layer_m = factor(n_layer_m)) %>% 
#     filter(n_layer_m != 3)
#     # mutate(n_layer = fct_rev(n_layer))
#   return(data_permute)
# })
# permute_res <- do.call(rbind,permute_list) %>%
#   mutate(source = "permutation")
# saveRDS(permute_res,"20240415_germ_layer_permute1000_results_merge_tissue_1to7and8+_gt2.rds")

df2 <- readRDS("20240407_germ_layer_permute1000_results_merge_tissue_1to7and8+.rds") %>% 
  mutate(percent = percent * 100)
df3 <- readRDS("20240415_germ_layer_permute1000_results_merge_tissue_1to7and8+_gt2.rds")

ggplot()+
  geom_bar(data = df1,aes(x = quartile,y = percent,fill = n_layer),stat = "identity")+
  geom_violin(data = df2,aes(x = quartile,y = percent,width = .75),color = "black",fill = "gray",scale = "width")+
  geom_violin(data = df3,aes(x = quartile,y = percent,width = .75),color = "black",fill = "gray",scale = "width")+
  xlab("# Tissues")+
  ylab("% Variants from 1/2/3 Germ Layer")+
  labs(fill = "# Germ layers")+
  scale_fill_manual(values = rev(str_c("#",c("f7e1d7","dedbd2","b0c4b1"))))+
  # scale_fill_brewer(palette = "Set2")+
  guides(fill = guide_legend(reverse = T))+
  theme_pubr()+
  theme(legend.position = "right")

# Figure S4 ---------------------------------------------------------------

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
ids <- Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,ID,POS,basesub_type) %>% 
  count() %>% 
  group_by(tissue) %>% 
  arrange(desc(n)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  filter(POS != 10380) %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," ")) %>% 
  mutate(ID = str_c(ID,var_type,sep = "\n")) %>% 
  filter(basesub_type %in% c("C>A","C>G")) %>% 
  select(ID,POS) %>% 
  unique() %>% 
  mutate(hotspot = T) %>% 
  arrange(POS)
p2 <- Vmtrna_Igtex_T_het_var %>% 
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
  scale_fill_manual(values = str_c("#",c("89c2d9","468faf","a4ac86","656d4a","a68a64","7f4f24","b3b3b3")))+
  ylab("Number of mutations")+
  labs(fill = "ROS Hotspots")+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

# row 3
p3 <- Igtex_T_CPN %>% 
  inner_join(gtex_color) %>% 
  ggplot(aes(tissue_site_detail,CPN,color = tissue_site_detail))+
  geom_boxplot()+
  scale_color_manual(values = gtex_color$tissue_color_hex)+
  ylab("Relative Copy Number")+
  theme_pubr()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
        axis.title.x = element_blank())+
  guides(color = "none")

p1 %>% insert_bottom(p2) %>% insert_bottom(p3)


# Figure S5 ---------------------------------------------------------------

data_plot <- Vmtrna_Igtex_T %>% 
  group_by(ID) %>% 
  summarise(max_hl = max(gt_AF)) %>% 
  left_join(Vmtrna_Gmt) %>% 
  mutate(var_type = str_remove_all(var_type,"\nconserved|\nnon-conserved")) %>% 
  mutate(people_freq = cut(gnomAD_AF_hom,c(-1,1e-10,1e-3,1e-2,1e-1,1),
                           labels = c("Not reported","<0.1%","0.1~1%","1~10%",">10%")),
         heteroplasmy = cut(max_hl,c(0,.005,.01,.05,.1,.9,1),
                            labels = c("<0.005","0.005~0.01","0.01~0.05","0.05~0.1","0.1~0.9",">0.9"),ordered_result = T)) %>% 
  mutate(heteroplasmy = fct_relevel(heteroplasmy,rev)) %>% 
  group_by(people_freq,heteroplasmy,var_type) %>% 
  count() %>% 
  group_by(var_type) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(var_type = factor(var_type,levels = c("D-loop","rRNA","Synonymous","Missense","LoF")))

ggbarplot(data_plot,"people_freq","percent",fill = "heteroplasmy",color = "heteroplasmy",
          facet.by = "var_type",legend = "right",
          palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E","#F3722C")))+
  coord_flip()+
  xlab("Population allele frequency")+
  ylab("Proportion of variants")+
  facet_wrap(~ var_type,ncol = 2)+
  labs(fill = "Clonal abundance\nBody-wide maximum",
       color = "Clonal abundance\nBody-wide maximum")+
  guides(fill = guide_legend(reverse = F),
         color = guide_legend(reverse = F))+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(),
        axis.text = element_text(size = 10)) #5*5

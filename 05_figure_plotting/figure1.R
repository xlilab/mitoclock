
# Figure 1B ----------------------------------------------------------------------

example <- tibble(Indiv = c("GTEX-11ZTT","GTEX-WY7C","GTEX-13VXU","GTEX-1C64N",
                            "GTEX-1GZ4H","GTEX-QV44","GTEX-1F6I4","GTEX-15CHQ","GTEX-1F52S"),
                  POS = c(16428,12389,9209,228,
                          3392,204,76,11778,7445))
tissues <- Vmtrna_Igtex_T %>% 
  filter(Indiv %in% example$Indiv) %>% 
  select(Indiv,tissue) %>% 
  unique()
ids <- Vmtrna_Igtex_T %>% 
  inner_join(example) %>% 
  select(Indiv,ID,POS) %>% 
  unique()

data_plot <- Vmtrna_Igtex_T %>% 
  inner_join(example) %>% 
  right_join(tissues) %>% 
  select(tissue,Indiv,gt_AF) %>% 
  mutate(gt_AF = ifelse(is.na(gt_AF),0,gt_AF)) %>% 
  left_join(ids) %>% 
  arrange(POS,Indiv) %>% 
  left_join(Igtex) %>% 
  # mutate(indiv_ID = str_c(str_pad(Indiv,width = 10,side = "right",pad = "_"),str_pad(ID,width = 14,side = "left",pad = "_"),sep = "_"))
  mutate(indiv_ID = str_c(Indiv,ID,sep = ": "))

data_color <- gtex_color %>% 
  right_join(data_plot) %>% 
  select(tissue,tissue_color_hex) %>% 
  unique()
var_color <- data_plot %>% 
  select(ID,indiv_ID) %>% 
  unique() %>% 
  left_join(Vmtrna_Gmt) %>% 
  select(indiv_ID,var_type,var_color) %>% 
  unique() %>% 
  mutate(var_type = as.character(var_type)) %>% 
  mutate(var_type = ifelse(str_detect(indiv_ID,"7445"),"Stop retained, SNHL",
                           ifelse(str_detect(indiv_ID,"11778"),"Missense\nconserved, LHON",var_type))) %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," "))

p1 <- ggplot(data_plot,aes(indiv_ID,gt_AF,color = tissue))+
  annotate("rect",xmin = -Inf,xmax = 1.5,ymin = -Inf,ymax = Inf,alpha = 0.3,fill = "gray")+
  annotate("rect",xmin = 2.5,xmax = 3.5,ymin = -Inf,ymax = Inf,alpha = 0.3,fill = "gray")+
  annotate("rect",xmin = 4.5,xmax = 5.5,ymin = -Inf,ymax = Inf,alpha = 0.3,fill = "gray")+
  annotate("rect",xmin = 6.5,xmax = 7.5,ymin = -Inf,ymax = Inf,alpha = 0.3,fill = "gray")+
  annotate("rect",xmin = 8.5,xmax = 9.6,ymin = -Inf,ymax = Inf,alpha = 0.3,fill = "gray")+
  geom_beeswarm(corral = "wrap",method = "center")+
  guides(color = F)+
  scale_color_manual(values = data_color$tissue_color_hex)+
  scale_x_discrete(limits = unique(data_plot$indiv_ID)[c(4,6,7,1,2,5,3,8,9)])+
  xlab("")+
  ylab("Clonal abundance")+
  theme_few()+
  theme(axis.text.y = element_text(color = var_color$var_color[c(4,6,7,1,2,5,3,8,9)]))+
  # axis.title.x = element_text(size = 14))+
  coord_flip()
p2 <- ggplot(var_color,aes(y = indiv_ID,x = 0))+
  geom_text(label = var_color$var_type,color = var_color$var_color,hjust = 0)+
  xlim(0,.2)+
  xlab("")+
  theme_few()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
p1 %>% insert_right(p2,width = .5)

# Figure 1D ---------------------------------------------------------------

ids <- c("chrM_9412_G_A","chrM_4061_C_A","chrM_4175_G_A","chrM_15280_C_T","chrM_2456_T_A","chrM_6269_A_G","chrM_2702_G_A")
data <- Vmtrna_Ipancreas8_C %>% 
  filter(AGE == 54,
         ID %in% ids)
indivs <- Vmtrna_Ipancreas8_C %>% 
  filter(AGE == 54) %>% 
  pull(Indiv) %>% 
  unique()
data_plus <- expand_grid(ID = ids,Indiv = indivs[!indivs %in% unique(data$Indiv)]) %>% 
  mutate(gt_AF = 0)
data <- data %>% 
  select(ID,Indiv,gt_AF) %>% 
  rbind(data_plus)

y_sort <- data %>% 
  group_by(ID) %>% 
  summarise(HF = sum(gt_AF) / nrow(data)) %>% 
  arrange(desc(HF)) %>% 
  mutate(POS = str_extract(ID,"\\d+") %>% 
           as.numeric())
x_sort <- data %>% 
  filter(gt_AF != 0) %>% 
  arrange(factor(ID,levels = y_sort$ID),desc(gt_AF)) %>% 
  select(Indiv) %>% 
  unique()
x_sort_plus <- data %>% 
  filter(gt_AF == 0) %>% 
  select(Indiv) %>% 
  unique()
x_sort <- rbind(x_sort,x_sort_plus)
y_col <- y_sort %>% 
  left_join(Vmtrna_Gmt) %>% 
  select(ID,var_type,var_color) %>% 
  unique() %>% 
  mutate(var_type = str_replace_all(var_type,"\n"," "))
data_cluster <- data %>% 
  select(ID,Indiv,gt_AF) %>% 
  # mutate(gt_AF = ifelse(gt_AF >= .01,gt_AF,0)) %>% 
  pivot_wider(names_from = Indiv,values_from = gt_AF) %>% 
  mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
  as.data.frame() %>% 
  {rownames(.) <- .[,1];.} %>% 
  .[,-1]
ha <- rowAnnotation(Bulk = anno_numeric(y_sort$HF %>% round(4)),
                    annotation_label = "",annotation_name_rot = 0,
                    var_type = anno_text(y_col$var_type,gp = gpar(col = y_col$var_color)))

p <- data_cluster[y_sort$ID,x_sort$Indiv] %>%
  Heatmap(cluster_rows = F,cluster_columns = F,border = T,
          show_column_names = F,name = "Clonal abundance",
          column_title = "Individual cells",column_title_side = "bottom",
          show_heatmap_legend = F,
          col = colorRamp2(c(0,.2), c("white", "red")),
          right_annotation = ha,
          row_names_gp = gpar(col = y_col$var_color),row_names_side = "left")
lgd <- Legend(col_fun = colorRamp2(c(0,.2), c("white", "red")),
              at = c(0,0.2),title = "Clonal abundance",
              direction = "horizontal",
              legend_width = unit(2,"cm"),
              title_position = "topcenter")
print(p)
draw(lgd,x = unit(1, "npc") - unit(0.4, "npc"), y = unit(1, "npc") - unit(0.15, "npc"))

# Figure 1E ---------------------------------------------------------------

# left
pos_exd_n <- c(blacklist,pos_tRNA,filt_pos_10bp,NUMT_FP$POS,pos_nc) %>% 
  sort() %>% 
  unique() %>% 
  length()
pos_var_n <- Vmtrna_Igtex_T %>% 
  pull(POS) %>% 
  unique() %>% 
  length()
# pos_var_n / (16569 - pos_exd_n)
tibble(type = c("Variant","No variant","Excluded"),
       n = c(pos_var_n,16569 - pos_var_n - pos_exd_n,pos_exd_n)) %>% 
  # mutate(percent = round(n /sum(n) * 100,1)) %>% 
  # mutate(percent_label = str_c(percent,"%")) %>% 
  ggpie("n",fill = "type",label = "n",color = "white",
        legend = "right",lab.pos = "in",legend.title = "16569 mtDNA sites",
        lab.font = c(4,"plain","white"),palette = c("#989898","#457b9d","#1d3557"))+
  font("legend.title",size = 14)

# right
exd_pos <- c(blacklist,pos_tRNA,filt_pos_10bp,NUMT_FP$POS,pos_nc) %>% 
  sort() %>% 
  unique()
p1 <- Pmtrna_T_coverage %>% 
  filter(tissue %in% c("WholeBlood","Liver","Esophagus-Mucosa")) %>% 
  left_join(gtex_color) %>% 
  mutate(tissue_site_detail = fct_relevel(tissue_site_detail,"Liver",after = 0)) %>% 
  ggplot(aes(POS,mean_coverage,color = tissue_site_detail))+
  geom_line()+
  scale_y_log10()+
  theme_pubr()+
  xlim(0,16569)+
  ylab("RNA coverage")+
  scale_color_manual(values = rev(c("#FF00BB","#552200","#AABB66")))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")

p2 <- chrM_gtf %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = gene_type))+
  geom_rect()+
  scale_fill_manual(values = unique(chrM_gtf$color)[c(4,3,2,1)])+
  xlab("POS")+
  theme_pubr()+
  xlim(1,16569)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")

discontinuity <- which(diff(exd_pos) != 1)
p3 <- data.frame(start = exd_pos[c(1,discontinuity + 1)], 
                 end = exd_pos[c(discontinuity,length(exd_pos))]) %>% 
  mutate(col = "Excluded position") %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1,fill = col))+
  geom_rect()+
  scale_fill_manual(values = "gray")+
  xlab("chrM coordinate")+
  theme_pubr()+
  xlim(1,16569)+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")

# aplot::plot_list(p1,p2,p3,ncol = 1,heights = c(8,1,1))
p1 %>% insert_bottom(p2,height = .1) %>% insert_bottom(p3,height = .1)

# Figure 1F ---------------------------------------------------------------

# left
Vmtrna_Igtex_T %>% 
  select(ID,basesub_type) %>% 
  unique() %>% 
  group_by(basesub_type) %>%
  count() %>%
  ungroup() %>%
  mutate(percent = n / sum(n)) %>%
  mutate(percent = round(n /sum(n) * 100,1)) %>%
  mutate(percent_label = str_c(percent,"%")) %>%
  ggpie("percent",fill = "basesub_type",label = "percent_label",color = "white",
        legend = "right",lab.pos = "in",legend.title = "Base substitution",
        lab.font = c(4,"plain","white"),palette = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"))+
  font("legend.title",size = 14)

# right
Vmtrna_Igtex_T %>% 
  left_join(Vmtrna_Gmt) %>% 
  select(ID,basesub_type,var_type) %>% 
  unique() %>% 
  group_by(basesub_type,var_type) %>%
  count() %>%
  # ungroup() %>%
  group_by(basesub_type) %>% 
  mutate(percent = n / sum(n) * 100) %>%
  mutate(basesub_type = fct_relevel(basesub_type,rev)) %>% 
  ggbarplot("basesub_type","percent",fill = "var_type",color = "var_type",
            ylab = "Variants %",xlab = "Base substitution",legend.title = "Variant type",
            legend = "right",
            palette = unique(Vmtrna_Gmt$var_color)[c(1,7,8,9,5,3,4)])+
  coord_flip()

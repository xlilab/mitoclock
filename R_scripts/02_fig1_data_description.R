
# Functions ---------------------------------------------------------------

plot_coverage_coordinate <- function(){
  # exd_pos <- c(blacklist,pos_tRNA,filt_pos_10bp,NUMT_FP$POS,pos_nc) %>% 
  exd_pos <- c(blacklist,pos_tRNA,filt_pos_10bp,pos_nc) %>% 
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
    ylab("Mean Coverage")+
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
    mutate(col = "Excluded Position") %>% 
    ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1,fill = col))+
    geom_rect()+
    scale_fill_manual(values = "gray")+
    xlab("chrM Coordinate")+
    theme_pubr()+
    xlim(1,16569)+
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "right")
  
  p <- p1 %>% insert_bottom(p2,height = .1) %>% insert_bottom(p3,height = .1)
  return(p)
}

plot_all_sites_pie <- function(){
  pos_exd_n <- c(blacklist,pos_tRNA,filt_pos_10bp,NUMT_FP$POS,pos_nc) %>% 
    sort() %>% 
    unique() %>% 
    length()
  pos_var_n <- Vmtrna_Igtex_T %>% 
    pull(POS) %>% 
    unique() %>% 
    length()
  p <- tibble(type = c("Variant","No variant","Excluded"),
              n = c(pos_var_n,16569 - pos_var_n - pos_exd_n,pos_exd_n)) %>% 
    mutate(percent = round(n /sum(n) * 100,1)) %>% 
    mutate(percent_label = str_c(percent,"%")) %>% 
    ggpie("percent",fill = "type",label = "percent_label",color = "white",
          legend = "right",lab.pos = "in",legend.title = "16569 mtDNA bases",
          lab.font = c(4,"plain","white"),palette = c("#989898","#457b9d","#1d3557"))+
    font("legend.title",size = 14)
  return(p)
}

plot_basesub_pie <- function(){
  p <- Vmtrna_Igtex_T %>% 
    select(ID,basesub_type) %>% 
    unique() %>% 
    group_by(basesub_type) %>%
    count() %>%
    ungroup() %>%
    mutate(percent = n / sum(n)) %>%
    mutate(percent = round(n /sum(n) * 100,1)) %>%
    mutate(percent_label = str_c(percent,"%")) %>%
    ggpie("percent",fill = "basesub_type",label = "percent_label",color = "white",
          legend = "right",lab.pos = "in",legend.title = "5851 unique variants",
          lab.font = c(4,"plain","white"),palette = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"))+
    font("legend.title",size = 14)
  return(p)
}

plot_basesub_percent_stack <- function(){
  p <- Vmtrna_Igtex_T %>% 
    left_join(Vmtrna_Gmt) %>% 
    select(ID,basesub_type,var_type) %>% 
    unique() %>% 
    # mutate(type = ifelse(basesub_type %in% c("T>C","C>T"),"SNV transition","SNV transversion")) %>% 
    # mutate(type = basesub_type) %>% 
    group_by(basesub_type,var_type) %>%
    count() %>%
    # ungroup() %>%
    group_by(basesub_type) %>% 
    mutate(percent = n / sum(n) * 100) %>%
    mutate(basesub_type = fct_relevel(basesub_type,rev)) %>% 
    ggbarplot("basesub_type","percent",fill = "var_type",color = "var_type",
              ylab = "Variants %",xlab = "Base Substitution",legend.title = "Variant type",
              legend = "right",
              palette = unique(Vmtrna_Gmt$var_color)[c(1,7,8,9,5,3,4)])+
    coord_flip()
  return(p)
}

plot_vartype_pie <- function(){
  p <- Vmtrna_Igtex_T %>% 
    left_join(Vmtrna_Gmt) %>% 
    group_by(var_type) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(percent = round(n / sum(n) * 100,1)) %>% 
    mutate(percent_label = str_c(percent,"%")) %>% 
    ggpie("percent",fill = "var_type",label = "percent_label",color = "white",
          legend = "right",lab.pos = "in",legend.title = "15 varaints per sample",
          lab.font = c(4,"plain","white"),palette = unique(Vmtrna_Gmt$var_color)[c(1,7,8,9,5,3,4)])+
    font("legend.title",size = 14)
  return(p)
}

plot_vartype_percent_stack <- function(){
  data_plot <- Vmtrna_Igtex_T %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.005,.01,.05,.1,.9,1),labels = c("<0.005","0.005~0.01","0.01~0.05","0.05~0.1","0.1~0.9","0.9~1"),ordered_result = T)) %>% 
    as_tibble() %>% 
    select(Indiv,tissue,ID,HF_interval,var_type) %>% 
    unique() %>% 
    group_by(HF_interval,var_type) %>% 
    count() %>% 
    # ungroup() %>% 
    group_by(var_type) %>% 
    mutate(percent = n / sum(n) * 100) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.005","0.005~0.01","0.01~0.05","0.05~0.1","0.1~0.9","0.9~1")))) %>% 
    mutate(var_type = fct_relevel(var_type,rev))
  p <- data_plot %>% 
    ggplot(aes(var_type,percent,fill = HF_interval))+
    geom_bar(stat = "identity")+
    ylab("Variants %")+
    xlab("Variant Type")+
    scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E","#F3722C")))+
    labs(fill = "Heteroplasmy")+
    theme_pubr()+
    theme(legend.position = "right")+
    coord_flip()
  return(p)
}

plot_clonal_expansion_example <- function(){
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
    ylab("Heteroplasmy")+
    theme_few()+
    theme(axis.text.y = element_text(color = var_color$var_color[c(4,6,7,1,2,5,3,8,9)]))+
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
  p <- p1 %>% insert_right(p2,width = .5)
  return(p)
}

plot_single_cell_example <- function(){
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
                      annotation_label = "Bulk",annotation_name_rot = 0,
                      var_type = anno_text(y_col$var_type,gp = gpar(col = y_col$var_color)))
  
  p <- data_cluster[y_sort$ID,x_sort$Indiv] %>%
    Heatmap(cluster_rows = F,cluster_columns = F,border = T,
            show_column_names = F,name = "Heteroplasmy",
            column_title = "Single cells",column_title_side = "bottom",
            show_heatmap_legend = F,
            col = colorRamp2(c(0,.2), c("white", "red")),
            right_annotation = ha,
            row_names_gp = gpar(col = y_col$var_color),row_names_side = "left")
  lgd <- Legend(col_fun = colorRamp2(c(0,.2), c("white", "red")),
                at = c(0,0.2),title = "Heteroplasmy",
                direction = "horizontal",
                legend_width = unit(2,"cm"),
                title_position = "topcenter")
  print(p)
  draw(lgd,x = unit(1, "npc") - unit(0.4, "npc"), y = unit(1, "npc") - unit(0.15, "npc"))
}

plot_3germ_layer_pie_prop_violin <- function(){
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
  
  # 饼图
  Vmtrna_Igtex_T %>% 
    as_tibble() %>% 
    inner_join(vars) %>% 
    filter(Indiv %in% indivs) %>% 
    anti_join(n_1tissue) %>% 
    left_join(gtex_color) %>%
    select(Indiv,ID,germ_layers) %>%
    unique() %>% 
    group_by(Indiv,ID) %>% 
    count() %>% 
    group_by(n) %>% 
    summarise(n_germ = n())
  
  data_plot <- tibble(name = c("Only 1 tissue","1 germ layer","2 germ layer","3 germ layer"),
                      n = c(10560,732,1549,782)) %>% 
    mutate(percent = round(n / sum(n) * 100,1)) %>% 
    # mutate(percent_label = str_c(n,"\n(",percent,"%)")) %>% 
    mutate(percent_label = str_c(percent,"%")) %>% 
    mutate(name = factor(name,levels = c("Only 1 tissue","1 germ layer","2 germ layer","3 germ layer")))
  p1 <- data_plot %>% 
    ggpie("percent",label = "percent_label",lab.pos = "in",lab.font = c(4,"plain","black"),
          fill = "name",palette = str_c("#",c("f4e9cd","f7e1d7","dedbd2","b0c4b1")),
          legend = "right",legend.title = "# Germ layers",color = "white")+
    font("legend.title",size = 14)
  
  # 比例堆积图
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
  
  ## 真实情况
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
           n_layer = factor(n_layer,levels = c(1,2,3),labels = c("1 germ layer","2 germ layer","3 germ layer"))) %>% 
    mutate(n_layer = fct_rev(n_layer)) %>%
    mutate(source = "GTEx")
  df2 <- readRDS("20240407_germ_layer_permute1000_results_merge_tissue_1to7and8+.rds") %>% 
    mutate(percent = percent * 100)
  df3 <- readRDS("20240415_germ_layer_permute1000_results_merge_tissue_1to7and8+_gt2.rds")
  
  p2 <- ggplot()+
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
  
  # 小提琴图
  data_plot <- Vmtrna_Igtex_T %>% 
    inner_join(Vmtrna_Igtex_N) %>% 
    group_by(Indiv,ID) %>% 
    summarise(gt_AF_diff = diff(range(gt_AF))) %>% 
    left_join(Vmtrna_Igtex_N) %>% 
    filter(n != 1)
  p3 <- data_plot %>% 
    ggviolin("quartile","gt_AF_diff",fill = "quartile",legend = "none",
             outlier.shape = 18,add = "boxplot",
             xlab = "# Tissues",ylab = "Heteroplasmy Range",
             palette = str_c("#",c("a9d6e5","89c2d9","61a5c2","468faf","2c7da0","2a6f97","014f86","01497c")))+
    scale_y_log10()+
    # stat_compare_means(comparisons = list(c("2","3"),c("3","4"),c("4","5"),c("5","6"),c("6","7"),c("7",">=8")),
    #                    label = "p.signif",hide.ns = T)
    stat_compare_means(comparisons = list(c("2","3"),c("3","4"),c("6","7"),c("7",">=8")),
                       label = "p.signif",hide.ns = T,step.increase = 0.03,tip.length = 0)
  
  # ROS与非ROS变异比例
  data_plot <- Vmtrna_Igtex_T %>% 
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
                              n == 2 ~ "2 germ layer",
                              n == 3 ~ "3 germ layer")) %>% 
    mutate(n_germ = factor(n_germ,levels = rev(c("Only 1 tissue","1 germ layer","2 germ layer","3 germ layer")))) %>% 
    left_join(Vmtrna_Igtex_T %>% 
                select(ID,basesub_type) %>% 
                unique()) %>% 
    # mutate(ROS = ifelse(basesub_type %in% c("C>A","C>G"),"Yes","No")) %>% 
    mutate(ROS = fct_relevel(basesub_type,"T:A>G:C","T:A>A:T",after = 0)) %>% 
    group_by(n_germ,ROS) %>%
    count() %>% 
    group_by(n_germ) %>%
    mutate(prop = n / sum(n) * 100) %>% 
    ungroup()
  
  # p4 <- ggbarplot(data_plot,"n_germ","prop",fill = "ROS",color = "ROS",palette = c("#5b8e7d","#f4a259"),
                  # xlab = "# Germ Layers",ylab = "Variants %",legend = "right",legend.title = "ROS variant")+
  p4 <- ggbarplot(data_plot,"n_germ","prop",fill = "ROS",color = "ROS",
                  palette = c("#014f86","#2a6f97","#2f8e44","#66bc4f","#f48c06","#ffba08"),
                  xlab = "# Germ Layers",ylab = "Variants %",legend = "right",legend.title = "")+
    guides(fill = guide_legend(reverse = T),
           color = guide_legend(reverse = T))+
    coord_flip()
  
  return(list(p1,p2,p3,p4))
}

# Main --------------------------------------------------------------------

plot_coverage_coordinate()
plot_all_sites_pie()
plot_basesub_pie()
plot_vartype_pie()
plot_basesub_percent_stack()
plot_vartype_percent_stack()
plot_clonal_expansion_example()
plot_single_cell_example()
plot_3germ_layer_pie_prop_violin()


# Functions ---------------------------------------------------------------

plot_var_density_circle <- function(){
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    group_by(POS) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(name = "chrM",
           start = POS,
           end = POS) %>% 
    select(name:end,n) %>% 
    mutate(n = log10(n))
  data_plot2 <- Vmtrna_Igtex_T_het_var %>% 
    filter(gt_AF > .05) %>% 
    group_by(POS) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(name = "chrM",
           start = POS,
           end = POS) %>% 
    select(name:end,n) %>% 
    mutate(n = log10(n))
  
  chrM_coord <- data.frame(name = "chrM",start = 1,end = 16569)
  
  circos.par(start.degree = 90,clock.wise = T,gap.degree = 0,track.margin = c(.005,.005))
  circos.genomicInitialize(chrM_coord,plotType = "labels",sector.names = " ")
  circos.genomicDensity(data_plot2,bg.border = NA,col = "#6c757d",track.height = .3)
  circos.genomicAxis(h = "top",major.by = 2000,col = "gray",labels.cex = .8,minor.ticks = 0,labels.col = "gray")
  circos.genomicDensity(data_plot,bg.border = NA,col = "#adb5bd",track.height = .3,area = T)
  circos.axis(h = "bottom",col = "gray",labels = F)
  circos.genomicTrack(chrM_gtf,panel.fun = function(region,value,...){
    circos.genomicRect(region,value,col = chrM_gtf$color,border = NA,...)
  },bg.border = NA,track.height = .1)
  circos.clear()
  
  lgd_rec <- Legend(at = c("D-loop","Protein coding","rRNA","tRNA"),
                    legend_gp = gpar(fill = c("#ffdc57","#9dbfdc","#c9a2cf","#fcac55")),
                    labels_gp = gpar(fontsize = 12),nrow = 1)
  draw(lgd_rec,x = unit(.5,"npc"),y = unit(0,"npc"))
}

plot_no.het_age_vartype <- function(){
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    left_join(Vmtrna_Gmt) %>% 
    group_by(Indiv,tissue,var_type) %>% 
    summarise(HFN = n()) %>% 
    ungroup() %>% 
    left_join(Igtex) %>% 
    select(Indiv,tissue,var_type,HFN,AGE)
  
  var_types <- c("D-loop","LoF","Missense\nconserved","Missense\nnon-conserved","rRNA\nconserved","rRNA\nnon-conserved","Synonymous")
  res_lm <- lapply(var_types,function(x){
    data_sub <- data_plot %>% 
      filter(var_type == x)
    res <- lm(HFN ~ AGE,data_sub) %>% 
      summary()
    return(tibble(var_type = x,
                  intercept = res$coefficients[1,1],
                  beta = res$coefficients[2,1],
                  p = res$coefficients[2,4]))
  })
  res_lm <- do.call(rbind,res_lm) %>% 
    arrange(var_type)
  print(res_lm)
  
  formula <- res_lm %>% 
    mutate(across(where(is.numeric),~ round(.,2))) %>% 
    mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
    pull(formula)
  label <- str_c(var_types,formula)[c(1,7,4,3,2,6,5)]
  color <- c("#ADB5BD","#E5383B","#FF7D00","#EDDEA4","#5A189A","#E4C1F9","#2D898B")[c(1,7,4,3,2,6,5)]
  
  p <- data_plot %>%
    ggplot(aes(AGE,HFN,color = var_type,fill = var_type))+
    geom_smooth(method = "lm")+
    scale_fill_manual(values = color,labels = label)+
    scale_color_manual(values = color,labels = label)+
    xlab("Age")+
    ylab("# Heteroplasmies")+
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL))+
    theme_pubr()+
    theme(legend.position =c(.3,.835),)
  return(p)
}

plot_population_individual_frequency <- function(){
  data_plot <- Vmtrna_Igtex_T %>% 
    group_by(ID) %>% 
    summarise(max_hl = max(gt_AF)) %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(var_type = str_replace(var_type,"\n"," ")) %>% 
    mutate(people_freq = cut(gnomAD_AF_hom,c(-1,1e-10,1e-3,1e-2,1e-1,1),
                             labels = c("Not reported","<0.1%","0.1~1%","1~10%",">10%")),
           heteroplasmy = cut(max_hl,c(0,.005,.01,.05,.1,.9,1),
                              labels = c("<0.005","0.005~0.01","0.01~0.05","0.05~0.1","0.1~0.9",">0.9"),ordered_result = T)) %>% 
    mutate(heteroplasmy = fct_relevel(heteroplasmy,rev)) %>% 
    group_by(people_freq,heteroplasmy,var_type) %>% 
    count() %>% 
    group_by(var_type) %>% 
    mutate(percent = n / sum(n) * 100) %>% 
    mutate(var_type = factor(var_type,levels = c("D-loop","rRNA conserved","rRNA non-conserved",
                                                 "Synonymous","Missense conserved","Missense non-conserved",
                                                 "LoF")))
  
  p <- ggbarplot(data_plot,"people_freq","percent",fill = "heteroplasmy",color = "heteroplasmy",
                 facet.by = "var_type",legend = "none",
                 palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E","#F3722C")))+
    coord_flip()+
    xlab("gnomAD Allele Frequency")+
    ylab("Proportion of Unique Variants")+
    labs(fill = "Max heteroplasmy")+
    guides(fill = guide_legend(reverse = F),
           color = guide_legend(reverse = F))+
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(),
          axis.text = element_text(size = 10))+
    guides(color = "none")
  return(p)
}

plot_population_individual_frequency_zoom <- function(){
  data_plot <- Vmtrna_Igtex_T %>% 
    group_by(ID) %>% 
    summarise(max_hl = max(gt_AF)) %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(var_type = str_replace(var_type,"\n"," ")) %>% 
    filter(str_detect(SYMBOL,"^MT-[^T|R]"),
           var_type == "Missense non-conserved") %>% 
    mutate(SYMBOL = ifelse(!SYMBOL %in% c("MT-ATP6","MT-ATP8"),"Other 11 genes",SYMBOL)) %>% 
    mutate(people_freq = cut(gnomAD_AF_hom,c(-1,1e-10,1e-3,1e-2,1e-1,1),
                             labels = c("Not reported","<0.1%","0.1~1%","1~10%",">10%")),
           heteroplasmy = cut(max_hl,c(0,.005,.01,.05,.1,.9,1),
                              labels = c("<0.005","0.005~0.01","0.01~0.05","0.05~0.1","0.1~0.9",">0.9"),ordered_result = T)) %>% 
    mutate(heteroplasmy = fct_relevel(heteroplasmy,rev)) %>% 
    group_by(SYMBOL,people_freq,heteroplasmy,var_type) %>% 
    count() %>% 
    group_by(SYMBOL,var_type) %>% 
    mutate(percent = n / sum(n) * 100)
  
  p <- ggbarplot(data_plot,"people_freq","percent",fill = "heteroplasmy",color = "heteroplasmy",
                 facet.by = c("SYMBOL","var_type"),legend = "right",
                 palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E","#F3722C")))+
    coord_flip()+
    # xlab("gnomAD Allele Frequency")+
    xlab("")+
    ylab("Proportion of Unique Variants")+
    labs(fill = "Max heteroplasmy")+
    guides(fill = guide_legend(reverse = F),
           color = "none")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(),
          axis.text = element_text(size = 10))
  return(p)
}

plot_no.het_tissue_cons_basesubtype <- function(){
  p <- Vmtrna_Igtex_T_het_var %>% 
    left_join(Igtex %>% 
                select(Indiv,AGE)) %>% 
    filter(between(AGE,50,70)) %>% 
    left_join(Vmtrna_Gmt) %>% 
    mutate(HFN = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),T,F)) %>%
    group_by(tissue,Indiv,basesub_type) %>%
    summarise(HFN = sum(HFN)) %>% 
    pivot_wider(id_cols = c("tissue","Indiv"),names_from = "basesub_type",values_from = "HFN") %>% 
    mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
    pivot_longer(names_to = "basesub_type",values_to = "HFN",cols = where(is.numeric)) %>% 
    group_by(tissue,basesub_type) %>% 
    summarise(HFN = mean(HFN)) %>% 
    left_join(gtex_color) %>% 
    mutate(basesub_type = factor(basesub_type,levels = c("T>C","C>T","C>A","C>G","T>G","T>A")[c(1,2,5,6,3,4)])) %>% 
    ggplot(aes(tissue_site_detail,HFN,fill = basesub_type))+
    geom_bar(stat = "identity")+
    theme_pubr()+
    scale_fill_manual(values = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97")[c(1,2,5,6,3,4)])+
    ylab("Number of mutations")+
    labs(fill = "")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          legend.position = "right")
}

plot_HF_tissue_cons <- function(){
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
    ylab("Total clonal abundance")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(fill = "none")
}

plot_no.het_tissue_cons <- function(){
  data_plot <- Vmtrna_Igtex_T_het_var %>%
    left_join(Igtex %>%
                select(Indiv,AGE)) %>%
    filter(between(AGE,50,70)) %>%
    left_join(Vmtrna_Gmt) %>%
    mutate(conserve = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),T,F)) %>%
    group_by(tissue,Indiv) %>%
    summarise(HFN = sum(conserve))

  data_simple <- bootstrap_multi_tissue(data_plot,gtex_49_tissues,"HFN","mean")
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

  p <- ggplot(data_simple,aes(tissue_site_detail,V1))+
    geom_bar(aes(fill = tissue_site_detail),stat = "identity")+
    geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
    theme_pubr()+
    scale_fill_manual(values = xlab_color)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("Number of mutations\n(Conserved)")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text())+
    guides(fill = "none")
  
  # data <- Vmtrna_Igtex_T_het_var %>%
  #   left_join(Igtex %>% 
  #               select(Indiv,AGE)) %>%
  #   filter(between(AGE,50,70)) %>% 
  #   left_join(Vmtrna_Gmt) %>% 
  #   mutate(gt_AF = ifelse(var_type %in% c("Missense\nconserved","rRNA\nconserved","LoF"),gt_AF,0)) %>% 
  #   mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),
  #                            ordered_result = T,include.lowest = T)) %>%
  #   mutate(HFN = ifelse(gt_AF != 0,1,0)) %>% 
  #   group_by(Indiv,tissue,HF_interval) %>%
  #   summarise(n = sum(HFN)) %>% 
  #   pivot_wider(id_cols = c(Indiv,tissue),names_from = HF_interval,values_from = n) %>% 
  #   ungroup() %>% 
  #   mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
  #   pivot_longer(cols = where(is.numeric),names_to = "HF_interval",values_to = "n") %>% 
  #   group_by(tissue,HF_interval) %>% 
  #   summarise(n = mean(n)) %>% 
  #   mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1")))) %>% 
  #   left_join(gtex_color)
  # 
  # p <- ggplot(data,aes(tissue_site_detail,n,fill = HF_interval))+
  #   geom_bar(stat = "identity")+
  #   ylab("Number of mutations")+
  #   scale_x_discrete(limits = xlab_sort)+
  #   scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
  #   labs(fill = "Clonal abundance")+
  #   theme_pubr()+
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         legend.position = "right")
  return(p)
}

plot_slope_tissue_cons <- function(){
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

  p <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    geom_hline(yintercept = 0,linetype = 2)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("# mutations per decade")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(color = "none")
  return(p)
}

plot_no.het_tissue_cons_prop <- function(){
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
  p <- ggplot(data_plot,aes(tissue_site_detail,n))+
    geom_bar(aes(fill = var_type),stat = "identity")+
    theme_pubr()+
    # scale_fill_manual(values = fill_color$var_color[c(1:3,6,4,5,7)])+
    scale_fill_manual(values = fill_color$var_color)+
    ylab("Number of mutations")+
    labs(fill = "Variant type")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          legend.position = "right")
  return(p)
}

plot_tissue_CPN <- function(){
  p <- Igtex_T_CPN %>% 
    inner_join(gtex_color) %>% 
    ggplot(aes(tissue_site_detail,CPN,color = tissue_site_detail))+
    geom_boxplot()+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    ylab("Relative Copy Number")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 11))+
    guides(color = "none")
  return(p)
}

plot_tissue_xCell_score <- function(){
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
  
  p <- ggplot(data,aes(tissue_site_detail,Celltype,fill = score))+
    geom_tile()+
    scale_fill_gradient(low = "white",high = "red")+
    ylab("Cell type")+
    labs(fill = "Enrichment\nscore")+
    theme_pubr()+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          legend.position = "right")
  return(p)
}

# Main --------------------------------------------------------------------

plot_var_density_circle()
plot_no.het_age_vartype()
plot_population_individual_frequency()
plot_population_individual_frequency_zoom()

p1 <- plot_no.het_tissue_cons()
p5 <- plot_slope_tissue_cons()
p4 <- plot_no.het_tissue_cons_prop()
p2 <- plot_tissue_CPN()
p3 <- plot_tissue_xCell_score()
p6 <- plot_no.het_tissue_cons_basesubtype()
p7 <- plot_HF_tissue_cons()
p7 %>% insert_bottom(p1) %>% insert_bottom(p5) %>% insert_bottom(p6) %>% insert_bottom(p4) %>% insert_bottom(p2) %>% insert_bottom(p3) #14*16

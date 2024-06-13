
# Functions ---------------------------------------------------------------

plot_slope_tissue <- function(Vmtrna_Igtex_T_HFN){
  # 克隆扩张与年龄斜率森林图
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

plot_slope_tissue_hotspot <- function(){
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
  
  p <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    geom_hline(yintercept = 0,linetype = 2)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("# mutations per decade\n(hotspot)")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(color = "none")
  return(p)
}

plot_slope_tissue_non_hotspot <- function(){
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
  
  p <- ggplot(data_plot,aes(tissue_site_detail,slope,color = tissue))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    geom_hline(yintercept = 0,linetype = 2)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("# mutations per decade\n(non-hotspot)")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(color = "none")
  return(p)
}

plot_tissue_difference_HF <- function(Vmtrna_Igtex_T_HF){
  data_simple <- bootstrap_multi_tissue(Vmtrna_Igtex_T_HF,gtex_49_tissues,"HF","median")
  # data_simple <- bootstrap_multi_tissue(Vmtrna_Igtex_T_HF,gtex_49_tissues,"HF","mean")
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
    ylab("Total clonal abundance")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 11))+
    guides(fill = "none")
  return(list(p = p,
              xlab_sort = xlab_sort))
}

plot_tissue_difference_HFN_stack <- function(Vmtrna_Igtex_T_het_var,xlab_sort){
  data <- Vmtrna_Igtex_T_het_var %>%
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

  p <- ggplot(data,aes(tissue_site_detail,n,fill = HF_interval))+
    geom_bar(stat = "identity")+
    ylab("Number of mutations")+
    scale_x_discrete(limits = xlab_sort)+
    scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
    labs(fill = "Clonal abundance")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          legend.position = "right")
  return(p)
}

plot_tissue_difference_hotspot <- function(Vmtrna_Igtex_T_het_var){
  data <- Vmtrna_Igtex_T_het_var %>%
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),
                             ordered_result = T,include.lowest = T)) %>%
    mutate(HFN = ifelse(gt_AF != 0,1,0))
  ids <- data %>% 
    group_by(tissue,ID) %>% 
    summarise(n = sum(HFN)) %>% 
    group_by(tissue) %>% 
    arrange(desc(n)) %>% 
    slice_head(n = 10) %>% 
    mutate(hotspot = "Yes") %>% 
    select(-n)
  
  data_plot <- data %>% 
    left_join(ids) %>% 
    mutate(hotspot = ifelse(is.na(hotspot),"No",hotspot)) %>% 
    group_by(Indiv,tissue,hotspot) %>%
    summarise(n = sum(HFN)) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = c(Indiv,tissue),names_from = hotspot,values_from = n) %>% 
    mutate(across(where(is.numeric),~replace_na(.,0))) %>% 
    pivot_longer(cols = where(is.numeric),names_to = "hotspot",values_to = "n") %>% 
    group_by(tissue,hotspot) %>% 
    summarise(n = mean(n)) %>% 
    group_by(tissue) %>% 
    mutate(prop = n / sum(n) * 100) %>% 
    left_join(gtex_color) %>% 
    mutate(hotspot = factor(hotspot,levels = c("No","Yes")))
  
  p <- ggplot(data_plot,aes(tissue_site_detail,prop,fill = hotspot))+
    geom_bar(stat = "identity")+
    ylab("Proportion of mutations")+
    scale_fill_manual(values = c("#cae9ff","#5fa8d3"))+
    labs(fill = "Hotspot")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          legend.position = "right")
  return(p)
}

plot_tissue_xCell_score <- function(Igtex_T_xCell_score,xlab_sort){
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
  return(p)
}

plot_brain_cibersort_score <- function(){
  load("estgtex.rda")
  p <- est.gtex$CIB$MB %>% 
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
    scale_x_discrete(limits = res1$xlab_sort[str_starts(res1$xlab_sort,"Brain")])+
    ylab("")+
    labs(fill = "Fraction")+
    theme_pubr()+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60,size = 11),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_blank(),
          legend.position = "right")
  return(p)
}

plot_tissue_CPN <- function(Igtex_T_CPN){
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

plot_no.het_age_basesubtype <- function(){
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    group_by(Indiv,tissue,basesub_type) %>% 
    summarise(HFN = n()) %>% 
    ungroup() %>% 
    left_join(Igtex) %>% 
    select(Indiv,tissue,basesub_type,HFN,AGE) %>% 
    mutate(basesub_type = factor(basesub_type,levels = c("T>C","C>T","C>A","C>G","T>G","T>A")))
  
  basesub_types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  res_lm <- lapply(basesub_types,function(x){
    data_sub <- data_plot %>% 
      filter(basesub_type == x)
    res <- lm(HFN ~ AGE,data_sub) %>% 
      summary()
    return(tibble(basesub_type = x,
                  intercept = res$coefficients[1,1],
                  beta = res$coefficients[2,1],
                  p = res$coefficients[2,4]))
  })
  res_lm <- do.call(rbind,res_lm) %>% 
    arrange(basesub_type)
  print(res_lm)
  
  formula <- res_lm %>% 
    mutate(across(where(is.numeric),~ round(.,2))) %>% 
    mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
    pull(formula)
  label <- str_c(basesub_types,formula)
  
  p <- data_plot %>%
    ggplot(aes(AGE,HFN,color = basesub_type,fill = basesub_type))+
    geom_smooth(method = "lm")+
    scale_fill_manual(values = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"),labels = label[c(5,3,1,2,6,4)])+
    scale_color_manual(values = c("#2f8e44","#66bc4f","#f48c06","#ffba08","#014f86","#2a6f97"),labels = label[c(5,3,1,2,6,4)])+
    xlab("Age")+
    ylab("# Heteroplasmies")+
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL))+
    theme_pubr()+
    theme(legend.position = c(.3,.84),)
  return(p)
}

plot_no.het_tissue_ROS <- function(){
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    mutate(ROS = ifelse(basesub_type %in% c("C>A","C>G"),T,F)) %>% 
    select(Indiv,tissue,ID,ROS) %>% 
    group_by(Indiv,tissue) %>% 
    summarise(HFN = sum(ROS))
  
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
  
  p1 <- ggplot(data_simple,aes(tissue_site_detail,V1))+
    geom_bar(aes(fill = tissue_site_detail),stat = "identity")+
    geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
    theme_pubr()+
    scale_fill_manual(values = xlab_color)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("# Heteroplasmies (ROS)")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text())+
    guides(fill = "none")
  
  p2 <- Igtex_T_CPN %>% 
    inner_join(gtex_color) %>% 
    ggplot(aes(tissue_site_detail,CPN,color = tissue_site_detail))+
    geom_boxplot()+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    ylab("Relative Copy Number")+
    theme_pubr()+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11))+
    guides(color = "none")
  
  p <- p1 %>% insert_bottom(p2)
  return(p)
}

plot_slope_no.het <- function(){
  hfn_at_60 <- Vmtrna_Igtex_T_HFN %>% 
    left_join(Igtex) %>% 
    filter(between(AGE,50,70)) %>% 
    group_by(tissue) %>% 
    summarise(HFN = mean(HFN),n_tissue = n())
  
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
    left_join(res_lm) %>% 
    left_join(gtex_color)
  p <- ggplot(data_plot,aes(beta,HFN,))+
    geom_point(aes(size = n_tissue,color = tissue))+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    scale_size_continuous(breaks = c(100,200,300))+
    ylim(0,11.5)+
    guides(color = "none")+
    ylab("Mean # Heteroplasmies")+
    xlab("# Heteroplasmies per 10 Years")+
    labs(size = "Sample Size (n)")+
    theme_pubr()+
    stat_cor()
  return(p)
}

plot_het_age_heart4 <- function(x){
  ## 变异数量图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == "Heart-LeftVentricle") %>%
    # filter(tissue == "Esophagus-Mucosa") %>%
    # filter(tissue == "Brain-Cortex") %>%
    left_join(Vmtrna_Gmt) %>% 
    mutate(HFN = ifelse(basesub_type %in% x,T,F)) %>%
    left_join(Igtex) %>% 
    mutate(HF_interval = cut(gt_AF,c(0,.01,.05,.1,1),labels = c("<0.01","0.01~0.05","0.05~0.1",">0.1"),ordered_result = T)) %>%
    mutate(HF_interval = fct_relevel(HF_interval,rev),
           AGE = cut(AGE,breaks = seq(20,70,5),
                     # labels = c("20~25","25~30","30~35","35~40","40~45","45~50","50~55","55~60","60~65","65~70"),include.lowest = T)) %>% 
                     labels = c("20","25","30","35","40","45","50","55","60","65"),include.lowest = T)) %>% 
    group_by(Indiv,AGE,HF_interval) %>%
    summarise(n = sum(HFN)) %>% 
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  if(all(x == c("T>C","C>T"))){
    p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                    palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                    # xlab = "Age",ylab = "# Heteroplasmies",legend = "none")+
                    xlab = "",ylab = "",legend = "none")+
      ylim(0,12)+
      theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  }else{
    p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                    palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                    # xlab = "Age",ylab = "# Heteroplasmies",legend = "none")+
                    xlab = "",ylab = "",legend = "none")+
      ylim(0,3)+
      theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  }
  
  ## 总克隆丰度图
  data_plot <- Vmtrna_Igtex_T_het_var %>% 
    filter(tissue == "Heart-LeftVentricle") %>%
    # filter(tissue == "Esophagus-Mucosa") %>%
    # filter(tissue == "Brain-Cortex") %>%
    left_join(Vmtrna_Gmt) %>% 
    mutate(gt_AF = ifelse(basesub_type %in% x,gt_AF,0)) %>%
    left_join(Igtex) %>% 
    group_by(Indiv,AGE) %>%
    summarise(HF = sum(gt_AF)) %>% 
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
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      ylim(0,.18)+
      # ylim(0,.15)+
      # ylim(0,.05)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())+
      guides(fill = "none")
  }else{
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      # ylim(0,.18)+
      # ylim(0,.15)+
      ylim(0,.03)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())+
      guides(fill = "none")
  }
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

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
    group_by(AGE,HF_interval) %>% 
    summarise(n = mean(n)) %>% 
    mutate(HF_interval = factor(HF_interval,levels = rev(c("<0.01","0.01~0.05","0.05~0.1",">0.1"))))
  if(all(x == c("T>C","C>T"))){
    p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                    palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                    xlab = "Age",ylab = "# Heteroplasmies",legend = "none")+
      # xlab = "",ylab = "",legend = "none")+
      ylim(0,12)+
      theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  }else{
    p1 <- ggbarplot(data_plot,"AGE","n",fill = "HF_interval",color = "HF_interval",
                    palette = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")),
                    xlab = "Age",ylab = "# Heteroplasmies",legend = "none")+
      # xlab = "",ylab = "",legend = "none")+
      ylim(0,3)+
      theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
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
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      ylim(0,.18)+
      # ylim(0,.15)+
      # ylim(0,.05)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())+
      guides(fill = "none")
  }else{
    p2 <- ggplot(data_simple,aes(AGE,V1))+
      geom_bar(aes(fill = AGE),stat = "identity")+
      geom_errorbar(aes(ymin = V2,ymax = V3),width = .5)+
      theme_pubr()+
      scale_fill_manual(values = str_c("#",c("99e2b4","88d4ab","78c6a3","67b99a","56ab91",
                                             "469d89","358f80","248277","14746f","036666")))+
      xlab("Age")+
      ylab("Sum of Heteroplasmy")+
      # ylim(0,.18)+
      # ylim(0,.15)+
      ylim(0,.03)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())+
      guides(fill = "none")
  }
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

# Main --------------------------------------------------------------------

# 同年龄段个体C>T/T>C变异
I_T_het_var <- Vmtrna_Igtex_T_het_var %>% 
  as_tibble() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70)) %>% 
  select(Indiv,tissue) %>% 
  unique()
Vmtrna_Igtex_T_het_var_nonROS_OldAge <- Vmtrna_Igtex_T_het_var %>% 
  filter(basesub_type %in% c("T>C","C>T")) %>% 
  right_join(I_T_het_var) %>% 
  mutate(gt_AF = ifelse(is.na(gt_AF),0,gt_AF))
Vmtrna_Igtex_T_HF_nonROS_OldAge <- Vmtrna_Igtex_T_het_var %>% 
  mutate(nonROS = ifelse(basesub_type %in% c("T>C","C>T"),gt_AF,0)) %>% 
  group_by(tissue,Indiv) %>% 
  summarise(HF = sum(nonROS)) %>% 
  ungroup() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70))
# Vmtrna_Igtex_T_HFN_nonROS_OldAge <- Vmtrna_Igtex_T_het_var %>%
Vmtrna_Igtex_T_HFN_nonROS <- Vmtrna_Igtex_T_het_var %>%
  mutate(nonROS = ifelse(basesub_type %in% c("T>C","C>T"),T,F)) %>%
  group_by(tissue,Indiv) %>%
  summarise(HFN = sum(nonROS)) %>%
  ungroup()
  # left_join(Igtex %>%
  #             select(Indiv,AGE)) %>%
  # filter(between(AGE,50,70))
Igtex_T_CPN_OldAge <- Igtex_T_CPN %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70))

p_res <- readRDS("20240507_p_res_linear.rds")
sig_ids <- p_res %>% 
  mutate(q = p.adjust(p,method = "fdr")) %>% 
  filter(q < .05) %>% 
  select(tissue,ID,beta)
# ids <- sig_ids[c(4:5,7:20,22:30,32:33,35,41:42,44:52),] %>% 
ids <- sig_ids[c(4:5,7:20,22:30,32:33,35,41:42,44:52),] %>% 
  mutate(hotspot = T) %>% 
  select(ID,hotspot) %>% 
  unique() %>% 
  filter(ID != "chrM_64_C_A")

p6 <- plot_slope_tissue(Vmtrna_Igtex_T_HFN_nonROS)
p8 <- plot_slope_tissue_hotspot()
p9 <- plot_slope_tissue_non_hotspot()
res1 <- plot_tissue_difference_HF(Vmtrna_Igtex_T_HF_nonROS_OldAge)
p2 <- plot_tissue_difference_HFN_stack(Vmtrna_Igtex_T_het_var_nonROS_OldAge,xlab_sort = res1$xlab_sort)
# p7 <- plot_tissue_difference_hotspot(Vmtrna_Igtex_T_het_var_nonROS_OldAge)
# p3 <- plot_tissue_CPN(Igtex_T_CPN_OldAge)
p4 <- plot_tissue_xCell_score(Igtex_T_xCell_score,xlab_sort = res1$xlab_sort)
p5 <- plot_brain_cibersort_score()
# p6 %>% insert_bottom(res1$p) %>% insert_bottom(p2) %>% insert_bottom(p7) %>% insert_bottom(p4,height = .8) %>% insert_bottom(p5,height = .25)
p6 %>% insert_bottom(p8) %>% insert_bottom(p9) %>% insert_bottom(res1$p) %>% insert_bottom(p2) %>% insert_bottom(p4,height = .8) %>% insert_bottom(p5,height = .25) #13*14

plot_no.het_tissue_ROS()
plot_no.het_age_basesubtype()
plot_slope_no.het()

basesub_types <- list(c("T>C","C>T"),c("C>A","C>G"))
ps <- lapply(basesub_types,plot_het_age_heart4)
aplot::plot_list(gglist = ps,nrow = 1,widths = c(1,1))
ps <- lapply(basesub_types,plot_het_age_brain)
aplot::plot_list(gglist = ps,nrow = 1,widths = c(1,1))

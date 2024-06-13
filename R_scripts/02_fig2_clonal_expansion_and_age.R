
# Functions ---------------------------------------------------------------

plot_no.het_age_liver <- function(){
  # 肝脏20~50岁分10个bin
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
                  xlab = "Age",ylab = "# Heteroplasmies",legend = "right")+
    labs(fill = "Heteroplasmy",color = "Heteroplasmy")+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),)
  
  # 肝脏增加总克隆丰度图
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
    ylab("Sum of Heteroplasmy")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(fill = "none")
  
  p <- p2 %>% insert_bottom(p1)
  return(p)
}

plot_no.het_age_tissue_difference <- function(){
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
  
  tissues <- c("ESPMCS","ESPMSL","HRTLV","LIVER","SKINNS","SKINS")
  data_plot <- data_plot %>% 
    filter(tissue_site_detail_abbr %in% tissues) %>%
    mutate(tissue_site_detail_abbr = factor(tissue_site_detail_abbr,levels = c("LIVER","HRTLV","ESPMCS","SKINS","SKINNS","ESPMSL")))
  
  formula <- res_lm %>% 
    left_join(gtex_color) %>% 
    filter(tissue_site_detail_abbr %in% tissues) %>% 
    arrange(tissue_site_detail_abbr) %>% 
    group_by(tissue) %>% 
    mutate(across(where(is.numeric),~ round(.,2))) %>% 
    mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
    pull(formula)
  label <- str_c(tissues,formula)
  
  p <- ggplot(data_plot,aes(AGE,HFN,color = tissue_site_detail_abbr,fill = tissue_site_detail_abbr))+
    geom_smooth(method = "lm")+
    scale_fill_manual(values = unique(data_plot$tissue_color_hex)[c(4,3,1,6,5,2)],labels = label[c(4,3,1,6,5,2)])+
    scale_color_manual(values = unique(data_plot$tissue_color_hex)[c(4,3,1,6,5,2)],labels = label[c(4,3,1,6,5,2)])+
    xlab("Age")+
    ylab("# Heteroplasmies")+
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL))+
    theme_pubr()+
    # theme(legend.position = "right")
    theme(legend.position =c(.3,.85),)
  return(p)
}

plot_correlation_tissue <- function(){
  # 克隆扩张与年龄相关性森林图
  data_corr <- Vmtrna_Igtex_T_HFN %>% 
    left_join(Igtex) %>% 
    select(Indiv,tissue,AGE,HFN)
  
  res_corr <- lapply(gtex_49_tissues,function(x){
    res_corr <- data_corr %>% 
      filter(tissue == x) %>% 
      select(AGE,HFN) %>% 
      bootstrap_corr()
    return(tibble(tissue = x,
                  corr = res_corr[1],
                  lower = res_corr[2],
                  upper = res_corr[3]))
  })
  res_corr <- do.call(rbind,res_corr)
  data_plot <- res_corr %>% 
    left_join(gtex_color)
  xlab_sort <- data_plot %>% 
    group_by(tissue_site_detail) %>% 
    summarise(corr = median(corr)) %>% 
    arrange(desc(corr)) %>% 
    pull(tissue_site_detail)
  
  p <- ggplot(data_plot,aes(tissue_site_detail,corr,color = tissue))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lower,ymax = upper),width = .5,size = 1)+
    scale_color_manual(values = gtex_color$tissue_color_hex)+
    geom_hline(yintercept = 0,linetype = 2)+
    scale_x_discrete(limits = xlab_sort)+
    ylab("Correlation with Age")+
    theme_pubr()+
    # theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
    #       axis.title.x = element_blank(),)+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(color = "none")
  return(p)
}

plot_slope_tissue <- function(){
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
    ylab("# Heteroplasmies\nper 10 Years")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(color = "none")
  return(p)
}

plot_HF_tissue <- function(Vmtrna_Igtex_T_HF){
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
    ylab("Sum of Heteroplasmy")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 11))+
    guides(fill = "none")
  return(p)
}

plot_HFN_tissue <- function(Vmtrna_Igtex_T_het_var){
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
    ylab("# Heteroplasmies")+
    # scale_x_discrete(limits = xlab_sort)+
    scale_fill_manual(values = rev(c("#577590","#43AA8B","#90BE6D","#F9C74F","#F8961E")))+
    labs(fill = "Heteroplasmy")+
    ylim(0,11.5)+
    theme_pubr()+
    theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60),
          axis.title.x = element_blank(),
          legend.position = "right")
  return(p)
}

# Main --------------------------------------------------------------------

Vmtrna_Igtex_T_het_var_OldAge <- Vmtrna_Igtex_T_het_var %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70))
Vmtrna_Igtex_T_HF_OldAge <- Vmtrna_Igtex_T_het_var %>% 
  group_by(tissue,Indiv) %>% 
  summarise(HF = sum(gt_AF)) %>% 
  ungroup() %>% 
  left_join(Igtex %>% 
              select(Indiv,AGE)) %>% 
  filter(between(AGE,50,70))

plot_no.het_age_liver()
plot_no.het_age_tissue_difference()

p1 <- plot_slope_tissue()
p2 <- plot_correlation_tissue()
p3 <- plot_HF_tissue(Vmtrna_Igtex_T_HF_OldAge)
p4 <- plot_HFN_tissue(Vmtrna_Igtex_T_het_var_OldAge)
p1 %>% insert_bottom(p2) %>% insert_bottom(p3) %>% insert_bottom(p4)

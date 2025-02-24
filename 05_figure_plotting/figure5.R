
# Figure 5A ---------------------------------------------------------------

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
formula <- c(": + 0.02 /y",": + 0.09 /y")
label <- str_c(c("L","Non L"),formula)
ggplot(data,aes(AGE,HFN,color = Macrogroup,fill = Macrogroup))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = Macrogroup,fill = Macrogroup,shape = Macrogroup),
             size = 1.5, alpha = 0.7)+
  scale_fill_manual(values = c("#1f77b4","#ff7f0e")[c(2,1)],labels = label)+
  scale_color_manual(values = c("#1f77b4","#ff7f0e")[c(2,1)],labels = label)+
  xlab("Age")+
  ylab("Number of mutations")+
  ggtitle("Skin - Sun Exposed")+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_none())+
  theme_pubr()+
  theme(legend.position =c(.25,.9),
        plot.title = element_text(hjust = .5)) # 4*5

# Figure 5B ---------------------------------------------------------------

Igtex_macrogrp <- Igtex %>% 
  filter(!is.na(Haplogroup)) %>% 
  # filter(between(AGE,50,70)) %>% 
  mutate(Macrogroup = ifelse(Macrogroup %in% c("B","R","J","T"),"R*",
                             ifelse(Macrogroup %in% c("X","A","N","W","I"),"N*",
                                    ifelse(Macrogroup %in% c("D","M","C","Z","L3"),"L3*",Macrogroup)))) %>%
  mutate(Macrogroup = ifelse(Macrogroup == "V","H",
                             ifelse(Macrogroup == "K","U",
                                    ifelse(Macrogroup %in% c("L1","L2"),"L0/1/2",Macrogroup))))

data_slope <- Vmtrna_Igtex_T_HFN %>% 
  filter(str_detect(tissue,"Skin")) %>%
  inner_join(Igtex_macrogrp)

conditions <- crossing(Macrogroup = unique(Igtex_macrogrp$Macrogroup),
                       tissue = c("Skin-NotSunExposed_Suprapubic","Skin-SunExposed_Lowerleg"))
res_slope <- apply(conditions,1,function(x){
  res_slope <- data_slope %>% 
    filter(Macrogroup == x[1],
           tissue == x[2]) %>% 
    select(AGE,HFN) %>% 
    bootstrap_slope()
  return(tibble(Macrogroup = x[1],
                tissue = x[2],
                slope = res_slope[1] * 10,
                lower = res_slope[2] * 10,
                upper = res_slope[3] * 10))
})
res_slope <- do.call(rbind,res_slope)
data_plot <- res_slope %>% 
  left_join(gtex_color)
# mutate(Macrogroup = factor(Macrogroup,levels = c("L0")))

p1 <- ggplot(data_plot,aes(Macrogroup,slope,color = tissue_site_detail))+
  geom_point(size = 3,position = position_dodge(width = .5))+
  geom_errorbar(aes(ymin = lower,ymax = upper),width = .3,size = 1,position = position_dodge(width = .5))+
  scale_color_manual(values = c("#0000FF","#7777FF"))+
  geom_hline(yintercept = 0,linetype = 2)+
  ylab("Mutations gained\nper 10y")+
  ylim(-1,2)+
  theme_pubr()+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

nwk <- "(L1:4,(L3:3,(N:2,(R:1,U:1,H:1):1[&&NHX:H=R]):1[&&NHX:H=N]):1[&&NHX:H=L3]);"
# tree <- read.tree(text = nwk)
tree <- treeio::read.nhx(textConnection(nwk))
tree@phylo$tip.label <- c("L0/1/2","L3*","N*","R*","U","H")
p3 <- ggtree(tree)+
  geom_label(aes(label=H),label.padding = unit(0.15,"lines"))+
  # geom_tiplab()+
  coord_flip()

## 每个类别里面RACE的比例用stacked bar展示
p2 <- Igtex_macrogrp %>% 
  mutate(RACE = case_when(RACE == 1 ~ "Asian",
                          RACE == 2 ~ "African American",
                          RACE == 3 ~ "White",
                          RACE == 4 ~ "American Indian",
                          RACE == 99 ~ "Unknown")) %>% 
  group_by(Macrogroup,RACE) %>% 
  count() %>% 
  group_by(Macrogroup) %>% 
  mutate(prop = n / sum(n) * 100) %>% 
  mutate(RACE = factor(RACE,levels = c("White","African American","Asian","American Indian","Unknown"))) %>% 
  ggbarplot("Macrogroup","prop",fill = "RACE",color = "RACE",palette = c("#6c8bba","#36b0fb","#24e077","#6866c9","#fd906e"),
            xlab = "Mitochondrial haplogroup",ylab = "Race %",legend = "right",legend.title = "Race")

p1 %>% insert_bottom(p2) %>% insert_bottom(p3,height = .4) # 7*5

# Figure 5C ---------------------------------------------------------------

# left
data <- Vmtrna_Igtex_T_HFN %>% 
  filter(tissue == "Breast-MammaryTissue") %>%
  left_join(Igtex) %>% 
  mutate(Macrogroup = ifelse(SEX == "M","Male","Female"))

res_lm <- lapply(c("Male","Female"),function(x){
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
formula <- c(": + 0.02 /y",": + 0.02 /y")
label <- str_c(c("Male","Female"),formula)
ggplot(data,aes(AGE,HFN,color = Macrogroup,fill = Macrogroup))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = Macrogroup,fill = Macrogroup,shape = Macrogroup),
             size = 1.5, alpha = 0.7)+
  scale_fill_manual(values = c("#b5d1ff","#657691")[c(2,1)],labels = label[c(2,1)])+
  scale_color_manual(values = c("#b5d1ff","#657691")[c(2,1)],labels = label[c(2,1)])+
  xlab("Age")+
  ylab("Number of mutations")+
  ggtitle("Breast - Mammary tissue")+
  coord_cartesian(ylim = c(0,12))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_none())+
  theme_pubr()+
  theme(legend.position =c(.25,.9),
        plot.title = element_text(hjust = .5)) # 4*5

# right
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

tissues <- c("Ovary","Testis")
data_plot <- data_plot %>% 
  filter(tissue_site_detail %in% tissues) %>% 
  mutate(tissue_site_detail = factor(tissue_site_detail,levels = c("Testis","Ovary")))

formula <- res_lm %>% 
  left_join(gtex_color) %>% 
  filter(tissue_site_detail %in% tissues) %>% 
  arrange(tissue_site_detail) %>% 
  group_by(tissue) %>% 
  mutate(across(where(is.numeric),~ round(.,2))) %>% 
  mutate(formula = str_c(" = ",intercept," + ",beta," × Age")) %>% 
  pull(formula)
formula <- c(": + 0 /y",": + 0.14 /y")
label <- str_c(tissues,formula)

ggplot(data_plot,aes(AGE,HFN,color = tissue_site_detail,fill = tissue_site_detail))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = tissue_site_detail,fill = tissue_site_detail,shape = tissue_site_detail),
             size = 1.5, alpha = 0.7)+
  scale_fill_manual(values = c("#FFAAFF","#AAAAAA")[c(2,1)],labels = label[c(2,1)])+
  scale_color_manual(values = c("#FFAAFF","#AAAAAA")[c(2,1)],labels = label[c(2,1)])+
  xlab("Age")+
  ylab("Number of mutations")+
  ggtitle("Reproductive organs")+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_none())+
  theme_pubr()+
  # theme(legend.position = "right")
  theme(legend.position =c(.25,.9),
        plot.title = element_text(hjust = .5)) # 4*5

# Figure 5D ---------------------------------------------------------------

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

p1 <- ggbarplot(data_up.df %>% 
                  filter(type_class == "Mitochondrial central dogma"),"LOEUF_group","n",fill = "Inheritance",color = "Inheritance",
                facet.by = "type_class",scales = "free_y",legend = "none",
                palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                xlab = "LOEUF Decile",ylab = "# Genes")+
  coord_cartesian(ylim = c(0,40))+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(),
        axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
p2 <- ggbarplot(data_up.df %>% 
                  filter(type_class == "Nuclear cellular machinery"),"LOEUF_group","n",fill = "Inheritance",color = "Inheritance",
                facet.by = "type_class",scales = "free_y",legend = "none",
                palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                xlab = "LOEUF Decile",ylab = "# Genes")+
  coord_cartesian(ylim = c(0,600))+
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

p3 <- ggbarplot(data_up.df %>% 
                  filter(type_class == "Mitochondrial central dogma"),"max_aFC_group","n",fill = "Inheritance",color = "Inheritance",
                facet.by = "type_class",scales = "free_y",legend = "",
                palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                xlab = "log2(aFC) Decile",ylab = "# Genes")+
  coord_cartesian(ylim = c(0,40))+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(),
        axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
p4 <- ggbarplot(data_up.df %>% 
                  filter(type_class == "Nuclear cellular machinery"),"max_aFC_group","n",fill = "Inheritance",color = "Inheritance",
                facet.by = "type_class",scales = "free_y",legend = "",
                palette = rev(c("#BB8588","#A5A58D","#669BBC","#CED4DA")),
                xlab = "log2(aFC) Decile",ylab = "# Genes")+
  coord_cartesian(ylim = c(0,600))+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(),
        axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))

aplot::plot_list(gglist = list(p1,p3,p2,p4),ncol = 2) # 8*6.5

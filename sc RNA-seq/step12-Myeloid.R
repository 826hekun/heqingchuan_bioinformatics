rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')
dir.create("12-T")
setwd("12-T/") 
getwd()
set.seed(12345)
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
table(sce.all$celltype)
sce1 = sce.all[, sce.all$celltype %in% c( 'Myeloid' ,'Mast')]
LayerData(sce1, assay = "RNA", layer = "counts")
sce1 <- JoinLayers(sce1)

as.data.frame(sce1@assays$RNA$counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident) 

sce = sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",  
                     scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000)  
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce) 

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.1)
table(sce@meta.data$RNA_snn_res.0.1)  

set.seed(321)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)

sce <- RunUMAP(object = sce, dims = 1:5, do.fast = TRUE)

p = DimPlot(sce,label=T,cols = mycolors) 
p

###先根据文中注释看看情况
genes_to_check = c( 
                   'CCL3L1' ,   #M2
                   'FABP4',  #M1
                   'CPA3' ,'CST3', 'KIT', 'TPSAB1','TPSB2',#肥大
                   'G0S2', 'S100A9','S100A8','CXCL8', #中性粒细胞
                   'S100B' ,  # CD8
                   'S100A9', # CD8
                   'CD68', 'CD163', 
                   'CD14',
                   'MKI67','PCNA', #proliferating
                   "VCAN","FCN1","CD300E", #
                   'TXN') 

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('myeloid_check_markers.pdf',height = 10,width = 7)

##左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)

axes <- result$axes
label <- result$label

umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
umap
ggsave('myeloid_RNA_snn_res.0.1_umap.pdf',width = 9,height = 7)



#####细胞生物学命名
celltype=data.frame(ClusterID=0:9,
                    celltype= 0:9) 


# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c( 1,6 ),2]='Mast'
celltype[celltype$ClusterID %in% c(0),2]='M2'
celltype[celltype$ClusterID %in% c(2,7  ),2]='M1'
celltype[celltype$ClusterID %in% c(4 ),2]='Granulocyte'
celltype[celltype$ClusterID %in% c( 3 ),2]='proliferating'
celltype[celltype$ClusterID %in% c(8,9 ),2]='TXN_DC'
celltype[celltype$ClusterID %in% c(5 ),2]='S100B_DC'

table(sce@meta.data$RNA_snn_res.0.1)
table(celltype$celltype)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)

result <- left_axes(sce)
axes <- result$axes
label <- result$label
celltype_umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.5,
                       group.by = "celltype",label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
celltype_umap
ggsave('myeloid_umap_by_celltype.pdf',width = 9,height = 7)

saveRDS(sce, "myeloidsce_celltype.rds")





###两组细胞比例图
myeloidsce = readRDS( "myeloidsce_celltype.rds")
myeloidsce$group = 'myeloid'
myeloidmeta = myeloidsce@meta.data
TNKsce=readRDS( "Tsce_celltype.rds")
TNKsce$group = 'T&NK'
TNKmeta = TNKsce@meta.data
colnames(myeloidmeta)
colnames(TNKmeta)

meta = rbind(myeloidmeta,TNKmeta)

##堆积柱状图
library(tidyr)
library(reshape2)
library(ggplot2)
tb=table(meta$sample, meta$celltype,meta$group)
head(tb)
library (gplots) 
library(dplyr)
balloonplot(tb)
bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(Var3) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
#write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#BBBE00","#41BED1")
colnames(bar_per)


library(ggthemes)  
library(ggpubr)
dat = bar_per[bar_per$Freq != 0 ,] 
str(dat)
dat$Var1 = factor(dat$Var1,levels = c('AIS' ,'MIA','IAC' ))
dat$Var2 = factor(dat$Var2,levels = c(   "TXN_DC",   "S100B_DC"  , "M1"  ,         
                                        "M2"      ,      "Mast"    , "proliferating",
                                        "Granulocyte" , "CCL4L2_CD8"  , "GZMK_CD8", "memory_CD4" ,  "GNLY_NK"    , 
                                         "Treg"   ,  "Fibro-like T"     ,"exhausted_CD4" , "FGFBP2_NK"   ))
levels(dat$Var2)
dat$Var3
p1 = ggplot(dat, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip()+
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=col)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  facet_grid(~Var3, scales = "free_x", space ="free_x") +
  theme_pubr(base_size = 10) +
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,"pt"),
        legend.position = "right",
        axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
p1

ggsave(filename="prop.pdf",width = 12,height = 7)


setwd('../')

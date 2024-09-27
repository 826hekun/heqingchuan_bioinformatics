rm(list=ls())
dir.create("9-monocle")
setwd('9-monocle/')
library(tidyverse)
library(tinyarray)
library(data.table) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
sce.all = readRDS('../7-epi/epi_sce_celltype.rds')
sce = sce.all
table(sce$celltype)

Idents(sce) = sce$celltype

# 如果用自己电脑，细胞量太大，可以每个细胞亚群抽样 
allCells=names(Idents(sce))
allType = levels(Idents(sce))

# choose_Cells = unlist(lapply(allType, function(x){
#   cgCells = allCells[Idents(sce)== x ]
#   cg=sample(cgCells,10)
#   cg
# }))

# cg_sce = sce[, allCells %in% choose_Cells]
cg_sce = sce
table(Idents(cg_sce))

library(monocle)
packageVersion("monocle")
#monocle构建CDS需要3个矩阵：expr.matrix、pd、fd
# 将Seurat中的对象转换为monocle识别的对象
#cds <- importCDS(GetAssayData(seurat.object))
#选择做拟时序的亚群
Mono_tj<-cg_sce

Mono_matrix<-as(as.matrix(GetAssayData(Mono_tj,slot = "counts")), 'sparseMatrix')
#构建featuredata，一般featuredata需要两个col，一个是gene_id,一个是gene_short_name,row对应counts的rownames
feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann)<-rownames(Mono_matrix)
Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
#Seurat object中的@meta.data一般会存放表型相关的信息如cluster、sample的来源、group等，所以选择将metadata转换为phenodata
sample_ann<-Mono_tj@meta.data
#rownames(sample_ann)<-colnames(Mono_matrix)

Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
#build new cell data set
Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())

#查看phenodata、featuredata
head(pData(Mono.cds))
head(fData(Mono.cds))
#预处理
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
#筛选基因,这里可以根据自己的需要筛选特定的基因
disp_table <- dispersionTable(Mono.cds)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)
#用DDRtree 进行降维分析
Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 2,
  method = 'DDRTree')
#计算psudotime值
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))

plot_cell_trajectory(Mono.cds,cell_size = 1)

p1 = plot_cell_trajectory(Mono.cds,color_by="celltype", size=1,show_backbone=TRUE)

p2 = plot_cell_trajectory(Mono.cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

p1+p2
ggsave('monocle.pdf',width = 12,height = 8)

###分面显示
plot_cell_trajectory(Mono.cds, color_by = "celltype") + facet_wrap("~sample", nrow = 1)
ggsave('monocle_sample.pdf',width = 14,height = 8)


#⚠️使用root_state参数可以设置拟时间轴的根
cds <- orderCells(Mono.cds, root_state = 3) #把State5设成拟时间轴的起始点
head(pData(cds))
pData(cds)$TGFBR2 = log2( exprs(cds) ['TGFBR2',]+1)
plot_cell_trajectory(cds,cell_size = 1)

p1 = plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE)

p2 = plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

p1+p2



setwd('../')

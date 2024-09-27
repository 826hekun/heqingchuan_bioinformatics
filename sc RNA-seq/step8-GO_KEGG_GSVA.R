rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(clusterProfiler)

source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')
setwd("7-epi/") 
sce=readRDS( "epi_sce_celltype.rds")
sce
Idents(sce)

###热图表示GO富集分析
library(grid)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GOplot)


if (!file.exists('epimarkers.csv')) {
  markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                            min.pct = 0.25, 
                            thresh.use = 0.25)
  write.csv(markers,file='epimarkers.csv')
} else {
  
  markers = read.csv('epimarkers.csv',row.names = 1)
}

#markers = read.csv('markers.csv',row.names = 1)
table(markers$cluster)
library(dplyr) 

gene = markers[markers$p_val <0.01 & markers$avg_log2FC >2,]$gene

entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
gene<- entrezIDs$ENTREZID

marker_new = markers[markers$gene %in% entrezIDs$SYMBOL,]
marker_new = marker_new[!duplicated(marker_new$gene),]
identical(marker_new$gene, entrezIDs$SYMBOL)

p = identical(marker_new$gene, entrezIDs$SYMBOL);p
if(!p) entrezIDs = entrezIDs[match(marker_new$gene,entrezIDs$SYMBOL),]
marker_new$ID = entrezIDs$ENTREZID


## GO
gcSample=split(marker_new$ID, marker_new$cluster) 
###参数可以更改，看看?compareCluster
#One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)

p <- dotplot(xx)
library(ggthemes)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))+theme_few()


###GSVA分析
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

#genesets <- msigdbr(species = "Homo sapiens") 
genesets <- msigdbr(species = "Homo sapiens", category = "C2") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
Idents(sce) <- sce$celltype
expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="gsva") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
gsva_d = gsva.res[sample(nrow(gsva.res),30),]
pheatmap::pheatmap(gsva_d, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


# 气泡图
gsva_long <- melt(gsva_d, id.vars = "Genesets")

# 创建气泡图
ggplot(gsva_long, aes(x = Var2, y = Var1, size = value, color = value)) +
  geom_point(alpha = 0.7) +  # 使用散点图层绘制气泡，alpha设置点的透明度
  scale_size_continuous(range = c(1, 6)) +  # 设置气泡大小的范围
  theme_minimal() + 
  scale_color_gradient(low='#427183',high='#D2D470') +
  labs(x = "Gene Set", y = "Sample", size = "GSVA Score")+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
  


setwd('../')

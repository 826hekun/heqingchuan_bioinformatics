rm(list = ls())  
load(file = 'step4output.Rdata')
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)

#(1)输入数据
gene_diff = deg$ENTREZID[deg$change != "stable"] 

#(2)富集
ekk <- enrichKEGG(gene = gene_diff,organism = 'hsa')
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
ego <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,
                ont = "ALL",readable = TRUE)
#setReadable和readable = TRUE都是把富集结果表格里的基因名称转为symbol
class(ekk)

#(3)可视化
barplot(ego, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
barplot(ekk)
#或者是dotplot

# 更多资料---
# GSEA：https://www.yuque.com/docs/share/a67a180f-dd2b-4f6f-96c2-68a4b86fe862?#
# Y叔的书：http://yulab-smu.top/clusterProfiler-book/index.html
# GOplot：https://mp.weixin.qq.com/s/LonwdDhDn8iFUfxqSJ2Wew
# 网上的资料和宝藏无穷无尽，学好R语言慢慢发掘~

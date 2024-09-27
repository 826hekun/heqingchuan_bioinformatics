rm(list = ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)

 
## Error in download.KEGG.Path(species)
# https://github.com/YuLab-SMU/clusterProfiler/pull/471
getOption("clusterProfiler.download.method")
#R.utils::setOption("clusterProfiler.download.method",'auto')
options(clusterProfiler.download.method = "wininet")
#options(clusterProfiler.download.method = "wget")
getOption("clusterProfiler.download.method")



# 读取差异分析结果
load(file = "data/Step03-edgeR_nrDEG.Rdata")
ls()

# 提取所有差异表达的基因名
DEG <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",2]
head(DEG)

## ===GO数据库, 输出所有结果，后续可根据pvalue挑选结果
ego_CC <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="CC", pvalueCutoff= 1,qvalueCutoff= 1)
ego_MF <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="MF", pvalueCutoff= 1,qvalueCutoff= 1)
ego_BP <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 1,qvalueCutoff= 1)

p_BP <- barplot(ego_BP,showCategory = 10,label_format = 100) + ggtitle("Biological process")
p_CC <- barplot(ego_CC,showCategory = 10,label_format = 100) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10,label_format = 100) + ggtitle("Molecular function")

ego_BP <- data.frame(ego_BP)
ego_CC <- data.frame(ego_CC)
ego_MF <- data.frame(ego_MF)
write.csv(ego_BP,'result/6.enrichGO_BP.csv')
write.csv(ego_CC,'result/6.enrichGO_CC.csv')
write.csv(ego_MF,'result/6.enrichGO_MF.csv')



## === KEGG
genelist <- bitr(gene=DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1)
p1 <- barplot(ekegg, showCategory=10,label_format = 100) + ggtitle("KEGG Pathway")
p2 <- dotplot(ekegg, showCategory=10,label_format=100)
plotc = p1/p2
plotc
ggsave('result/6.enrichKEGG.png', plot = plotc, width = 8, height = 10)

ekegg <- data.frame(ekegg)
write.csv(ekegg,'result/6.enrichKEGG.csv')

library(patchwork)
p <- (p_BP/p_CC) |(p_MF/p1) 
p
ggsave('result/6.enrichKEGG_GO.png', plot = p, width = 16, height = 9)


## === 其他数据库通路
geneset <- read.gmt("data/MsigDB/v7.4/h.all.v7.4.symbols.gmt")
table(geneset$term)
geneset$term <- gsub(pattern = "HALLMARK_","", geneset$term)
geneset$term <- str_to_title(geneset$term)

my_path <- enricher(gene=DEG, pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE=geneset)
p1 <- barplot(my_path, showCategory=15,color = "pvalue")
p1
ggsave("result/6.enrich_HALLMARK.png", plot = p1, width = 8, height = 7)
  
my_path <- data.frame(my_path)
write.csv(my_path,"result/6.enrich_HALLMARK.csv")
  


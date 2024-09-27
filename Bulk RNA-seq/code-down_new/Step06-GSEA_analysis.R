# 清空当前环境变量
rm(list = ls())
options(stringsAsFactors = F)
 
# 加载包
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stats)

# 加载数据
load("data/Step03-edgeR_nrDEG.Rdata")
DEG <- DEG_edgeR_symbol

## 构造GSEA分析数据
# 去掉没有配对上symbol的行
DEG <- DEG[!is.na(DEG$SYMBOL),]

# 去掉重复行
DEG <- DEG[!duplicated(DEG$SYMBOL),]

geneList <- DEG$logFC
names(geneList) <- DEG$SYMBOL
head(geneList)
geneList <- sort(geneList,decreasing = T)
head(geneList)
tail(geneList)


# 选择gmt文件
geneset <- read.gmt("data/MsigDB/h.all.v2023.2.Hs.symbols.gmt")

# 运行,输出全部结果
egmt <- GSEA(geneList, TERM2GENE=geneset, pvalueCutoff = 1)
colnames(egmt@result)

#出点图
dotplot(egmt, label_format = 100)

#按p值出点图
dotplot(egmt,color="qvalue")   

# 单个通路图
# 按照通路名
gseaplot2(egmt, "HALLMARK_ADIPOGENESIS",  
          title = "HALLMARK_ADIPOGENESIS")
# 按照行数
gseaplot2(egmt, 5, color="red", pvalue_table = T)

#按第一到第十个出图，不显示p值
gseaplot2(egmt, 1:6, color="red") 


# 保存结果
go_gsea <- as.data.frame(egmt@result)
write.csv(go_gsea,"result/6.gsea_go_fc.csv",row.names = F)



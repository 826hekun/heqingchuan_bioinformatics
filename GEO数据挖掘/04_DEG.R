rm(list = ls()) 
load(file = "step2output.Rdata")
#差异分析
library(limma)
design = model.matrix(~Group)
fit = lmFit(exp,design)
fit = eBayes(fit)
deg = topTable(fit,coef = 2,number = Inf)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg = mutate(deg,probe_id = rownames(deg))
#2.加上探针注释
ids = distinct(ids,symbol,.keep_all = T)
#其他去重方式在zz.去重方式.R
deg = inner_join(deg,ids,by="probe_id")
nrow(deg) #如果行数为0就是你找的探针注释是错的。

#3.加change列,标记上下调基因
logFC_t = 1
p_t = 0.05
#思考，如何使用padj而非p值
k1 = (deg$P.Value < p_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < p_t)&(deg$logFC > logFC_t)
deg = mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)
#火山图
library(ggplot2)
ggplot(data = deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",linewidth=0.8) +
  geom_hline(yintercept = -log10(p_t),lty=4,col="black",linewidth=0.8) +
  theme_bw()

# 差异基因热图----
# 表达矩阵行名替换为基因名
exp = exp[deg$probe_id,]
rownames(exp) = deg$symbol
diff_gene = deg$symbol[deg$change !="stable"]
n = exp[diff_gene,]
library(pheatmap)
annotation_col = data.frame(group = Group)
rownames(annotation_col) = colnames(n) 
pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         #cluster_cols = F, 
         annotation_col=annotation_col,
         breaks = seq(-3,3,length.out = 100)
) 

#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e = bitr(deg$symbol, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)#人类,注意物种
#一部分基因没匹配上是正常的。<30%的失败都没事。
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
nrow(deg)
deg = inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
#多了几行少了几行都正常，SYMBOL与ENTREZID不是一对一的。
nrow(deg)
save(exp,Group,deg,logFC_t,p_t,file = "step4output.Rdata")

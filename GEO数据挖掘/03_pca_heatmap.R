rm(list = ls())  
load(file = "step2output.Rdata")
#输入数据：exp和Group
#Principal Component Analysis
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials

# 1.PCA 图----
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Group, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
# 2.top 1000 sd 热图---- 
g = names(tail(sort(apply(exp,1,sd)),1000)) #day7-apply的思考题
n = exp[g,]
library(pheatmap)
annotation_col = data.frame(row.names = colnames(n),
                            Group = Group)
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row", #按行标准化，只保留行内差别，不保留行间差别，会把数据范围缩放到大概-5~5之间
         breaks = seq(-3,3,length.out = 100) #设置色带分布范围为-3~3之间，超出此范围的数字显示极限颜色
         ) 

# 关于scale的进一步学习：zz.scale.R


rm(list = ls())
load("step4output.Rdata")
gene_up= deg[deg$change == 'up','symbol'] 
gene_down=deg[deg$change == 'down','symbol']
gene_diff = c(gene_up,gene_down)

# 1.制作string的输入数据
write.table(gene_diff,
            file="diffgene.txt",
            row.names = F,
            col.names = F,
            quote = F)
# 从string网页获得string_interactions.tsv
# 2.准备cytoscape的输入文件
p = deg[deg$change != "stable",
        c("symbol","logFC")]
head(p)
write.table(p,
            file = "deg.txt",
            sep = "\t",
            quote = F,
            row.names = F)
# string_interactions.tsv是网络文件
# deg.txt是属性表格
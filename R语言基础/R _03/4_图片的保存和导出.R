#图片保存的三种方法

#1.基础包作图的保存
pdf("iris_box_ggpubr.pdf")
boxplot(iris[,1]~iris[,5])
text(6.5,4, labels = 'hello')
dev.off()

#2.ggplot系列图(包括ggpubr)通用的简便保存 ggsave
p <- ggboxplot(iris, x = "Species", 
               y = "Sepal.Length",
               color = "Species", 
               shape = "Species",
               add = "jitter")
ggsave(p,filename = "iris_box_ggpubr.png")

#3.eoffice包 导出为ppt,全部元素都是可编辑模式
library(eoffice)
topptx(p,"iris_box_ggpubr.pptx")

#https://mp.weixin.qq.com/s/p7LLLvzR5LPgHhuRGhYQBQ

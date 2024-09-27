# ggpubr 搜代码直接用，基本不需要系统学习

# sthda上有大量ggpubr出的图
library(ggpubr)
ggscatter(iris,x="Sepal.Length",
          y="Petal.Length",
          color="Species")

p <- ggboxplot(iris, x = "Species", 
               y = "Sepal.Length",
               color = "Species", 
               shape = "Species",
               add = "jitter")
p
my_comparisons <- list( c("setosa", "versicolor"), 
                        c("setosa", "virginica"), 
                        c("versicolor", "virginica") )
p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 9) 


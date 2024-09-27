#作图分三类

#1.基础包 略显陈旧 了解一下
plot(iris[,1],iris[,3],col = iris[,5]) 
text(6.5,4, labels = 'hello')

dev.off() #关闭画板

#2.ggplot2 中坚力量，语法有个性
library(ggplot2)
ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species))

#3.ggpubr 新手友好型 ggplot2简化和美化 褒贬不一
library(ggpubr)
ggscatter(iris,
          x="Sepal.Length",
          y="Petal.Length",
          color="Species")

library(ggplot2)
#1.入门级绘图模板：作图数据，横纵坐标
ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length))
#2.属性设置（颜色、大小、透明度、点的形状，线型等）

#2.1 手动设置，需要设置为有意义的值

ggplot(data = iris) + 
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length), 
             color = "blue")

ggplot(data = iris) + 
  geom_point(mapping = aes(x = Sepal.Length, y = Petal.Length), 
             size = 5,     # 点的大小5mm
             alpha = 0.5,  # 透明度 50%
             shape = 8)  # 点的形状

#2.2 映射：按照数据框的某一列来定义图的某个属性
ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species))

## Q1 能不能自行指定映射的具体颜色？

ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species))+
  scale_color_manual(values = c("blue","grey","red"))

## Q2 区分color和fill两个属性
### Q2-1 空心形状和实心形状都用color设置颜色
ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species),
             shape = 17) #17号，实心的例子

ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species),
             shape = 2) #2号，空心的例子
### Q2-2 既有边框又有内心的，才需要color和fill两个参数

ggplot(data = iris)+
  geom_point(mapping = aes(x = Sepal.Length,
                           y = Petal.Length,
                           color = Species),
             shape = 24,
             fill = "black") #24号，双色的例子

#3.分面
ggplot(data = iris) + 
  geom_point(mapping = aes(x = Sepal.Length, y = Petal.Length)) + 
  facet_wrap(~ Species) 
#双分面
dat = iris
dat$Group = sample(letters[1:5],150,replace = T)
ggplot(data = dat) + 
  geom_point(mapping = aes(x = Sepal.Length, y = Petal.Length)) + 
  facet_grid(Group ~ Species) 

#4.几何对象

#局部设置和全局设置

ggplot(data = iris) + 
  geom_smooth(mapping = aes(x = Sepal.Length, 
                          y = Petal.Length))+
  geom_point(mapping = aes(x = Sepal.Length, 
                           y = Petal.Length))

ggplot(data = iris,mapping = aes(x = Sepal.Length, y = Petal.Length))+
  geom_smooth()+
  geom_point()

#5.统计变换-直方图
View(diamonds)
table(diamonds$cut)

ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut))

ggplot(data = diamonds) + 
  stat_count(mapping = aes(x = cut))

#统计变换使用场景
#5.1.不统计，数据直接做图
fre = as.data.frame(table(diamonds$cut))
fre

ggplot(data = fre) +
  geom_bar(mapping = aes(x = Var1, y = Freq), stat = "identity")
#5.2count改为prop
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, y = ..prop.., group = 1))


#6.位置关系

# 6.1抖动的点图
ggplot(data = iris,mapping = aes(x = Species, 
                                 y = Sepal.Width,
                                 fill = Species)) + 
  geom_boxplot()+
  geom_point()

ggplot(data = iris,mapping = aes(x = Species, 
                                 y = Sepal.Width,
                                 fill = Species)) + 
  geom_boxplot()+
  geom_jitter()

# 6.2堆叠直方图
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut,fill=clarity))

# 6.3 并列直方图
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "dodge")

#7.坐标系

#翻转coord_flip()

ggplot(data = mpg, mapping = aes(x = class, y = hwy)) + 
  geom_boxplot() +
  coord_flip()
#极坐标系coord_polar()
bar <- ggplot(data = diamonds) + 
  geom_bar(
    mapping = aes(x = cut, fill = cut), 
    width = 1
  ) + 
  theme(aspect.ratio = 1) +
  labs(x = NULL, y = NULL)
bar
bar + coord_flip()
bar + coord_polar()



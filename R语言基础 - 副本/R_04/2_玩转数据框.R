test <- iris[c(1:2,51:52,101:102),]
rownames(test) =NULL # 去掉行名，NULL是“什么都没有”
test

# arrange，数据框按照某一列排序

library(dplyr)
arrange(test, Sepal.Length) #从小到大
arrange(test, desc(Sepal.Length)) #从大到小

# distinct，数据框按照某一列去重复
distinct(test,Species,.keep_all = T)

# mutate，数据框新增一列
mutate(test, new = Sepal.Length * Sepal.Width)

# 连续的步骤

# 1.多次赋值，产生多个变量

x1 = filter(iris,Sepal.Width>3)
x2 = select(x1, Sepal.Length,Sepal.Width)
x3 = arrange(x2,Sepal.Length)

# 2.管道符号传递，简洁明了
x = iris %>% 
  filter(Sepal.Width>3) %>% 
  select(Sepal.Length,Sepal.Width)%>%
  arrange(Sepal.Length)

# 3. 嵌套，代码不易读
arrange(select(filter(iris,Sepal.Width>3),
               Sepal.Length,Sepal.Width),
        Sepal.Length)

rm(list = ls())
## apply()族函数

### 1.apply 处理矩阵或数据框

#apply(X, MARGIN, FUN, …) 
#其中X是数据框/矩阵名；
#MARGIN为1表示行，为2表示列，FUN是函数

test<- iris[1:6,1:4]

apply(test, 2, mean)

apply(test, 1, sum)

### 2.lapply(list, FUN, …) 
# 对列表/向量中的每个元素（向量）实施相同的操作

test <- list(x = 36:33,y = 32:35,z = 30:27);test

#返回值是列表，对列表中的每个元素（向量）求均值(试试方差var,分位数quantile)

lapply(test,mean)
lapply(test,fivenum)
### 3.sapply 简化结果，返回矩阵或向量

sapply(test,mean)
sapply(test,fivenum)

class(sapply(test,fivenum))


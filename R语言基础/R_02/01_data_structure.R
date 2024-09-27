#重点：数据框
#1.数据框来源
# （1）用代码新建
# （2）由已有数据转换或处理得到
# （3）读取表格文件
# （4）R语言内置数据

#2.新建和读取数据框
df1 <- data.frame(gene   = paste0("gene",1:4),
                 change  = rep(c("up","down"),each = 2),
                 score   = c(5,3,-2,-4))
df1

df2 <- read.csv("gene.csv")
df2

#3.数据框属性
#
dim(df1)
nrow(df1)
ncol(df1)
#
rownames(df1)
colnames(df1)

#4.数据框取子集
df1$gene  #删掉score，按tab键试试
mean(df1$score)

## 按坐标
df1[2,2]
df1[2,]
df1[,2]
df1[c(1,3),1:2]
## 按名字
df1[,"gene"]
df1[,c('gene','change')]
## 按条件（逻辑值）
df1[df1$score>0,]

## 代码思维
#如何取数据框的最后一列？
df1[,3]
df1[,ncol(df1)]
#如何取数据框除了最后一列以外的其他列？
df1[,-ncol(df1)]

#筛选score > 0的基因
df1[df1$score > 0,1]
df1$gene[df1$score > 0]

#5.数据框修改

#改一个格
df1[3,3] <- 5
df1
#改一整列
df1$score <- c(12,23,50,2)     
df1
#？
df1$p.value <- c(0.01,0.02,0.07,0.05) 
df1

#改行名和列名
rownames(df1) <- c("r1","r2","r3","r4")
#只修改某一行/列的名
colnames(df1)[2] <- "CHANGE"

#6.两个数据框的连接
test1 <- data.frame(name = c('jimmy','nicker','Damon','Sophie'), 
                    blood_type = c("A","B","O","AB"))
test1
test2 <- data.frame(name = c('Damon','jimmy','nicker','tony'),
                    group = c("group1","group1","group2","group2"),
                    vision = c(4.2,4.3,4.9,4.5))
test2

test3 <- data.frame(NAME = c('Damon','jimmy','nicker','tony'),
                    weight = c(140,145,110,138))
test3
merge(test1,test2,by="name")
merge(test1,test3,by.x = "name",by.y = "NAME")

##### 矩阵和列表
m <- matrix(1:9, nrow = 3)
colnames(m) <- c("a","b","c") #加列名
m
m[2,]
m[,1]
m[2,3]
m[2:3,1:2]
m
t(m)
as.data.frame(m)
#列表
l <- list(m1 = matrix(1:9, nrow = 3),
          m2 = matrix(2:9, nrow = 2))
l

l[[2]]
l$m1

# 补充：元素的名字

scores = c(100,59,73,95,45)
names(scores) = c("jimmy","nicker","Damon","Sophie","tony")
scores
scores["jimmy"]
scores[c("jimmy","nicker")]

names(scores)[scores>60]

# 删除 
rm(l)
rm(df1,df2)
rm(list = ls()) 

# match练习题
load("matchtest.Rdata")

# 1.生成1到15之间所有偶数
seq(from = 1,to = 15,by = 2)
seq(from = 2,to = 15,by = 2)
# 2.生成向量，内容为："student2"  "student4"  "student6"  "student8"  "student10" "student12"
# "student14" 
# 提示：paste0
paste0(rep("student",times = 7),seq(from = 2, to = 15,by = 2))
# 3.将两种不同类型的数据用c()组合在一起，看输出结果
c(1,"a")
c(TRUE,"a")
c(1,TRUE)

# 4.用函数计算向量g的长度
load("gands.Rdata")
length(g)
# 5.筛选出向量g中下标为偶数的基因名。
seq(2,100,2)
g[seq(2,100,2)]
# 6.向量g中有多少个元素在向量s中存在(要求用函数计算出具体个数)？将这些元素筛选出来
# 提示：%in%
table(g %in% s)
g[g %in% s]
# 7.生成10个随机数: rnorm(n=10,mean=0,sd=18)，用向量取子集的方法，取出其中小于-2的值
z = rnorm(n=10,mean=0,sd=18)
z
z[z<-2]
z

z = rnorm(n=10,mean=0,sd=18)
z
z[z< -2]
z[z<(-2)]

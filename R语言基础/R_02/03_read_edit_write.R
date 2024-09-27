#文件读写部分
#1.读取ex1.txt
ex1 <- read.table("ex1.txt")
ex1 <- read.table("ex1.txt",header = T)
#2.读取ex2.csv
ex2 <- read.csv("ex2.csv")
ex2 <- read.csv("ex2.csv",row.names = 1,check.names = F)

#注意：数据框不允许重复的行名
rod = read.csv("rod.csv",row.names = 1)
rod = read.csv("rod.csv")

#3.读取soft.txt
soft <- read.table("soft.txt")
soft <- read.table("soft.txt",header = T,fill = T) #其实不对
soft2 <- read.table("soft.txt",header = T,sep = "\t")

#4.soft 的行数列数是多少？列名是什么
dim(soft)
colnames(soft)
#5.将soft导出为csv
write.csv(soft,file = "soft.csv")
#6.将soft保存为Rdata并加载。
save(soft,file = "soft.Rdata")
rm(list = ls())
load(file = "soft.Rdata")


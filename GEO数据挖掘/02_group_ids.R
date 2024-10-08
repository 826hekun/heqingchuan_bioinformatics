# Group(实验分组)和ids(探针注释)
rm(list = ls())  
load(file = "step1output.Rdata")
# 1.Group----
library(stringr)
# 标准流程代码是二分组，多分组数据的分析后面另讲
# 生成Group向量的三种常规方法，三选一，选谁就把第几个逻辑值写成T，另外两个为F。如果三种办法都不适用，可以继续往后写else if
if(F){
  # 第一种方法，有现成的可以用来分组的列
  
}else if(F){
  # 第二种方法，眼睛数，自己生成
  Group = rep(c("Disease","Normal"),each = 10)
  # rep函数的其他用法？相间、两组的数量不同？
}else if(T){
  # 第三种方法，使用字符串处理的函数获取分组
  k = str_detect(pd$title,"Normal");table(k)
  Group = ifelse(k,"Normal","Disease")
}

# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,levels = c("Normal","Disease"))
Group

#2.探针注释的获取-----------------
#捷径
library(tinyarray)
find_anno(gpl_number) #辅助写出找注释的代码
library(hgu133plus2.db);ids <- toTable(hgu133plus2SYMBOL) 
#这是从28行运行结果里复制下来的代码👆
#如果能打出代码就不需要再管其他方法。
#如果使用复制下来的AnnoProbe::idmap('xxx')代码发现报错了，请注意尝试不同的type参数
#如果显示no annotation avliable in Bioconductor and AnnoProbe则要去GEO网页上看GPL表格里找啦。

#四种方法，方法1里找不到就从方法2找，以此类推。
#捷径里面包含了全部的R包、一部分表格、一部分自主注释
#方法1 BioconductorR包(最常用，已全部收入find_anno里面，不用看啦)
if(F){
  gpl_number #看看编号是多少
  #http://www.bio-info-trainee.com/1399.html #在这里搜索，找到对应的R包
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db") #列出R包里都有啥
  ids <- toTable(hgu133plus2SYMBOL) #把R包里的注释表格变成数据框
}
# 方法2 读取GPL网页的表格文件，按列取子集
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
# 方法3 官网下载注释文件并读取
# 方法4 自主注释，了解一下
#https://mp.weixin.qq.com/s/mrtjpN8yDKUdCSvSUuUwcA
save(exp,Group,ids,file = "step2output.Rdata")

# 1.match-----
load("matchtest.Rdata")
x
y
## 如何把y的列名正确替换为x里面的ID？

## (1)分步解法
a = colnames(y)
b = x$file_name
k = match(a,b);k
#match(a,b)的意思是a里的每个元素在b的第几个位置上。
#是b的下标，可以给b取子集，也可以给与b对应的其他向量取子集。
colnames(y) = x$ID[k]

## (2)一步解法
load("matchtest.Rdata")
colnames(y) = x$ID[match(colnames(y),x$file_name)]

## (3)放弃match的解法
load("matchtest.Rdata")
rownames(x) = x$file_name
x = x[colnames(y),]
colnames(y) = x$ID

# 2.一些搞文件的函数----
dir() # 列出工作目录下的文件
dir(pattern = ".R$") #列出工作目录下以.R结尾的文件

file.create("douhua.txt") #用代码创建文件
file.exists("douhua.txt") #某文件在工作目录下是否存在
file.remove("douhua.txt") #用代码删除文件
file.exists("douhua.txt") #删掉了就不存在啦

## 可以批量的新建和删除
f = paste0("douhua",1:100,".txt")
file.create(f)
file.remove(f)

rm(list = ls())

## 一.条件语句

###1.if(){ }

#### (1)只有if没有else，那么条件是FALSE时就什么都不做

i = -1
if (i<0) print('up')
if (i>0) print('up')

#理解下面代码
if(!require(tidyr)) install.packages('tidyr')

#### (2)有else
i =1
if (i>0){
  print('+')
} else {
  print("-")
}
i = 1
ifelse(i>0,"+","-")

x = rnorm(3)
x
ifelse(x>0,"+","-")

#ifelse()+str_detect(),王炸
samples = c("tumor1","tumor2","tumor3","normal1","normal2","normal3")
k1 = str_detect(samples,"tumor");k1
ifelse(k1,"tumor","normal")
k2 = str_detect(samples,"normal");k2
ifelse(k2,"normal","tumor")

#### (3)多个条件
i = 0
if (i>0){
  print('+')
} else if (i==0) {
  print('0')
} else if (i< 0){
  print('-')
}

ifelse(i>0,"+",ifelse(i<0,"-","0"))

## 二、循环语句

### 1.for循环
x <- c(5,6,0,3)
s=0
for (i in x){
  s=s+i
  print(c(i,s))
}

x <- c(5,6,0,3)
s = 0
for (i in 1:length(x)){
  s=s+x[[i]]
  print(c(x[[i]],s))
}

#如何将结果存下来?
s = 0
result = list()
for(i in 1:length(x)){
  s=s+x[[i]]
  result[[i]] = c(x[[i]],s)
}
result
do.call(cbind,result)


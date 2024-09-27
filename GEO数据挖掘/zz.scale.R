require(stats)
# 1.示例数据
x <- matrix(sample(1:30,30), ncol = 6)
rownames(x) = paste0("gene",1:5)
colnames(x) = paste0("sample",1:6)

# 2.标准化
scale(x) #函数只能按列标准化，但是我们需要按行
x
y = t(scale(t(x)))

# 3.标准化前后，某gene的表达量点图比较，大小趋势不变。

par(mfrow = c(2,2))
plot(x[1,])
plot(y[1,])
plot(x[2,])
plot(y[2,])

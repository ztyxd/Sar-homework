# 引用相关的包，就下面几个，source的几个包是老师自己写的
source("myread.ENVI.R")
source("imagematrix.R")
require(ggplot2)
require(reshape2)
require(ggthemes)
# 导入数据，只有一个HH波段的数据，数据比较大，群里百度网盘里有数据包，最后上传github时候记得拿出来别一起传上去
# 数据在这 https://pan.baidu.com/s/1yIOIsA1-tNR_KMyN1kqnSA
imagepath <- "../Statistics-SAR-Intensity-master/Data/Images/ESAR/"
HH_Complex <- myread.ENVI(paste(imagepath,
                                "ESAR97HH.DAT", sep = ""), 
                          paste(imagepath, "ESAR97HH.hdr", sep = ""))
# 数据有虚数部分不好用所以给他平方一下
HH_Intensity <- (Mod(HH_Complex))^2

# 数据是4000X1599大小的图片 太大了，取一小块，我这里取的是一块城市的图片
example <- HH_Intensity[1500:1599,1500:1599]


# 把数据拉直成向量好处理
vexample <- data.frame(HH=as.vector(example))


## 显示一下刚选的一小块图片，并且保存到当前路径下面，命名为example.png
plot(imagematrix(equalize(example)))
imagematrixPNG(name = "./example.png", imagematrix(equalize(example)))


# 这个是包装成一个data.frame 方便后面显示，顺便显示了个summary，可以去掉summary
vexample <- data.frame(HH=as.vector(example))
summary(vexample)

## 后面就是显示统计的直方图
binwidth_complete <- 2*IQR(vexample$HH)*length(vexample$HH)^(-1/3)
ggplot(data=vexample, aes(x=HH)) + 
  geom_histogram(aes(y=..density..), 
                 binwidth = binwidth_complete) + 
  xlab("Intensities") +
  ylab("Proportions") +
  ggtitle("Complete Histogram") +
  theme_few()
ggsave(filename = "./HistogramExample.pdf")


## 下面是两种估计方法
require(maxLik)

## 第一种是矩估计
GI0.Estimator.m1m2 <- function(z, L) {
  m1 <- mean(z)
  m2 <- mean(z^2)
  m212 <- m2/m1^2
  
  a <- -2 - (L+1) / (L * m212)
  g <- m1 * (2 + (L+1) / (L * m212))
  
  return(list("alpha"=a, "gamma"=g))
}
# 矩估计结果
estim.example <- GI0.Estimator.m1m2(example, 1)

# 这是最大似然估计
LogLikelihoodLknown <- function(params) {
  
  p_alpha <- -abs(params[1])
  p_gamma <- abs(params[2])
  p_L <- abs(params[3])
  
  n <- length(vexample$HH)
  
  return(
    n*(lgamma(p_L-p_alpha) - p_alpha*log(p_gamma) - lgamma(-p_alpha)) + 
      (p_alpha-p_L)*sum(log(p_gamma + z*p_L)) 
  )
}
# 最大似然估计结果
estim.exampleML <- maxNR(LogLikelihoodLknown, 
                         start=c(estim.example$alpha, estim.example$gamma,1), 
                         activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]

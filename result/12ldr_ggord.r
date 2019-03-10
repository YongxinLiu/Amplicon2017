## LDA排序和添加椭圆置信区间



### ggord绘制LDA排序


# Installation
devtools::install_github('fawda123/ggord')

# 加载lda包
library(MASS)

# 查看测试数据
head(iris)

# 按物种分组排序
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

# 展示LDA分析
library(ggord)
p <- ggord(ord, iris$Species)
p


# 计算置信椭圆函数

get_lda_ell <- function(ord_in, grp_in, ellipse_pro = 0.97){
  ## adapted from https://github.com/fawda123/ggord/blob/master/R/ggord.R
  require(plyr)
  axes = c('LD1', 'LD2')
  obs <- data.frame(predict(ord_in)$x[, axes])
  obs$Groups <- grp_in
  names(obs)[1:2] <- c('one', 'two')
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(obs, 'Groups', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$one, x$two))
    mu <- c(mean(x$one), mean(x$two))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
  })
  names(ell)[2:3] <- c('one', 'two')
  ell <- ddply(ell, .(Groups), function(x) x[chull(x$one, x$two), ])
  ell
}

# 计算置信椭圆，并添加至原图
anotherEll <- get_lda_ell(ord, iris$Species, 0.97)
## Loading required package: plyr
p + geom_polygon(data = anotherEll, 
                 aes_string(color = 'Groups', group = 'Groups'),
                 lty=2, fill = NA)

# 菌群实战数据

# 读入实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 
# 读取OTU表
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")
# 转换原始数据为百分比
norm = t(t(otu_table)/colSums(otu_table,na=T)) * 100 # normalization to total 100
# 按mad值排序取前6波动最大的OTUs
mad.5 = head(norm[order(apply(norm,1,mad), decreasing=T),],n=6)
row.names(mad.5)=c("Streptophyta","Rubrivivax","Methylibium","Streptosporangiaceae","Streptomyces","Niastella")
data=as.data.frame(t(mad.5))
# 添加分组信息
data$group=design[row.names(data),]$genotype

# 按实验基因组分组排序
ord <- lda(group ~ ., data)

# 展示LDA分析
library(ggord)
p <- ggord(ord, data$group, ellipse_pro = 0.68)
p

anotherEll <- get_lda_ell(ord, data$group, 0.90)
p + geom_polygon(data = anotherEll, 
                 aes_string(color = 'Groups', group = 'Groups'),
                 lty=2, fill = NA)


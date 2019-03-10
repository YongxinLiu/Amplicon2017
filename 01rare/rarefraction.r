# 安装和载入相关包和数据文件
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ggplot2","reshape2")) # 没安装过此类包的请手动运行这两行代码安装包
library("ggplot2") # load related packages
library("reshape2")
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 
rare = read.table("alpha_rare.txt", header=T, row.names= 1, sep="\t") 

# 提取样品组信息
sampFile = as.data.frame(design$genotype,row.names = row.names(design))
colnames(sampFile)[1] = "group"

# 转换宽表格为ggplot通用长表格格式
rare$x = rownames(rare) # 添加x轴列
rare_melt = melt(rare, id.vars=c("x")) # 转换为长表格
rare_melt$x = factor(rare_melt$x, levels=1:100) # 设置x轴顺序

# 添加分组信息
rare_melt3 = merge(sampFile,rare_melt, by.x="row.names", by.y="variable")
rare_melt3$variable=rare_melt3$Row.names

# 按样品分组，按组上色
p = ggplot(rare_melt3, aes(x = x, y = value, group = variable, color = group )) + 
  geom_line()+
  xlab("Rarefraction Percentage")+ylab("Richness (Observed OTUs)")+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) + theme_classic()
p
ggsave(paste("alpha_rare_samples.pdf", sep=""), p, width = 8, height = 5)



# 求各组均值
# 读取usearch rarefraction文件，上面己经修改，必须重新读入
rare = read.table("alpha_rare.txt", header=T, row.names= 1, sep="\t") 
mat_t = merge(sampFile, t(rare), by="row.names")[,-1]
# 按第一列合并求均值
mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
# 修正行名
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

rare=as.data.frame(round(mat_mean_final))
rare$x = rownames(rare)
rare_melt = melt(rare, id.vars=c("x"))

# 求各组标准误
# 转置rare表格与实验设计合并，并去除第一列样品名
se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se) # se 为什么全是NA
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

rare_se=as.data.frame(round(mat_se_final))
rare_se$x = rownames(rare_se)
rare_se_melt = melt(rare_se, id.vars=c("x"))

# 添加标准误到均值中se列
rare_melt$se=rare_se_melt$value

rare_melt$x = factor(rare_melt$x, levels=c(1:100))

p = ggplot(rare_melt, aes(x = x, y = value, group = variable, color = variable )) + 
  geom_line()+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
  xlab("Percentage")+ylab("Richness (Observed OTUs)")+
  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+
  theme_classic()
p
ggsave(paste("alpha_rare_groups.pdf", sep=""), p, width = 8, height = 5)

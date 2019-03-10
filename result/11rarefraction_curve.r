if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","reshape2"))
}

# 安装和载入相关包和数据文件
source("https://bioconductor.org/biocLite.R")
biocLite(c("ggplot2","reshape2"))
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
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
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
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
p
ggsave(paste("alpha_rare_groups.pdf", sep=""), p, width = 8, height = 5)





# Install related packages
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","grid","scales","vegan","agricolae","dplyr"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result") # set work directory
library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")
library("agricolae")
library("reshape2")
library("dplyr")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=7),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=7),
                    text=element_text(family="sans", size=7))

# Public file 1. "design.txt"  Design of experiment
design = read.table("/mnt/bai/yongxin/wheat/NP/doc/design.txt", header=T, row.names= 1, sep="\t") 

# setting subset design
if (TRUE){
	sub_design = subset(design,groupID %in% c("BDHNp","JMNp","KN2011Np","RHTANp","XY54Np","RsBDHNp","RsJMNp","RsKN2011Np","RsRHTANp","RsXY54Np","BSNp") ) # select group1
}else{
	sub_design = design
}
if (TRUE){
	sub_design = subset(sub_design,compartment %in% c("rhizosphere","root","soil") ) # select group2
}

# Set group style, single or combine
if (FALSE){
	sub_design$group=paste(sub_design$groupID,sub_design$compartment,sep = "")
}else{
	sub_design$group=sub_design$groupID
}

# Set group order
if ("TRUE" == "TRUE") {
    sub_design$group  = factor(sub_design$group, levels=c("BDHNp","JMNp","KN2011Np","RHTANp","XY54Np","RsBDHNp","RsJMNp","RsKN2011Np","RsRHTANp","RsXY54Np","BSNp"))   # set group order
}
sub_design$compartment  = factor(sub_design$compartment)   # set group order
print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers



# 读取usearch rarefraction文件
rare = read.table("alpha_rare.txt", header=T, row.names= 1, sep="\t") 
# 提取样品组信息
sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"

# 直接展示样品
rare$x = rownames(rare) # 添加x轴列
rare_melt = melt(rare, id.vars=c("x")) # 转换为长表格
rare_melt$x = factor(rare_melt$x, levels=1:100) # 设置x轴顺序

rare_melt3 = merge(sampFile,rare_melt, by.x="row.names", by.y="variable")
rare_melt3$variable=rare_melt3$Row.names

# 按样品分组，按组上色
p = ggplot(rare_melt3, aes(x = x, y = value, group = variable, color = group )) + 
  geom_line()+
  xlab("Rarefraction Percentage")+ylab("Richness (Observed OTUs)")+main_theme+
#  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  #stat_smooth(method="auto", se=FALSE)+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
ggsave(paste("alpha_rare_samples.pdf", sep=""), p, width = 8, height = 5)
ggsave(paste("alpha_rare_samples.png", sep=""), p, width = 8, height = 5)



# 求各组均值
# 读取usearch rarefraction文件，上面己经修改，必须重新读入
rare = read.table("alpha_rare.txt", header=T, row.names= 1, sep="\t") 
# 文件的行名为纯数字，转置将整个矩阵num变为字符chr，需要提前修改行名；是因为上方修改了格式，重读解决
# rownames(rare)=c(paste("a",1:100,sep=""))
# 转置rare表格与实验设计合并，并去除第一列样品名
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
#rare_melt$x = factor(rare_melt$x, levels=1:100)

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
# 去除之前添加的x轴字符，转换为纯数字，第一步样品绘图修改了原始数据，第二步要重新读入
# rare_melt$x1=gsub("a","",rare_melt ,perl=TRUE)
# 添加levels顺序，否则
rare_melt$x = factor(rare_melt$x, levels=c(1:100))

p = ggplot(rare_melt, aes(x = x, y = value, group = variable, color = variable )) + 
  geom_line()+
  #geom_line(stat = "identity",position="identity")+ 
  #stat_smooth(method="auto", se=FALSE)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
  xlab("Percentage")+ylab("Richness (Observed OTUs)")+main_theme+
  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 

ggsave(paste("alpha_rare_groups.pdf", sep=""), p, width = 8, height = 5)
ggsave(paste("alpha_rare_groups.png", sep=""), p, width = 8, height = 5)


print("All done!!!")


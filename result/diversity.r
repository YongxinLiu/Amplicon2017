# 运行前，请在Rstudio中菜单栏选择“Session - Set work directory -- Choose directory”，弹窗选择之前分析目录中的result文件夹

# 安装相关软件包，如果末安装改为TRUE运行即可安装
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","grid","scales","vegan"))
}

# 加载相关软件包
library("ggplot2") # load related packages

# 读入实验设计和Alpha多样性值
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 
alpha = read.table("alpha.txt", header=T, row.names= 1, sep="\t")

# 以Observed OTU为例进行可视化和统计分析，其它指数将observed_otus替换为shannon, chao1, PD_whole_tree即可计算

# 合并Alpha指数与实验设计
index = cbind(alpha, design[match(rownames(alpha), rownames(design)), ]) 
# 绘图代码、预览、保存PDF
p = ggplot(index, aes(x=genotype, y=observed_otus, color=genotype))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="observed_otus index")
p
ggsave(paste("alpha_observed_otus.png", sep=""), p, width = 5, height = 3)

# 统计组间是否显著差异
# anova对指数与分组统计
observed_otus_stats <- aov(observed_otus ~ genotype, data = index)
# 使用TukeyHSD对组间进行检验，效正pvalue
Tukey_HSD_observed_otus <- TukeyHSD(observed_otus_stats, ordered = FALSE, conf.level = 0.95)
# 结果中提取需要的结果
Tukey_HSD_observed_otus_table <- as.data.frame(Tukey_HSD_observed_otus$genotype)
# 预览结果
Tukey_HSD_observed_otus_table
# 保存结果到文件
write.table(Tukey_HSD_observed_otus_table[order(Tukey_HSD_observed_otus_table$p, decreasing=FALSE), ], file="alpha_observed_otus_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)


#############################################################
# Title: Beta diversity - PCoA
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 03/14/2017
# Description: Script to draw beta diversity - scatterplot PCoA
# Version 1.1
# Run enviroment: R3.3.1, ggplot2
# Input File required in the blow list: 
# 1. "beta/bray_curtis_otu_table_css.txt" output by beta_diversity.py from qiime
# 2. "beta/unweighted_unifrac_otu_table_css.txt" 
# 3. "beta/weighted_unifrac_otu_table_css.txt"
#############################################################
# 运行前，请在Rstudio中菜单栏选择“Session - Set work directory -- Choose directory”，弹窗选择之前分析目录中的result文件夹

# 安装相关软件包，如果末安装改为TRUE运行即可安装
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","vegan"))
}

# 加载相关软件包
library("ggplot2")
library("vegan")

# 读入实验设计和Alpha多样性值
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 

# 以bray_curtis为例，其它距离只需替换为weighted_unifrac,unweighted_unifrac
bray_curtis = read.table("beta/bray_curtis_otu_table_css.txt", sep="\t", header=T, check.names=F)

# 过滤数据并排序
idx = rownames(design) %in% colnames(bray_curtis) 
sub_design = design[idx,]
bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# 将距离矩阵进行主坐标轴分析
pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])

# 绘制主坐标准轴的第1，2轴
p = ggplot(points, aes(x=x, y=y, color=genotype)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA")
p

p + geom_text(aes(x, y, label=rownames(points)))+ theme_classic()

library(ggrepel)

p + geom_text_repel(aes(x, y, label=rownames(points)))+ theme_classic()


ggplot(points, aes(x=x, y=y)) +geom_point(alpha=.7, size=2) + geom_label_repel(aes(x, y, fill=factor(genotype), label=rownames(points)))+ theme_classic()


ggsave("beta_pcoa_bray_curtis.pdf", p, width = 5, height = 3)
ggsave("beta_pcoa_bray_curtis.png", p, width = 5, height = 3)

# 统计WT与KO间是否有差异显著
# 取两组实验设计子集
design2 = subset(sub_design, genotype %in% c("WT","KO"))
# 获取对应的子距离矩阵并排序
sub_dis_table = bray_curtis[rownames(design2),rownames(design2)]
# 计算距离矩阵
sub_dis_table <- as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
# 统计按genotype分组下，组间差异的显著性水平；检验10000次
adonis_table = adonis(sub_dis_table~genotype, data=design2, permutations = 10000) 
# 获得pvalue值
adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
# 显示组间的pvalue值
adonis_pvalue



# Function for analysis CPCoA/CCA result

# CCA分析功能函数
variability_table = function(cca){
  chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# 读入CSS标准化的OTU表，并与实验设计比对筛选和数据重排
otu_table = read.table("otu_table_css.txt", sep="\t", header=T, row.names= 1) # CSS norm otu table
idx = rownames(sub_design) %in% colnames(otu_table) 
sub_design = sub_design[idx,]
sub_otu_table = otu_table[, rownames(sub_design)] 

# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(sub_otu_table) ~ genotype, data=sub_design, add=F, sqrt.dist=T, distance="bray") 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen$CCA$wa[, 1:2]
points = as.data.frame(points)
colnames(points) = c("x", "y")
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=genotype)) +
  geom_point(alpha=.7, size=1.5) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) 
p
ggsave(paste( "CPCoA.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste( "CPCoA.png", sep=""), p, width = 5, height = 3)



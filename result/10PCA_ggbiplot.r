# https://github.com/vqv/ggbiplot/blob/master/README.md

# 安装包，安装过请跳过
install.packages("devtools", repo="http://cran.us.r-project.org")
library(devtools)
install_github("vqv/ggbiplot")


# 最简单帅气的例子
data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
# 演示样式
ggbiplot(wine.pca, obs.scale = 1, var.scale = 1,
         groups = wine.class, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


# 菌群数据实战
# 读入实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 

# 读取OTU表
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# 过滤数据并排序
idx = rownames(design) %in% colnames(otu_table) 
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]

# 基于OTU表PCA分析
otu.pca <- prcomp(t(count), scale. = TRUE)

# 绘制PCA图，并按组添加椭圆
ggbiplot(otu.pca, obs.scale = 1, var.scale = 1,
         groups = sub_design$genotype, ellipse = TRUE,var.axes = F)

# 显著高丰度菌的影响

# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100

# 筛选mad值大于0.5的OTU
mad.5 = norm[apply(norm,1,mad)>0.5,]
# 按mad值排序取前6波动最大的OTUs
mad.5 = head(norm[order(apply(norm,1,mad), decreasing=T),],n=6)
row.names(mad.5)
# 手动从result/rep_seqs_tax.txt文件中对应注释相应物种，不同菌分类程度不同，需要手动选择有准确注释的名称
row.names(mad.5)=c("Streptophyta","Rubrivivax","Methylibium","Streptosporangiaceae","Streptomyces","Niastella")
# 计算PCA和菌与菌轴的相关性
otu.pca <- prcomp(t(mad.5))
ggbiplot(otu.pca, obs.scale = 1, var.scale = 1,
         groups = sub_design$genotype, ellipse = TRUE,var.axes = T)



if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2"))
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
ggsave(paste("alpha_observed_otus.pdf", sep=""), p, width = 5, height = 3)

# 统计组间是否显著差异
# anova对指数与分组统计
observed_otus_stats <- aov(observed_otus ~ genotype, data = index)
# 使用TukeyHSD对组间进行检验，效正pvalue
Tukey_HSD_observed_otus <- TukeyHSD(observed_otus_stats, ordered = FALSE, conf.level = 0.95)
# 结果中提取需要的结果
Tukey_HSD_observed_otus_table <- as.data.frame(Tukey_HSD_observed_otus$genotype)
# 预览结果
Tukey_HSD_observed_otus_table
# 保存结果到文件，按Pvaule值由小到大排序
write.table(Tukey_HSD_observed_otus_table[order(Tukey_HSD_observed_otus_table$p, decreasing=FALSE), ], file="alpha_observed_otus_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

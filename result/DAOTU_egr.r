# 运行前，请在Rstudio中菜单栏选择“Session - Set work directory -- Choose directory”，弹窗选择之前分析目录中的result文件夹

# # 安装相关软件包，如果末安装改为TRUE运行即可安装
# if (FALSE){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite(c("ggplot2","vegan"))
# }
# 
# # 加载相关软件包
# library("ggplot2")
# library("vegan")

# 读入实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 

# 读取OTU表
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# 过滤数据并排序
idx = rownames(design) %in% colnames(otu_table) 
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]

# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100

# 计算所有样品间相关系数
sim=cor(norm,method="pearson")

# 使用热图可视化，并保存为8x8英寸的PDF
library("gplots")
library("RColorBrewer")
png(file=paste("heat_cor_samples.png", sep=""), height = 8, width = 8, units = "in", res = 300)

#pdf(file=paste("heat_cor_samples.pdf", sep=""), height = 8, width = 8)
heatmap.2(sim, Rowv=TRUE, Colv=TRUE, dendrogram='both', trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()



# 使用edgeR统计组间差异OTU，以OE vs WT为例

library(edgeR)
# create DGE list
d = DGEList(counts=count, group=sub_design$genotype)
d = calcNormFactors(d)

# 生成实验设计矩阵
design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(design$genotype)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# # function DA edgeR
# sampleA = "OE"
# sampleB = "WT"
# #	design2 = subset(sub_design, genotype %in% c(sampleA,sampleB))
# print(paste("Start DA OTU", sampleA, sampleB, ":",sep=" "))
# SampAvsB=paste(sampleA,"-", sampleB, sep="")
# print(SampAvsB)
# 设置比较组
BvsA <- makeContrasts(contrasts = "OE-WT", levels=design.mat)
# 组间比较,统计Fold change, Pvalue
lrt = glmLRT(fit,contrast=BvsA)
# FDR检验，控制假阳性率小于5%
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

# 导出计算结果
x=lrt$table
x$sig=de_lrt
OE_enriched = row.names(subset(x,sig==1))
length(OE_enriched)
OE_depleted = row.names(subset(x,sig==-1))
length(OE_depleted)

# 同理计算KO vs WT间的差异OTU
BvsA <- makeContrasts(contrasts = "KO-WT", levels=design.mat)
lrt = glmLRT(fit,contrast=BvsA)
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
x=lrt$table
x$sig=de_lrt
KO_enriched = row.names(subset(x,sig==1))
length(KO_enriched)
KO_depleted = row.names(subset(x,sig==-1))
length(KO_depleted)

# 加载包和颜色方案
library(VennDiagram)
color_v <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
# 绘制两组比较Venn
p=venn.diagram( x = list(OE_depleted=OE_depleted, KO_depleted=KO_depleted), filename=NULL,fill = color_v[1:2])
grid.draw(p)
dev.off()
# 绘制三组比较Venn
p=venn.diagram( x = list(OE_enriched=OE_enriched, OE_depleted=OE_depleted, KO_depleted=KO_depleted), filename=NULL,fill = color_v[1:3])
grid.draw(p)
dev.off()
# 绘制四组比较Venn
p=venn.diagram( x = list(OE_enriched=OE_enriched, OE_depleted=OE_depleted, KO_enriched=KO_enriched, KO_depleted=KO_depleted), filename=NULL,fill = color_v[1:4])
grid.draw(p)
dev.off()
# 绘制五组比较Venn，并保存PDF
pdf(file="venn5.pdf", onefile=FALSE, paper="special", width=4, height=4, pointsize=8)
p=venn.diagram( x = list(OE_enriched=OE_enriched, OE_depleted=OE_depleted, Fifth=OE_depleted, KO_enriched=KO_enriched, KO_depleted=KO_depleted), filename=NULL,fill = color_v)
grid.draw(p)
dev.off()

## 三元图
# merge group to mean
## 按样品名合并实验组与转置的OTU
mat_t2 = merge(sub_design[c("genotype")], t(norm), by="row.names")[,-1]
## 按实验设计求组平均值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean 
# 重新转载并去除组名
per3=t(mat_mean[,-1])
colnames(per3) = mat_mean$genotype
per3=as.data.frame(per3[rowSums(per3)>0,]) # remove all 0 OTU
#per3=per3[,tern] # reorder per3 as input
color=c(c_green,c_orange,c_red,c_grey) 


# 两底角相对于顶点显著富集的OTU，分共有和特有，类似维恩图
per3$color=color[4] # set all default # 设置点默认颜色为灰
AvC = KO_enriched
BvC = OE_enriched
C = intersect(row.names(AvC), row.names(BvC))
A = setdiff(AvC, C) 
B = setdiff(BvC, C) 
if (length(A)>0){per3[A,]$color=color[1]} 
if (length(B)>0){per3[B,]$color=color[2]} 
if (length(C)>0){per3[C,]$color=color[3]}
## output pdf and png in 8x8 inches
per3lg=log2(per3[,1:3]*100+1) # 对数变换，差OTU千分比的差距，点大小更均匀
pdf(file=paste("ter_",tern[1],tern[2],tern[3],"venn.pdf", sep=""), height = 8, width = 8)
tern_e(per3lg[,1:3], prop=T, col=per3$color, grid_color="black", labels_color="transparent", pch=19, main="Tenary Plot")
dev.off()


## 热图展示差异OTU
pair_group = subset(sub_design, genotype %in% c("OE", "WT"))
# Sig OTU in two genotype
DE=c(enriched,depleted)
sub_norm = as.matrix(norm[DE, rownames(pair_group)])
#colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
pdf(file=paste("heat_otu_OEvsWT_sig.pdf", sep=""), height = 8, width = 8)
# scale in row, dendrogram only in row, not cluster in column
heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
dev.off()


## Manhattan图展示差异OTU和Taxonomy

# 读取taxonomy，并添加各列名称
taxonomy = read.delim("rep_seqs_tax.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")

# 标记差异OTU类型
x$level = as.factor(ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")))
x$otu = rownames(x)
# 转换Pvalue为负对数
x$neglogp = -log(x$PValue)

# Taxonomy排序，并筛选OTU表中存在的
library(dplyr)
taxonomy$id=rownames(taxonomy)
taxonomy = arrange(taxonomy, phylum, class, order, family, genus, species)
rownames(taxonomy) = taxonomy$id
idx = rownames(taxonomy) %in% x$otu
tax = taxonomy[idx, ] # subset taxonomy from used OTU

# 手动筛选显著的组
x = x[rownames(tax), ] # reorder according to tax
x$tax = gsub("p__","",tax$phylum,perl=TRUE) 
top_phylum=c("Bacteroidetes","Firmicutes","Planctomycetes","Proteobacteria","Verrucomicrobia")
x[!(x$tax %in% top_phylum),]$tax = "Low Abundance" # no level can get value
# 设置各类的level对应顺序
x$otu = factor(x$otu, levels=x$otu)   # set x order
x$level = factor(x$level, levels=c("enriched","depleted","nosig"))
levels(x$tax)=c(top_phylum,"Low Abundance")
# 调整Y轴范围更美观
x[x$neglogp>15,]$neglogp  = 15

# Manhattan plot
## 添加显著阈值线
FDR = min(x$neglogp[x$level=="depleted"])
library(ggplot2)
p = ggplot(x, aes(x=otu, y=neglogp, color=tax, size=logCPM, shape=level)) +
  geom_point(alpha=.7) + 
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(17, 25, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="OTU", y="-loge(P)") +
  theme(axis.line.x=element_line(size=.5, colour="black"),legend.key=element_blank(),axis.line.y=element_line(size=.5, colour="black"),panel.background=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="top")
p
ggsave(file=paste("man_otu.pdf", sep=""), p, width = 10, height = 3, useDingbats=F)
ggsave(file=paste("man_otu.png", sep=""), p, width = 10, height = 3)



theme(panel.background=element_blank(),
      panel.grid=element_blank(),
      
      axis.ticks=element_line(color="black"),
      axis.text=element_text(color="black", size=7),
      legend.position="right",
      legend.background=element_blank(),
      
      legend.text= element_text(size=7),
      text=element_text(family="sans", size=7))

# 绘制火山图
if (max(x$logFC)>4){x[x$logFC>4,]$logFC = 4} # norm x axis
if (min(x$logFC)< -4){x[x$logFC< -4,]$logFC = -4} # norm x axis
x$level = as.factor(ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")))

# Volcanol plot of fold change vs abundance plot
p = ggplot(x, aes(x=logFC, y=logCPM, color=level, size=logCPM, shape=tax)) + geom_point()  + 
  scale_colour_manual(values=c("red","green","grey"))+ xlim(-4, 4)+
  labs(x="log2(fold change)",y="log2(count per million)", title=paste("OE vs WT", sep=" "))
p
ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = 8, height = 8)
ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".png", sep=""), p, width = 8, height = 5)







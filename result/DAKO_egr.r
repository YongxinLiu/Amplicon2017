# 运行前，请在Rstudio中菜单栏选择“Session - Set work directory -- Choose directory”，弹窗选择之前分析目录中的result文件夹

# 安装相关软件包，如果末安装改为TRUE运行即可安装
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","vegan"))
}

# 加载相关软件包
library("ggplot2")
library("vegan")

# 读入实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 

# 读取ko表, 以L3级为例
ko_table = read.delim("metagenome_predictions.L3.txt", row.names= 1,  header=T, sep="\t")

# 过滤数据并排序
idx = rownames(design) %in% colnames(ko_table) 
sub_design = design[idx,]
count = ko_table[, rownames(sub_design)]

# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100

# 计算所有样品间相关系数
sim=cor(norm,method="pearson")

# 使用热图可视化，并保存为8x8英寸的PDF
library("gplots")
library("RColorBrewer")
png(file=paste("KO3_heat_cor_samples.png", sep=""), height = 8, width = 8, units = "in", res = 300)
#pdf(file=paste("KO3_heat_cor_samples.pdf", sep=""), height = 8, width = 8)
heatmap.2(sim, Rowv=TRUE, Colv=TRUE, dendrogram='both', trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()
# 从KO角度，WT和KO区分不明显，而且还分为两类


# 使用edgeR统计组间差异ko，以OE vs WT为例

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

# 设置比较组
BvsA <- makeContrasts(contrasts = "OE-WT", levels=design.mat)
# 组间比较,统计Fold change, Pvalue
lrt = glmLRT(fit,contrast=BvsA)
# FDR检验，控制假阳性率小于5%
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.0000005)

# 导出计算结果
x=lrt$table
x$sig=de_lrt
OE_enriched = row.names(subset(x,sig==1))
length(OE_enriched)
OE_depleted = row.names(subset(x,sig==-1))
length(OE_depleted)

## 热图展示差异ko
pair_group = subset(sub_design, genotype %in% c("OE", "WT"))
# Sig ko in two genotype, 数据太多，修改p 0.05至0.0000005下降到50个以内方便展示
DE=c(OE_enriched,OE_depleted)
sub_norm = as.matrix(norm[DE, rownames(pair_group)])
#colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
pdf(file=paste("heat_ko_OEvsWT_sig.pdf", sep=""), height = 8, width = 8)
# scale in row, dendrogram only in row, not cluster in column
# 文字太多，调整画面、字体大小，让文字显示清楚；
heatmap.2(sub_norm, scale="row", Colv=TRUE, Rowv=TRUE,dendrogram="both", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none",margins = c(3, 17))
dev.off()


# 同理计算KO vs WT间的差异ko
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
dev.off()
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
# 定义常用颜色 Defined color with transparent
alpha = .7
c_yellow =          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
c_blue =            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
c_orange =          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
c_green =           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
c_dark_green =      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_very_dark_green = rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_sea_green =       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
c_black =           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
c_grey =            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
c_dark_brown =      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_red =             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_red =        rgb(255 / 255, 130 / 255,   0 / 255, alpha)

# 三元图函数，无须理解直接调用即可 Function of ternary plot
tern_e=function (x, scale = 1, dimnames = NULL, dimnames_position = c("corner",
                                                                      "edge", "none"), dimnames_color = "black", id = NULL, id_color = "black",
                 coordinates = FALSE, grid = TRUE, grid_color = "gray", labels = c("inside",
                                                                                   "outside", "none"), labels_color = "darkgray", border = "black",
                 bg = "white", pch = 19, cex = 1, prop_size = FALSE, col = "red",
                 main = "ternary plot", newpage = TRUE, pop = TRUE, ...)
{
  labels = match.arg(labels)
  if (grid == TRUE)
    grid = "dotted"
  if (coordinates)
    id = paste("(", round(x[, 1] * scale, 1), ",", round(x[,
                                                           2] * scale, 1), ",", round(x[, 3] * scale, 1), ")",
               sep = "")
  dimnames_position = match.arg(dimnames_position)
  if (is.null(dimnames) && dimnames_position != "none")
    dimnames = colnames(x)
  if (is.logical(prop_size) && prop_size)
    prop_size = 3
  if (ncol(x) != 3)
    stop("Need a matrix with 3 columns")
  if (any(x < 0))
    stop("X must be non-negative")
  s = rowSums(x)
  if (any(s <= 0))
    stop("each row of X must have a positive sum")
  x = x/s
  top = sqrt(3)/2
  if (newpage)
    grid.newpage()
  xlim = c(-0.03, 1.03)
  ylim = c(-1, top)
  pushViewport(viewport(width = unit(1, "snpc")))
  if (!is.null(main))
    grid.text(main, y = 0.9, gp = gpar(fontsize = 18, fontstyle = 1))
  pushViewport(viewport(width = 0.8, height = 0.8, xscale = xlim,
                        yscale = ylim, name = "plot"))
  eps = 0.01
  grid.polygon(c(0, 0.5, 1), c(0, top, 0), gp = gpar(fill = bg,
                                                     col = border), ...)
  if (dimnames_position == "corner") {
    grid.text(x = c(0, 1, 0.5), y = c(-0.02, -0.02, top +
                                        0.02), label = dimnames, gp = gpar(fontsize = 12))
  }
  if (dimnames_position == "edge") {
    shift = eps * if (labels == "outside")
      8
    else 0
    grid.text(x = 0.25 - 2 * eps - shift, y = 0.5 * top +
                shift, label = dimnames[2], rot = 60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.75 + 3 * eps + shift, y = 0.5 * top +
                shift, label = dimnames[1], rot = -60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.5, y = -0.02 - shift, label = dimnames[3],
              gp = gpar(col = dimnames_color))
  }
  if (is.character(grid))
    for (i in 1:4 * 0.2) {
      grid.lines(c(1 - i, (1 - i)/2), c(0, 1 - i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(1 - i, 1 - i + i/2), c(0, i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(i/2, 1 - i + i/2), c(i, i) * top, gp = gpar(lty = grid,
                                                               col = grid_color))
      if (labels == "inside") {
        grid.text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 *
                    top, label = i * scale, gp = gpar(col = labels_color),
                  rot = 120)
        grid.text(x = 1 - i + i/4 + eps, y = i/2 * top -
                    eps, label = (1 - i) * scale, gp = gpar(col = labels_color),
                  rot = -120)
        grid.text(x = 0.5, y = i * top + eps, label = i *
                    scale, gp = gpar(col = labels_color))
      }
      if (labels == "outside") {
        grid.text(x = (1 - i)/2 - 6 * eps, y = (1 - i) *
                    top, label = (1 - i) * scale, gp = gpar(col = labels_color))
        grid.text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 -
                                                      i) * top + 5 * eps, label = i * scale, rot = -120,
                  gp = gpar(col = labels_color))
        grid.text(x = i + eps, y = -0.05, label = (1 -
                                                     i) * scale, vjust = 1, rot = 120, gp = gpar(col = labels_color))
      }
    }
  xp = x[, 2] + x[, 3]/2
  yp = x[, 3] * top
  size = unit(if (prop_size)
    #emiel inserted this code. x are proportions per row.  x*s is original data matrix. s = rowsums of original data matrix (x*s)
    prop_size * rowSums(x*x*s) / max(  rowSums(x*x*s) )
    #prop_size * rowSums(    (x*s) * ((x*s)/s)) / max(  rowSums(    (x*s) * ((x*s)/s)) )
    else cex, "lines")
  grid.points(xp, yp, pch = pch, gp = gpar(col = col), default.units = "snpc",
              size = size, ...)
  if (!is.null(id))
    grid.text(x = xp, y = unit(yp - 0.015, "snpc") - 0.5 *
                size, label = as.character(id), gp = gpar(col = id_color,
                                                          cex = cex))
  if (pop)
    popViewport(2)
  else upViewport(2)
}
# merge group to mean
## 按样品名合并实验组与转置的ko
mat_t2 = merge(sub_design[c("genotype")], t(norm), by="row.names")[,-1]
## 按实验设计求组平均值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean 
# 重新转载并去除组名
per3=t(mat_mean[,-1])
colnames(per3) = mat_mean$genotype
per3=as.data.frame(per3[rowSums(per3)>0,]) # remove all 0 ko
#per3=per3[,tern] # reorder per3 as input
color=c(c_green,c_orange,c_red,c_grey) 


# 两底角相对于顶点显著富集的ko，分共有和特有，类似维恩图
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
per3lg=log2(per3[,1:3]*10000+1) # 对数变换，差ko千分比的差距，点大小更均匀
pdf(file=paste("KO3_ter_",tern[1],tern[2],tern[3],"venn.pdf", sep=""), height = 8, width = 8)
tern_e(per3lg[,1:3], prop=T, col=per3$color, grid_color="black", labels_color="transparent", pch=19, main="Tenary Plot")
dev.off()
# KO差异小，三元图分不开


## Manhattan图展示差异ko和Taxonomy

# 读取taxonomy，并添加各列名称
taxonomy = read.delim("rep_seqs_tax.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")

# 标记差异ko类型
x$level = as.factor(ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")))
x$ko = rownames(x)
# 转换Pvalue为负对数
x$neglogp = -log(x$PValue)

# Taxonomy排序，并筛选ko表中存在的
library(dplyr)
taxonomy$id=rownames(taxonomy)
taxonomy = arrange(taxonomy, phylum, class, order, family, genus, species)
rownames(taxonomy) = taxonomy$id
idx = rownames(taxonomy) %in% x$ko
tax = taxonomy[idx, ] # subset taxonomy from used ko

# 手动筛选显著的组
x = x[rownames(tax), ] # reorder according to tax
x$tax = gsub("p__","",tax$phylum,perl=TRUE) 
top_phylum=c("Bacteroidetes","Firmicutes","Planctomycetes","Proteobacteria","Verrucomicrobia")
x[!(x$tax %in% top_phylum),]$tax = "Low Abundance" # no level can get value
# 设置各类的level对应顺序
x$ko = factor(x$ko, levels=x$ko)   # set x order
x$level = factor(x$level, levels=c("enriched","depleted","nosig"))
levels(x$tax)=c(top_phylum,"Low Abundance")
# 调整Y轴范围更美观
x[x$neglogp>15,]$neglogp  = 15

# Manhattan plot
## 添加显著阈值线
FDR = min(x$neglogp[x$level=="depleted"])
library(ggplot2)
p = ggplot(x, aes(x=ko, y=neglogp, color=tax, size=logCPM, shape=level)) +
  geom_point(alpha=.7) + 
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(17, 25, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="ko", y="-loge(P)") +
  theme(axis.line.x=element_line(size=.5, colour="black"),legend.key=element_blank(),axis.line.y=element_line(size=.5, colour="black"),panel.background=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="top")
p
ggsave(file=paste("man_ko.pdf", sep=""), p, width = 10, height = 3, useDingbats=F)
ggsave(file=paste("man_ko.png", sep=""), p, width = 10, height = 3)



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
ggsave(file=paste("vol_ko_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = 8, height = 8)
ggsave(file=paste("vol_ko_", sampleA, "vs", sampleB, ".png", sep=""), p, width = 8, height = 5)







1. 首先整合各个数据的FPKM表达文件到一个文件里
WGCNA代码：
```console
goodSample <- c("FPKM.H10m.1","FPKM.H10m.2","FPKM.H10m.3","FPKM.H1h.1","FPKM.H1h.2","FPKM.H1h.3","FPKM.H30m.1","FPKM.H30m.2","FPKM.H30m.3","FPKM.H4h.1","FPKM.H4h.2","FPKM.H4h.3","FPKM.S10m.1","FPKM.S10m.2","FPKM.S10m.3","FPKM.S1h.1","FPKM.S1h.2","FPKM.S1h.3","FPKM.S30m.1","FPKM.S30m.2","FPKM.S30m.3","FPKM.S4h.1","FPKM.S4h.2","FPKM.S4h.3")
cot_goodSample <- total_FPKM[,goodSample]
cot_goodSample_FPKM <- total_FPKM[,goodSample]
cot_goodSample_phone <- cot_pheno[goodSample,]
datExpr <- cot_goodSample_FPKM
```

# 筛选方法：mad>1且top5000
```console
WGCNA_matrix = t(datExpr[order(apply(datExpr,1,mad),decreasing = T)[1:5000],])
datExpr_filted <- WGCNA_matrix
```

检查离群样品
如果你的是TCGA上下载的数据，几百个临床样品，这一步是必须的，可能会有离群值，需要把它们去掉。
```console
library(WGCNA)
gsg = goodSamplesGenes(datExpr_filted, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr_filted), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
save(datExpr_filted, file = "datExpr_filted_hclust.RData")
```

去除离群值代码：
#离群值height大于8000
```console
clust <- cutreeStatic(sampleTree, cutHeight = 8000, minSize = 10) 
table(clust)
keepSamples <- (clust==1)
datExpr_remove <- datExpr_filted[keepSamples, ]
datExpr_filted <- datExpr_remove
```

确定软阈值
# Constructing a weighted gene network entails the choice of the soft thresholding power to which coexpression similarity is raised to calculate adjacency.
# Set up a bunch of power gradients(设定一些列power梯度)
```console
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr_filted, powerVector = powers, verbose = 5) #this step will take some time # The "sft" object contains the network characteristics that calculated for each power value(在sft这个对象中保存了每个power值计算出来的值)
str(sft) # make figures to show the soft thresholding power
par(mfrow = c(1,2))
cex1 = 0.9  # x axis is the Soft threshold (power)，y axis is evaluation parameters for scale-free networks(纵轴是无标度网络的评估参数)
# The higher R^2 ,the more scale-free the network is.(数值越高，网络越符合无标度特征 (non-scale))
# "FitIndices" stores the characteristics of each network corresponding to its power. (fitIndices保存了每个power对应的网络的特征)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# R^2 value(h) usually around 0.9. But if there's some big changes between your samples, R^2 value will be lower. 
# In my case, I use 0.7 to set the criteria.
abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# System will recommend the best power value to you. 
# But you can not totally trust it, depend on your data.
sft$powerEstimate
```

网络构建
# power: 上面得到的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)（4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个）
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
```console
net = blockwiseModules(datExpr_filted, power = sft$powerEstimate, maxBlockSize = 5000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "RNA_seq_datExpr_filted_TOM",
                       verbose = 3)
# You can check how many modules and how many genes in each module
# "0" represent "grey" module
table(net$colors)
 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
  40 1014  959  858  430  265  263  248  230  195  159  147   77   69   46 
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

可视化基因网络
上面我们构建了网络，还可以将基因共表达相似度矩阵绘制成热图，和模块的聚类一起展示：
# shown in a heat map according to 1-TOM value (distance between genes)
# This step you could use 5000 genes to make heatmap figure, or you could use randomly some genes from each module
```console
geneTree = net$dendrograms[[1]]
moduleColors = labels2colors(net$colors)
dissTOM = 1 - TOMsimilarityFromExpr(
  datExpr_filted,
  power = 14)
plotTOM = dissTOM ^ 7
diag(plotTOM) = NA
# The more genes you choose, the longer time this step will take:
# Actually, this step take a long time. This is why we need to filt genes before WGCNA analysis.
#可视化这一步，你可以选择用5000个基因来画，但是很费时间，我这个破电脑大概１０多分钟吧
#当然你也可以选择只画一部分基因，具体如何操作，请见：https://www.jianshu.com/p/3618f5ff3eb0
TOMplot(
  plotTOM,
  geneTree,
  moduleColors,
  main = "Network heatmap plot, all genes"
)
#不要急，可以用下面的代码把颜色转过来：
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)
```











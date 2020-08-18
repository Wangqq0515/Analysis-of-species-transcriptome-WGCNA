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
  power = sft$powerEstimate)
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

保存你的结果：
```console
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors) #查看所有模块里基因的数量
moduleColors
      black        blue       brown        cyan       green greenyellow        grey 
        248         959         858          46         265         147          40 
    magenta        pink      purple         red      salmon         tan   turquoise 
        195         230         159         263          69          77        1014 
     yellow 
        430 
MEs = net$MEs
MEs_col = MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs <- MEs_col
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "RNA_seq_FPKM_filted_networkConstruction.RData")
```

想查看每一个模块里的基因吗？代码：
```console
unique(moduleColors)
 [1] "blue"      "green"     "brown"     "grey"      "turquoise"
 [6] "red"       "black"     "yellow"    "pink"      "magenta"  
head(colnames(datExpr_filted)[moduleColors=="black"])

#　如果不想看前几个，想调出来全部的，可以这样：
brown_module <- colnames(datExpr_filted)[moduleColors=="brown"]
write.table(data.frame(brown_module),file = "brown_module.txt")
```

选一个你感兴趣的模块，看看基因的表达情况
```console
# You could extract the FPKM of the genes in the intested module, make heatmap of gene expression
which.module="black"
plotMat(t(scale(datExpr_filted[,moduleColors==which.module ]) ),
    nrgcols=30,
    rlabels=F,
    rcols=which.module,
    main=which.module,
    cex.main=2
  )
# You can also displaying module heatmap and the eigengene
sizeGrWindow(8,7)
which.module="black"
ME=MEs[,paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr_filted[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, 
        col=which.module, 
        main="", 
        cex.main=2, 
        ylab="eigengene expression",
        xlab="samples")
```

模块之间的相关性
```console
MEs = net$MEs
MEs_col = MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
MEs <- MEs_col
# The correlation map of each module was obtained by clustering according to the gene expression
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
```

模块与性状/实验条件/临床特征的相关性
对于我来说，这一步是最关键的，我想知道我的实验处理条件和哪一个或哪几个模块最相关：
#先处理一下样品属性的matrix，或者你有临床属性也可以
```console
datTraits <- as.data.frame(datTraits)
design = model.matrix(~0+ datTraits$sample_time)
c = as.factor(datTraits$sample_time)
levels(c)
colnames(design) = levels(c)
moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr_filted, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
nSamples = nrow(datExpr_filted)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot (correlation between module and conditions)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
```

你可以指定一个感兴趣的表型/condition，可以得到与这个形状相关性最高的模块：
# You can get the most relevant module by specifying a condition of interest
```console
which.trait <- "S1h"
moduleTraitCor[, which.trait]
    MEblack      MEpink      MEblue    MEyellow     MEbrown MEturquoise 
  0.6348693   0.6386766  -0.2192051   0.2993171   0.1302572  -0.2362041 
  MEmagenta     MEgreen       MEred      MEgrey 
 -0.1391687   0.0269631  -0.1231228   0.4544833 
```

对某一个性状/条件/临床特征相对应的某个模块进行具体分析
# 对p/brown(条件/模块)具体分析,模块内基因与表型数据关联:
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达值算出相关系数
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要

# Firstly, calculate the correlation matrix between module and genes.(首先计算模块与基因的相关性矩)
# names (colors) of the modules
```console
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr_filted, MEs, use = "p"))
# Calculate the Pearson correlation coefficient matrix of each module and its genes(算出每个模块跟基因的皮尔森相关系数矩)
# MEs is the value of each module in each sample.(MEs是每个模块在每个样本里面的)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Secondly, calculate the correlation matrix between conditions and genes. (再计算性状与基因的相关性矩)
# Only continuous properties can be computed. If the variables are discrete, the matrix is converted to 0-1 when the sample table is constructed.(只有连续型性状才能只有计算,如果是离散变量，在构建样品表时就转为0-1矩阵)
# Here, the variable whether or not it belongs to the P condition is numeralized with 0 and 1.(这里把是否属于 P 实验条件这个变量进行数值化，０和１表示)
P = as.data.frame(design[,1]) # choose an interested condition!!
names(P) = "proliferating"
geneTraitSignificance = as.data.frame(cor(datExpr_filted, P, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(P), sep="")
names(GSPvalue) = paste("p.GS.", names(P), sep="")

# Then, combine aboved two correlation matrixes(最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析)
module = "brown" # choose interested module
column = match(module, modNames)
# get the genes in the interested module(获取模块内的基因)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# Genes that are highly correlated with conditions are also highly associated with modules (与性状高度相关的基因，也是与性状相关的模型的关键基因)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for proliferating",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 
```

导出网络到Cytoscape
可以导出指定module对应的基因共表达网络：
```console
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr_filted, power = sft$powerEstimate)
# Select modules
module = "brown"
# Select module probes
probes = colnames(datExpr_filted) #基因名
inModule = (moduleColors == module)
modProbes = probes[inModule]
# modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#export to VisANT
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,                               
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
```

筛选hub基因
```console
#（1）Calculate Intramodular connectivity
moduleColors <- labels2colors(net$colors)
connet=abs(cor(datExpr_filted,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)
# (2) calculate relationship between gene significance and intramodular connectivity
#which.module="brown"
P= as.data.frame(design[,1]) # change specific 
names(P) = "Proliferating"
GS1 = as.numeric(cor(P,datExpr_filted, use="p"))
GeneSignificance=abs(GS1)
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]]
  restrict1 = (moduleColors == whichmodule)
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
```

下面的hub基因筛选比较重要：
这里注意筛选标准：abs(GS1)>0.9 可以根据实际情况调整;　abs(datKME$MM.black)>0.8 (至少大于 >0.8)
```console
#(3）Calculate the connectivity of all genes in the module.Screen the hub genes.
# abs(GS1)> 0.9 (could be adjusted depend on your data)
# abs(datKME$MM.black)>0.8 (at least more than 0.8)
# Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr_filted, MEs, outputColumnName="MM.")
head(datKME)
write.csv(datKME,file = "datKME.csv")
# Finding genes with high gene significance and high intramodular connectivity in interesting modules 
FilterGenes= abs(GS1)> 0.7 & abs(datKME$MM.brown) > 0.7
table(FilterGenes)
x = datKME[(GS1>0.8 & datKME$MM.brown >0.8),]
View(x)
write.csv(x, file="module_brown_hub_gene_GS0.8.csv")
```
